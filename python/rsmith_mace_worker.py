"""Line-oriented MACE worker for rsmith.

The parent process sends one JSON request per line and receives one JSON
response per line. This keeps the Rust/Python boundary stable while leaving
PyTorch, MACE, and ASE in the Python environment where they are normally used.
"""

import json
import math
import sys
import traceback
from contextlib import redirect_stdout


atoms = None
calc = None
np = None
torch = None
r_max = None
num_interactions = None
symbols = None
positions = None
box = None
delta_mode = "full"
accepted_node_energies = None
accepted_total_energy = None
pending_local_update = None


def respond(payload):
    print(json.dumps(payload), flush=True)


def require_atoms():
    if atoms is None:
        raise RuntimeError("MACE worker has not been initialized")
    return atoms


def handle_init(request):
    global atoms, calc, np, torch, r_max, num_interactions, symbols, positions, box
    global delta_mode, accepted_node_energies, accepted_total_energy
    global pending_local_update

    with redirect_stdout(sys.stderr):
        try:
            import numpy
            import torch
            from ase import Atoms
            from mace.calculators import MACECalculator
        except Exception as exc:
            raise RuntimeError(
                "failed to import MACE dependencies; install mace-torch, torch, and ase "
                "in the selected Python environment"
            ) from exc

        torch_threads = request.get("torch_threads")
        if torch_threads is not None:
            torch.set_num_threads(int(torch_threads))
            try:
                torch.set_num_interop_threads(int(torch_threads))
            except RuntimeError:
                pass

        calc = MACECalculator(
            model_paths=request["model"],
            device=request.get("device", "cpu"),
        )

        np = numpy
        symbols = list(request["species"])
        positions = np.asarray(request["positions"], dtype=float)
        box = np.asarray(request["box"], dtype=float)
        r_max = float(calc.r_max)
        raw_num_interactions = getattr(calc.models[0], "num_interactions", None)
        if raw_num_interactions is not None:
            try:
                num_interactions = int(raw_num_interactions.item())
            except AttributeError:
                num_interactions = int(raw_num_interactions)
        else:
            num_interactions = len(getattr(calc.models[0], "interactions"))

        atoms = Atoms(
            symbols=symbols,
            positions=positions,
            cell=box,
            pbc=True,
        )
        atoms.calc = calc
        delta_mode = request.get("delta", "full")
        pending_local_update = None
        if delta_mode == "local":
            accepted_total_energy, accepted_node_energies = energy_only(atoms)
        else:
            accepted_total_energy = None
            accepted_node_energies = None


def handle_energy():
    if delta_mode == "local" and accepted_total_energy is not None:
        if pending_local_update is None:
            return float(accepted_total_energy)
        return float(accepted_total_energy + pending_local_update[2])
    with redirect_stdout(sys.stderr):
        energy, _ = energy_only(require_atoms())
        return energy


def handle_move(request):
    global positions
    current_atoms = require_atoms()
    atom = int(request["atom"])
    position = np.asarray(request["position"], dtype=float)
    current_atoms.positions[atom] = position
    positions[atom] = position


def handle_metadata():
    require_atoms()
    return {
        "r_max": r_max,
        "num_interactions": num_interactions,
        "local_central_radius": local_radius(),
        "local_context_radius": local_radius(),
    }


def handle_local_delta(request):
    global pending_local_update
    atom = int(request["atom"])
    old_position = np.asarray(request["old_position"], dtype=float)
    new_position = np.asarray(request["new_position"], dtype=float)
    central = request.get("central_atoms")
    before_cluster = request.get("before_cluster")
    after_cluster = request.get("after_cluster")
    if central is None or after_cluster is None:
        # Compatibility path for older parents and custom site integrations.
        central = central_atom_indices(atom, old_position, new_position)
        before, _ = local_energy_sum(central, atom, old_position)
    else:
        central = [int(index) for index in central]
        if accepted_node_energies is None:
            if before_cluster is None:
                raise RuntimeError(
                    "precomputed local MACE request omitted the before cluster "
                    "without an accepted energy cache"
                )
            before, _ = local_energy_sum_precomputed(
                atom, old_position, before_cluster
            )
        else:
            before = float(np.sum(accepted_node_energies[central]))
    handle_move({"atom": atom, "position": new_position})
    if after_cluster is None:
        after, _ = local_energy_sum(central, atom, new_position)
        new_central_energies = None
    else:
        after, new_central_energies = local_energy_sum_precomputed(
            atom, new_position, after_cluster
        )
    delta = after - before
    if new_central_energies is not None and accepted_node_energies is not None:
        pending_local_update = (central, new_central_energies, delta)
    return {
        "delta": delta,
        "central_atoms": len(central),
        "cluster_atoms": len(after_cluster["atoms"])
        if after_cluster is not None
        else None,
    }


def handle_accept_local():
    global accepted_total_energy, pending_local_update
    if pending_local_update is None:
        raise RuntimeError("local MACE accept requested without a pending trial")
    central, new_central_energies, delta = pending_local_update
    accepted_node_energies[central] = new_central_energies
    accepted_total_energy += delta
    pending_local_update = None


def handle_reject_local(request):
    global pending_local_update
    handle_move({"atom": request["atom"], "position": request["old_position"]})
    pending_local_update = None


def local_radius():
    return float(num_interactions) * float(r_max)


def minimum_image_distance(a, b):
    delta = np.asarray(a, dtype=float) - np.asarray(b, dtype=float)
    delta -= box * np.rint(delta / box)
    return float(np.linalg.norm(delta))


def central_atom_indices(atom, old_position, new_position):
    radius = local_radius()
    central = [atom]
    for idx, position in enumerate(positions):
        if idx == atom:
            continue
        if (
            minimum_image_distance(position, old_position) <= radius
            or minimum_image_distance(position, new_position) <= radius
        ):
            central.append(idx)
    return central


def local_energy_sum(central, moved_atom, moved_position):
    from ase import Atoms

    radius = local_radius()
    eval_positions = positions.copy()
    eval_positions[moved_atom] = moved_position
    central_set = set(central)
    central_positions = {idx: eval_positions[idx].copy() for idx in central}

    shift_ranges = [
        range(-int(math.ceil(radius / length)) - 1, int(math.ceil(radius / length)) + 2)
        for length in box
    ]
    cluster_symbols = []
    cluster_positions = []
    central_cluster_indices = []

    for atom_idx, position in enumerate(eval_positions):
        for sx in shift_ranges[0]:
            for sy in shift_ranges[1]:
                for sz in shift_ranges[2]:
                    shift = np.asarray([sx, sy, sz], dtype=float) * box
                    image_position = position + shift
                    keep = False
                    for central_position in central_positions.values():
                        if np.linalg.norm(image_position - central_position) <= radius + 1.0e-8:
                            keep = True
                            break
                    if not keep:
                        continue
                    cluster_index = len(cluster_symbols)
                    cluster_symbols.append(symbols[atom_idx])
                    cluster_positions.append(image_position)
                    if atom_idx in central_set and sx == 0 and sy == 0 and sz == 0:
                        central_cluster_indices.append(cluster_index)

    if len(central_cluster_indices) != len(central):
        raise RuntimeError(
            f"local MACE cluster missed central atoms: expected {len(central)}, "
            f"got {len(central_cluster_indices)}"
        )

    cluster = Atoms(
        symbols=cluster_symbols,
        positions=cluster_positions,
        pbc=False,
    )
    with redirect_stdout(sys.stderr):
        _, energies = energy_only(cluster)
    central_energies = energies[central_cluster_indices]
    return float(np.sum(central_energies)), central_energies


def local_energy_sum_precomputed(moved_atom, moved_position, cluster_spec):
    from ase import Atoms

    atom_indices = np.asarray(cluster_spec["atoms"], dtype=int)
    image_shifts = np.asarray(cluster_spec["shifts"], dtype=float)
    if atom_indices.ndim != 1 or image_shifts.shape != (len(atom_indices), 3):
        raise RuntimeError("invalid precomputed local MACE cluster shape")
    if np.any(atom_indices < 0) or np.any(atom_indices >= len(positions)):
        raise RuntimeError("precomputed local MACE cluster atom index is out of range")

    eval_positions = positions.copy()
    eval_positions[moved_atom] = moved_position
    cluster_positions = eval_positions[atom_indices] + image_shifts * box
    cluster_symbols = [symbols[index] for index in atom_indices]
    central_cluster_indices = np.asarray(cluster_spec["central_indices"], dtype=int)
    if (
        central_cluster_indices.ndim != 1
        or np.any(central_cluster_indices < 0)
        or np.any(central_cluster_indices >= len(atom_indices))
    ):
        raise RuntimeError("invalid precomputed local MACE central indices")

    cluster = Atoms(
        symbols=cluster_symbols,
        positions=cluster_positions,
        pbc=False,
    )
    with redirect_stdout(sys.stderr):
        _, energies = energy_only(cluster)
    central_energies = energies[central_cluster_indices]
    return float(np.sum(central_energies)), central_energies


def energy_only(eval_atoms):
    """Evaluate energies without the force/stress autograd used by ASE calculate()."""
    batch_base = calc._atoms_to_batch(eval_atoms)
    node_energies = []
    total_energies = []
    real_atom_count = len(eval_atoms)

    with torch.no_grad():
        for model in calc.models:
            batch = calc._clone_batch(batch_base)
            model_dtype = next(model.parameters()).dtype
            for key in batch.keys:
                value = batch[key]
                if torch.is_tensor(value) and torch.is_floating_point(value):
                    batch[key] = value.to(dtype=model_dtype)
            output = model(
                batch.to_dict(),
                compute_force=False,
                compute_stress=False,
                compute_virials=False,
                compute_displacement=False,
                compute_hessian=False,
            )
            node_energies.append(output["node_energy"][:real_atom_count].detach())
            total_energies.append(output["energy"].reshape(-1)[0].detach())

    unit_conversion = float(calc.energy_units_to_eV)
    mean_node_energies = (
        torch.stack(node_energies, dim=0).mean(dim=0).cpu().numpy()
        * unit_conversion
    )
    mean_total_energy = (
        torch.stack(total_energies, dim=0).mean().cpu().item() * unit_conversion
    )
    return float(mean_total_energy), np.asarray(mean_node_energies, dtype=float)


def main():
    for line in sys.stdin:
        try:
            request = json.loads(line)
            cmd = request.get("cmd")
            if cmd == "init":
                handle_init(request)
                respond({"ok": True})
            elif cmd == "energy":
                respond({"ok": True, "energy": handle_energy()})
            elif cmd == "metadata":
                response = {"ok": True}
                response.update(handle_metadata())
                respond(response)
            elif cmd == "move":
                handle_move(request)
                respond({"ok": True})
            elif cmd == "local_delta":
                response = {"ok": True}
                response.update(handle_local_delta(request))
                respond(response)
            elif cmd == "accept_local":
                handle_accept_local()
                respond({"ok": True})
            elif cmd == "reject_local":
                handle_reject_local(request)
                respond({"ok": True})
            elif cmd == "shutdown":
                respond({"ok": True})
                break
            else:
                respond({"ok": False, "error": f"unknown command: {cmd}"})
        except Exception:
            respond({"ok": False, "error": traceback.format_exc()})


if __name__ == "__main__":
    main()
