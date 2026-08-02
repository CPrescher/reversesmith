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
r_max = None
num_interactions = None
symbols = None
positions = None
box = None


def respond(payload):
    print(json.dumps(payload), flush=True)


def require_atoms():
    if atoms is None:
        raise RuntimeError("MACE worker has not been initialized")
    return atoms


def handle_init(request):
    global atoms, calc, np, r_max, num_interactions, symbols, positions, box

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


def handle_energy():
    with redirect_stdout(sys.stderr):
        return float(require_atoms().get_potential_energy())


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
    atom = int(request["atom"])
    old_position = np.asarray(request["old_position"], dtype=float)
    new_position = np.asarray(request["new_position"], dtype=float)
    central = central_atom_indices(atom, old_position, new_position)
    before = local_energy_sum(central, atom, old_position)
    handle_move({"atom": atom, "position": new_position})
    after = local_energy_sum(central, atom, new_position)
    return {
        "delta": after - before,
        "central_atoms": len(central),
    }


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
    cluster.calc = calc
    with redirect_stdout(sys.stderr):
        cluster.get_potential_energy()
    energies = np.asarray(calc.results["energies"], dtype=float)
    return float(np.sum(energies[central_cluster_indices]))


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
            elif cmd == "shutdown":
                respond({"ok": True})
                break
            else:
                respond({"ok": False, "error": f"unknown command: {cmd}"})
        except Exception:
            respond({"ok": False, "error": traceback.format_exc()})


if __name__ == "__main__":
    main()
