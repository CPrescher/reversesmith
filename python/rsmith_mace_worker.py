"""Line-oriented MACE worker for rsmith.

The parent process sends one JSON request per line and receives one JSON
response per line. This keeps the Rust/Python boundary stable while leaving
PyTorch, MACE, and ASE in the Python environment where they are normally used.
"""

import json
import math
import sys
import time
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
accepted_model_layers = None
pending_local_update = None


def respond(payload):
    print(json.dumps(payload), flush=True)


def require_atoms():
    if atoms is None:
        raise RuntimeError("MACE worker has not been initialized")
    return atoms


def validate_incremental_models():
    from mace.modules.models import MACE, ScaleShiftMACE

    for model in calc.models:
        if not isinstance(model, (MACE, ScaleShiftMACE)):
            raise RuntimeError(
                "incremental mode supports only standard MACE and ScaleShiftMACE models"
            )
        if hasattr(model, "pair_repulsion") or hasattr(model, "joint_embedding"):
            raise RuntimeError(
                "incremental mode does not yet support pair repulsion or joint embeddings"
            )
        if any(hasattr(interaction, "conv_fusion") for interaction in model.interactions):
            raise RuntimeError(
                "incremental mode does not yet support fused interaction kernels"
            )


def handle_init(request):
    global atoms, calc, np, torch, r_max, num_interactions, symbols, positions, box
    global delta_mode, accepted_node_energies, accepted_total_energy
    global accepted_model_layers, pending_local_update

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
            default_dtype=request.get("dtype") or "",
            compile_mode=request.get("compile_mode"),
        )
        if request.get("delta") == "incremental":
            validate_incremental_models()

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
        if delta_mode in ("local", "incremental"):
            (
                accepted_total_energy,
                accepted_node_energies,
                accepted_model_layers,
            ) = energy_only(atoms, capture_layers=delta_mode == "incremental")
        else:
            accepted_total_energy = None
            accepted_node_energies = None
            accepted_model_layers = None


def handle_energy():
    if delta_mode in ("local", "incremental") and accepted_total_energy is not None:
        if pending_local_update is None:
            return float(accepted_total_energy)
        return float(accepted_total_energy + pending_local_update[2])
    with redirect_stdout(sys.stderr):
        energy, _, _ = energy_only(require_atoms())
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
        "local_context_radius": r_max
        if delta_mode == "incremental"
        else local_radius(),
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
        incremental_profile = None
    else:
        if delta_mode == "incremental":
            after, new_central_energies, layer_updates, incremental_profile = (
                incremental_energy_sum_precomputed(
                    atom, old_position, new_position, after_cluster
                )
            )
        else:
            after, new_central_energies = local_energy_sum_precomputed(
                atom, new_position, after_cluster
            )
            layer_updates = None
            incremental_profile = None
    delta = after - before
    if new_central_energies is not None and accepted_node_energies is not None:
        pending_local_update = (
            central,
            new_central_energies,
            delta,
            layer_updates,
        )
    response = {
        "delta": delta,
        "central_atoms": len(central),
        "cluster_atoms": len(after_cluster["atoms"])
        if after_cluster is not None
        else None,
    }
    if incremental_profile is not None:
        response.update(incremental_profile)
    return response


def handle_accept_local():
    global accepted_total_energy, accepted_model_layers, pending_local_update
    if pending_local_update is None:
        raise RuntimeError("local MACE accept requested without a pending trial")
    central, new_central_energies, delta, layer_updates = pending_local_update
    accepted_node_energies[central] = new_central_energies
    accepted_total_energy += delta
    if layer_updates is not None:
        for model_layers, model_updates in zip(
            accepted_model_layers, layer_updates
        ):
            for layer_cache, (base_atoms, new_features) in zip(
                model_layers, model_updates
            ):
                layer_cache[base_atoms] = new_features
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
        _, energies, _ = energy_only(cluster)
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
        _, energies, _ = energy_only(cluster)
    central_energies = energies[central_cluster_indices]
    return float(np.sum(central_energies)), central_energies


def incremental_energy_sum_precomputed(
    moved_atom, old_position, new_position, cluster_spec
):
    """Recompute only the message-passing causal cone of a moved atom."""
    from ase import Atoms
    from mace.modules.models import MACE, ScaleShiftMACE, prepare_graph

    if accepted_model_layers is None:
        raise RuntimeError("incremental MACE layer cache has not been initialized")

    started = time.perf_counter()
    atom_indices = np.asarray(cluster_spec["atoms"], dtype=int)
    image_shifts = np.asarray(cluster_spec["shifts"], dtype=float)
    central_cluster_indices = np.asarray(
        cluster_spec["central_indices"], dtype=int
    )
    if atom_indices.ndim != 1 or image_shifts.shape != (len(atom_indices), 3):
        raise RuntimeError("invalid incremental MACE cluster shape")
    if (
        central_cluster_indices.ndim != 1
        or np.any(central_cluster_indices < 0)
        or np.any(central_cluster_indices >= len(atom_indices))
    ):
        raise RuntimeError("invalid incremental MACE central indices")

    central_base_atoms = atom_indices[central_cluster_indices]
    if len(np.unique(central_base_atoms)) != len(central_base_atoms):
        raise RuntimeError("incremental MACE central atom images are not unique")
    central_lookup = {
        int(base_atom): int(cluster_atom)
        for base_atom, cluster_atom in zip(
            central_base_atoms, central_cluster_indices
        )
    }

    eval_positions = positions.copy()
    eval_positions[moved_atom] = new_position
    cluster_positions = eval_positions[atom_indices] + image_shifts * box
    cluster = Atoms(
        symbols=[symbols[index] for index in atom_indices],
        positions=cluster_positions,
        pbc=False,
    )
    assembled = time.perf_counter()
    batch_base = calc._atoms_to_batch(cluster)
    graph_built = time.perf_counter()
    new_node_energies = []
    all_model_updates = []

    def dirty_base_atoms(layer):
        radius = float(layer + 1) * float(r_max)
        return np.asarray(
            [
                int(base_atom)
                for base_atom in central_base_atoms
                if base_atom == moved_atom
                or minimum_image_distance(positions[base_atom], old_position)
                <= radius + 1.0e-8
                or minimum_image_distance(positions[base_atom], new_position)
                <= radius + 1.0e-8
            ],
            dtype=int,
        )

    with torch.no_grad():
        for model_index, model in enumerate(calc.models):
            if not isinstance(model, (MACE, ScaleShiftMACE)):
                raise RuntimeError(
                    "incremental mode supports only standard MACE and ScaleShiftMACE models"
                )
            if hasattr(model, "pair_repulsion") or hasattr(model, "joint_embedding"):
                raise RuntimeError(
                    "incremental mode does not yet support pair repulsion or joint embeddings"
                )
            if hasattr(model.interactions[0], "conv_fusion"):
                raise RuntimeError(
                    "incremental mode does not yet support fused interaction kernels"
                )

            batch = calc._clone_batch(batch_base)
            model_dtype = next(model.parameters()).dtype
            for key in batch.keys:
                value = batch[key]
                if torch.is_tensor(value) and torch.is_floating_point(value):
                    batch[key] = value.to(dtype=model_dtype)
            data = batch.to_dict()
            ctx = prepare_graph(
                data,
                compute_virials=False,
                compute_stress=False,
                compute_displacement=False,
                lammps_mliap=False,
            )
            edge_index = data["edge_index"]
            node_attrs = data["node_attrs"]
            node_heads = ctx.node_heads.to(torch.int64)
            cluster_atoms_tensor = torch.as_tensor(
                atom_indices, dtype=torch.long, device=edge_index.device
            )
            central_tensor = torch.as_tensor(
                central_cluster_indices,
                dtype=torch.long,
                device=edge_index.device,
            )

            node_feats = model.node_embedding(node_attrs)
            trial_features = []
            model_updates = []
            previous_dirty_bases = None
            previous_dirty_features = None

            for layer, (interaction, product) in enumerate(
                zip(model.interactions, model.products)
            ):
                dirty_bases = dirty_base_atoms(layer)
                dirty_cluster = torch.as_tensor(
                    [central_lookup[int(base)] for base in dirty_bases],
                    dtype=torch.long,
                    device=edge_index.device,
                )

                if layer > 0:
                    node_feats = accepted_model_layers[model_index][layer - 1][
                        cluster_atoms_tensor
                    ].clone()
                    for base, feature in zip(
                        previous_dirty_bases, previous_dirty_features
                    ):
                        node_feats[cluster_atoms_tensor == int(base)] = feature

                edge_mask = torch.isin(edge_index[1], dirty_cluster)
                selected_edges = edge_index[:, edge_mask]
                vectors = ctx.vectors[edge_mask]
                lengths = ctx.lengths[edge_mask]
                edge_attrs = model.spherical_harmonics(vectors)
                edge_feats, cutoff = model.radial_embedding(
                    lengths,
                    node_attrs,
                    selected_edges,
                    model.atomic_numbers,
                )
                messages, sc = interaction(
                    node_attrs=node_attrs,
                    node_feats=node_feats,
                    edge_attrs=edge_attrs,
                    edge_feats=edge_feats,
                    edge_index=selected_edges,
                    cutoff=cutoff,
                    first_layer=(layer == 0),
                    lammps_class=None,
                    lammps_natoms=(0, 0),
                )
                dirty_features = product(
                    node_feats=messages[dirty_cluster],
                    sc=sc[dirty_cluster] if sc is not None else None,
                    node_attrs=node_attrs[dirty_cluster],
                )
                trial_features.append(
                    (dirty_bases, dirty_cluster, dirty_features)
                )
                model_updates.append(
                    (dirty_bases, dirty_features.detach().clone())
                )
                previous_dirty_bases = dirty_bases
                previous_dirty_features = dirty_features

            central_heads = node_heads[central_tensor]
            central_arange = torch.arange(
                len(central_tensor), device=edge_index.device
            )
            node_e0 = model.atomic_energies_fn(node_attrs[central_tensor])[
                central_arange, central_heads
            ]
            readout_energies = []
            for readout_index, readout in enumerate(model.readouts):
                feature_layer = (
                    len(model.interactions) - 1
                    if len(model.readouts) == 1
                    else readout_index
                )
                central_features = accepted_model_layers[model_index][
                    feature_layer
                ][torch.as_tensor(central_base_atoms, device=edge_index.device)].clone()
                dirty_bases, _, dirty_features = trial_features[feature_layer]
                for base, feature in zip(dirty_bases, dirty_features):
                    central_features[
                        torch.as_tensor(
                            central_base_atoms == int(base),
                            device=edge_index.device,
                        )
                    ] = feature
                readout_energies.append(
                    readout(central_features, central_heads)[
                        central_arange, central_heads
                    ]
                )

            node_inter_energy = torch.sum(
                torch.stack(readout_energies, dim=0), dim=0
            )
            if isinstance(model, ScaleShiftMACE):
                node_inter_energy = model.scale_shift(
                    node_inter_energy, central_heads
                )
            model_node_energy = node_e0.double() + node_inter_energy.double()
            new_node_energies.append(model_node_energy.detach())
            all_model_updates.append(model_updates)

    unit_conversion = float(calc.energy_units_to_eV)
    mean_node_energies = (
        torch.stack(new_node_energies, dim=0).mean(dim=0).cpu().numpy()
        * unit_conversion
    )
    converted_updates = [
        [
            (
                torch.as_tensor(
                    base_atoms,
                    dtype=torch.long,
                    device=new_features.device,
                ),
                new_features,
            )
            for base_atoms, new_features in model_updates
        ]
        for model_updates in all_model_updates
    ]
    finished = time.perf_counter()
    return (
        float(np.sum(mean_node_energies)),
        np.asarray(mean_node_energies, dtype=float),
        converted_updates,
        {
            "incremental_assemble_ms": (assembled - started) * 1.0e3,
            "incremental_graph_ms": (graph_built - assembled) * 1.0e3,
            "incremental_model_ms": (finished - graph_built) * 1.0e3,
        },
    )


def energy_only(eval_atoms, capture_layers=False):
    """Evaluate energies without the force/stress autograd used by ASE calculate()."""
    batch_base = calc._atoms_to_batch(eval_atoms)
    node_energies = []
    total_energies = []
    model_layers = []
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
            if capture_layers:
                widths = [
                    int(product.linear.irreps_out.dim)
                    for product in model.products
                ]
                model_layers.append(
                    [
                        features.detach().clone()
                        for features in torch.split(
                            output["node_feats"][:real_atom_count], widths, dim=-1
                        )
                    ]
                )

    unit_conversion = float(calc.energy_units_to_eV)
    mean_node_energies = (
        torch.stack(node_energies, dim=0).mean(dim=0).cpu().numpy()
        * unit_conversion
    )
    mean_total_energy = (
        torch.stack(total_energies, dim=0).mean().cpu().item() * unit_conversion
    )
    return (
        float(mean_total_energy),
        np.asarray(mean_node_energies, dtype=float),
        model_layers if capture_layers else None,
    )


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
