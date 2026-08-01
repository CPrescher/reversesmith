"""Line-oriented MACE worker for rsmith.

The parent process sends one JSON request per line and receives one JSON
response per line. This keeps the Rust/Python boundary stable while leaving
PyTorch, MACE, and ASE in the Python environment where they are normally used.
"""

import json
import sys
import traceback
from contextlib import redirect_stdout


atoms = None


def respond(payload):
    print(json.dumps(payload), flush=True)


def require_atoms():
    if atoms is None:
        raise RuntimeError("MACE worker has not been initialized")
    return atoms


def handle_init(request):
    global atoms

    with redirect_stdout(sys.stderr):
        try:
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

        atoms = Atoms(
            symbols=request["species"],
            positions=request["positions"],
            cell=request["box"],
            pbc=True,
        )
        atoms.calc = calc


def handle_energy():
    with redirect_stdout(sys.stderr):
        return float(require_atoms().get_potential_energy())


def handle_move(request):
    current_atoms = require_atoms()
    current_atoms.positions[int(request["atom"])] = request["position"]


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
            elif cmd == "move":
                handle_move(request)
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
