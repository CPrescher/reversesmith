#!/usr/bin/env python3
"""Download a small MIT-licensed MACE-MP checkpoint for integration tests."""

import inspect
import sys


def main() -> int:
    try:
        from mace.calculators import mace_mp
    except Exception as exc:
        print(
            "failed to import mace.calculators; install mace-torch first",
            file=sys.stderr,
        )
        print(str(exc), file=sys.stderr)
        return 1

    module = inspect.getmodule(mace_mp)
    path = module.download_mace_mp_checkpoint("small")
    print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
