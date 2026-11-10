#!/usr/bin/env python3
"""Run the unchanged Batch 15 CG-enabled Batch 13 module with a path alias.

The current SourceBoundary dataclass exposes `.path`, while the Batch 15
full-tensor code reads `.seed_path`. This entrypoint supplies an equivalent
read-only property for the lifetime of the process. No physics inputs,
coefficients, existing source files, or repository outputs are modified.

Usage (from project root):
    python -m scripts.RunT3MasterWilsonBatch13CGFix --all --debug
All command-line arguments are forwarded to the original Batch 13 runner.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from Numerical.joint_scalar_matching.fermion_boundary_bridge import SourceBoundary

if not hasattr(SourceBoundary, 'seed_path'):
    if 'path' not in SourceBoundary.__dataclass_fields__:
        raise RuntimeError(
            'Incompatible SourceBoundary: neither seed_path nor the verified path field exists'
        )
    SourceBoundary.seed_path = property(lambda self: self.path)

from scripts.RunT3MasterWilsonBatch13 import main

if __name__ == '__main__':
    raise SystemExit(main())
