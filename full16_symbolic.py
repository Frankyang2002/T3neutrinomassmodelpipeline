"""Canonical general-interaction 16-model *symbolic* study entry point.

    python -X utf8 full16_symbolic.py
    python -X utf8 full16_symbolic.py --reports-only
    python -X utf8 full16_symbolic.py --model A_-2 --study general_symbolic_test

This intentionally does not call the restricted sequential-matching pipeline.
"""
from scripts.RunT3GeneralSymbolicStudy import main

if __name__ == "__main__":
    raise SystemExit(main())
