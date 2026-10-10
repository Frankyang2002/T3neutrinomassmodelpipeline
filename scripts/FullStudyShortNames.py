"""Full study using short *study directory* names instead of dimensional keys.

The underlying UV and EFT physical run names, representation objects and
numerical JSON contents are unchanged. A full study is expensive; use only
when new full-run calculations are actually desired.

Usage: python scripts/FullStudyShortNames.py
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from Numerical.orchestration.ModelNames import model_label
import studies.FullT3Study as study
import pipeline


def main() -> int:
    # FullT3Study._model_key is used for study subdirectories and generated
    # benchmark configuration filenames, not for defining representations.
    study._model_key = model_label
    sys.argv = [str(ROOT / 'pipeline.py'), '--full', *sys.argv[1:]]
    return pipeline.main()


if __name__ == '__main__':
    raise SystemExit(main())
