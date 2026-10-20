"""Execute an existing Python plot module with presentation-only model labels.

Project path: tools/PlotDisplayWrapper.py

Usage: python -m tools.PlotDisplayWrapper Numerical.plotting.PlotFittedSensitivityFigures --input ...
Used by RebuildCurrentReports.py to handle its separate plotting process.
"""
from __future__ import annotations

import runpy
import sys

from Numerical.plotting.DisplayModelLabels import presentation_labels_on_save


def main() -> int:
    if len(sys.argv) < 2:
        raise SystemExit('Usage: python -m tools.PlotDisplayWrapper <module> [module options]')
    module = sys.argv[1]
    if module != 'Numerical.plotting.PlotFittedSensitivityFigures':
        raise SystemExit(f'Unrecognised presentation-only plotting module: {module}')
    sys.argv = [module] + sys.argv[2:]
    with presentation_labels_on_save():
        runpy.run_module(module, run_name='__main__', alter_sys=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
