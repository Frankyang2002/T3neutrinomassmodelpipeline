"""Presentation-only Matplotlib save redirection: full -> display/full_studies.

Never changes input paths or numerical JSON/CSV. Use only while generating plots.
"""
from __future__ import annotations
from contextlib import contextmanager
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FULL = ROOT / 'Reports' / 'output' / 'full'
DISPLAY = ROOT / 'Reports' / 'output' / 'display' / 'full_studies'


def display_path(path):
    if not isinstance(path, (str, Path)):
        return path
    source = Path(path)
    try:
        relative = source.resolve().relative_to(FULL.resolve())
    except ValueError:
        return path
    # Keep original study/scenario/subtree, so labels and provenance stay clear.
    return DISPLAY / relative


@contextmanager
def redirect_matplotlib_saves():
    from matplotlib.figure import Figure
    original = Figure.savefig
    def relocated(self, fname, *args, **kwargs):
        output = display_path(fname)
        if isinstance(output, Path):
            output.parent.mkdir(parents=True, exist_ok=True)
        return original(self, output, *args, **kwargs)
    Figure.savefig = relocated
    try:
        yield
    finally:
        Figure.savefig = original
