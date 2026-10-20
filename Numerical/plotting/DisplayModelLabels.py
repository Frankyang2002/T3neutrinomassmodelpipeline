"""Presentation-only T3 model labels for the current saved-report workflow.

Project path: Numerical/plotting/DisplayModelLabels.py

The context manager temporarily normalises Matplotlib Text artists immediately
before saving a figure. Data arrays, model keys, source paths, and computations
are untouched. This also handles legacy generators which do not use ModelNames
for their figure legends. Not installed globally outside the context manager.
"""
from __future__ import annotations

from contextlib import contextmanager
import re

from Numerical.orchestration.ModelNames import model_label, shorten

_LONG_ID = re.compile(r'(?<![\w])T3_dS1_\d+_dS2_\d+_dF_\d+_alpha_[mp]\d+(?!\w)')
_SHORT_ID = re.compile(r'(?<![\w])T3-[A-E]-[mp]\d+(?!\w)')
_REPRESENTATION = re.compile(
    r'\(\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*\)\s*,\s*'
    r'(?:alpha\s*=\s*)?([+-])\s*(\d+)'
)


def display_label(value: str) -> str:
    """Convert known model IDs/representation labels within a display string."""
    def from_id(match):
        return shorten(match.group()).replace('T3-', 'T3', 1)

    def from_rep(match):
        d1, d2, df, sign, magnitude = match.groups()
        alpha = int(magnitude) * (-1 if sign == '-' else 1)
        try:
            return model_label(int(d1), int(d2), int(df), alpha).replace('T3-', 'T3', 1)
        except ValueError:
            return match.group()

    value = _LONG_ID.sub(from_id, value)
    value = _SHORT_ID.sub(from_id, value)
    return _REPRESENTATION.sub(from_rep, value)


def standardize_figure_labels(fig) -> None:
    """Relabel texts, tick marks, and legend entries on one rendered figure."""
    from matplotlib.text import Text

    for artist in fig.findobj(match=Text):
        old = artist.get_text()
        if old:
            new = display_label(old)
            if new != old:
                artist.set_text(new)

    legends = list(fig.legends)
    legends.extend(ax.get_legend() for ax in fig.axes if ax.get_legend() is not None)
    for legend in legends:
        names = legend.get_texts()
        # Legacy 16-model figures used 5pt labels; make them legible once shortened.
        if len(names) >= 8:
            for name in names:
                if name.get_fontsize() < 7:
                    name.set_fontsize(7)


@contextmanager
def presentation_labels_on_save():
    """Apply to the active report rebuild without affecting numerical routines."""
    from matplotlib.figure import Figure

    original = Figure.savefig

    def save_with_readable_names(self, *args, **kwargs):
        standardize_figure_labels(self)
        return original(self, *args, **kwargs)

    Figure.savefig = save_with_readable_names
    try:
        yield
    finally:
        Figure.savefig = original


@contextmanager
def presentation_envelope_labels():
    """Use the same short spellings in fitted/common-input HTML plot selectors.

    Only the generated `label` field is affected; internal `model` IDs in the
    JavaScript plot payload remain untouched. No scans or summary data change.
    """
    import Numerical.plotting.BuildFittedModelEnvelopes as builder

    original = builder.short_name

    def short_display_name(model):
        return display_label(original(model))

    builder.short_name = short_display_name
    try:
        yield
    finally:
        builder.short_name = original
