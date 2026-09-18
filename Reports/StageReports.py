from __future__ import annotations

"""Stage-aware report path helpers."""

from pathlib import Path
from typing import Iterable

from common.Paths import REPORT_OUTPUT_DIR


def final_eft_stage_label(records: Iterable[object]) -> str:
    """Return the final EFT stage label present in the supplied run records."""
    for record in records:
        stages = getattr(record, "eft_stages", None)
        if stages:
            return stages[-1].label
    return "EFT_2_after_S1_S2"


def lagrangian_report_path(
    stage_label: str,
    report_root: Path | None = None,
) -> Path:
    """Return the study-aware .tex path for one Lagrangian-stage report."""
    root = report_root or REPORT_OUTPUT_DIR
    return root / "Lagrangian" / f"{stage_label}.tex"


def rge_report_path(
    stage_label: str,
    report_root: Path | None = None,
) -> Path:
    """Return the study-aware .tex path for one RGE-stage report."""
    root = report_root or REPORT_OUTPUT_DIR
    return root / "RGE" / f"{stage_label}.tex"


def group_factor_report_path(
    stage_label: str,
    report_root: Path | None = None,
) -> Path:
    """Return the .tex path for one stage-aware analytic group-factor report."""
    root = report_root or REPORT_OUTPUT_DIR
    return root / "GroupFactors" / f"GF_{stage_label}.tex"


def c5_report_path(
    report_root: Path | None = None,
) -> Path:
    """Return the study-aware .tex path for the standalone C5 report."""
    root = report_root or REPORT_OUTPUT_DIR
    return root / "Lagrangian" / "C5.tex"
