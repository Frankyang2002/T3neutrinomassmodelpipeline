from __future__ import annotations

"""Stage-aware report path helpers."""

from pathlib import Path
from typing import Iterable

from common.Paths import REPORT_OUTPUT_DIR


def final_eft_stage_label(records: Iterable[object]) -> str:
    """Return the final recorded physical EFT-stage label.

    Stage labels come from each run's physical threshold metadata. If no
    stage metadata is available, use the neutral fallback ``EFT`` rather than
    assuming the ordinary-T3 ``S1,S2`` field content.
    """
    for record in records:
        stages = getattr(record, "eft_stages", None)
        if stages:
            return stages[-1].label

    return "EFT"


def lagrangian_report_path(
    stage_label: str,
    report_root: Path | None = None,
) -> Path:
    root = report_root or REPORT_OUTPUT_DIR
    return root / "Lagrangian" / f"{stage_label}.tex"


def rge_report_path(
    stage_label: str,
    report_root: Path | None = None,
) -> Path:
    root = report_root or REPORT_OUTPUT_DIR
    return root / "RGE" / f"{stage_label}.tex"


def group_factor_report_path(
    stage_label: str,
    report_root: Path | None = None,
) -> Path:
    root = report_root or REPORT_OUTPUT_DIR
    return root / "GroupFactors" / f"GF_{stage_label}.tex"


def c5_report_path(report_root: Path | None = None) -> Path:
    root = report_root or REPORT_OUTPUT_DIR
    return root / "Lagrangian" / "C5.tex"
