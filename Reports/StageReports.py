from __future__ import annotations

from pathlib import Path

from common.Paths import REPORT_OUTPUT_DIR
from common.Records import RunRecord


def final_eft_stage_label(records: list[RunRecord]) -> str:
    """Return the common final EFT stage label for one comparison run."""

    for record in records:
        if record.eft_stages:
            return record.eft_stages[-1].label

        stages = record.summary.get("EFTStages", [])
        if stages:
            label = stages[-1].get("Label")
            if label:
                return str(label)

    return "EFT"


def lagrangian_report_path(stage_label: str, report_root: Path | None = None) -> Path:
    """Return the .tex path for one stage-aware Lagrangian report."""

    root = report_root or REPORT_OUTPUT_DIR
    return root / "Lagrangian" / f"{stage_label}.tex"


def rge_report_path(stage_label: str, report_root: Path | None = None) -> Path:
    """Return the .tex path for one stage-aware RGE report."""

    root = report_root or REPORT_OUTPUT_DIR
    return root / "RGE" / f"{stage_label}.tex"
