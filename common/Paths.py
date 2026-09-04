from __future__ import annotations

from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parent.parent

LAGRANGIAN_DIR = PROJECT_ROOT / "Lagrangian"
RGE_DIR = PROJECT_ROOT / "RGE"

# Machine-readable calculation output remains here.
OUTPUT_DIR = PROJECT_ROOT / "output"

# Human-readable report output lives separately.
REPORTS_DIR = PROJECT_ROOT / "Reports"
REPORT_OUTPUT_DIR = REPORTS_DIR / "output"

RUN_MODEL_SCRIPT = LAGRANGIAN_DIR / "RunModel.wl"

# Weinberg 5D at 1 loop
EFT_ORDER = 5
LOOP_ORDER = 1
