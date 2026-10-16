"""Read-only inventory of generated report artifacts and numerical dependencies.

Run from the repository root:
    python -m tools.AuditReportArtifacts
    python -m tools.AuditReportArtifacts --report-root Reports/output/full

No files are deleted, moved, renamed or overwritten. Output manifest defaults
outside Reports/output to avoid confusing it with numerical result inputs.
"""
from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_REPORT_ROOT = PROJECT_ROOT / "Reports" / "output" / "full"
DEFAULT_MANIFEST = PROJECT_ROOT / "report_artifact_audit.csv"

# These directories are confirmed consumers of saved numerical scan data.
NUMERICAL_DATA_STUDIES = {
    "fitted_sensitivity", "fit_feasibility", "fit_feasibility_refined",
    "fit_feasibility_refined2", "staged_rg_comparison",
    "weinberg_running_comparison", "fermion_scale",
}
PRESENTATION_SUFFIXES = {".png", ".jpg", ".jpeg", ".webp", ".svg", ".html"}
DOCUMENT_SUFFIXES = {".md", ".txt", ".tex", ".pdf"}
NUMERICAL_SUFFIXES = {".json", ".csv", ".npz", ".npy", ".h5", ".hdf5"}
BUILD_SUFFIXES = {".aux", ".log", ".out", ".toc", ".synctex.gz"}


def classify(path: Path, root: Path) -> tuple[str, str]:
    relative = path.relative_to(root)
    suffix = ".synctex.gz" if path.name.endswith(".synctex.gz") else path.suffix.lower()
    if suffix in NUMERICAL_SUFFIXES:
        return "KEEP_NUMERICAL", "Saved results/configurations may be read by other studies"
    if suffix in BUILD_SUFFIXES:
        return "REVIEW_BUILD", "Likely generated compilation artifact; inspect before cleanup"
    if suffix in PRESENTATION_SUFFIXES:
        return "REVIEW_REGENERABLE", "Generated presentation candidate; verify exact source and rebuild"
    if suffix in DOCUMENT_SUFFIXES:
        return "REVIEW_DOCUMENT", "May contain unique derivations, context or thesis material"
    if relative.parts and relative.parts[0] in NUMERICAL_DATA_STUDIES:
        return "REVIEW_STUDY_FILE", "Inside study directory used by numerical/reporting workflows"
    return "REVIEW_UNKNOWN", "Unclassified: do not delete without inspection"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report-root", type=Path, default=DEFAULT_REPORT_ROOT)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    args = parser.parse_args(argv)
    root = args.report_root.resolve()
    manifest = args.manifest.resolve()
    if not root.is_dir():
        parser.error(f"Report root not found: {root}")
    if manifest == root or root in manifest.parents:
        parser.error("Manifest must be outside the report tree")
    files = sorted(p for p in root.rglob("*") if p.is_file())
    counts: Counter[str] = Counter()
    sizes: Counter[str] = Counter()
    manifest.parent.mkdir(parents=True, exist_ok=True)
    with manifest.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(("path", "bytes", "classification", "reason"))
        for path in files:
            cls, reason = classify(path, root)
            size = path.stat().st_size
            writer.writerow((path.relative_to(PROJECT_ROOT).as_posix()
                             if path.is_relative_to(PROJECT_ROOT) else str(path), size, cls, reason))
            counts[cls] += 1
            sizes[cls] += size
    print(f"Scanned {len(files)} files; manifest: {manifest}")
    for cls in sorted(counts):
        print(f"{cls:22} {counts[cls]:7d} files  {sizes[cls] / (1024 ** 2):9.2f} MiB")
    print("READ-ONLY: no report or numerical files were changed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
