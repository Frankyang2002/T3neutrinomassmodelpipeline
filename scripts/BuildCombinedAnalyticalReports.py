"""Rebuild genuine across-model analytical tables from previously saved full-study runs.

No Matchete, RGBeta, ODE integration or changes to individual model reports.
Run from project root:
    python scripts/BuildCombinedAnalyticalReports.py --check
    python scripts/BuildCombinedAnalyticalReports.py

Expected saved structure, for one complete scenario:
    output/full/comparison/<scenario>/<model-key>/t3_model_comparison.json
    output/full/comparison/<scenario>/<model-key>/<physical-run>/data/...
Model keys can be long (T3_dS1_...) or short (T3-A-m4).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from common.RunRecords import RunRecord

MODELS = (
    ("A", 1, 3, 2, -4), ("A", 1, 3, 2, -2), ("A", 1, 3, 2, 0),
    ("B", 2, 2, 1, -3), ("B", 2, 2, 1, -1), ("B", 2, 2, 1, 1),
    ("C", 2, 2, 3, -3), ("C", 2, 2, 3, -1), ("C", 2, 2, 3, 1),
    ("D", 3, 1, 2, -2), ("D", 3, 1, 2, 0), ("D", 3, 1, 2, 2),
    ("E", 3, 3, 2, -4), ("E", 3, 3, 2, -2), ("E", 3, 3, 2, 0), ("E", 3, 3, 2, 2),
)


def model_keys(item: tuple[str, int, int, int, int]) -> tuple[str, str]:
    cls, s1, s2, f, a = item
    suffix = f"m{-a}" if a < 0 else f"p{a}"
    return f"T3-{cls}-{suffix}", f"T3_dS1_{s1}_dS2_{s2}_dF_{f}_alpha_{suffix}"


def _load_one_summary(model_directory: Path) -> dict:
    file = model_directory / "t3_model_comparison.json"
    if not file.is_file():
        raise FileNotFoundError(f"Missing saved run summary: {file}")
    payload = json.loads(file.read_text(encoding="utf-8-sig"))
    if not isinstance(payload, list) or len(payload) != 1 or not isinstance(payload[0], dict):
        raise ValueError(f"Expected one model summary in {file}; got {type(payload).__name__}")
    return payload[0]


def reconstruct_records(source: Path) -> list[RunRecord]:
    """Reconstruct exactly the 16 ordinary models, before writing reports."""
    records: list[RunRecord] = []
    missing: list[str] = []
    for model in MODELS:
        short, long = model_keys(model)
        options = [source / short, source / long]
        directories = [p for p in options if p.is_dir()]
        if len(directories) > 1:
            raise ValueError(f"Ambiguous short and long output folders for {short}: {directories}")
        if not directories:
            missing.append(f"{short}: expected {options[0]} or {options[1]}")
            continue
        folder = directories[0]
        summary = _load_one_summary(folder)
        # A physical run is nested immediately under the study folder. Avoid
        # relying on the report-folder presentation name; its metadata lives
        # within the single-model run and is needed for all existing readers.
        uv_relative = summary.get("UVRGEFile")
        if not isinstance(uv_relative, str) or not uv_relative:
            raise ValueError(f"Missing UVRGEFile in saved summary: {folder}")
        run_dirs = []
        for p in folder.iterdir():
            if p.is_dir() and (p / uv_relative).is_file():
                run_dirs.append(p)
        if len(run_dirs) != 1:
            raise ValueError(f"Expected exactly one physical run containing {uv_relative} below {folder}; found {len(run_dirs)}")
        run_dir = run_dirs[0]
        label, ds1, ds2, df, alpha = model
        if summary.get("UVRGEStatus") != "Success":
            raise ValueError(f"{short}: saved UV RGE was not successful")
        if not isinstance(summary.get("EFTStages"), list):
            raise ValueError(f"{short}: EFTStages metadata missing; cannot compare all EFT reports")
        records.append(RunRecord(
            name=short, alpha=alpha, d_s1=ds1, d_s2=ds2, d_f=df,
            return_code=0, summary=summary, output_dir=run_dir,
            shared_scalar=False,
        ))
    if missing:
        raise FileNotFoundError("Missing models:\n" + "\n".join(missing))
    if len(records) != 16:
        raise RuntimeError(f"Expected 16 reconstructed records, got {len(records)}")
    return records


def generate(records: list[RunRecord], destination: Path) -> dict:
    from Reports.ReportGeneration import (
        write_bsm_uv_field_table, write_bsm_matched_field_table,
        write_c5_coefficient_report, compile_latex_document,
    )
    from Reports.RGEComparison import (
        write_and_compile_rge_comparison,
        write_and_compile_eft1_rge_comparison,
        write_and_compile_final_eft_rge_comparison,
    )
    from Reports.GroupFactorReports import write_and_compile_stage_group_factor_reports

    destination.mkdir(parents=True, exist_ok=True)
    created = []
    uv_lagrangian = write_bsm_uv_field_table(records, report_root=destination)
    compile_latex_document(uv_lagrangian)
    created.append(str(uv_lagrangian))
    created.append(str(write_bsm_matched_field_table(records, report_root=destination)))
    c5 = write_c5_coefficient_report(records, report_root=destination)
    compile_latex_document(c5)
    created.append(str(c5))
    created.extend(str(p) for p in (
        write_and_compile_rge_comparison(records, report_root=destination),
        write_and_compile_eft1_rge_comparison(records, report_root=destination),
        write_and_compile_final_eft_rge_comparison(records, report_root=destination),
    ))
    created.extend(str(p) for p in write_and_compile_stage_group_factor_reports(
        records, report_root=destination,
    ))
    result = {"status": "Complete", "models": [r.name for r in records],
              "source_runs": [str(r.output_dir) for r in records],
              "report_tex_files": created}
    (destination / "combined_analytical_manifest.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scenario", default="smallY_smallL",
                        choices=("smallY_smallL", "smallY_largeL", "largeY_smallL", "largeY_largeL"))
    parser.add_argument("--source", type=Path,
                        help="Override input scenario folder")
    parser.add_argument("--output", type=Path,
                        default=PROJECT_ROOT / "Reports" / "output" / "analytical")
    parser.add_argument("--check", action="store_true",
                        help="Validate 16 saved records without generating reports")
    args = parser.parse_args(argv)
    source = args.source or PROJECT_ROOT / "output" / "full" / "comparison" / args.scenario
    records = reconstruct_records(source)
    print(f"Validated {len(records)} saved analytical models in {source}")
    for record in records:
        print(f"  {record.name}: {record.output_dir}")
    if args.check:
        return 0
    result = generate(records, args.output)
    print(f"Combined analytical report files: {len(result['report_tex_files'])} in {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
