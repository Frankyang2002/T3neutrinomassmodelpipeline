"""Build model-centric four-scenario reports from existing numerical JSON only.

Usage: python -m Numerical.orchestration.BuildComparisonReports [--strict]
No symbolic matching, RGE integration, or numerical resume is performed.
"""
from __future__ import annotations

import argparse
import html
import json
from pathlib import Path

from Numerical.plotting.RunningTrajectoryData import SCENARIOS, discover_saved_runs
from Numerical.plotting.WithinModelComparison import generate_model_figures, LABELS
from Numerical.plotting.WithinModelInteractive import write_interactive_comparison

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = PROJECT_ROOT / "output" / "full" / "comparison"
DEFAULT_OUTPUT = PROJECT_ROOT / "Reports" / "output" / "full" / "within_model_comparison"


def build_reports(input_root: Path = DEFAULT_INPUT, output_root: Path = DEFAULT_OUTPUT,
                  *, strict: bool = False) -> dict:
    models = discover_saved_runs(input_root)
    if not models:
        raise FileNotFoundError(f"No saved full-comparison diagnostics under {input_root}")
    summary = {"status": "Complete", "models": {}, "input_root": str(input_root),
               "output_root": str(output_root)}
    for model_key, runs in sorted(models.items()):
        missing = [scenario for scenario in SCENARIOS if scenario not in runs]
        if missing and strict:
            raise ValueError(f"{model_key}: missing scenarios {missing}")
        directory = Path(output_root) / model_key
        figures, issues = generate_model_figures(runs, directory / "figures")
        interactive_path, interactive_warnings = write_interactive_comparison(
            runs, directory / "interactive_comparison.html")
        issues = list(dict.fromkeys([*issues, *interactive_warnings]))
        if missing or issues:
            summary["status"] = "Partial"
        source_rows = "\n".join(
            f'<tr><td>{html.escape(scenario)}</td><td>{html.escape(LABELS[scenario])}</td>'
            f'<td>{html.escape(str(runs[scenario].path)) if scenario in runs else "MISSING"}</td></tr>'
            for scenario in SCENARIOS
        )
        image_blocks = "\n".join(
            f'<section><h2>{html.escape(name.replace("_", " ").removesuffix(".png"))}</h2>'
            f'<img src="figures/{html.escape(name)}" alt="{html.escape(name)}"></section>'
            for name in figures
        )
        notes = "\n".join(f"<li>{html.escape(note)}</li>" for note in issues)
        page = f'''<!doctype html><html lang="en"><head><meta charset="utf-8">
<title>T3 within-model comparison: {html.escape(model_key)}</title>
<style>body{{font-family:system-ui,sans-serif;margin:2rem auto;max-width:1100px;padding:0 1rem;color:#222}}
img{{max-width:100%;height:auto}}section{{margin:2rem 0;border-top:1px solid #ddd}}
table{{border-collapse:collapse;width:100%;font-size:0.9rem}}td,th{{text-align:left;border:1px solid #ddd;padding:0.5rem;overflow-wrap:anywhere}}
code{{overflow-wrap:anywhere}}</style></head><body>
<h1>{html.escape(model_key)}</h1>
<p>Four fixed benchmark scenarios from previously computed diagnostics. No physics recomputation.</p>
<p>Intermediate dotted Weinberg curves, when present, represent a validated
threshold-anchored hard-plus-direct display reconstruction. They are <strong>not</strong>
an intermediate-EFT Wilson coefficient. Actual direct LLSS-to-Weinberg contributions
are shown separately. Class-E lambda values use the saved benchmark prescription,
whose multiplicity justification remains under validation.</p>
<table><thead><tr><th>Scenario</th><th>Reference couplings</th><th>Source diagnostics</th></tr></thead>
<tbody>{source_rows}</tbody></table>
<p><a href="interactive_comparison.html"><strong>Open interactive four-scenario graph</strong></a> (quantity selector, scenario toggles, zoom, hover)</p>
<h2>Missing scenarios</h2><p>{html.escape(', '.join(missing) if missing else 'None')}</p>
<h2>Reconstruction warnings</h2><ul>{notes or '<li>None</li>'}</ul>
{image_blocks}</body></html>'''
        directory.mkdir(parents=True, exist_ok=True)
        (directory / "comparison_report.html").write_text(page, encoding="utf-8")
        summary["models"][model_key] = {"scenarios": list(runs), "missing": missing,
                                          "reconstruction_warnings": issues,
                                          "figures": figures, "interactive_report": interactive_path.name}
    if len(models) != 16:
        summary["status"] = "Partial"
    Path(output_root).mkdir(parents=True, exist_ok=True)
    (Path(output_root) / "report_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    p.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    p.add_argument("--strict", action="store_true", help="require all four scenarios per model")
    args = p.parse_args(argv)
    result = build_reports(args.input, args.output, strict=args.strict)
    print(f"Status: {result['status']}; models: {len(result['models'])}; output: {args.output}")
    return 0 if result["status"] == "Complete" else 1


if __name__ == "__main__":
    raise SystemExit(main())
