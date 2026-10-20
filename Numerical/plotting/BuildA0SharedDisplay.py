"""Build standalone A(0) shared-Yukawa displays from already completed diagnostics.

Presentation-only. No RGEs, fitting, matching, or changes to older displays.
Run: python -m Numerical.plotting.BuildA0SharedDisplay
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
from pathlib import Path
import shutil

import numpy as np

from Numerical.plotting.InteractiveModelComparison import write_dashboard
from Numerical.plotting.RunningResultFigures import generate_running_result_figures

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "output/studies/small_lambda_1e-7/a0_shared_benchmark"
DEST = ROOT / "Reports/output/display_a0_shared"
CHI2_STRICT = 5.86
CHI2_NOMINAL = 11.070497693516351


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def short(name):
    parts = name.split("_")
    try:
        dims = tuple(int(parts[parts.index(k)+1]) for k in ("dS1", "dS2", "dF"))
        cls = {(1,3,2):"A", (2,2,1):"B", (2,2,3):"C", (3,1,2):"D", (3,3,2):"E"}[dims]
        a = parts[parts.index("alpha")+1]
        sign = "-" if a.startswith("m") else "+"
        number = a[1:]
        return f"{cls}({sign}{number})"
    except (ValueError, KeyError, IndexError):
        return name


def get_chi2(comparison: Path, diagnostics: list[Path]):
    rows = {}
    if comparison.is_file():
        with comparison.open(newline="", encoding="utf-8-sig") as f:
            for row in csv.DictReader(f):
                name = row.get("model") or row.get("model_key")
                if not name:
                    continue
                for key in ("chi2", "common_chi2", "shared_chi2", "fit_chi2"):
                    if row.get(key) not in (None, ""):
                        try:
                            rows[name] = float(row[key])
                        except ValueError:
                            pass
                        break
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--output", type=Path, default=DEST)
    parser.add_argument("--resume", action="store_true", help="Skip existing PNG regeneration; HTML and CSV still rebuild")
    args = parser.parse_args()
    source, output = args.source.resolve(), args.output.resolve()
    files = sorted((source / "models").glob("*/running_diagnostics.json"))
    files = [p for p in files if read_json(p).get("status") == "Success"]
    if not files:
        parser.error(f"No successful running diagnostics found under {source / 'models'}")
    if source == output or source in output.parents or output in source.parents:
        parser.error("Output must be separate from source directory")
    output.mkdir(parents=True, exist_ok=True)
    figures = output / "figures"
    figures.mkdir(exist_ok=True)
    chi2_by_model = get_chi2(source / "comparison.csv", files)
    records = []
    hashes = {}
    for i, path in enumerate(files, 1):
        model = path.parent.name
        payload = read_json(path)
        dst = figures / model
        if not args.resume or not (dst / "uv_coupling_running.png").is_file():
            generate_running_result_figures(path, dst)
        record = {"model":model, "label":short(model), "chi2":chi2_by_model.get(model),
                  "diagnostics":str(path.relative_to(source))}
        if record["chi2"] is not None:
            record["strict_pass"] = record["chi2"] < CHI2_STRICT
            record["nominal_pass"] = record["chi2"] <= CHI2_NOMINAL
        else:
            record["strict_pass"] = None
            record["nominal_pass"] = None
        records.append(record)
        hashes[str(path.relative_to(source))] = hashlib.sha256(path.read_bytes()).hexdigest()
        print(f"[{i:02d}/{len(files):02d}] {record['label']}: figure display ready", flush=True)

    dashboard = output / "interactive_comparison.html"
    write_dashboard(source / "models", dashboard, "T3: A(0) shared Yukawa, lambdaT3=1e-7 (A-D), 5e-8 (E)")
    with (output / "model_summary.csv").open("w", encoding="utf-8", newline="") as f:
        w=csv.DictWriter(f, fieldnames=list(records[0]))
        w.writeheader(); w.writerows(records)
    meta = {"benchmark":"A(0) shared UV Yukawas", "source":str(source),
            "strict_screen_cutoff":CHI2_STRICT,"nominal_95_screen_cutoff":CHI2_NOMINAL,
            "scope":"saved diagnostics only; screening thresholds are not fitted-parameter confidence regions",
            "source_sha256":hashes,"models":records}
    (output / "display_manifest.json").write_text(json.dumps(meta, indent=2)+"\n",encoding="utf-8")

    cards = []
    for r in sorted(records, key=lambda r: float('inf') if r['chi2'] is None else r['chi2']):
        chi = "not in CSV" if r['chi2'] is None else f"{r['chi2']:.6g}"
        status = "not classified" if r['chi2'] is None else ("strict pass" if r['strict_pass'] else "95% screening pass" if r['nominal_pass'] else "fails screening")
        model = html.escape(r['model'])
        image_links = []
        for file in sorted((figures / r['model']).glob("*.png")):
            rel = file.relative_to(output).as_posix()
            image_links.append(f'<a href="{html.escape(rel)}">{html.escape(file.stem.replace("_"," "))}</a>')
        cards.append(f'<tr><td title="{model}">{html.escape(r["label"])}</td><td>{chi}</td><td>{status}</td><td>{" · ".join(image_links)}</td></tr>')
    num_strict = sum(x['strict_pass'] is True for x in records)
    num_nominal = sum(x['nominal_pass'] is True for x in records)
    summary = f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>T3 A(0) shared benchmark</title><style>
body{{font:15px system-ui,Segoe UI,sans-serif;margin:2rem auto;padding:0 1.5rem;max-width:1200px;background:#f8fafc;color:#172438}}h1{{font-size:1.7rem}}.sub{{color:#536578}}a{{color:#245db0}}table{{border-collapse:collapse;width:100%;background:white}}td,th{{padding:.7rem;border-bottom:1px solid #dbe2e9;text-align:left;vertical-align:top}}tr:hover{{background:#f3f7ff}}td:last-child{{font-size:.82rem;line-height:1.8}}.note{{padding:1rem 0;color:#4f5e70}}
</style></head><body><h1>A(0) shared-Yukawa benchmark</h1><p class="sub">One identical UV Yukawa pair evaluated in 16 T3 representations; lambdaT3=1e-7 for A–D and 5e-8 for E.</p><p><strong>{num_strict}/{len(records)}</strong> strict screening passes (chi² &lt; 5.86); <strong>{num_nominal}/{len(records)}</strong> nominal-95% screening passes (chi² ≤ 11.0705).</p><p><a href="interactive_comparison.html"><strong>Open interactive cross-model running dashboard</strong></a> · <a href="model_summary.csv">CSV table</a> · <a href="display_manifest.json">Input provenance</a></p><table><thead><tr><th>Model</th><th>chi²</th><th>Classification</th><th>Running figures</th></tr></thead><tbody>{''.join(cards)}</tbody></table><p class="note">Screening against five NuFIT observables, not a calibrated confidence region for 18 fitted Yukawa parameters. The interactive Weinberg continuation is a display reconstruction, not an independently evolved LLSS Wilson coefficient. Diagnostics and fit values are reused without further optimisation.</p></body></html>'''
    (output / "index.html").write_text(summary, encoding="utf-8")
    print(f"Built display: {output}")
    print(f"Successful diagnostics: {len(files)}, strict pass: {num_strict}, nominal95 pass: {num_nominal}")
    if len(chi2_by_model)<len(files):
        print("WARNING: Some model chi2 values were absent from comparison.csv; see model_summary.csv")

if __name__ == "__main__":
    main()
