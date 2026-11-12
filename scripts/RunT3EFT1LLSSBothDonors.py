#!/usr/bin/env python3
"""Replot completed all-16 A(-4) and D(-2) LLSS ablation results.

Project path: scripts/RunT3EFT1LLSSBothDonors.py

Uses existing local CSVs when present. Otherwise downloads the two published
CSVs and provenance files from the repository's GitHub main branch. NO RGE,
matching, or oscillation physics is rerun. The unchanged underlying plotting
runner renders title-free heatmaps with Matplotlib default fonts. Font sizes
are fully configurable, with larger defaults for thesis readability.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
import sys
from pathlib import Path
from urllib.error import URLError
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]
RUNNER = ROOT / "scripts" / "RunT3EFT1LLSSIsolation.py"
OUTPUT_ROOT = ROOT / "output" / "parallel_t3"
PUBLIC_RESULTS = Path("output/parallel_t3/eft1_llss_both_donors_end")
REPO_RAW = "https://raw.githubusercontent.com/Frankyang2002/T3neutrinomassmodelpipeline/main/"
DONORS = (("A_-4", "A_m4", "a4_source"), ("D_-2", "D_m2", "d2_source"))
MODELS = (
    "A_-4", "A_-2", "A_+0", "B_-3", "B_-1", "B_+1",
    "C_-3", "C_-1", "C_+1", "D_-2", "D_+0", "D_+2",
    "E_-4", "E_-2", "E_+0", "E_+2",
)
# Multiplier relative to --font-size; explicit flags override these.
FONT_FACTORS = {
    "axis_label_fontsize": 1.15,
    "bottom_tick_fontsize": 1.0,
    "top_tick_fontsize": 1.25,
    "y_tick_fontsize": 1.10,
    "cell_fontsize": 0.90,
    "colorbar_fontsize": 1.0,
}


def _validate_saved(path: Path, donor: str) -> dict:
    csv_path = path / "all16_three_way.csv"
    provenance_path = path / "provenance.json"
    if not csv_path.is_file() or not provenance_path.is_file():
        raise FileNotFoundError(f"{donor}: required results not found in {path}")
    with csv_path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != len(MODELS) or set(r.get("model") for r in rows) != set(MODELS):
        raise ValueError(f"{donor}: source does not contain 16 distinct T3 models")
    if any(r.get("reference") != donor for r in rows):
        raise ValueError(f"{donor}: CSV contains wrong donor references")
    if any(r.get("status") != "VALID_EXPLORATORY" for r in rows):
        raise ValueError(f"{donor}: not every saved model has a valid completed result")
    for row in rows:
        for name in ("dm21_sq_eV2", "dm31_sq_eV2", "sin2_theta12", "sin2_theta13", "sin2_theta23"):
            for case in ("W_FULL", "W_FROZEN", "LLSS_OFF"):
                key = f"{name}_{case}"
                if not row.get(key):
                    raise ValueError(f"{donor}/{row['model']}: missing {key}")
                try:
                    number = float(row[key])
                except ValueError as exc:
                    raise ValueError(f"{donor}/{row['model']}: invalid {key}") from exc
                if not (-float("inf") < number < float("inf")):
                    raise ValueError(f"{donor}/{row['model']}: non-finite {key}")
    meta = json.loads(provenance_path.read_text(encoding="utf-8-sig"))
    if meta.get("reference") != donor:
        raise ValueError(f"{donor}: provenance identifies {meta.get('reference')}")
    return {
        "csv_sha256": hashlib.sha256(csv_path.read_bytes()).hexdigest(),
        "provenance_sha256": hashlib.sha256(provenance_path.read_bytes()).hexdigest(),
        "models_valid": len(rows),
    }


def _download_bytes(url: str) -> bytes:
    req = Request(url, headers={"User-Agent": "T3-EFT1-LLSS-replot/1.0"})
    with urlopen(req, timeout=30) as response:
        return response.read()


def _resolve_source(args: argparse.Namespace, donor: str, tag: str, attribute: str,
                    output: Path) -> tuple[Path, str, dict]:
    user_path = getattr(args, attribute)
    if user_path is not None:
        path = (ROOT / user_path).resolve()
        return path, "explicit_local", _validate_saved(path, donor)
    local = ROOT / PUBLIC_RESULTS / tag
    if not args.github_only and local.is_dir() and (local / "all16_three_way.csv").is_file():
        return local, "local_saved", _validate_saved(local, donor)
    source = output / "github_inputs" / tag
    source.mkdir(parents=True, exist_ok=False)
    try:
        for file in ("all16_three_way.csv", "provenance.json"):
            relative = (PUBLIC_RESULTS / tag / file).as_posix()
            content = _download_bytes(REPO_RAW + relative)
            (source / file).write_bytes(content)
    except (OSError, URLError) as exc:
        raise RuntimeError(
            f"Unable to download {donor} saved inputs from GitHub: {exc}. "
            "Alternatively pull the latest repository or set --a4-source/--d2-source."
        ) from exc
    return source, "github_main", _validate_saved(source, donor)


def _style_options(args: argparse.Namespace) -> list[str]:
    if not 0 < args.font_size <= 72:
        raise ValueError("--font-size must be between 0 and 72")
    opts = []
    for name, factor in FONT_FACTORS.items():
        explicit = getattr(args, name)
        size = float(explicit if explicit is not None else args.font_size * factor)
        if not 0 < size <= 72:
            raise ValueError(f"invalid {name}={size}")
        opts.extend(("--" + name.replace("_", "-"), str(size)))
    if args.figure_width <= 0 or args.row_height <= 0 or not 72 <= args.dpi <= 1200:
        raise ValueError("invalid figure width/row height/DPI")
    opts.extend(("--figure-width", str(args.figure_width), "--row-height", str(args.row_height),
                 "--dpi", str(args.dpi)))
    return opts


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--a4-source", type=Path, help="Optional saved A(-4) results directory")
    p.add_argument("--d2-source", type=Path, help="Optional saved D(-2) results directory")
    p.add_argument("--github-only", action="store_true", help="Download both donors from public GitHub instead of local saved results")
    p.add_argument("--output", type=Path, default=Path("output/parallel_t3/eft1_llss_both_donors_large_font"))
    p.add_argument("--font-size", type=float, default=14.0, help="Base text font size, default 14 pt")
    for name in FONT_FACTORS:
        p.add_argument("--" + name.replace("_", "-"), type=float, help="Override a specific font size in points")
    p.add_argument("--figure-width", type=float, default=21.0)
    p.add_argument("--row-height", type=float, default=0.67)
    p.add_argument("--dpi", type=int, default=300)
    args = p.parse_args(argv)
    try:
        if not RUNNER.is_file():
            raise FileNotFoundError(f"Missing plotting runner: {RUNNER}")
        out = (ROOT / args.output).resolve()
        if out == OUTPUT_ROOT.resolve() or OUTPUT_ROOT.resolve() not in out.parents:
            raise ValueError("--output must be a new folder inside output/parallel_t3")
        if out.exists():
            raise ValueError(f"Output already exists: {out}; choose a new output directory")
        style = _style_options(args)
    except (ValueError, FileNotFoundError) as exc:
        p.error(str(exc))

    out.mkdir(parents=True)
    report = {
        "schema": "t3-eft1-llss-two-donor-replot-large-font-v1",
        "font_family": "Matplotlib default",
        "font_size_base_points": args.font_size,
        "title": None,
        "physics_recomputed": False,
        "reference": "W_FROZEN",
        "running": "100*(O_W_FULL-O_W_FROZEN)/O_W_FROZEN",
        "mixing_off": "100*(O_LLSS_OFF-O_W_FROZEN)/O_W_FROZEN",
        "source_repository": "Frankyang2002/T3neutrinomassmodelpipeline",
        "source_directory": str(PUBLIC_RESULTS),
        "donors": [],
    }
    succeeded = True
    for donor, tag, attribute in DONORS:
        record = {"donor": donor, "output": str(out / tag)}
        try:
            source, mode, details = _resolve_source(args, donor, tag, attribute, out)
            record.update({"input_source": mode, "input_directory": str(source), **details})
            cmd = [sys.executable, str(RUNNER), "--plot-only", str(source),
                   "--output", str(out / tag), *style]
            print(f"{donor}: plotting {details['models_valid']} saved models from {mode}", flush=True)
            result = subprocess.run(cmd, cwd=ROOT, check=False)
            record["exit_code"] = int(result.returncode)
            if result.returncode:
                succeeded = False
        except Exception as exc:
            record["error"] = f"{type(exc).__name__}: {exc}"
            succeeded = False
            print(f"{donor}: BLOCKED: {record['error']}", flush=True)
        report["donors"].append(record)
        (out / "both_donors_summary.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return 0 if succeeded else 1


if __name__ == "__main__":
    raise SystemExit(main())
