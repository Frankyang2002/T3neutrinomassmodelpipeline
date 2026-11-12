#!/usr/bin/env python3
"""Thesis-consistent EFT1 LLSS isolation for A(-4) and D(-2) Yukawa donors.

Running  = 100*(O_1111 - O_1001)/O_1001
Mixing off = 100*(O_LLSS_OFF - O_1001)/O_1001

The 1111/1001 observable values are read *unaltered* from the original
thesis stage-freeze results (EFT1 stage frozen = renormalisable + Wilson).
Only the missing LLSS_OFF counterfactual is evaluated. This guarantees that
Running reproduces the source thesis EFT1 stage-freeze data exactly.

The LLSS_OFF evaluation evolves UV, both-scalar EFT1 background, all other
Wilson terms, matches both scalars together, and evolves the SMEFT; only
[lambda_HHSS C_LLSS]_{LLHH} is removed in each Wilson RHS evaluation.

Scope: fixed Yukawa, conditional exploratory basis/ADM/matching, NOT a
certified physical MSbar calculation. Does not modify GitHub or production code.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
import traceback
from contextlib import ExitStack
from pathlib import Path
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

OBSERVABLES = ("dm21_sq_eV2", "dm31_sq_eV2", "sin2_theta12",
               "sin2_theta13", "sin2_theta23")
DONORS = ("A_-4", "D_-2")
ORIGINAL_STAGE = Path("output/parallel_t3/thesis_two_benchmarks_all16_stagefreeze_v2")
ORIGINAL_TABLES = Path("output/parallel_t3/thesis_generalised_tables_v1")
ORIGINAL_FITS = Path("output/parallel_t3/A4_Dm2_targeted_fit01")
SELECTED_SIX = ("A_+0", "B_+1", "C_+1", "D_-2", "E_+0", "E_-4")


def percentage(num: float, baseline: float):
    """Return signed percent; None for a zero/nonfinite frozen denominator."""
    try:
        a, b = float(num), float(baseline)
    except (ValueError, TypeError):
        return None
    return 100.0 * (a - b) / b if math.isfinite(a) and math.isfinite(b) and b != 0 else None


def stagefreeze_lookup(path: Path) -> dict[tuple[str, str], dict]:
    """Read only EFT1-off (1001) rows from the actual thesis running CSV."""
    with path.open(encoding="utf-8-sig", newline="") as fh:
        records = list(csv.DictReader(fh))
    found = {}
    for row in records:
        if row.get("stage_frozen") != "EFT1":
            continue
        key = (row["reference"], row["model"])
        if key in found:
            raise ValueError(f"Duplicate thesis EFT1 stage-freeze result {key}")
        if row.get("status") != "VALID_EXPLORATORY":
            continue
        numbers = {}
        for obs in OBSERVABLES:
            full = float(row[f"{obs}_full"])
            frozen = float(row[f"{obs}_frozen"])
            reported = float(row[f"{obs}_pct"])
            expected = percentage(full, frozen)
            if expected is None or not math.isclose(reported, expected, rel_tol=2e-9, abs_tol=2e-10):
                raise ValueError(f"Corrupt thesis percentage: {key}, {obs}")
            numbers[obs] = (full, frozen, reported)
        found[key] = numbers
    return found


def style_from_args(args) -> dict:
    names = ("axis_label_fontsize", "bottom_tick_fontsize", "top_tick_fontsize",
             "y_tick_fontsize", "cell_fontsize", "colorbar_fontsize")
    base = float(args.font_size)
    if not math.isfinite(base) or base <= 0:
        raise ValueError("--font-size must be positive")
    offsets = dict(axis_label_fontsize=1.0, bottom_tick_fontsize=0.0,
                   top_tick_fontsize=2.0, y_tick_fontsize=0.0,
                   cell_fontsize=-1.0, colorbar_fontsize=0.0)
    style = {name: float(getattr(args, name)) if getattr(args, name) is not None
             else max(1.0, base + offsets[name]) for name in names}
    style.update(figure_width=float(args.figure_width),
                 row_height=float(args.row_height), dpi=int(args.dpi))
    if any(not math.isfinite(v) or v <= 0 for v in style.values()):
        raise ValueError("All font sizes, figure dimensions, and DPI must be positive")
    return style


def make_row(reference: str, model: str, baseline: dict,
             mixing_off: dict, *, source: str) -> dict:
    result = {"reference": reference, "model": model,
              "status": "VALID_EXPLORATORY", "error": "",
              "thesis_source": source}
    for obs in OBSERVABLES:
        full, frozen, thesis_pct = baseline[obs]
        removed = mixing_off.get(obs)
        if percentage(full, frozen) is None or percentage(removed, frozen) is None:
            raise ValueError(f"Unusable or missing observable {obs}: {reference}/{model}")
        result[f"{obs}_FULL_1111"] = full
        result[f"{obs}_FROZEN_1001"] = frozen
        result[f"{obs}_LLSS_OFF"] = float(removed)
        # Preserve original reported running number exactly for thesis equality.
        result[f"{obs}_Running_pct"] = thesis_pct
        result[f"{obs}_Mixing_off_pct"] = percentage(removed, frozen)
        result[f"{obs}_LLSS_effect_pct"] = thesis_pct - result[f"{obs}_Mixing_off_pct"]
    return result


def plot_three_way(rows: list[dict], inventory: tuple[str, ...], out: Path, style: dict):
    # Reuse the exact title-free, default-Matplotlib-font heatmap drawing routine
    # of the previous project figure, not its Wilson-only numerical baseline.
    from scripts.RunT3EFT1LLSSIsolation import _signed_log_heatmap, _write_matrix_csv
    all_rows = {r["model"]: r for r in rows if r.get("status") == "VALID_EXPLORATORY"}
    matrix = np.full((len(inventory), 10), np.nan)
    for i, model in enumerate(inventory):
        row = all_rows.get(model)
        if row is None:
            continue
        for j, obs in enumerate(OBSERVABLES):
            for k, name in enumerate(("Running", "Mixing_off")):
                val = row.get(f"{obs}_{name}_pct")
                if val is not None and math.isfinite(float(val)):
                    matrix[i, 2*j+k] = float(val)
    if not np.isfinite(matrix).any():
        raise RuntimeError("No valid numerical results available for heatmap")
    heads = [r"$\Delta m^2_{21}$", r"$\Delta m^2_{3\ell}$",
             r"$\sin^2\theta_{12}$", r"$\sin^2\theta_{13}$", r"$\sin^2\theta_{23}$"]
    labels = [x for _ in OBSERVABLES for x in ("Running", "Mixing off")]
    model_labels = [m.replace("_", "(") + ")" for m in inventory]
    (out / "plots").mkdir(parents=True, exist_ok=True)
    (out / "csv").mkdir(parents=True, exist_ok=True)
    _write_matrix_csv(out / "csv/all16_heatmap_table.csv", matrix, list(inventory), labels)
    _signed_log_heatmap(matrix, model_labels, labels, heads,
                        out / "plots/all16_heatmap_percent.png", 2, style=style)
    for name, selection, subheads in (("masssq", list(range(4)), heads[:2]),
                                      ("mixing", list(range(4, 10)), heads[2:])):
        sub = matrix[:, selection]
        sublabels = [labels[i] for i in selection]
        _write_matrix_csv(out / f"csv/all16_{name}_heatmap_table.csv",
                          sub, list(inventory), sublabels)
        _signed_log_heatmap(sub, model_labels, sublabels, subheads,
                            out / f"plots/all16_{name}_heatmap_percent.png", 2, style=style)


def run_llss_off(root: Path, model: str, donor: dict, cfg: dict,
                 provenance: dict, *, batch1: str, reuse_cached: bool,
                 keep_trial_files: bool):
    """One unchanged thesis evaluator, with just the LLSS->C5 source removed."""
    from Numerical.parallel_t3.complete.all16_rough_yukawa_fit import FreshAll16Evaluator
    from Numerical.parallel_t3.background import runner as background_runner
    from Numerical.parallel_t3.rough_inputs import builder as builder_module
    from Numerical.parallel_t3.threshold import preflight as preflight_module
    from Numerical.parallel_t3.threshold import rough_numerics as numerics_module
    from Numerical.parallel_t3.wilson import workflow, rough_stream2
    from Numerical.parallel_t3.wilson.basis import build_basis
    from scripts.RunT3All16ExploratoryFit import (
        corrected_quartic_unit_tensors, consistent_project_c5,
        exploratory_preflight, exploratory_rough_case)
    from scripts.RunT3EFT1LLSSIsolation import (
        make_llss_off_beta, make_cached_batch1_verifier)
    from scripts.RunT3GeneralisedRunningComparisons import _cleanup

    evaluator = FreshAll16Evaluator(
        root, model, cfg, batch1=batch1,
        gamma_factor=float(provenance["gamma_factor_assumed"]),
        mu_low=float(provenance["mu_low_gev"]), vev=float(cfg["vev_GeV"]),
        reduced_potential=True, mass_bilinears=provenance["mass_bilinears"],
        wilson_solver=provenance["wilson_solver"],
        wilson_steps=int(provenance["wilson_steps"]), keep_trial_files=True,
    )
    trial_token = None
    capture: dict = {}
    with ExitStack() as stack:
        stack.enter_context(patch.object(builder_module, "quartic_unit_tensors", corrected_quartic_unit_tensors))
        stack.enter_context(patch.object(workflow, "project_c5", consistent_project_c5))
        stack.enter_context(patch.object(rough_stream2, "project_c5", consistent_project_c5))
        stack.enter_context(patch.object(preflight_module, "preflight_pair", exploratory_preflight))
        stack.enter_context(patch.object(numerics_module, "run_rough_case", exploratory_rough_case))
        if reuse_cached:
            verifier = make_cached_batch1_verifier(
                background_runner.verify_batch1_provenance, capture)
            stack.enter_context(patch.object(background_runner, "verify_batch1_provenance", verifier))
        stack.enter_context(patch.object(
            workflow, "beta", make_llss_off_beta(workflow.beta, build_basis(model),
                                                  capture=capture,
                                                  projector=consistent_project_c5)))
        try:
            trial = evaluator(np.asarray(donor["x"], dtype=float).copy())
            trial_token = trial.meta.get("scratch_token")
            if not trial_token:
                raise RuntimeError("Evaluator returned no scratch_token; cannot audit trial")
            if capture.get("rhs_calls", 0) <= 0:
                raise RuntimeError("LLSS-off Wilson beta was never called")
            return dict(trial.observables), float(trial.chi2), capture
        finally:
            if trial_token and not keep_trial_files:
                _cleanup(root, trial_token)


def write_csv(path: Path, rows: list[dict]):
    fields = ["reference", "model", "status", "error", "thesis_source", "chi2_LLSS_OFF"]
    for obs in OBSERVABLES:
        fields.extend(f"{obs}_{x}" for x in
                      ("FULL_1111", "FROZEN_1001", "LLSS_OFF", "Running_pct",
                       "Mixing_off_pct", "LLSS_effect_pct"))
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_reports(out: Path, reference: str, rows: list[dict], inventory: tuple[str, ...]):
    rep = out / "reports"
    rep.mkdir(exist_ok=True)
    heads = ["Δm²₂₁", "Δm²₃ℓ", "sin²θ₁₂", "sin²θ₁₃", "sin²θ₂₃"]
    md = [f"# Thesis-consistent EFT1 LLSS isolation — {reference}", "",
          "Same frozen denominator for both comparisons: original EFT1-off (1001).",
          "", "Running = 100 × (O₁₁₁₁ − O₁₀₀₁) / O₁₀₀₁.",
          "Mixing off = 100 × (O_LLSS-off − O₁₀₀₁) / O₁₀₀₁.", "",
          "The Running numbers are copied from the source thesis stage-freeze CSV and audited against its raw observables. The LLSS-off numbers are freshly calculated.",
          "", "| Model | " + " | ".join(x + " Running / Mixing off (%)" for x in heads) + " |",
          "|---|" + "---:|"*5]
    look = {row["model"]: row for row in rows}
    for model in inventory:
        row = look.get(model, {})
        cells = []
        for obs in OBSERVABLES:
            a = row.get(f"{obs}_Running_pct")
            b = row.get(f"{obs}_Mixing_off_pct")
            cells.append("N/A" if a is None or b is None else f"{a:+.6g} / {b:+.6g}")
        md.append("| " + model + " | " + " | ".join(cells) + " |")
    md.extend(["", "## Scope and interpretation", "",
               "Frozen: UV on; both renormalisable and Wilson EFT1 running off; scalar matching and SMEFT on.",
               "Running: UV, renormalisable EFT1, Wilson EFT1, scalar matching and SMEFT all on.",
               "LLSS off: identical to Running except the quartic LLSS→LLHH Wilson beta source is removed at each RHS evaluation.",
               "Both scalar multiplets are integrated out together at the scalar threshold.",
               "Exploratory only: spinor normalization, full basis map, ADM and independent finite hard matching are not certified.",
               "Thesis-vs-new source mismatches invalidate claims of exact reproducibility; input provenance is logged."])
    (rep / "thesis_consistent_LLSS_report.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    table = [r"\begin{table}[htbp]", r"\centering", r"\scriptsize",
             r"\begin{tabular}{lrrrrrrrrrr}", r"\hline",
             "Model & " + " & ".join(x for _ in OBSERVABLES for x in ("Running", "Mixing off")) + r" \\",
             r"\hline"]
    for model in inventory:
        row = look.get(model, {})
        vals = []
        for obs in OBSERVABLES:
            for nm in ("Running", "Mixing_off"):
                v = row.get(f"{obs}_{nm}_pct")
                vals.append("--" if v is None else f"{float(v):+.4g}")
        table.append(model.replace("_", r"\_") + " & " + " & ".join(vals) + r" \\")
    table += [r"\hline", r"\end{tabular}",
              r"\caption{Thesis-consistent EFT1 stage-freeze and LLSS-mixing ablation. Percent changes relative to the identical original EFT1-off baseline. Exploratory calculation.}",
              r"\end{table}"]
    (rep / "thesis_consistent_LLSS_table.tex").write_text("\n".join(table) + "\n", encoding="utf-8")


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--reference", choices=("A_-4", "D_-2", "both"), default="both")
    p.add_argument("--models", nargs="+", help="Optional model subset for smoke testing")
    p.add_argument("--stagefreeze", type=Path, default=ORIGINAL_STAGE)
    p.add_argument("--tables", type=Path, default=ORIGINAL_TABLES)
    p.add_argument("--reference-results", type=Path, default=ORIGINAL_FITS)
    p.add_argument("--batch1", default="output/batch1_general_contract/batch1_native_v1")
    p.add_argument("--reuse-cached-batch1", action="store_true", default=True,
                   help="Use existing verified beta JSON exports despite the one known changed Wolfram source (default)")
    p.add_argument("--strict-batch1-sources", action="store_true",
                   help="Disable the known source-hash exception")
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--keep-trial-files", action="store_true")
    p.add_argument("--debug", action="store_true")
    p.add_argument("--plot-only", action="store_true", help="Regenerate figures from previous thesis-consistent CSV in --output; no physics")
    p.add_argument("--font-size", type=float, default=14.0)
    for name in ("axis-label", "top-tick", "bottom-tick", "y-tick", "cell", "colorbar"):
        p.add_argument("--" + name + "-fontsize", type=float, default=None,
                       dest=name.replace("-", "_") + "_fontsize")
    p.add_argument("--figure-width", type=float, default=21.0)
    p.add_argument("--row-height", type=float, default=0.57)
    p.add_argument("--dpi", type=int, default=300)
    args = p.parse_args(argv)
    try:
        style = style_from_args(args)
    except ValueError as exc:
        p.error(str(exc))
    from Numerical.parallel_t3.complete.all16_rough_yukawa_fit import (
        BENCHMARK, all_models, read_json)
    from scripts.BuildT3GeneralisedThesisTables import load_reference
    inventory = tuple(all_models())
    models = tuple(args.models or inventory)
    if not models or len(set(models)) != len(models) or any(x not in inventory for x in models):
        p.error("Invalid --models; must be distinct production T3 models")
    refs = DONORS if args.reference == "both" else (args.reference,)
    root = ROOT.resolve()
    out = (root / args.output).resolve()
    required = (root / "output/parallel_t3").resolve()
    if required not in out.parents:
        p.error("--output must be under output/parallel_t3")
    if out.exists() and not args.plot_only:
        p.error("--output already exists; choose a new output path or use --plot-only")
    if args.plot_only and not out.exists():
        p.error("--plot-only requires the existing --output directory")
    source = (root / args.stagefreeze).resolve() / "running_comparisons.csv"
    if not source.is_file():
        p.error(f"Missing original stage-freeze CSV: {source}")
    baselines = stagefreeze_lookup(source)
    if args.plot_only:
        for ref in refs:
            donor_out = out / ref.replace("-", "m").replace("+", "p")
            file = donor_out / "thesis_consistent_three_way.csv"
            if not file.is_file():
                p.error(f"No completed thesis-consistent CSV: {file}")
            with file.open(encoding="utf-8-sig", newline="") as fh:
                rows = list(csv.DictReader(fh))
            for row in rows:
                if row["reference"] != ref:
                    p.error(f"Donor mismatch in {file}")
            plot_three_way(rows, inventory, donor_out, style)
            write_reports(donor_out, ref, rows, inventory)
            print(f"REPLOTTED {ref}: {donor_out}", flush=True)
        return 0
    tables = (root / args.tables).resolve()
    original_provenance = read_json(tables / "provenance.json")
    # Verify that saved original stage-freeze FULL values really share the
    # reference full results used for the subsequent new LLSS-off trial.
    table_results = read_json(tables / "results.json")
    cfg = read_json(root / BENCHMARK)
    donors = {r: load_reference(r, [(root / args.reference_results).resolve()]) for r in refs}
    for r in refs:
        expected_source = original_provenance["reference_benchmarks"][r]["result_json"]
        if str(donors[r]["source"]) != expected_source:
            p.error(f"Reference donor provenance mismatch: {r}")
        for model in models:
            key = (r, model)
            if key not in baselines:
                p.error(f"Thesis EFT1 freeze baseline missing/blocked: {key}")
            old = table_results[r].get(model, {})
            if old.get("status") != "VALID_EXPLORATORY":
                p.error(f"Thesis FULL baseline missing: {key}")
            for obs in OBSERVABLES:
                original_full = old["observables"].get(obs)
                if original_full is None or not math.isclose(
                    float(original_full), baselines[key][obs][0], rel_tol=1e-9, abs_tol=1e-16
                ):
                    p.error(f"Thesis FULL/stage-freeze disagreement: {r}/{model}/{obs}")
    out.mkdir(parents=True)
    failed = 0
    for ref in refs:
        donor_out = out / ref.replace("-", "m").replace("+", "p")
        donor_out.mkdir(parents=True)
        rows = []
        for model in models:
            record = {"reference": ref, "model": model,
                      "status": "BLOCKED", "error": ""}
            print(f"{ref} | {model} | LLSS_OFF (thesis FULL/FROZEN reused)", flush=True)
            try:
                off, chi2, capture = run_llss_off(
                    root, model, donors[ref], cfg, original_provenance,
                    batch1=args.batch1,
                    reuse_cached=bool(args.reuse_cached_batch1) and not args.strict_batch1_sources,
                    keep_trial_files=args.keep_trial_files,
                )
                record = make_row(ref, model, baselines[(ref, model)], off,
                                  source=str(source))
                record["chi2_LLSS_OFF"] = chi2
                (donor_out / f"{model.replace('-', 'm').replace('+', 'p')}_diagnostic.json").write_text(
                    json.dumps({"reference": ref, "model": model,
                                "rhs_calls": capture.get("rhs_calls"),
                                "accepted_stale_source": capture.get("batch1_source_mismatches", {}),
                                "chi2_LLSS_OFF": chi2}, indent=2) + "\n", encoding="utf-8")
            except Exception as exc:
                failed += 1
                record["error"] = f"{type(exc).__name__}: {exc}"
                if args.debug:
                    record["traceback"] = traceback.format_exc(limit=12)
                print(f"BLOCKED {ref}/{model}: {record['error']}", flush=True)
            rows.append(record)
            # Save incrementally so an interrupted run does not destroy results.
            write_csv(donor_out / "thesis_consistent_three_way.csv", rows)
        if any(r["status"] == "VALID_EXPLORATORY" for r in rows):
            plot_three_way(rows, inventory, donor_out, style)
        write_reports(donor_out, ref, rows, inventory)
        (donor_out / "provenance.json").write_text(json.dumps({
            "reference": ref, "stagefreeze_source": str(source),
            "thesis_full": "1111: from original stage-freeze CSV",
            "thesis_frozen": "1001: both renormalisable and Wilson EFT1 frozen, from original stage-freeze CSV",
            "LLSS_off": "same as FULL, subtract only quartic LLSS->LLHH source per Wilson beta RHS call",
            "running_formula": "100*(1111-1001)/1001",
            "mixing_off_formula": "100*(LLSS_OFF-1001)/1001",
            "original_running_reproduced_without_recomputing": True,
            "cached_batch1_exception": not args.strict_batch1_sources,
            "plot_style": style,
            "conditional_exploratory": True,
            "models_completed": sum(r["status"] == "VALID_EXPLORATORY" for r in rows),
            "models_requested": list(models),
        }, indent=2) + "\n", encoding="utf-8")
        print(f"DONE {ref}: {sum(r['status']=='VALID_EXPLORATORY' for r in rows)}/{len(rows)} valid; {donor_out}", flush=True)
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
