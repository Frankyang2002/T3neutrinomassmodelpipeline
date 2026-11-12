#!/usr/bin/env python3
"""Three-way T3 EFT1 WILSON ablation, on an identical evolved scalar background.

Project path: scripts/RunT3EFT1LLSSIsolation.py

Cases (all retain UV, renormalisable EFT1, both scalar matching, SMEFT RG):
    W_FROZEN   seven Wilson beta channels set to zero;
    LLSS_OFF   only quartic LLSS -> LLHH term removed from scalar_pair beta;
    W_FULL     unmodified full Wilson evolution.

Main thesis output:
    a 16x10 heatmap (or selected subset), with five observables and two columns
    per observable:
        1) full running compared with W_FROZEN,
        2) LLSS_OFF compared with W_FROZEN.

This is NOT the earlier stage-frozen EFT1-off heatmap (which also froze the
renormalisable EFT1 background), and is not a certified physical MSbar result.
It uses the existing conditional 45-Weyl, 7-slot candidate ADM and identical
fixed-Yukawa thesis donor inputs in all three cases. No GitHub writes.
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

CASES = ("W_FROZEN", "LLSS_OFF", "W_FULL")
SELECTED_SIX = ("A_+0", "B_+1", "C_+1", "D_-2", "E_+0", "E_-4")
DEFAULT_TABLES = Path("output/parallel_t3/thesis_generalised_tables_v1")
DEFAULT_REFERENCES = Path("output/parallel_t3/A4_Dm2_targeted_fit01")
DEFAULT_OUTPUT = Path("output/parallel_t3/eft1_llss_isolation_v1")
OBSERVABLES = ("dm21_sq_eV2", "dm31_sq_eV2", "sin2_theta12", "sin2_theta13", "sin2_theta23")
HEATMAP_COLUMNS = tuple(
    (observable, comparison)
    for observable in OBSERVABLES
    for comparison in ("pct_vs_W_FROZEN", "pct_LLSS_OFF_vs_W_FROZEN")
)
OBS_LABEL = {
    "dm21_sq_eV2": r"Δm²21",
    "dm31_sq_eV2": r"Δm²3ℓ",
    "sin2_theta12": r"sin²θ12",
    "sin2_theta13": r"sin²θ13",
    "sin2_theta23": r"sin²θ23",
}
COMP_LABEL = {
    "pct_vs_W_FROZEN": "Running",
    "pct_LLSS_OFF_vs_W_FROZEN": "Mixing off",
}
LOOP = 1.0 / (16.0 * math.pi**2)

# The one known stale-source mismatch seen in the existing Batch 1 inputs.
# This is NOT permission to accept arbitrary source or cached-beta changes.
ALLOWED_CACHED_BATCH1_SOURCE_MISMATCH = (
    "RGE/running/rgbeta/RunGeneralT3QuadraticRGBeta.wl",
)


def make_cached_batch1_verifier(original_verifier, audit: dict):
    """Use existing Batch 1 exported beta JSONs, with an explicit source-only exception.

    This invokes the repository's original verifier with strict_sources=False.
    It still verifies the cached UV and EFT1 export SHA256, mandatory portal
    checks, model identities, and that every recorded source file exists.
    Only the specific changed Mathematica source above may differ. All
    other required-for-B-reuse source mismatches remain hard errors.
    """
    def verified(root, batch_root, model, record, uv_path, eft_path, *, strict_sources=True):
        result = original_verifier(
            root, batch_root, model, record, uv_path, eft_path, strict_sources=False
        )
        observed = result["local_sources"]
        mismatched = {
            path: item for path, item in observed.items()
            if item["required_for_B_reuse"] and not item["same_as_batch1"]
        }
        disallowed = sorted(set(mismatched) - set(ALLOWED_CACHED_BATCH1_SOURCE_MISMATCH))
        if disallowed:
            raise ValueError(
                f"{model}: unapproved Batch 1 source mismatches: {disallowed}. "
                "Only the explicitly identified RGBeta Wolfram driver may differ."
            )
        audit["cached_batch1_export_hashes_verified"] = True
        audit["batch1_source_mismatches"] = {
            path: {"current_sha256": item["sha256"],
                   "batch1_sha256": item["batch1_sha256"]}
            for path, item in mismatched.items()
        }
        audit["batch1_source_mismatch_accepted_for_exploratory_use"] = bool(mismatched)
        return result
    return verified


def llss_to_llhh_quartic_source(c: np.ndarray, lam: np.ndarray, basis) -> np.ndarray:
    """Return ONLY [lambda_HHSS C_LLSS] in 16pi^2 beta_C's LLHH block.

    The heavy scalars are 4:ns, Higgs real components 0:4, and SM leptons
    are the three explicit SU(2) doublets in basis.leptons.  The original
    scalar_pair equation contains 0.5*(lambda_abcd*C_ijcd + a<->b).
    We remove just its LLSS -> LLHH quartic source; all other tensor entries
    and all seven other contributions remain unmodified.
    """
    c = np.asarray(c, dtype=np.complex128)
    lam = np.asarray(lam, dtype=float)
    expected = (basis.nf, basis.nf, basis.ns, basis.ns)
    if c.shape != expected or lam.shape != (basis.ns,) * 4:
        raise ValueError("LLSS ablation: tensor/background basis mismatch")
    if basis.ns <= 4 or tuple(basis.higgs) != (0, 1, 2, 3):
        raise ValueError("LLSS ablation: expected first four real scalars to be Higgs")
    leptons = [f for pair in basis.leptons for f in pair]
    if len(leptons) != 6 or len(set(leptons)) != 6:
        raise ValueError("LLSS ablation: expected three complete lepton doublets")
    H = list(basis.higgs)
    S = list(range(4, basis.ns))
    cab = c[np.ix_(leptons, leptons, S, S)]
    q = lam[np.ix_(H, H, S, S)]
    block = np.einsum("abAB,ijAB->ijab", q, cab, optimize=True)
    block = 0.5 * (block + block.swapaxes(2, 3))
    result = np.zeros_like(c)
    result[np.ix_(leptons, leptons, H, H)] = block
    return result


def make_llss_off_beta(original_beta, basis, *, capture=None, projector=None):
    """Remove only the LLSS quartic source of the LLHH Wilson beta.

    This subtraction is made inside the ODE right-hand side at EACH time,
    using the evolving C_LLSS and lambda_HHSS.  It does not subtract a fixed
    post-hoc logarithm and does not freeze C_LLSS self-running.
    """
    calls = 0

    def beta_isolated(c, y, lam, gauge=(), **kwargs):
        nonlocal calls
        channels, tag = original_beta(c, y, lam, gauge, **kwargs)
        if "scalar_pair" not in channels:
            raise RuntimeError("Required scalar_pair channel missing")
        source = llss_to_llhh_quartic_source(c, lam, basis)
        if channels["scalar_pair"].shape != source.shape:
            raise RuntimeError("Unexpected scalar_pair beta tensor shape")
        channels = dict(channels)
        channels["scalar_pair"] = channels["scalar_pair"] - source
        calls += 1
        if capture is not None:
            capture["rhs_calls"] = calls
            if "beta_LLSS_to_C5_initial" not in capture and projector is not None:
                c5, residual = projector(source, basis)
                capture["beta_LLSS_to_C5_initial"] = np.asarray(c5, complex).copy()
                capture["source_projection_absolute_residual_initial"] = float(residual)
                capture["source_full_tensor_norm_initial"] = float(np.linalg.norm(source))
        return channels, tag + "+LLSS_to_LLHH_quartic_source_removed"

    return beta_isolated


def norm(matrix) -> float:
    return float(np.linalg.norm(np.asarray(matrix, dtype=complex)))


def _read_c5_arrays(root: Path, token: str, model: str) -> dict:
    from Numerical.parallel_t3.background.contracts import tag

    slug = tag(model)
    wpath = root / "output/parallel_t3/rough_wilson" / token / slug / "wilson_at_MS.npz"
    mpath = root / "output/parallel_t3/rough_results" / token / slug / "C5_smeft_at_MS.npz"
    lpath = root / "output/parallel_t3/rough_results" / token / slug / "C5_smeft_low.npz"
    with np.load(wpath, allow_pickle=False) as data:
        c5_eft1 = np.asarray(data["C5_EFT1_ms"], complex).copy()
    with np.load(mpath, allow_pickle=False) as data:
        c5_match = np.asarray(data["C5_scalar_bubble_MS"], complex).copy()
        c5_smeft_ms = np.asarray(data["C5_SMEFT_MS"], complex).copy()
    with np.load(lpath, allow_pickle=False) as data:
        c5_low = np.asarray(data["C5_SMEFT_low"], complex).copy()
    matrices = {"eft1": c5_eft1, "match": c5_match,
                "smeft_ms": c5_smeft_ms, "low": c5_low}
    for name, value in matrices.items():
        if value.shape != (3, 3) or not np.isfinite(value).all():
            raise RuntimeError(f"Bad {name} C5 matrix in {model}/{token}")
    if not np.allclose(c5_eft1 + c5_match, c5_smeft_ms, rtol=1e-7, atol=1e-23):
        raise RuntimeError("C5 threshold bookkeeping mismatch; additional hard source? ")
    return matrices


def _finite_ratio(num: float, den: float):
    return (float(num) / float(den)) if den > 0 else None


def _cosine(a: np.ndarray, b: np.ndarray):
    den = norm(a) * norm(b)
    if den == 0:
        return None
    return float(np.vdot(a, b).real / den)


def _percentage(full, frozen):
    return 100.0 * (full - frozen) / frozen if frozen is not None and frozen != 0 else None


def summarize_model(case_data: dict, initial_beta=None, mf=1e5, ms=1e3) -> tuple[dict, dict]:
    frozen, removed, full = (case_data[k]["matrices"] for k in CASES)
    d_llss = full["eft1"] - removed["eft1"]
    d_all = full["eft1"] - frozen["eft1"]
    log_estimate = (math.log(ms / mf) * LOOP * initial_beta
                    if initial_beta is not None else None)
    result = {
        "C5_LLSS_effect_norm_GeVinv": norm(d_llss),
        "C5_all_Wilson_effect_norm_GeVinv": norm(d_all),
        "LLSS_fraction_norm": _finite_ratio(norm(d_llss), norm(d_all)),
        "LLSS_all_Wilson_alignment_cosine": _cosine(d_llss, d_all),
        "C5_LLSS_over_scalar_matching": _finite_ratio(norm(d_llss), norm(full["match"])),
        "C5_running_over_matching": _finite_ratio(norm(full["eft1"]), norm(full["match"])),
        "C5_LLSS_leading_log_estimate_norm_GeVinv":
            norm(log_estimate) if log_estimate is not None else None,
        "leading_log_vs_exact_LLSS_relative_matrix_residual":
            _finite_ratio(norm(log_estimate - d_llss), norm(d_llss))
            if log_estimate is not None else None,
        "matching_change_LLSS_off_relative":
            _finite_ratio(norm(full["match"] - removed["match"]), norm(full["match"])),
        "matching_change_W_frozen_relative":
            _finite_ratio(norm(full["match"] - frozen["match"]), norm(full["match"])),
    }
    for label in CASES:
        result[f"C5_EFT1_{label}_norm_GeVinv"] = norm(case_data[label]["matrices"]["eft1"])
        result[f"C5_matching_{label}_norm_GeVinv"] = norm(case_data[label]["matrices"]["match"])
        result[f"C5_SMEFT_MS_{label}_norm_GeVinv"] = norm(case_data[label]["matrices"]["smeft_ms"])
        result[f"chi2_{label}"] = case_data[label]["chi2"]
    for observable in OBSERVABLES:
        val_full = case_data["W_FULL"]["observables"].get(observable)
        val_frozen = case_data["W_FROZEN"]["observables"].get(observable)
        val_off = case_data["LLSS_OFF"]["observables"].get(observable)
        result[f"{observable}_W_FULL"] = val_full
        result[f"{observable}_W_FROZEN"] = val_frozen
        result[f"{observable}_LLSS_OFF"] = val_off
        result[f"{observable}_pct_vs_W_FROZEN"] = _percentage(val_full, val_frozen)
        result[f"{observable}_pct_LLSS_OFF_vs_W_FROZEN"] = _percentage(val_off, val_frozen)
        # Same denominator makes the difference equal to the isolated LLSS
        # effect on this observable, measured relative to frozen Wilson running.
        full_pct = result[f"{observable}_pct_vs_W_FROZEN"]
        off_pct = result[f"{observable}_pct_LLSS_OFF_vs_W_FROZEN"]
        result[f"{observable}_pct_LLSS_effect_vs_W_FROZEN"] = (
            full_pct - off_pct if full_pct is not None and off_pct is not None else None
        )
    matrices = {
        "C5_effect_LLSS": d_llss,
        "C5_effect_all_Wilson": d_all,
    }
    for case in CASES:
        for name, value in case_data[case]["matrices"].items():
            matrices[f"C5_{name}_{case}"] = value
    if log_estimate is not None:
        matrices["C5_LLSS_leading_log_estimate"] = log_estimate
    return result, matrices


def _run_one_model(root, model, reference, x, cfg, provenance, *, solver, steps,
                   batch1, keep_trial_files, debug, reuse_cached_batch1=False):
    from Numerical.parallel_t3.complete.all16_rough_yukawa_fit import FreshAll16Evaluator
    from Numerical.parallel_t3.background import runner as background_runner
    from Numerical.parallel_t3.rough_inputs import builder as builder_module
    from Numerical.parallel_t3.threshold import preflight as preflight_module
    from Numerical.parallel_t3.threshold import rough_numerics as numerics_module
    from Numerical.parallel_t3.wilson import workflow, rough_stream2
    from Numerical.parallel_t3.wilson.basis import build_basis
    from scripts import RunT3GeneralisedRunningComparisons as general
    from scripts.RunT3All16ExploratoryFit import (
        corrected_quartic_unit_tensors, consistent_project_c5,
        exploratory_preflight, exploratory_rough_case)

    basis = build_basis(model)
    outcomes = {}
    capture = {}
    for case in CASES:
        evaluator = FreshAll16Evaluator(
            root, model, cfg, batch1=batch1,
            gamma_factor=float(provenance["gamma_factor_assumed"]),
            mu_low=float(provenance["mu_low_gev"]), vev=float(cfg["vev_GeV"]),
            reduced_potential=True, mass_bilinears=provenance["mass_bilinears"],
            wilson_solver=solver, wilson_steps=steps, keep_trial_files=True)
        trial_token = None
        print(f"{reference} | {model} | {case}", flush=True)
        try:
            with ExitStack() as stack:
                stack.enter_context(patch.object(builder_module, "quartic_unit_tensors", corrected_quartic_unit_tensors))
                stack.enter_context(patch.object(workflow, "project_c5", consistent_project_c5))
                stack.enter_context(patch.object(rough_stream2, "project_c5", consistent_project_c5))
                stack.enter_context(patch.object(preflight_module, "preflight_pair", exploratory_preflight))
                stack.enter_context(patch.object(numerics_module, "run_rough_case", exploratory_rough_case))
                if reuse_cached_batch1:
                    if case == CASES[0]:
                        print("  EXPLORATORY CACHED BATCH1: source hash exception for "
                              + ALLOWED_CACHED_BATCH1_SOURCE_MISMATCH[0] + "; "
                              "UV/EFT1 cached-beta SHA256 checks remain active", flush=True)
                    verifier = make_cached_batch1_verifier(
                        background_runner.verify_batch1_provenance, capture
                    )
                    stack.enter_context(patch.object(
                        background_runner, "verify_batch1_provenance", verifier))
                if case == "W_FROZEN":
                    stack.enter_context(patch.object(workflow, "beta", general._frozen_wilson_beta))
                elif case == "LLSS_OFF":
                    original_beta = workflow.beta
                    isolated = make_llss_off_beta(original_beta, basis, capture=capture,
                                                 projector=consistent_project_c5)
                    stack.enter_context(patch.object(workflow, "beta", isolated))
                trial = evaluator(np.asarray(x, float).copy())
                trial_token = trial.meta.get("scratch_token")
            if not trial_token:
                raise RuntimeError("Expected trial-owned files for C5 extraction")
            matrices = _read_c5_arrays(root, trial_token, model)
            outcomes[case] = {
                "matrices": matrices, "chi2": float(trial.chi2),
                "observables": trial.observables, "meta": trial.meta,
            }
        finally:
            if trial_token and not keep_trial_files:
                general._cleanup(root, trial_token)
    return outcomes, capture


def _write_csv(path, rows):
    fields = ["reference", "model", "status", "error"]
    fields.extend(k for row in rows for k in row if k not in fields)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _write_matrix_csv(path: Path, matrix: np.ndarray, rows: list[str], cols: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["model", *cols])
        for model, values in zip(rows, matrix, strict=True):
            writer.writerow([model, *[("" if not np.isfinite(v) else repr(float(v))) for v in values]])


def _plot_style(args) -> dict:
    style = {
        "axis_label_fontsize": float(args.axis_label_fontsize),
        "bottom_tick_fontsize": float(args.bottom_tick_fontsize),
        "top_tick_fontsize": float(args.top_tick_fontsize),
        "y_tick_fontsize": float(args.y_tick_fontsize),
        "cell_fontsize": float(args.cell_fontsize),
        "colorbar_fontsize": float(args.colorbar_fontsize),
        "figure_width": float(args.figure_width),
        "row_height": float(args.row_height),
        "dpi": int(args.dpi),
    }
    for name, value in style.items():
        if name == "dpi":
            if not 72 <= value <= 1200:
                raise ValueError("--dpi must be 72..1200")
        elif not (math.isfinite(value) and value > 0):
            raise ValueError(f"Invalid positive plot style option {name}={value}")
    return style


def _build_heatmap_matrix(rows: list[dict], model_order: tuple[str, ...]) -> tuple[np.ndarray, list[str], list[str]]:
    lookup = {row["model"]: row for row in rows if row.get("status") == "VALID_EXPLORATORY"}
    valid_models = list(model_order)  # Preserve all 16 rows; blocked entries stay blank.
    col_labels = [COMP_LABEL[comp] for _, comp in HEATMAP_COLUMNS]
    matrix = np.full((len(valid_models), len(HEATMAP_COLUMNS)), np.nan, float)
    for i, model in enumerate(valid_models):
        row = lookup.get(model, {})
        for j, (obs, comp) in enumerate(HEATMAP_COLUMNS):
            key = f"{obs}_{comp}"
            value = row.get(key)
            try:
                v = float(value)
            except (TypeError, ValueError):
                continue
            if math.isfinite(v):
                matrix[i, j] = v
    return matrix, valid_models, col_labels


def _signed_log_heatmap(matrix: np.ndarray, row_labels: list[str], column_labels: list[str],
                        group_labels: list[str], output_path: Path, group_size: int, *,
                        linthresh: float = 1.0, style: dict | None = None):
    """Title-free thesis heatmap: observable headings on top, case labels below."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import SymLogNorm
    from matplotlib.ticker import FixedLocator, FuncFormatter

    style = dict(style or {})
    axis_fontsize = style.get("axis_label_fontsize", 13)
    bottom_fontsize = style.get("bottom_tick_fontsize", 10)
    top_fontsize = style.get("top_tick_fontsize", 13)
    y_fontsize = style.get("y_tick_fontsize", 11)
    cell_fontsize = style.get("cell_fontsize", 8)
    cb_fontsize = style.get("colorbar_fontsize", 10)
    width = style.get("figure_width", 15.0)
    height = max(5.5, style.get("row_height", 0.49) * len(row_labels) + 2.0)
    dpi = style.get("dpi", 300)

    a = np.asarray(matrix, float)
    valid = np.isfinite(a)
    if not valid.any():
        raise ValueError(f"No finite heatmap cells for {output_path}")
    if a.shape[1] != len(column_labels) or a.shape[0] != len(row_labels):
        raise ValueError("Heatmap axis length does not match matrix")
    if a.shape[1] != group_size * len(group_labels):
        raise ValueError("Heatmap observable heading count differs from paired columns")
    vmax = max(float(np.max(np.abs(a[valid]))), linthresh)
    norm = SymLogNorm(linthresh=linthresh, linscale=1.0,
                      vmin=-vmax, vmax=vmax, base=10)

    # Respect Matplotlib's default text and math fonts; only set export and
    # axis-formatting options required for consistent thesis figures.
    with plt.rc_context({"axes.formatter.use_mathtext": True,
                         "axes.unicode_minus": False,
                         "pdf.fonttype": 42, "ps.fonttype": 42,
                         "svg.fonttype": "none"}):
        cmap = plt.get_cmap("RdBu_r").copy()
        cmap.set_bad("white")
        fig = plt.figure(figsize=(width, height))
        # A reserved colourbar column prevents the colourbar/twin x-axis from
        # overlaying the rightmost "Mixing off" data column.
        gs = fig.add_gridspec(1, 2, width_ratios=(33, 1),
                              left=0.075, right=0.945,
                              bottom=0.090, top=0.925, wspace=0.055)
        ax = fig.add_subplot(gs[0, 0])
        cax = fig.add_subplot(gs[0, 1])
        image = ax.imshow(np.ma.masked_invalid(a), cmap=cmap, norm=norm, aspect="auto")
        # No figure title and no bottom x-axis title; the caption owns the explanation.
        ax.set_xticks(np.arange(a.shape[1]))
        ax.set_xticklabels(column_labels, fontsize=bottom_fontsize)
        ax.set_yticks(np.arange(a.shape[0]))
        ax.set_yticklabels(row_labels, fontsize=y_fontsize)
        ax.set_ylabel("T3 model", fontsize=axis_fontsize)
        ax.tick_params(axis="x", length=0, pad=9)
        ax.tick_params(axis="y", length=0, pad=7)

        centres = [i*group_size + 0.5*(group_size-1) for i in range(len(group_labels))]
        upper = ax.twiny()
        upper.set_xlim(ax.get_xlim())
        upper.set_xticks(centres)
        upper.set_xticklabels(group_labels, fontsize=top_fontsize)
        upper.tick_params(axis="x", length=0, pad=10)
        upper.spines["top"].set_visible(False)

        for i in range(a.shape[0]):
            for j in range(a.shape[1]):
                value = a[i, j]
                if not np.isfinite(value):
                    continue
                magnitude = abs(value)
                if magnitude >= 1e4:
                    label = f"{value:+.1e}"
                elif magnitude >= 10:
                    label = f"{value:+.1f}"
                elif magnitude >= 1:
                    label = f"{value:+.2f}"
                elif magnitude >= 0.1:
                    label = f"{value:+.3f}"
                elif magnitude >= 0.01:
                    label = f"{value:+.4f}"
                else:
                    label = f"{value:+.2e}"
                r, g, b, _ = cmap(norm(value))
                luminance = 0.2126*r + 0.7152*g + 0.0722*b
                ax.text(j, i, label, ha="center", va="center", fontsize=cell_fontsize,
                        color="white" if luminance < 0.48 else "black")
        ax.set_xticks(np.arange(-0.5, a.shape[1], 1), minor=True)
        ax.set_yticks(np.arange(-0.5, a.shape[0], 1), minor=True)
        ax.grid(which="minor", color="white", linewidth=0.7)
        ax.tick_params(which="minor", bottom=False, left=False)
        for x in range(group_size, a.shape[1], group_size):
            ax.axvline(x-0.5, color="black", linewidth=1.3)
        top = int(math.floor(math.log10(vmax))) if vmax >= 10 else 0
        decades = [10.0**i for i in range(top+1)]
        ticks = sorted(set([-d for d in decades] + [0.0] + decades))
        ticks = [v for v in ticks if abs(v) <= vmax]
        cbar = fig.colorbar(image, cax=cax)
        cbar.locator = FixedLocator(ticks)
        cbar.formatter = FuncFormatter(
            lambda v, pos: "0" if v == 0 else (f"{v:.0e}" if abs(v) >= 1e4 else f"{v:g}"))
        cbar.update_ticks()
        cbar.set_label("Signed change (%)", fontsize=cb_fontsize)
        cbar.ax.tick_params(labelsize=cb_fontsize)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        fig.savefig(output_path.with_suffix(".pdf"), bbox_inches="tight")
        plt.close(fig)


def _render_heatmaps(out: Path, rows: list[dict], model_order: tuple[str, ...],
                     stem: str, style: dict | None = None) -> None:
    matrix, models, col_labels = _build_heatmap_matrix(rows, model_order)
    if not models or not np.isfinite(matrix).any():
        return
    short_rows = [model.replace("_", "(") + ")" for model in models]
    plots, tables = out / "plots", out / "csv"
    observable_headers = [
        r"$\Delta m^2_{21}$", r"$\Delta m^2_{3\ell}$",
        r"$\sin^2\theta_{12}$", r"$\sin^2\theta_{13}$", r"$\sin^2\theta_{23}$",
    ]
    _write_matrix_csv(tables / f"{stem}_heatmap_table.csv", matrix, models, col_labels)
    _signed_log_heatmap(
        matrix, short_rows, col_labels, observable_headers,
        plots / f"{stem}_heatmap_percent.png", 2, style=style)

    # Split plots remain available for readability in thesis pages.
    for name, ix, headings in (
        ("masssq", [0, 1, 2, 3], observable_headers[:2]),
        ("mixing", [4, 5, 6, 7, 8, 9], observable_headers[2:]),
    ):
        subset = matrix[:, ix]
        labels = [col_labels[i] for i in ix]
        _write_matrix_csv(tables / f"{stem}_{name}_heatmap_table.csv", subset, models, labels)
        _signed_log_heatmap(subset, short_rows, labels, headings,
                            plots / f"{stem}_{name}_heatmap_percent.png", 2, style=style)


def _write_long_heatmap_csv(path: Path, rows: list[dict], model_order: tuple[str, ...]) -> None:
    valid = {row["model"]: row for row in rows if row.get("status") == "VALID_EXPLORATORY"}
    fieldnames = ("model", "observable", "comparison", "relative_change_percent")
    output_rows = []
    for model in model_order:
        row = valid.get(model)
        if row is None:
            continue
        for observable, comparison in HEATMAP_COLUMNS:
            key = f"{observable}_{comparison}"
            value = row.get(key)
            output_rows.append({
                "model": model,
                "observable": observable,
                "comparison": comparison,
                "relative_change_percent": "" if value is None else repr(float(value)),
            })
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(output_rows)


def _rederive_heatmap_percentages(rows: list[dict]) -> list[dict]:
    """Regenerate both plotted comparisons from the THREE raw observables.

    Supports the earlier CSV with full-vs-off percentages and ensures replot
    does not depend on those obsolete fields. A missing/zero frozen denominator
    produces a blank cell rather than an artificial zero.
    """
    result = []
    for original in rows:
        row = dict(original)
        if row.get("status") != "VALID_EXPLORATORY":
            result.append(row)
            continue
        for obs in OBSERVABLES:
            values = []
            for case in CASES:
                raw = row.get(f"{obs}_{case}")
                try:
                    value = float(raw)
                except (TypeError, ValueError):
                    value = float("nan")
                values.append(value)
            frozen, off, full = values
            if not all(math.isfinite(x) for x in values) or frozen == 0:
                row[f"{obs}_pct_vs_W_FROZEN"] = None
                row[f"{obs}_pct_LLSS_OFF_vs_W_FROZEN"] = None
                row[f"{obs}_pct_LLSS_effect_vs_W_FROZEN"] = None
            else:
                row[f"{obs}_pct_vs_W_FROZEN"] = _percentage(full, frozen)
                row[f"{obs}_pct_LLSS_OFF_vs_W_FROZEN"] = _percentage(off, frozen)
                row[f"{obs}_pct_LLSS_effect_vs_W_FROZEN"] = 100.0 * (full - off) / frozen
        result.append(row)
    return result


def _fmt_pct(value) -> str:
    try:
        v = float(value)
    except (TypeError, ValueError):
        return "--"
    if not math.isfinite(v):
        return "--"
    return f"{v:+.3e}" if abs(v) >= 10000 or (0 < abs(v) < 0.001) else f"{v:+.3f}"


def _write_donor_report(out: Path, rows: list[dict], inventory: tuple[str, ...],
                        reference: str, style: dict) -> None:
    """Generate a numerical report from existing three-way observables only."""
    tag = reference.replace("-", "m").replace("+", "p")
    display_ref = reference.replace("_", "(") + ")"
    report_dir = out / "reports"
    report_dir.mkdir(parents=True, exist_ok=True)
    indexed = {row.get("model"): row for row in rows}
    valid = [m for m in inventory if indexed.get(m, {}).get("status") == "VALID_EXPLORATORY"]
    blocked = [m for m in inventory if m not in valid]
    headings = {
        "dm21_sq_eV2": r"$\Delta m^2_{21}$",
        "dm31_sq_eV2": r"$\Delta m^2_{3\ell}$",
        "sin2_theta12": r"$\sin^2\theta_{12}$",
        "sin2_theta13": r"$\sin^2\theta_{13}$",
        "sin2_theta23": r"$\sin^2\theta_{23}$",
    }
    lines = [
        f"# EFT1 LLSS mixing isolation: {display_ref} donor", "",
        f"Completed exploratory models: **{len(valid)}/{len(inventory)}**.", "",
        "## Method and percentage definitions", "",
        "All three cases use the same UV-to-fermion matching, renormalisable EFT1",
        "background, simultaneous scalar-threshold matching, and SMEFT evolution.",
        "Only the dimension-five Wilson evolution is changed: frozen (all Wilson beta",
        "functions zero), mixing off (only the quartic LLSS-to-LLHH source removed),",
        "or running (complete candidate Wilson beta).", "",
        r"$\delta_{\rm Running}=100(O_{\rm Running}-O_{\rm Frozen})/O_{\rm Frozen}$.", "",
        r"$\delta_{\rm MixingOff}=100(O_{\rm MixingOff}-O_{\rm Frozen})/O_{\rm Frozen}$.", "",
        "Both are signed percentages with a **common frozen denominator**.",
        "The difference between these columns isolates the observable response to",
        "restoring the LLSS-to-Weinberg quartic mixing term.", "",
        "## Numerical comparison", "",
    ]
    report_data = {"reference": reference, "completed": len(valid), "requested": len(inventory),
                   "blocked": blocked, "observables": {}}
    tex = [
        r"% Generated numerical report; exploratory results only.",
        rf"\subsection*{{EFT1 LLSS mixing isolation: {display_ref} donor}}",
        r"All percentage differences use the same frozen-Wilson denominator:",
        r"\[\delta_{\mathrm{Running}}=100\frac{O_{\mathrm{Running}}-O_{\mathrm{Frozen}}}{O_{\mathrm{Frozen}}},\qquad"
        r"\delta_{\mathrm{MixingOff}}=100\frac{O_{\mathrm{MixingOff}}-O_{\mathrm{Frozen}}}{O_{\mathrm{Frozen}}}.\]",
    ]
    for group_title, obs_group in (
        ("Mass-squared differences", OBSERVABLES[:2]),
        ("Mixing-angle observables", OBSERVABLES[2:]),
    ):
        cols = ["Model"]
        for obs in obs_group:
            cols += [f"{headings[obs]} Running", f"{headings[obs]} Mixing off"]
        lines.extend([f"### {group_title}", "", "| " + " | ".join(cols) + " |",
                      "| " + " | ".join(["---"] + ["---:"] * (len(cols)-1)) + " |"])
        tex.extend([rf"\paragraph{{{group_title}}}", r"\begin{center}",
                    r"\begin{tabular}{l" + "r" * (len(cols)-1) + "}", r"\hline",
                    "Model & " + " & ".join(
                        [rf"{headings[o]} R" if kind == "R" else rf"{headings[o]} M-off"
                         for o in obs_group for kind in ("R", "M")]) + r" \\", r"\hline"])
        for model in inventory:
            row = indexed.get(model, {})
            obs_vals = []
            for obs in obs_group:
                obs_vals.extend([
                    _fmt_pct(row.get(f"{obs}_pct_vs_W_FROZEN")),
                    _fmt_pct(row.get(f"{obs}_pct_LLSS_OFF_vs_W_FROZEN")),
                ])
            label = model.replace("_", "(") + ")"
            lines.append("| " + " | ".join([label, *obs_vals]) + " |")
            tex.append(label + " & " + " & ".join(obs_vals) + r" \\")
        lines.append("")
        tex.extend([r"\hline", r"\end{tabular}", r"\end{center}"])
    lines.extend(["## Observable-level diagnostics", "",
                  "| Observable | Median abs(Running) (%) | Median abs(Mixing off) (%) |",
                  "| --- | ---: | ---: |"])
    for obs in OBSERVABLES:
        a, b = [], []
        for m in valid:
            row = indexed[m]
            try:
                full = float(row.get(f"{obs}_pct_vs_W_FROZEN"))
                off = float(row.get(f"{obs}_pct_LLSS_OFF_vs_W_FROZEN"))
            except (TypeError, ValueError):
                continue
            if math.isfinite(full) and math.isfinite(off):
                a.append(abs(full))
                b.append(abs(off))
        med_a = float(np.median(a)) if a else None
        med_b = float(np.median(b)) if b else None
        report_data["observables"][obs] = {
            "n_finite": len(a), "median_abs_running_percent": med_a,
            "median_abs_mixing_off_percent": med_b,
        }
        lines.append(f"| {headings[obs]} | {_fmt_pct(med_a)} | {_fmt_pct(med_b)} |")
    lines.extend(["", "## Figures", "",
                  "- `plots/all16_heatmap_percent.png` (and `.pdf`): 16 x 10 comparison.",
                  "- `plots/all16_masssq_heatmap_percent.png` (and `.pdf`).",
                  "- `plots/all16_mixing_heatmap_percent.png` (and `.pdf`).", "",
                  "## Scope and limitations", "",
                  "These are conditional exploratory numerical results. The chosen candidate",
                  "Wilson anomalous dimension, source-to-Weyl normalisation, scalar basis",
                  "maps and independent hard one-loop matching are not fully certified.",
                  "The `Frozen` case freezes the Wilson RGEs only: the renormalisable",
                  "EFT1 scalar evolution is active in every comparison. The resulting",
                  "plots are not numerically interchangeable with the original EFT1-off",
                  "figure that also freezes the scalar background.", "",
                  "Large mass-squared changes can coexist with very small mixing-angle",
                  "changes when the flavour matrix is approximately rescaled.", "",
                  "A conclusion that LLSS dominates is justified only when the Running",
                  "column is much larger than the Mixing-off column for the relevant",
                  "observables; exceptions must be reported rather than discarded.", "",
                  f"Blocked or missing models: {', '.join(blocked) if blocked else 'none'}.", ""])
    tex.extend([r"\noindent\emph{Exploratory, uncertified Wilson/matching normalisation.} ",
                r"All three conditions retain renormalisable EFT1 running; only the Wilson source is isolated."])
    (report_dir / f"eft1_llss_{tag}_report.md").write_text("\n".join(lines), encoding="utf-8")
    (report_dir / f"eft1_llss_{tag}_report.tex").write_text("\n".join(tex) + "\n", encoding="utf-8")
    report_data["plot_style"] = style
    (report_dir / f"eft1_llss_{tag}_summary.json").write_text(
        json.dumps(report_data, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--reference", choices=("A_-4", "D_-2"), default="A_-4")
    p.add_argument("--models", nargs="+", help="Subset, e.g. A_+0 B_+1 C_+1 D_-2 E_+0 E_-4")
    p.add_argument("--selected-six", action="store_true", help="Run exactly the six comparison models")
    p.add_argument("--tables", type=Path, default=DEFAULT_TABLES)
    p.add_argument("--reference-results", type=Path, default=DEFAULT_REFERENCES)
    p.add_argument("--batch1", default="output/batch1_general_contract/batch1_native_v1")
    p.add_argument("--reuse-cached-batch1", action="store_true",
                   help="Exploratory use of existing hashed UV/EFT1 beta exports despite "
                        "the known RunGeneralT3QuadraticRGBeta.wl source mismatch; "
                        "still verifies cached beta JSON SHA256 and other source hashes")
    p.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    p.add_argument("--keep-trial-files", action="store_true")
    p.add_argument("--debug", action="store_true")
    p.add_argument("--figure-width", type=float, default=16.0)
    p.add_argument("--row-height", type=float, default=0.48)
    p.add_argument("--dpi", type=int, default=300)
    p.add_argument("--axis-label-fontsize", type=float, default=12.0)
    p.add_argument("--bottom-tick-fontsize", type=float, default=10.0)
    p.add_argument("--top-tick-fontsize", type=float, default=13.0)
    p.add_argument("--y-tick-fontsize", type=float, default=11.0)
    p.add_argument("--cell-fontsize", type=float, default=8.0)
    p.add_argument("--colorbar-fontsize", type=float, default=10.0)
    p.add_argument("--plot-only", type=Path, metavar="EXISTING_RESULTS_DIR",
                   help="Replot saved three-way CSV with common frozen denominator; does not rerun physics")
    args = p.parse_args(argv)
    try:
        style = _plot_style(args)
    except ValueError as exc:
        p.error(str(exc))

    if args.plot_only is not None:
        source = (ROOT / args.plot_only).resolve()
        out = (ROOT / args.output).resolve()
        owned_root = (ROOT / "output/parallel_t3").resolve()
        if out == owned_root or owned_root not in out.parents or out.exists():
            p.error("--output must be a NEW directory under output/parallel_t3")
        src_csv = source / "all16_three_way.csv"
        if not src_csv.is_file():
            p.error(f"--plot-only requires an existing all16_three_way.csv: {src_csv}")
        with src_csv.open(newline="", encoding="utf-8-sig") as f:
            imported = list(csv.DictReader(f))
        from Numerical.parallel_t3.complete.all16_rough_yukawa_fit import all_models
        inventory = tuple(all_models())
        rows = _rederive_heatmap_percentages(imported)
        out.mkdir(parents=True)
        _write_csv(out / "all16_three_way.csv", rows)
        _write_csv(out / "selected6_three_way.csv", [r for r in rows if r["model"] in SELECTED_SIX])
        _write_long_heatmap_csv(out / "csv" / "all16_heatmap_long.csv", rows, inventory)
        _write_long_heatmap_csv(out / "csv" / "selected6_heatmap_long.csv", rows, SELECTED_SIX)
        source_meta = source / "provenance.json"
        provenance = json.loads(source_meta.read_text(encoding="utf-8-sig")) if source_meta.is_file() else {}
        ref = provenance.get("reference", "unknown")
        _render_heatmaps(out, rows, inventory, "all16", style=style)
        _render_heatmaps(out, rows, SELECTED_SIX, "selected6", style=style)
        _write_donor_report(out, rows, inventory, ref, style)
        (out / "provenance.json").write_text(json.dumps({
            "replot_of": str(source), "reference": ref,
            "physics_recomputed": False,
            "plot_style": style,
            "comparison_1": "100*(O_W_FULL-O_W_FROZEN)/O_W_FROZEN",
            "comparison_2": "100*(O_LLSS_OFF-O_W_FROZEN)/O_W_FROZEN",
            "original_physics_status": provenance.get("status"),
        }, indent=2) + "\n", encoding="utf-8")
        (out / "summary.json").write_text(json.dumps({
            "models_in_input": len(rows),
            "models_valid": sum(r.get("status") == "VALID_EXPLORATORY" for r in rows),
            "plot_only": True,
        }, indent=2) + "\n", encoding="utf-8")
        print(f"Replotted common-denominator EFT1 heatmaps from {source} to {out}", flush=True)
        return 0

    from Numerical.parallel_t3.complete.all16_rough_yukawa_fit import (
        BENCHMARK, all_models, read_json, save_json)
    from scripts.BuildT3GeneralisedThesisTables import load_reference

    if args.selected_six and args.models:
        p.error("Use --selected-six OR --models")
    inventory = tuple(all_models())
    models = SELECTED_SIX if args.selected_six else tuple(args.models or inventory)
    if not models or len(set(models)) != len(models) or any(m not in inventory for m in models):
        p.error("Unknown or duplicate model in --models")
    tables = (ROOT / args.tables).resolve()
    refs = (ROOT / args.reference_results).resolve()
    out = (ROOT / args.output).resolve()
    owned_root = (ROOT / "output/parallel_t3").resolve()
    if out == owned_root or owned_root not in out.parents or out.exists():
        p.error("--output must be a NEW directory under output/parallel_t3")
    provenance_path = tables / "provenance.json"
    if not provenance_path.is_file() or not (tables / "results.json").is_file():
        p.error("Missing existing 32-case tables and provenance")
    provenance = read_json(provenance_path)
    if args.reference not in provenance.get("reference_benchmarks", {}):
        p.error("Reference donor not present in baseline provenance")
    donor = load_reference(args.reference, [refs])
    if str(donor["source"]) != provenance["reference_benchmarks"][args.reference]["result_json"]:
        p.error("Donor Yukawas do not match baseline provenance")
    cfg = read_json(ROOT / BENCHMARK)
    mf = float(cfg["scales_GeV"]["F"])
    ms = float(cfg["scales_GeV"]["S"])
    solver = provenance["wilson_solver"]
    steps = int(provenance["wilson_steps"])
    if solver not in ("RK2", "DOP853") or steps < 2:
        p.error("Unsupported solver/steps in source baseline")

    out.mkdir(parents=True)
    meta = {
        "status": "CONDITIONAL_NOT_PHYSICAL_MSbar", "reference": args.reference,
        "models_requested": list(models), "modes": list(CASES),
        "scope": "one-loop candidate 45-Weyl seven-slot tensor Wilson beta",
        "all_three_keep_renormalisable_EFT1_running": True,
        "all_three_keep_UV_F_matching_scalar_matching_SMEFT_running": True,
        "W_FROZEN": "freeze all seven EFT1 Wilson beta terms only",
        "LLSS_OFF": "subtract lambda_HHSS C_LLSS source from scalar_pair beta LLHH block at every RHS call",
        "W_FULL": "unmodified full Wilson beta",
        "heatmap_definition_vs_frozen": "100*(W_FULL - W_FROZEN)/W_FROZEN",
        "heatmap_definition_mixing_off_vs_frozen": "100*(LLSS_OFF - W_FROZEN)/W_FROZEN",
        "heatmap_columns": [
            {"observable": observable, "comparisons": ["running vs frozen", "mixing off vs frozen"]}
            for observable in OBSERVABLES
        ],
        "not_same_as_original_EFT1_off_heatmap": True,
        "no_posthoc_leading_log_added": True,
        "all_heavy_fermions_out_before_EFT1": True,
        "all_scalars_out_together_after_EFT1": True,
        "benchmark_scales_GeV": {"F": mf, "S": ms},
        "conditional_approximations": ["source-to-Weyl normalization uncertified",
            "full scalar tensor and candidate ADM not independently certified",
            "fermion hard one-loop C5 assumed zero", "reduced scalar potential"],
        "source_tables": str(tables), "reference_donor": str(donor["source"]),
        "wilson_solver": solver, "wilson_steps": steps,
        "reuse_cached_batch1": bool(args.reuse_cached_batch1),
        "plot_style": style,
        "allowed_changed_source_for_cached_batch1": (
            list(ALLOWED_CACHED_BATCH1_SOURCE_MISMATCH) if args.reuse_cached_batch1 else []
        ),
        "cached_beta_hash_validation_retained": True,
        "cached_input_interpretation": (
            "Saved Batch 1 UV/EFT1 beta functions are reused and are NOT regenerated from current Wolfram source"
            if args.reuse_cached_batch1 else "Require full Batch 1 source identity"
        ),
    }
    save_json(out / "provenance.json", meta)
    rows = []
    for model in models:
        row = {"reference": args.reference, "model": model, "status": "BLOCKED", "error": ""}
        try:
            outcomes, capture = _run_one_model(ROOT, model, args.reference, donor["x"],
                                                cfg, provenance, solver=solver, steps=steps,
                                                batch1=args.batch1,
                                                keep_trial_files=args.keep_trial_files,
                                                debug=args.debug,
                                                reuse_cached_batch1=args.reuse_cached_batch1)
            estimate = capture.get("beta_LLSS_to_C5_initial")
            summary, matrices = summarize_model(outcomes, estimate, mf, ms)
            row.update(summary)
            row["LLSS_beta_projection_residual_initial"] = capture.get(
                "source_projection_absolute_residual_initial")
            row["LLSS_source_tensor_norm_initial"] = capture.get("source_full_tensor_norm_initial")
            row["LLSS_ablation_rhs_calls"] = capture.get("rhs_calls")
            row["cached_batch1_export_hashes_verified"] = capture.get(
                "cached_batch1_export_hashes_verified", None)
            row["batch1_source_mismatches_json"] = json.dumps(
                capture.get("batch1_source_mismatches", {}), sort_keys=True)
            row["batch1_source_mismatch_accepted_for_exploratory_use"] = capture.get(
                "batch1_source_mismatch_accepted_for_exploratory_use", False)
            row["status"] = "VALID_EXPLORATORY"
            np.savez_compressed(out / f"{model}_C5_three_cases.npz", **matrices)
            save_json(out / f"{model}_observables.json", {
                "model": model, "reference": args.reference,
                "status": row["status"],
                "cases": {k: {"observables": v["observables"], "chi2": v["chi2"],
                               "meta": {key: value for key, value in v["meta"].items()
                                        if key != "scratch_token"}}
                          for k, v in outcomes.items()},
            })
        except Exception as exc:
            row["error"] = f"{type(exc).__name__}: {exc}"
            if args.debug:
                (out / f"{model}_error.txt").write_text(traceback.format_exc(), encoding="utf-8")
            print(f"BLOCKED {model}: {row['error']}", flush=True)
        rows.append(row)
        _write_csv(out / "all16_three_way.csv", rows)
        _write_csv(out / "selected6_three_way.csv", [r for r in rows if r["model"] in SELECTED_SIX])
        save_json(out / "summary.json", {
            "models_requested": len(models),
            "models_completed": sum(r["status"] == "VALID_EXPLORATORY" for r in rows),
            "models_blocked": [r["model"] for r in rows if r["status"] != "VALID_EXPLORATORY"],
            "reference": args.reference,
            "heatmap_interpretation": {
                "running vs frozen": "all EFT1 Wilson-running effects relative to zero Wilson beta",
                "mixing off vs frozen": "all other EFT1 Wilson-running effects with LLSS-to-C5 quartic source disabled, relative to the same frozen baseline",
                "column difference": "LLSS contribution to observable, normalised by the frozen observable",
            },
            "interpretation": "LLSS_fraction_norm close to 1 AND alignment_cosine close to 1 indicates LLSS source dominance, within candidate ADM",
        })
        print(f"DONE {model}: {row['status']}", flush=True)
    _write_long_heatmap_csv(out / "csv" / "all16_heatmap_long.csv", rows, inventory)
    _write_long_heatmap_csv(out / "csv" / "selected6_heatmap_long.csv", rows, SELECTED_SIX)
    _render_heatmaps(out, rows, inventory, "all16", style=style)
    _render_heatmaps(out, rows, SELECTED_SIX, "selected6", style=style)
    _write_donor_report(out, rows, inventory, args.reference, style)
    completed = sum(r["status"] == "VALID_EXPLORATORY" for r in rows)
    print(f"Three-way EFT1 Wilson isolation {completed}/{len(models)}; outputs: {out}", flush=True)
    return 0 if completed == len(models) else 1


if __name__ == "__main__":
    raise SystemExit(main())
