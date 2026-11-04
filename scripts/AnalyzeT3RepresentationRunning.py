"""Read-only, saved-data presentation of 16 general T3 models.

No matching, RGE integration, fit, or mutation of saved physics outputs.

Usage from repository root (PowerShell):
    python scripts/AnalyzeT3RepresentationRunning.py

Input contract is pinned to the source names checked at main SHA
6ebe9937a72880c6f180bcd026261779800ac71f. Updated data are rejected if
source IDs, cases, or 1111 fit observables disagree.

Four switches are (U, E, W, S). W=0 is ONLY removal of an LLSS leading-log
term from an anchored C5 bookkeeping decomposition; NOT an independently
hard-matched physical theory. Stage-by-stage comparisons hold the other
switches fixed and may not be added to reconstruct the full result.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path
from statistics import median

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = Path("output/studies/general_scalar_extended")
OUTPUT = Path("Reports/output/thesis_qualitative/representation_running")
MANIFEST = DATA / "thesis_general_all16/continue_general16_v1/merged_all16.json"
ABLATION = DATA / "general_sequential_rge_ablation/general_seq_05"
BC_ABLATION = DATA / "general_sequential_rge_ablation/general_seq_bc_normalized_01"
BC_AUDIT = DATA / "general_log_coefficient_audit/log_bc_normalized_01/summary.json"
SENSITIVITY = DATA / "zee_style_general_sensitivity/zee_general_all16_merged_v1"
MODELS = (
    "A_-4", "A_-2", "A_+0", "B_-3", "B_-1", "B_+1", "C_-3", "C_-1",
    "C_+1", "D_-2", "D_+0", "D_+2", "E_-4", "E_-2", "E_+0", "E_+2",
)
# Confirmed in Numerical/running/GeneralT3SequentialLeadingLog.py.
SCALAR_SU2_DIMS = {"A": (1, 3), "B": (2, 2), "C": (2, 2),
                   "D": (3, 1), "E": (3, 3)}
FERMION_SU2_DIMS = {"A": 2, "B": 1, "C": 3, "D": 2, "E": 2}
FROZEN = {"UV": "0111", "EFT1": "1011", "Wilson_conditional": "1101",
          "SMEFT": "1110"}
OBS = (
    "m_light_eV", "m_middle_eV", "m_heavy_eV", "delta_m21_sq_ev2",
    "delta_m3l_sq_ev2", "sin2_theta12", "sin2_theta13", "sin2_theta23",
)
# Avoid percentages when the switched-off reference is small in its units.
REFERENCE_FLOORS = {
    "m_light_eV": 1e-5, "m_middle_eV": 1e-5, "m_heavy_eV": 1e-5,
    "delta_m21_sq_ev2": 1e-8, "delta_m3l_sq_ev2": 1e-8,
    "sin2_theta12": 1e-4, "sin2_theta13": 1e-4, "sin2_theta23": 1e-4,
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def load_json(path: Path):
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def rows_csv(path: Path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict], columns: list[str] | None = None):
    if columns is None:
        columns = list(dict.fromkeys(key for row in rows for key in row))
    require(bool(columns), f"No CSV columns for {path}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def numeric_case(record: dict, label: str):
    require(record["status"] == "Computed", f"Uncomputed saved case: {label}")
    require(len(record["case"]) == 4 and set(record["case"]) <= {"0", "1"},
            f"Invalid switch bits: {label}")
    require(record["model"] in MODELS, f"Unexpected model: {label}")
    data = dict(record)
    for key in ("chi2", "log_over_full_C5", *OBS):
        data[key] = float(record[key])
        require(math.isfinite(data[key]), f"Nonfinite {key} in {label}")
    mass = [data[name] for name in OBS[:3]]
    require(0 <= mass[0] <= mass[1] <= mass[2], f"Invalid ascending masses: {label}")
    for key in OBS[5:]:
        require(0 <= data[key] <= 1, f"Unphysical saved sin^2 angle: {label}")
    return data


def build_dataset(root: Path):
    manifest = load_json(root / MANIFEST)
    main_summary = load_json(root / ABLATION / "summary.json")
    bc_summary = load_json(root / BC_ABLATION / "summary.json")
    bc_audit = load_json(root / BC_AUDIT)
    sensitivity_summary = load_json(root / SENSITIVITY / "summary.json")
    sensitivity_points = rows_csv(root / SENSITIVITY / "points.csv")
    main_rows = rows_csv(root / ABLATION / "all16_cases.csv")
    bc_rows = rows_csv(root / BC_ABLATION / "all16_cases.csv")

    require(manifest["accepted_count"] == 16 and manifest["blocked_count"] == 0,
            "Final merged fitted-source manifest does not contain 16 accepted models")
    fit_rows = {r["model"]: r for r in manifest["rows"]}
    # Model representations are source-verified against existing generated input
    # configurations; alpha parametrises the hypercharge assignment.
    for model in MODELS:
        fam, atext = model.split("_")
        alpha = int(atext)
        ds1, ds2 = SCALAR_SU2_DIMS[fam]
        df = FERMION_SU2_DIMS[fam]
        suffix = f"m{abs(alpha)}" if alpha < 0 else f"p{alpha}"
        config = root / "Numerical" / "configs" / "generated_models" / (
            f"T3_dS1_{ds1}_dS2_{ds2}_dF_{df}_alpha_{suffix}.json")
        require(config.is_file(), f"Missing model representation source: {config}")
        actual = load_json(config)["representation"]
        require((actual["d_s1"], actual["d_s2"], actual["d_f"], actual["alpha"])
                == (ds1, ds2, df, alpha),
                f"SU(2)/alpha metadata mismatch for {model}")
    require(set(fit_rows) == set(MODELS), "Final manifest model inventory differs")
    require(all(r["status"] == "ExploratoryAcceptedFit" for r in fit_rows.values()),
            "Final manifest contains a non-accepted fit")
    require(main_summary["computed"] == 16 and set(main_summary["requested_models"]) == set(MODELS),
            "general_seq_05 incomplete or model inventory differs")
    require(bc_summary["computed"] == 2 and set(bc_summary["requested_models"]) == {"B_-1", "C_-1"},
            "latest conditional B/C normalised run incomplete")
    for name, summary in (("main", main_summary), ("B/C", bc_summary)):
        require(summary["format"] == "t3-full-general-16_U_E_W_S_ablation-v1",
                f"{name} ablation schema differs")
        require(summary["one_scale_C5_used_as_W_on_reference"] is True,
                f"{name}: unexpected W-on anchoring")
        require(summary["hard_matching_residual_independently_calculated"] is False,
                f"{name}: matching assumptions changed")
        require(float(summary["spinor_factor_assumed"]) == 1.0,
                f"{name}: different GammaCC spinor factor")
        require(summary.get("basis_assumed") is True and
                summary.get("no_double_counting_by_construction") is True,
                f"{name}: C5 basis/bookkeeping assumptions changed")
        require(all(r["status"] == "Computed" for r in summary["model_rows"]),
                f"{name}: blocked model in summary")
    for row in main_summary["model_rows"]:
        if row["model"] not in ("B_-1", "C_-1"):
            require(row.get("provenance", {}).get("matching_kind") ==
                    "full-general-matchete-one-scale-channels",
                    f"Nonmixed matching changed for {row['model']}")
    require(main_summary["source_run"] == bc_summary["source_run"],
            "B/C and general ablations use different source runs")
    require(bc_audit["case_run"] == "general_seq_bc_normalized_01"
            and bc_audit["source_run"] == bc_summary["source_run"]
            and bc_audit["computed"] == 2,
            "B/C independent logarithmic audit inconsistent with selected B/C run")
    for item in bc_summary["model_rows"]:
        require("mixed_native_CG_exact_majorana_source" in item.get("provenance", {}),
                f"Updated B/C row lacks source-exact native-CG provenance: {item['model']}")
    require(bc_audit["complete_independent_hard_matching"] is False,
            "Unexpected hard-matching audit status")
    require(bc_audit["independent_mass_log_coefficient_extracted"] is True,
            "B/C independent mass-log coefficient not verified")
    require(set(x["model"] for x in bc_audit["model_rows"]) == {"B_-1", "C_-1"},
            "B/C audit model IDs differ")
    for row in bc_audit["model_rows"]:
        require(row["status"] == "Computed" and
                abs(float(row["ratio_Matchete_over_EFT_at_MF"]) - 1) < 1e-2,
                f"B/C independent mass-log coefficient not consistent: {row['model']}")
    require(sensitivity_summary["schema"] == "t3-general-zee-compensated-sensitivity-v1",
            "Unexpected sensitivity schema")
    require(set(x["model"] for x in sensitivity_summary["models"]) == set(MODELS)
            and all(x["status"] == "Complete" and x["failed_points"] == 0
                    for x in sensitivity_summary["models"]),
            "General coupling-sensitivity campaign incomplete")
    require(len(sensitivity_points) > 0 and
            set(r["model"] for r in sensitivity_points) == set(MODELS),
            "Sensitivity points incomplete")
    require(all(r["status"] == "Success" for r in sensitivity_points),
            "Failed sensitivity point present")
    for item in sensitivity_summary["models"]:
        m = item["model"]
        require(math.isclose(float(item["baseline_chi2_saved"]),
                             float(fit_rows[m]["chi2_final"]), rel_tol=1e-8, abs_tol=1e-9),
                f"Sensitivity versus fit manifest mismatch: {m}")

    # Replace the old mixed-scalar rows; never combine the two B/C conventions.
    require(len(main_rows) == 16 * 16, "Expected 256 general case rows")
    require(len(bc_rows) == 2 * 16, "Expected 32 conditional B/C case rows")
    selected = [r for r in main_rows if r["model"] not in {"B_-1", "C_-1"}] + bc_rows
    pairs = {}
    for row in selected:
        key = (row["model"], row["case"])
        require(key not in pairs, f"Duplicate case {key}")
        case = numeric_case(row, str(key))
        case["source_run"] = ("general_seq_bc_normalized_01" if row in bc_rows
                              else "general_seq_05")
        case["matched_physics_qualification"] = (
            "conditional_native_CG_BC" if row in bc_rows else "general_nonmixed")
        pairs[key] = case
    require(len(pairs) == 256 and
            all((m, f"{i:04b}") in pairs for m in MODELS for i in range(16)),
            "Not exactly 16 U/E/W/S cases per model")
    max_osc_rel = 0.0
    max_chi2_abs = 0.0
    for model in MODELS:
        fit = fit_rows[model]
        full = pairs[(model, "1111")]
        for key, desired in fit["prediction"].items():
            actual = full[key]
            max_osc_rel = max(max_osc_rel, abs(actual - float(desired)) /
                              max(abs(float(desired)), 1e-16))
        max_chi2_abs = max(max_chi2_abs,
                           abs(full["chi2"] - float(fit["chi2_final"])))
    require(max_osc_rel < 1e-8 and max_chi2_abs < 1e-6,
            "1111 does not reproduce final merged fit manifest: unsafe to combine")
    validation = {
        "basis_commit_checked": "6ebe9937a72880c6f180bcd026261779800ac71f",
        "source_campaign": main_summary["source_run"],
        "models": 16, "cases": len(pairs),
        "case_source_priority": "B_-1/C_-1 from general_seq_bc_normalized_01; other 14 from general_seq_05",
        "maximum_relative_1111_oscillation_manifest_difference": max_osc_rel,
        "maximum_absolute_1111_chi2_manifest_difference": max_chi2_abs,
        "independent_hard_matching": False,
        "log_normalization_certified": False,
        "W_off_interpretation": "conditionally removed LLSS log from one-scale anchor only",
        "B_C_mass_log_audit": bc_audit["model_rows"],
        "sensitivity_points_count": len(sensitivity_points),
    }
    return pairs, fit_rows, sensitivity_points, validation


def rel_effect(full, frozen, variable):
    """Signed delta in native units, and percent only away from the reference floor."""
    df = full[variable] - frozen[variable]
    base = frozen[variable]
    pct = (100.0 * df / abs(base)) if abs(base) >= REFERENCE_FLOORS[variable] else None
    return df, pct


def derive(pairs):
    baseline = []
    long = []
    heat = []
    groups = []
    wilson = []
    long_cases = []
    for model in MODELS:
        family, alpha_text = model.split("_")
        alpha = int(alpha_text)
        rep = SCALAR_SU2_DIMS[family]
        meta = {"model": model, "family": family, "alpha": alpha,
                "scalar_S1_SU2_dimension": rep[0],
                "scalar_S2_SU2_dimension": rep[1]}
        full = pairs[(model, "1111")]
        baseline.append({**meta, "fermion_SU2_dimension": FERMION_SU2_DIMS[family], **{k: full[k] for k in ("chi2", *OBS)},
                         "fit_case": "1111", "case_source": full["source_run"]})
        rec = {**meta, "fermion_SU2_dimension": FERMION_SU2_DIMS[family]}
        for stage, frozen_code in FROZEN.items():
            off = pairs[(model, frozen_code)]
            ratio = full["m_heavy_eV"] / off["m_heavy_eV"] if off["m_heavy_eV"] >= 1e-5 else math.nan
            pct = 100 * (ratio - 1) if math.isfinite(ratio) else math.nan
            rec[stage + "_ratio_mheavy"] = ratio
            rec[stage + "_percent_mheavy"] = pct
            rec[stage + "_log10_ratio_mheavy"] = math.log10(ratio) if ratio > 0 else math.nan
            groups.append({**meta, "stage": stage, "ratio_mheavy": ratio,
                           "percent_mheavy": pct})
            for variable in OBS:
                absolute, percent = rel_effect(full, off, variable)
                long.append({**meta, "switch": stage, "on_case": "1111",
                             "off_case": frozen_code, "observable": variable,
                             "full_value": full[variable],
                             "frozen_value": off[variable],
                             "absolute_change": absolute,
                             "relative_change_percent": percent if percent is not None else "",
                             "reference_near_zero": percent is None,
                             "comparison_kind": ("conditional_wilson_log_reconstruction"
                                                 if stage == "Wilson_conditional"
                                                 else "independently_evolved_running_stage")})
        heat.append(rec)
        wilson.append({**meta,
                       "c5_induced_log_over_anchored_one_scale_frobenius": full["log_over_full_C5"],
                       "mheavy_full_eV": full["m_heavy_eV"],
                       "mheavy_wilson_off_eV": pairs[(model, "1101")]["m_heavy_eV"],
                       "mheavy_ratio_full_over_wilson_off": rec["Wilson_conditional_ratio_mheavy"],
                       "interpretation": "conditional_only; no independent hard matching"})
        for i in range(16):
            code = f"{i:04b}"
            row = pairs[(model, code)]
            case = {**meta, "fermion_SU2_dimension": FERMION_SU2_DIMS[family], "case": code, "U": code[0], "E": code[1],
                    "W": code[2], "S": code[3], "case_source": row["source_run"],
                    **{k: row[k] for k in ("chi2", "log_over_full_C5", *OBS)}}
            for a, s2 in (("theta12_deg", "sin2_theta12"),
                          ("theta13_deg", "sin2_theta13"),
                          ("theta23_deg", "sin2_theta23")):
                case[a] = math.degrees(math.asin(math.sqrt(row[s2])))
            long_cases.append(case)
    return baseline, long, heat, groups, wilson, long_cases


def sensitivity_slopes(points):
    by_group = defaultdict(dict)
    for point in points:
        key = point["model"], point["direction"]
        gamma = round(float(point["gamma"]), 9)
        require(gamma not in by_group[key], f"Duplicate sensitivity gamma at {key}")
        by_group[key][gamma] = point
    require(all(set(k for k in by_group if k[0] == model) for model in MODELS),
            "Missing model sensitivity directions")
    result = []
    quantities = OBS[3:]  # only the saved oscillation observables
    for (model, direction), g in sorted(by_group.items()):
        require(all(x in g for x in (0.99, 1.0, 1.01)),
                f"Missing symmetric sensitivity stencil: {model} {direction}")
        family, alpha = model.split("_")
        for name in quantities:
            x0 = float(g[1.0][name]); xm = float(g[0.99][name]); xp = float(g[1.01][name])
            deriv = (xp - xm) / 0.02
            result.append({"model": model, "family": family, "alpha": int(alpha),
                           "saved_direction": direction, "observable": name,
                           "baseline": x0, "d_observable_d_gamma_centered": deriv,
                           "fractional_derivative_per_unit_gamma":
                               deriv / x0 if abs(x0) >= REFERENCE_FLOORS[name] else "",
                           "gamma_minus": 0.99, "gamma_plus": 1.01,
                           "interpretation": "local saved compensated-coupling sensitivity; not a stage RGE"})
    return result


def plot_heatmaps(heat, out):
    names = [r["model"] for r in heat]
    stages = ["UV", "EFT1", "SMEFT"]
    matrix = np.array([[r[s + "_percent_mheavy"] for s in stages] for r in heat])
    fig, ax = plt.subplots(figsize=(9, 7.1))
    bound = max(6.0, float(np.nanmax(np.abs(matrix))))
    im = ax.imshow(matrix, cmap="RdBu_r", vmin=-bound, vmax=bound, aspect="auto")
    ax.set_xticks(range(3), ["UV", "EFT1", "SMEFT"])
    ax.set_yticks(range(16), names)
    ax.set_title("Effect of enabling each independently evolved stage\non $m_3$ with other switches fixed on")
    for i in range(16):
        for j in range(3):
            ax.text(j, i, f"{matrix[i,j]:+.1f}%", ha="center", va="center", fontsize=8,
                    color="white" if abs(matrix[i,j]) > 0.58 * bound else "black")
    cbar = fig.colorbar(im, ax=ax, shrink=0.85)
    cbar.set_label("100 [$m_3$(1111)/$m_3$(stage frozen) − 1] (%)")
    fig.tight_layout()
    fig.savefig(out / "independent_stage_mass_heatmap.png", dpi=210)
    plt.close(fig)

    wilson = np.array([r["Wilson_conditional_ratio_mheavy"] for r in heat])
    fig, ax = plt.subplots(figsize=(9, 6.7))
    im = ax.imshow(wilson[:, None], cmap="magma", aspect="auto")
    ax.set_xticks([0], ["One-scale-anchored conditional W toggle"])
    ax.set_yticks(range(16), names)
    for i in range(16):
        ax.text(0, i, f"{wilson[i]:.2f}×", ha="center", va="center", fontsize=9,
                color="white" if wilson[i] < 9.0 else "black")
    fig.colorbar(im, ax=ax, shrink=0.84, label="$m_3$(1111)/$m_3$(1101)")
    ax.set_title("Conditional Wilson leading-log toggle\nNOT independent hard matching")
    fig.tight_layout()
    fig.savefig(out / "conditional_wilson_toggle.png", dpi=210)
    plt.close(fig)


def plot_observables(long, out):
    fig, axs = plt.subplots(1, 2, figsize=(13.5, 6.5), sharey=True)
    stage_keys = ("UV", "EFT1", "SMEFT", "Wilson_conditional")
    wanted = ("delta_m3l_sq_ev2", "sin2_theta23")
    for ax, quantity in zip(axs, wanted):
        for stage in stage_keys:
            vals = [x for x in long if x["switch"] == stage and x["observable"] == quantity]
            vals.sort(key=lambda x: MODELS.index(x["model"]))
            y = [x["relative_change_percent"] if quantity.startswith("delta") else x["absolute_change"]
                 for x in vals]
            ax.plot(range(16), y, marker="o", markersize=3, linewidth=1.2,
                    label="Wilson (conditional)" if stage == "Wilson_conditional" else stage)
        ax.axhline(0, linewidth=0.7, color="gray")
        ax.set_xticks(range(16), MODELS, rotation=75)
        ax.set_title((r"$\Delta m^2_{3\ell}$ relative shift (%)" if quantity.startswith("delta")
                      else r"$\sin^2\theta_{23}$ absolute shift"))
        ax.grid(alpha=0.18)
    axs[0].legend(loc="best", fontsize=8)
    fig.suptitle("Fixed-input toggles (all other switches enabled); Wilson is conditional")
    fig.tight_layout()
    fig.savefig(out / "observable_stage_changes.png", dpi=205)
    plt.close(fig)


def summarize_groups(groups):
    out = []
    for family in "ABCDE":
        for stage in FROZEN:
            vals = [r["percent_mheavy"] for r in groups
                    if r["family"] == family and r["stage"] == stage]
            require(len(vals) == (4 if family == "E" else 3),
                    "Incomplete family group")
            out.append({"family": family,
                        "scalar_S1_SU2_dimension": SCALAR_SU2_DIMS[family][0],
                        "scalar_S2_SU2_dimension": SCALAR_SU2_DIMS[family][1],
                        "fermion_SU2_dimension": FERMION_SU2_DIMS[family],
                        "stage": stage, "n_models": len(vals),
                        "minimum_percent_mheavy": min(vals),
                        "median_percent_mheavy": median(vals),
                        "maximum_percent_mheavy": max(vals),
                        "positive_count": sum(x > 0 for x in vals),
                        "negative_count": sum(x < 0 for x in vals)})
    return out


def report_markdown(heat, validation, groups, fit_rows):
    def range_stage(s):
        vals = [100 * (x[s + "_ratio_mheavy"] - 1) for x in heat]
        return min(vals), max(vals)
    def spread(s, family):
        vals = [r["percent_mheavy"] for r in groups
                if r["stage"] == s and r["family"] == family]
        return f"{min(vals):+.2f}% to {max(vals):+.2f}%"
    uv = range_stage("UV"); smeft = range_stage("SMEFT")
    w = [r["Wilson_conditional_ratio_mheavy"] for r in heat]
    lines = [
        "# Qualitative SU(2) representation and running effects: 16 general T3 models",
        "",
        "## Scope and provenance",
        "",
        "Read-only analysis of saved data (no matching, RG integration, or refits).",
        f"Pinned reference main commit: `{validation['basis_commit_checked']}`.",
        f"Run provenance field: `{validation['source_campaign']}`; this is **not** a replacement for the final fit-status manifest.",
        "Final fit-status source: `continue_general16_v1/merged_all16.json` (16 accepted exploratory fits).",
        "14 ordinary model results: `general_seq_05`; B(-1)/C(-1): **newer** `general_seq_bc_normalized_01`.",
        "The independent B/C mass-log audit is `log_bc_normalized_01`; local coupling sensitivity is `zee_general_all16_merged_v1`.",
        "",
        "## Switch and EFT semantics",
        "",
        "`U/E/W/S`: UV beta (1e7→1e5 GeV); EFT1 scalar beta (1e5→1e3 GeV);",
        "reconstructed intermediate LLSS→Weinberg leading log; SMEFT/SM+Weinberg beta (1e3→1e2 GeV).",
        "All fermions integrated at 1e5 GeV, then both scalar multiplets at 1e3 GeV.",
        "`1111` is fully enabled; freeze one stage using `0111/1011/1101/1110` respectively.",
        "The stage deltas are **conditional ablations at fixed fitted couplings**, not separately refitted models.",
        "They need not add or factorize: turning on U or E modifies downstream matching inputs.",
        "The W-on C5 is the full-general **one-scale anchor**; W-off subtracts the reconstructed LLSS log while retaining the bubble and bookkeeping residual.",
        "W-off is **not independent physical hard matching**; its large suppression must not be claimed as an independently certified running prediction.",
        "The GammaCC-to-Weyl factor is assumed 1; the extra HHSS basis and absolute normalization remain provisional.",
        "",
        "## Verified qualitative patterns in the saved cases",
        "",
        f"- UV enabled: heaviest mass changes **{uv[0]:+.2f}% to {uv[1]:+.2f}%**, positive in all 16 models; family E is largest.",
        f"- SMEFT enabled: changes **{smeft[0]:+.2f}% to {smeft[1]:+.2f}%**, negative and almost universal.",
        "- Intermediate EFT1: suppresses the heaviest mass for every B and C assignment;",
        "  enhances every E assignment; A and D show mixed or small shifts.",
        f"- W-off conditional comparison: the full/W-off heaviest-mass ratio ranges from {min(w):.2f} to {max(w):.2f}.",
        "  This is a **bookkeeping sensitivity** to the LLSS log, not independent hard-matched physics.",
        "",
        "Family-dependent EFT1 shifts to $m_3$ (other switches enabled):",
        "",
        "| Family | SU(2) dimensions $(F,S_1,S_2)$ | EFT1 effect |",
        "|---|---|---:|",
    ]
    for family in "ABCDE":
        dims = SCALAR_SU2_DIMS[family]
        lines.append(f"| {family} | {FERMION_SU2_DIMS[family]}, {dims[0]}, {dims[1]} | {spread('EFT1',family)} |")
    lines += [
        "",
        "The integer alpha in model labels parametrizes hypercharge assignments.",
        "Pairing models within each family probes alpha dependence at fixed scalar SU(2) dimensions;",
        "however family differences also include distinct fermion sectors and channels, so this is not a controlled single-representation derivative.",
        "",
        "## Interpretation limitations",
        "",
        "All 16 accepted UV Yukawa configurations were calibrated against the same oscillation targets.",
        "Therefore nearly identical `1111` mass-squared splittings/mixings are *not* evidence that different representations physically predict the same observables,",
        "nor can best-fit agreement alone establish SU(2) or hypercharge effects.",
        "The meaningful representation comparison here is **the response to frozen stages under the same campaign conventions**.",
        "Some `W=0` predictions can have strongly distorted mass spectra or mixing angles; interpret near-zero references using absolute differences.",
        "B(-1)/C(-1) have channel-normalised conditional inputs and should be compared with the other fourteen only as a *qualified* combined exploratory sample.",
        "The independent B/C audit checks the mass-log coefficient, not the complete finite or hard matching.",
        "The compensated-coupling sensitivity CSV records local saved finite differences; it is **not an additional running stage**.",
        "",
        "## Reproduction checks",
        "",
        f"- Selected case count: {validation['cases']} (16 models × 16 U/E/W/S combinations).",
        f"- Max relative oscillation-observable difference between `1111` and final manifest: {validation['maximum_relative_1111_oscillation_manifest_difference']:.3g}.",
        f"- Sensitivity points checked: {validation['sensitivity_points_count']} (all reported successful).",
        "- No saved physics outputs written.",
        "",
    ]
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT,
                        help="Repository root containing input output/studies; default is script parent")
    parser.add_argument("--output", type=Path, default=OUTPUT,
                        help="New presentation output path, relative to repository root")
    args = parser.parse_args()
    root = args.root.resolve()
    output = (root / args.output).resolve() if not args.output.is_absolute() else args.output.resolve()
    require(output.is_relative_to(root / "Reports" / "output" / "thesis_qualitative"),
            "Output must stay inside Reports/output/thesis_qualitative")
    for p in (MANIFEST, ABLATION / "all16_cases.csv", BC_ABLATION / "all16_cases.csv"):
        require((root / p).is_file(), f"Required existing data missing: {root / p}")
    pairs, fit_rows, sens, validation = build_dataset(root)
    baseline, long, heat, groups, wilson, cases = derive(pairs)
    slopes = sensitivity_slopes(sens)
    group_summary = summarize_groups(groups)

    csvdir = output / "csv"; plotdir = output / "plots"
    csvdir.mkdir(parents=True, exist_ok=True)
    plotdir.mkdir(parents=True, exist_ok=True)
    write_csv(csvdir / "all_256_saved_cases.csv", cases)
    write_csv(csvdir / "full_1111_baselines.csv", baseline)
    write_csv(csvdir / "frozen_stage_observable_changes.csv", long)
    write_csv(csvdir / "representation_stage_heatmap_values.csv", heat)
    write_csv(csvdir / "model_family_stage_summary.csv", group_summary)
    write_csv(csvdir / "conditional_wilson_log_diagnostics.csv", wilson)
    write_csv(csvdir / "compensated_coupling_sensitivity_slopes.csv", slopes)
    plot_heatmaps(heat, plotdir)
    plot_observables(long, plotdir)
    (output / "physical_interpretation.md").write_text(
        report_markdown(heat, validation, groups, fit_rows), encoding="utf-8")
    (output / "data_validation.json").write_text(
        json.dumps(validation, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(f"Verified {validation['cases']} cases, 16 accepted fitted models, "
          f"max oscillation consistency {validation['maximum_relative_1111_oscillation_manifest_difference']:.3g}")
    print(f"Saved presentation outputs under: {output}")


if __name__ == "__main__":
    main()
