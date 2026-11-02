"""Batch 3: diagnostic-only UV-versus-conditional-EFT neutrino matrix comparison.

No EFT matching is run and no physical normalisation is inferred.  The UV
reference and EFT response use DIFFERENT coupling/scale choices; the best
complex scale is a linear-algebra diagnostic, not a measured GammaCC factor.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import re
from typing import Any

import numpy as np

MODELS = ("A_-4", "B_-1")
VEV_GEV = 246.22
EV_PER_GEV = 1.0e9
FORMAT = "t3-batch3-uv-vs-conditional-eft-diagnostic-v1"
SNAPSHOT_REL = "Numerical/lo_certification/benchmark/batch3_github_reference.json"
LOCAL_EFT_REL = "output/general_observable/lo_joint_scalar_matching/conditional_all16_batch2_v1/summary.json"
LOCAL_UV_DIR = "output/general_observable/lo_certification_direct_uv/native_uv_reference_01/models"
BOUNDARY_REL = "output/general_observable/lo_fermion_boundary/releases/source_021dcc1e8fe8ba0f"
BOUNDARY_SLUG = {"A_-4": "A_m4", "B_-1": "B_m1"}
OUTPUT_REL = "output/general_observable/lo_certification_benchmark"


class ComparisonError(ValueError):
    pass


def require(condition: bool, reason: str) -> None:
    if not condition:
        raise ComparisonError(reason)


def read_json(path: Path) -> dict:
    require(path.is_file(), f"Missing source file: {path}")
    data = json.loads(path.read_text(encoding="utf-8-sig"))
    require(isinstance(data, dict), f"Expected JSON object: {path}")
    return data


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def matrix(obj: Any, name: str) -> np.ndarray:
    require(isinstance(obj, dict) and set(obj) == {"real", "imag"},
            f"{name}: expected real/imag matrices")
    try:
        r, im = np.asarray(obj["real"], dtype=float), np.asarray(obj["imag"], dtype=float)
    except (TypeError, ValueError) as exc:
        raise ComparisonError(f"{name}: invalid numeric entry: {exc}") from exc
    require(r.shape == im.shape == (3, 3), f"{name}: expected 3 x 3")
    x = r + 1j * im
    require(bool(np.isfinite(x).all()), f"{name}: non-finite matrix")
    norm = float(np.linalg.norm(x))
    require(norm > 0, f"{name}: zero matrix cannot be compared")
    require(float(np.linalg.norm(x - x.T)) / norm < 1e-10,
            f"{name}: not complex symmetric")
    return x


def verify_singular_values(x: np.ndarray, values: Any, name: str) -> list[float]:
    try:
        recorded = np.asarray(values, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ComparisonError(f"{name}: invalid singular values: {exc}") from exc
    require(recorded.shape == (3,) and bool(np.isfinite(recorded).all()) and
            bool(np.all(recorded >= 0)), f"{name}: expected three nonnegative finite masses")
    observed = np.sort(np.linalg.svd(x, compute_uv=False))
    error = float(np.max(abs(np.sort(recorded) - observed)))
    scale = float(max(np.max(observed), 1e-12))
    require(error / scale < 2e-8, f"{name}: recorded singular masses disagree with matrix")
    return sorted(observed.tolist(), reverse=True)


def _record(uv: dict, eft: dict, *, origin: str) -> dict:
    """Normalize an existing GitHub source record or full local report."""
    uv_full = "complex_symmetric_matrix_gev" in uv
    eft_full = "conditional_mnu_per_N_eV" in eft
    if uv_full:
        require(uv.get("status") == "independent_exact_spectral_reference_from_conditional_input",
                "UV source is not the expected independent reference")
        require(uv.get("physical_source_normalisation_certified") is False,
                "UV source unexpectedly asserts physical normalisation")
        uv_x, uv_masses = uv["complex_symmetric_matrix_gev"], uv["takagi_masses_ev_from_independent_singular_values"]
        uv_scale = uv["input_scale_gev"]
    else:
        require(uv.get("status") == "independent_exact_spectral_reference_from_conditional_input" and
                uv.get("physical_source_normalisation_certified") is False,
                "Snapshot UV status invalid")
        uv_x, uv_masses = uv["matrix_gev"], uv["singular_masses_ev"]
        uv_scale = uv["input_scale_gev"]
    if eft_full:
        require(eft.get("status") == "conditional_C5_per_unknown_NGamma_exploratory" and
                eft.get("physical_C5_certified") is False and
                eft.get("absolute_GammaCC_to_Weyl_factor") is None and
                eft.get("complete_one_loop_prediction") is False,
                "EFT report is not explicitly conditional")
        eft_x, eft_masses = eft["conditional_mnu_per_N_eV"], eft["mass_singular_values_eV_per_abs_N"]
        mu_f, mu_s = eft["mu_f_gev"], eft["mu_s_gev"]
    else:
        require(eft.get("status") == "conditional_C5_per_unknown_NGamma_exploratory" and
                eft.get("physical_C5_certified") is False and
                eft.get("absolute_GammaCC_to_Weyl_factor") is None and
                eft.get("complete_one_loop_prediction") is False,
                "Snapshot EFT status invalid")
        eft_x, eft_masses = eft["matrix_ev_per_N"], eft["singular_masses_ev_per_abs_N"]
        mu_f, mu_s = eft["mu_f_gev"], eft["mu_s_gev"]
    require(all(isinstance(q, (int, float)) and math.isfinite(q) and q > 0
                for q in (uv_scale, mu_f, mu_s)) and uv_scale > mu_f > mu_s,
            "Invalid scales: expected UV input > all-F threshold > simultaneous scalar threshold")
    require(uv_scale == 1e7 and mu_f == 1e5 and mu_s == 1e3,
            "Different thresholds than pinned A_-4/B_-1 source comparison")
    return {"origin": origin, "uv_matrix_gev": uv_x, "uv_masses_ev": uv_masses,
            "eft_matrix_ev_per_N": eft_x, "eft_masses_ev_per_abs_N": eft_masses,
            "uv_input_scale_gev": uv_scale, "mu_f_gev": mu_f, "mu_s_gev": mu_s}


def flavour_support_from_boundary(d: dict, model: str) -> dict:
    require(d.get("model") == model and
            d.get("status") == "NumericalLLSSBySourceVerified_SpinorNormalizationOpen" and
            d.get("spinor_normalization", {}).get("certified") is False,
            f"{model}: boundary source format/status changed")
    rows = d.get("summed_independent_symmetric_coordinates")
    require(isinstance(rows, list) and len(rows) > 0, f"{model}: empty LLSS boundary")
    keys = [x["coordinate"] for x in rows]
    require(all(len(k) == 6 and k[1] in (1, 2, 3) and
                k[3] in (1, 2, 3) for k in keys),
            f"{model}: unexpected LLSS flavour coordinate")
    off = sum(k[1] != k[3] for k in keys)
    return {"total_independent_LLSS_coordinates": len(rows),
            "lepton_flavour_offdiagonal_coordinates": off,
            "lepton_flavour_diagonal_coordinates": len(rows)-off}


def load_snapshot(root: Path) -> tuple[dict[str, dict], dict]:
    path = root / SNAPSHOT_REL
    data = read_json(path)
    require(data.get("format") == "t3-batch3-source-snapshot-v1" and
            data.get("physical_results_certified") is False and
            set(data.get("models", {})) == set(MODELS),
            "Wrong GitHub reference snapshot or model inventory")
    result = {}
    flavour = data.get("fermion_boundary_flavour_support", {})
    require(set(flavour) == set(MODELS), "Snapshot LLSS flavour support inventory missing")
    for m in MODELS:
        row = data["models"][m]
        result[m] = _record(row["uv"], row["eft"], origin="github_source_snapshot")
        f = flavour[m]
        require(f["total_independent_LLSS_coordinates"] ==
                f["lepton_flavour_diagonal_coordinates"] +
                f["lepton_flavour_offdiagonal_coordinates"],
                f"{m}: inconsistent snapshot flavour support")
        result[m]["fermion_boundary_flavour_support"] = f
    return result, {"source_mode": "github-reference", "snapshot_file": SNAPSHOT_REL,
                    "snapshot_sha256": file_sha256(path),
                    "github_commit_inspected": data["github_commit_inspected"],
                    "source_git_blobs": data["source_git_blobs"]}


def load_local(root: Path) -> tuple[dict[str, dict], dict]:
    eft_path = root / LOCAL_EFT_REL
    summary = read_json(eft_path)
    require(summary.get("format") == "t3-conditional-all16-consolidated-batch2-v1" and
            summary.get("status") == "All16ConditionalC5Exploratory" and
            summary.get("physical_C5_certified") is False and
            summary.get("absolute_GammaCC_to_Weyl_factor") is None and
            summary.get("complete_one_loop_neutrino_mass_prediction") is False,
            "Batch 2 summary is not the expected conditional source")
    records = summary.get("models")
    require(isinstance(records, list), "Batch 2 model records missing")
    index = {x["model"]: x for x in records}
    require(len(index) == len(records) and all(m in index for m in MODELS),
            "Batch 2 is missing A_-4 or B_-1")
    out, hashes = {}, {LOCAL_EFT_REL: file_sha256(eft_path)}
    for model in MODELS:
        uv_rel = f"{LOCAL_UV_DIR}/{model}.json"
        uv_path = root / uv_rel
        uv = read_json(uv_path)
        require(uv.get("model") == model, f"{model}: UV source model mismatch")
        out[model] = _record(uv, index[model], origin="local_source_reports")
        bound_rel = f"{BOUNDARY_REL}/{BOUNDARY_SLUG[model]}/boundary.json"
        bound_path = root / bound_rel
        out[model]["fermion_boundary_flavour_support"] = flavour_support_from_boundary(
            read_json(bound_path), model)
        hashes[uv_rel] = file_sha256(uv_path)
        hashes[bound_rel] = file_sha256(bound_path)
    return out, {"source_mode": "local", "local_source_sha256": hashes,
                 "source_sha256_are_byte_level_not_semantic": True}


def compare_one(model: str, record: dict) -> dict:
    require(model in MODELS, f"Unsupported model: {model}")
    u = matrix(record["uv_matrix_gev"], f"{model} UV") * EV_PER_GEV
    e = matrix(record["eft_matrix_ev_per_N"], f"{model} EFT")
    uv_sing = verify_singular_values(u, record["uv_masses_ev"], f"{model} UV")
    eft_sing = verify_singular_values(e, record["eft_masses_ev_per_abs_N"], f"{model} EFT")
    # Diagnostic minimisation over a SINGLE complex overall scale.  In the
    # chosen COMMON REPORTED flavour axes, this cannot repair off-diagonal
    # disagreement. It is not a physical fitted GammaCC coefficient.
    z = np.vdot(e, u) / np.vdot(e, e)
    residual = float(np.linalg.norm(u - z * e) / np.linalg.norm(u))
    uv_offdiag = u - np.diag(np.diag(u))
    eft_offdiag = e - np.diag(np.diag(e))
    frac_uv = float(np.linalg.norm(uv_offdiag) / np.linalg.norm(u))
    frac_eft = float(np.linalg.norm(eft_offdiag) / np.linalg.norm(e))
    flavour = record["fermion_boundary_flavour_support"]
    require(flavour["total_independent_LLSS_coordinates"] ==
            flavour["lepton_flavour_diagonal_coordinates"] +
            flavour["lepton_flavour_offdiagonal_coordinates"],
            "Fermion boundary flavour-support counts inconsistent")
    require(-1e-12 <= residual <= 1+1e-12 and 0 <= frac_uv <= 1+1e-12,
            "Invalid relative matrix geometry")
    return {
        "model": model,
        "source_origin": record["origin"],
        "UV_input_scale_gev": record["uv_input_scale_gev"],
        "EFT_mu_F_gev": record["mu_f_gev"],
        "EFT_mu_S_gev": record["mu_s_gev"],
        "uv_neutrino_mass_singular_values_ev_desc": uv_sing,
        "eft_mass_singular_values_ev_per_abs_N_desc": eft_sing,
        "uv_offdiagonal_frobenius_fraction": frac_uv,
        "eft_offdiagonal_frobenius_fraction": frac_eft,
        "fermion_boundary_flavour_support": flavour,
        "fermion_boundary_contains_only_same_flavour_lepton_pairs":
            flavour["lepton_flavour_offdiagonal_coordinates"] == 0,
        "source_flavour_mismatch_is_a_proved_physics_bug": False,
        "best_complex_scale_reported_axes": {"real": float(z.real), "imag": float(z.imag)},
        "residual_after_best_complex_scale": residual,
        "best_complex_scale_is_physical_GammaCC_normalisation": False,
        "flavour_basis_and_scale_compatibility_established": False,
        "physical_equality_expected": False,
        "UV_matrix_units": "GeV converted to eV",
        "EFT_matrix_units": "eV per unknown signed N_GammaCC",
    }


def compare_all(root: Path, source_mode: str = "github-reference") -> dict:
    require(source_mode in ("github-reference", "local"), "Unknown source mode")
    root = root.resolve()
    recs, source = (load_snapshot(root) if source_mode == "github-reference"
                    else load_local(root))
    models = [compare_one(m, recs[m]) for m in MODELS]
    return {
        "format": FORMAT,
        "status": "COMPARISON_ONLY_NOT_PHYSICAL_MATCHING",
        "model_count": len(models),
        "sources": source,
        "models": models,
        "comparison_scope": "fixed-reported flavour axes; best single complex scale per model",
        "matching_sequence": "UV -> integrate out ALL heavy BSM fermions at muF -> EFT1 retaining BOTH BSM scalars -> integrate BOTH scalar multiplets at muS -> SMEFT",
        "known_incompatibilities": [
            "Direct UV uses frozen high-scale inputs; conditional EFT uses illustrative running endpoint",
            "Native Matchete GammaCC/Weyl prefactor and cross-package CG/phases not certified",
            "Same flavour axes and running Yukawa/mass boundary not independently certified",
            "Hard one-loop heavy-fermion terms and full EFT1 operator matching missing",
        ],
        "physical_GammaCC_normalization": None,
        "physical_C5_gev_inverse": None,
        "physical_mnu_ev": None,
        "physical_certification": False,
    }


def markdown_report(data: dict) -> str:
    lines = ["# Batch 3 — UV vs conditional EFT matrix comparison", "",
             "**DIAGNOSTIC ONLY — NOT PHYSICAL MATCHING.**", "",
             f"Source mode: `{data['sources']['source_mode']}`. The two calculations use different source scales/coupling states and unmatched operator content.",
             "", "| Model | UV masses (eV) | EFT masses (eV / abs N) | UV off-diagonal fraction | Best-scale residual |",
             "|---|---|---|---:|---:|"]
    for r in data["models"]:
        u = ", ".join(f"{v:.6g}" for v in r["uv_neutrino_mass_singular_values_ev_desc"])
        e = ", ".join(f"{v:.6g}" for v in r["eft_mass_singular_values_ev_per_abs_N_desc"])
        lines.append(f"| `{r['model']}` | {u} | {e} | {r['uv_offdiagonal_frobenius_fraction']:.6f} | {r['residual_after_best_complex_scale']:.6f} |")
    lines += ["", "The residual is `||M_UV-z M_EFT||_F / ||M_UV||_F`, where `z=vdot(M_EFT,M_UV)/vdot(M_EFT,M_EFT)` and both matrices are expressed in eV. It is evaluated in the **reported** lepton-flavour axes; there is no basis optimisation.",
              "", "The inspected source-resolved F-tree LLSS tensors contain only lepton flavour pairs `p=q` for both models, while the reported direct-UV matrices have nonzero off-diagonal flavour entries. This is an **unresolved source/basis/input alignment diagnostic**; it is not proof the tensors should agree before common-scale matching.",
             "", "An overall spinor phase or factor alone cannot align the matrices in these axes. This is NOT evidence of a physical matching failure until source couplings, scales, flavour basis, CG normalisations, and missing hard/soft matching are reconciled.",
              "", "Both heavy scalar multiplets are integrated together at `muS`. Both direct-UV and EFT masses are conditional; no physical `GammaCC`, `C5` or `mnu` is assigned.", ""]
    return "\n".join(lines)


def emit(root: Path, data: dict, run_id: str) -> Path:
    require(re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]{0,63}", run_id) is not None,
            "Unsafe run ID")
    dest = root.resolve() / OUTPUT_REL / f"batch3_uv_eft_{run_id}"
    require(not dest.exists(), f"Output already exists; choose a new --run-id: {dest}")
    dest.mkdir(parents=True)
    (dest / "summary.json").write_text(json.dumps(data, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    (dest / "comparison.md").write_text(markdown_report(data), encoding="utf-8")
    return dest
