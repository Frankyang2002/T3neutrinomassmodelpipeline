"""A_-4 dense-fitted common-input exploratory two-threshold calculation.

Reuses the source-pinned frozen benchmark, canonical RGBeta UV/EFT1 runner,
Matchete source Yukawa factorization + CG projector, established LLSS orbit
adapter and joint scalar kernel. No signed GammaCC or physical C5 is inferred.

Owns only output/general_observable/lo_certification_benchmark/batch4_*.
"""
from __future__ import annotations

from hashlib import sha256
from pathlib import Path
from types import SimpleNamespace
import json
import math
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
MODEL = "A_-4"
FORMAT = "t3-a-m4-fitted-common-input-conditional-two-threshold-v1"
BENCH = Path("Numerical/benchmarks/fitted_reference/T3_dS1_1_dS2_3_dF_2_alpha_m4.json")
BENCH_SHA = "3cb8dc3e94f2bd6f60571a7a85477638e8fc02109bccfcec0e04055852d81247"
FERMION_SEED = Path("output/general_observable/fermion_first_all16_v1/A_m4/data/eft1_after_F_wilson_seed.json")
QUARTIC_SEED = Path("output/general_observable/fermion_first_all16_v1/A_m4/data/eft1_after_F_scalar_quartic_seed.json")
OLD_BOUNDARY = Path("output/general_observable/lo_fermion_boundary/releases/source_021dcc1e8fe8ba0f/A_m4/boundary.json")
OLD_FLAVOUR = Path("output/general_observable/general_wilson_boundary_v1/A_m4/wilson_boundary.json")
UV_REFERENCE = Path("output/general_observable/lo_certification_direct_uv/native_uv_reference_01/models/A_-4.json")
OWNED = Path("output/general_observable/lo_certification_benchmark")
RUN_ID = "full16_general_v1"


class Batch4Error(ValueError):
    pass


def require(ok: bool, why: str) -> None:
    if not ok:
        raise Batch4Error(why)


def load(path: Path) -> dict:
    require(path.is_file(), f"Missing source {path}")
    d = json.loads(path.read_text(encoding="utf-8-sig"))
    require(isinstance(d, dict), f"Not a JSON object: {path}")
    return d


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def z_matrix(raw: Any, name: str) -> np.ndarray:
    require(isinstance(raw, dict) and "real" in raw and "imag" in raw,
            f"{name}: missing real/imag")
    a = np.asarray(raw["real"], float) + 1j * np.asarray(raw["imag"], float)
    require(a.shape == (3, 3) and bool(np.isfinite(a).all()),
            f"{name}: not a finite 3x3 complex matrix")
    return a


def z_json(a: np.ndarray) -> dict:
    a = np.asarray(a, complex)
    require(a.shape == (3, 3) and bool(np.isfinite(a).all()), "Invalid complex 3x3")
    return {"real": a.real.tolist(), "imag": a.imag.tolist()}


def overrides_from_frozen(frozen: dict, export_keys: set[str]) -> dict:
    """All 20 values from one frozen UV state, NO other source/default input."""
    require([frozen["representation"][x] for x in ("d_s1", "d_s2", "d_f", "alpha")]
            == [1, 3, 2, -4], "Frozen model representation mismatch")
    require(float(frozen["base_state"]["mu_gev"]) == 1e7,
            "Frozen reference not at 1e7 GeV")
    scales = frozen["scales"]
    require((scales["mu_fermion_threshold_gev"],
             scales["mu_scalar_threshold_gev"]) == (1e5, 1e3),
            "Changed F/S threshold scales")
    src = {**frozen["base_state"]["sm"], **frozen["base_state"]["ordinary"]}
    require(set(src) == export_keys,
            f"Not one complete frozen/export UV state; frozen only={sorted(set(src)-export_keys)}, "
            f"export only={sorted(export_keys-set(src))}")
    overrides = {}
    for name in sorted(export_keys):
        raw = src[name]
        if isinstance(raw, dict) and set(raw) == {"real", "imag"} and isinstance(raw["real"], list):
            mat = z_matrix(raw, name)
            overrides[name] = [[{"real": float(v.real), "imag": float(v.imag)}
                                for v in row] for row in mat]
        elif isinstance(raw, dict) and set(raw) == {"real", "imag"}:
            v = complex(float(raw["real"]), float(raw["imag"]))
            require(math.isfinite(v.real) and math.isfinite(v.imag), f"Nonfinite {name}")
            overrides[name] = {"real": v.real, "imag": v.imag}
        else:
            require(type(raw) in (int, float) and math.isfinite(raw), f"Invalid {name}")
            overrides[name] = float(raw)
    for name in ("y1", "y2", "MF"):
        require(name in overrides, f"Missing {name}")
    return overrides


def decode_state(raw: Any, name: str) -> np.ndarray:
    require(isinstance(raw, list) and len(raw) == 3, f"{name}: expected 3x3 lists")
    try:
        values = np.asarray([[complex(z["real"], z["imag"]) for z in row]
                             for row in raw], dtype=complex)
    except (TypeError, KeyError, ValueError) as exc:
        raise Batch4Error(f"{name}: incorrect state encoding") from exc
    require(values.shape == (3, 3) and bool(np.isfinite(values).all()),
            f"{name}: invalid matrix")
    return values


def invariant_flavour(a: np.ndarray, b: np.ndarray, masses: np.ndarray) -> np.ndarray:
    """Independent three-family reference for two *barred* native source Yukawas."""
    require(a.shape == b.shape == (3, 3) and masses.shape == (3,), "Invalid heavy basis")
    require(bool((masses > 0).all()), "Nonpositive heavy mass")
    return np.einsum("pr,qr,r->pq", a.conj(), b.conj(), 1/masses, optimize=True)


def check_source(root: Path, symbolic_run_id: str = RUN_ID) -> tuple[dict, dict, dict, dict]:
    """Preflight sources without executing ODE or modifying any file."""
    frozen_path = root / BENCH
    require(digest(frozen_path) == BENCH_SHA, "Frozen fitted reference SHA256 drift")
    frozen = load(frozen_path)
    tree = load(root / FERMION_SEED)
    quartic = load(root / QUARTIC_SEED)
    old = load(root / OLD_BOUNDARY)
    require(tree.get("IntegratedFields") == ["F"] and tree.get("ActiveHeavyFields") == ["S1", "S2"],
            "Tree source does not integrate all fermions and retain both scalars")
    pl = [x for x in tree["TreeWilsonTerms"] if x.get("Chirality") == "PL"]
    require(len(pl) == 1 and pl[0]["Index"] == 2 and tree["TreeWilsonTermCount"] == 2,
            "A_-4 original PL/PR source inventory changed")
    require("DiracProduct[GammaCC, Proj[-1]]" in pl[0]["TermInputForm"],
            "Native PL source changed")
    require(quartic.get("IntegratedFields") == ["F"] and
            quartic.get("ActiveHeavyFields") == ["S1", "S2"] and
            quartic["ScalarQuarticTermCount"] == len(quartic["ScalarQuarticTerms"]),
            "Scalar quartic source not an F-only EFT1 seed")
    require(old.get("model") == MODEL and len(old["sources"]) == 1 and
            old["sources"][0]["inputform_operator"] == pl[0]["TermInputForm"] and
            old["spinor_normalization"]["certified"] is False,
            "Old Matchete source/certification changed")
    return frozen, tree, quartic, old


def source_component_rows(tree: dict, old: dict, running: dict) -> tuple[dict, dict, dict]:
    """Project original native CG, with NEW Yukawas/MF and original source order."""
    from scripts.RunT3GeneralWilsonBoundary import _normalise_term, flavour_matrix
    from Numerical.matching.GeneralT3FlavourC5 import prepare_heavy_flavour
    from Numerical.running.GeneralT3OrderedTensor import project_ordered_term
    from output.general_observable.lo_fermion_boundary.build_boundary import operator_rows

    require(running.get("model") == MODEL and running.get("status") == "NumericalRenormalisableRunningOnly",
            "Wrong running snapshot")
    s = running["uv_at_fermion_threshold"]
    uv = SimpleNamespace(MF=decode_state(s["MF"], "MF(muF)"),
                         y1=decode_state(s["y1"], "y1(muF)"),
                         y2=decode_state(s["y2"], "y2(muF)"),
                         general_couplings={},
                         representation=SimpleNamespace(majorana_fermion=False))
    masses, yukawas, basis = prepare_heavy_flavour(uv, {"y1", "y2"})
    require(np.min(masses) > 0 and np.max(masses)/np.min(masses) < 1e6,
            "Unusable heavy eigenmass basis")
    src = [x for x in tree["TreeWilsonTerms"] if x["Chirality"] == "PL"][0]
    term = _normalise_term(src["TermInputForm"])
    term["source_index"] = src["Index"]
    require(term["heavy_mass_power"] == -1 and len(term["factors"]) == 2,
            "Wrong tree topology")
    require([(f["name"], f["conjugated"]) for f in term["factors"]] ==
            [("y1", True), ("y2", True)], "A_-4 source Yukawa conjugations changed")
    term["flavour_matrix_axes_external_indices"] = [f["external_index"] for f in term["factors"]]
    W = flavour_matrix(term, masses, yukawas)
    W_ref = invariant_flavour(yukawas["y1"], yukawas["y2"], masses)
    require(np.max(abs(W-W_ref)) <= max(np.max(abs(W)),1e-30)*1e-12,
            "Independent source Yukawa contraction mismatches")
    term["ordered_flavour_matrix_gev_inv"] = z_json(W)
    op = project_ordered_term(term, tree["CGRegistry"], MODEL)
    rows, perm = operator_rows(op, term, old["real_scalar_dimension"])
    require(rows and perm == old["sources"][0]["flavour_permutation_lepton_order"],
            "Native source orientation or projection changed")
    offdiag = [row for row in rows if row["coordinate"][1] != row["coordinate"][3]]
    require(bool(offdiag), "Dense-fitted input did not generate off-diagonal LLSS rows")
    diag = {"masses_at_muF_gev": masses.tolist(),
            "heavy_basis_diagnostics": basis,
            "ordered_flavour_matrix_axes": term["flavour_matrix_axes_external_indices"],
            "ordered_flavour_matrix_gev_inv": z_json(W),
            "new_canonical_coordinate_count": len(rows),
            "new_offdiagonal_flavour_coordinate_count": len(offdiag),
            "old_canonical_coordinate_count": len(old["summed_independent_symmetric_coordinates"]),
            "old_offdiagonal_flavour_coordinate_count": sum(
                x["coordinate"][1] != x["coordinate"][3]
                for x in old["summed_independent_symmetric_coordinates"]),
            "lepton_order_permutation": perm,
            "native_source_index": src["Index"]}
    return term, {"rows": rows, "operator": op}, diag


def make_boundary(old: dict, rows: list[dict], term: dict, running_hash: str,
                  source_hash: str, quartic_hash: str) -> dict:
    """Adapter-compatible NEW provisional release, never a replacement for old."""
    require(len(rows) > len(old["summed_independent_symmetric_coordinates"]),
            "Aligned matrix failed to add flavour support")
    import copy
    doc = copy.deepcopy(old)
    s = doc["sources"][0]
    s["independent_symmetric_coordinates"] = rows
    s["flavour_permutation_lepton_order"] = [1, 0]
    s["zero_at_this_benchmark"] = False
    doc["summed_independent_symmetric_coordinates"] = rows
    doc["checksums_sha256"] = {"source_tree": source_hash,
                                  "source_quartic_seed": quartic_hash,
                                  "fitted_running_batch4": running_hash}
    doc["batch4_provenance"] = {
        "status": "REGENERATED_COMMON_INPUT_CONDITIONAL_NOT_CERTIFIED",
        "source_parameter_point": "original frozen fitted benchmark evolved with canonical UV RGBeta",
        "native_source_template": term["source_term_inputform"],
        "all_original_source_reference_hashes_intentionally_not_reused": True,
        "hard_one_loop_at_muF": False,
        "gamma_cc_to_weyl": None}
    return doc


def calculate(root: Path, run_id: str, *, symbolic_run_id: str = RUN_ID) -> dict:
    """Calculate all stages in memory before committing to a unique output dir."""
    import re
    require(re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]{0,50}", run_id) is not None,
            "Unsafe run ID")
    root = root.resolve()
    destination = root / OWNED / ("batch4_common_input_" + run_id)
    require(not destination.exists(), "Output exists; choose a NEW run ID: " + str(destination))
    frozen, tree, quartic, old = check_source(root, symbolic_run_id)
    old_flavour = load(root / OLD_FLAVOUR)
    require(old_flavour.get("model") == MODEL and
            len(old_flavour.get("operator_terms", [])) == 1 and
            old_flavour["operator_terms"][0]["source_term_inputform"] ==
            tree["TreeWilsonTerms"][1]["TermInputForm"],
            "Pinned tree source differs from native source-factorization release")
    from Numerical.orchestration.GeneralRunningBenchmark import (
        load_canonical_general_exports, run as run_rg)
    uv, eft1, export_prov = load_canonical_general_exports(MODEL, symbolic_run_id)
    overrides = overrides_from_frozen(frozen, set(uv["report_betas"]))
    require(set(overrides) == set(uv["report_betas"]), "Not all exported UV couplings pinned")
    result = run_rg(uv, eft1, mu_uv=1e7, mu_f=1e5, mu_s=1e3, overrides=overrides)
    require(result["uv_solver"]["success"] and result["eft1_solver"]["success"],
            "A UV/EFT1 RGE solver failed")
    for name in ("y1", "y2", "MF"):
        check = decode_state(result["uv_initial"][name], name+"(UV)")
        original = z_matrix(frozen["base_state"]["ordinary"][name], name)
        require(np.array_equal(check, original), "Input state changed by runner: "+name)
    for name in ("gY", "g2", "g3", "lambdaH"):
        require(result["uv_initial"][name] == frozen["base_state"]["sm"][name],
                "SM UV input drift: " + name)
    for name in ("lambdaT3",):
        require(result["uv_initial"][name] == frozen["base_state"]["ordinary"][name],
                "UV scalar coupling drift: " + name)
    result.update(symbolic_provenance=export_prov,
                  physical_mass_parameters_available=False,
                  basis_mapping_required_for_physical_interpretation=True,
                  illustrative_benchmark=True,
                  eft1_wilson_running_performed=False,
                  source_benchmark_sha256=BENCH_SHA)
    term, component, flav = source_component_rows(tree, old, result)
    from Numerical.joint_scalar_matching.fermion_boundary_bridge import load_source
    from Numerical.joint_scalar_matching.source_quartic_v4 import build_source_channels
    from Numerical.joint_scalar_matching.numerical_benchmark_v5 import (
        load_running_snapshot, build_ordinary_background, evaluate_ordinary_response)
    from Numerical.joint_scalar_matching.large_batch import parse_matrix
    # The release will be written to a NEW path; do not validate with pinned
    # EXPECTED A_-4 counts (they were tied to the *old diagonal* Yukawas).
    new_boundary = make_boundary(old, component["rows"], term, "pending-run-hash",
                                 digest(root/FERMION_SEED), digest(root/QUARTIC_SEED))
    require(new_boundary["spinor_normalization"]["certified"] is False and
            new_boundary["spinor_normalization"]["gamma_cc_to_master_weyl_factor"] is None,
            "Physical spinor phase wrongly assigned")
    import tempfile
    import shutil
    destination.parent.mkdir(parents=True, exist_ok=True)
    # Keep all files isolated in a temporary directory and publish once.
    staging = Path(tempfile.mkdtemp(prefix=".batch4_staging_", dir=destination.parent))
    try:
        rpath = staging / "fitted_running.json"
        rpath.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n", encoding="utf-8")
        new_boundary["checksums_sha256"]["fitted_running_batch4"] = digest(rpath)
        bpath = staging / "A_m4_fitted_LLSS_provisional.json"
        bpath.write_text(json.dumps(new_boundary,indent=2,allow_nan=False)+"\n", encoding="utf-8")
        source = load_source(bpath, expected_model=MODEL)
        seed = build_source_channels(root/QUARTIC_SEED, MODEL)
        require(seed.source_basis_labels == source.labels, "Quartic/source component bases differ")
        snapshot = load_running_snapshot(rpath, MODEL)
        bg_dir = staging / "numeric_backgrounds"
        info = build_ordinary_background(seed, source, snapshot, bg_dir)
        conditional = evaluate_ordinary_response(source, Path(info["background_json"]))
        response = conditional["response"]
        C5 = parse_matrix(response["C5_LLSS_total_one_loop_per_N_gev_inverse"])
        require(np.max(abs(C5-C5.T)) <= 1e-10 * max(np.max(abs(C5)),1e-30),
                "Computed conditional C5 is not symmetric")
        masses = np.sort(np.linalg.svd(C5,compute_uv=False)*((246.22**2)/2)*1e9)[::-1]
        uvref = load(root/UV_REFERENCE) if (root/UV_REFERENCE).is_file() else None
        summary = {
            "format": FORMAT, "model": MODEL,
            "status": "COMMON_FITTED_INPUT_CONDITIONAL_LLSS_SCALAR_RESPONSE",
            "source_frozen_reference_sha256": BENCH_SHA,
            "source_files_sha256": {str(k):digest(root/k) for k in
                (FERMION_SEED,QUARTIC_SEED,OLD_BOUNDARY,OLD_FLAVOUR)},
            "canonical_symbolic_provenance": export_prov,
            "RG": {"UV": result["scales_gev"]["UV"],"mu_F": result["scales_gev"]["F"],
                   "mu_S": result["scales_gev"]["S"],
                   "UV_success":True,"EFT1_success":True},
            "F_tree_native_source_and_flavour": flav,
            "F_tree_GammaCC_signed_factor": None,
            "conditional_C5_per_N_GeV_inv": z_json(C5),
            "conditional_mnu_eV_per_N": z_json(C5*((246.22**2)/2)*1e9),
            "conditional_neutrino_masses_eV_per_abs_N": masses.tolist(),
            "scalar_response": {"muF_gev": response["mu_f_gev"],
              "muS_gev":response["mu_s_gev"],
              "weinberg_beta_residual":response["weinberg_beta_projection_residual"],
              "weinberg_bubble_residual":response["weinberg_bubble_projection_residual"],
              "matching_scale_cancellation_residual":response["relative_matching_scale_cancellation_residual"],
              "same_muS_for_both_scalars":True,
              "scalar_masses_gev":response["scalar_eigenstate_masses_gev"]},
            "direct_UV_frozen_reference": None if uvref is None else {
                "masses_ev": uvref.get("takagi_masses_ev_from_independent_singular_values"),
                "input_scale_gev": uvref.get("input_scale_gev"),
                "same_RGE_evolved_muF_inputs_as_new_EFT":False,
                "UV_masses_not_a_same_scale_matching_test":True},
            "physical_certification": {
                "physical_C5_certified":False,
                "full_one_loop_neutrino_mass_certified":False,
                "gammaCC_to_Weyl":None,
                "Matchete_to_RGBeta_quartic_alias_normalisation_certified":False,
                "muF_hard_fermion_one_loop_terms_complete":False,
                "EFT1_full_dim5_Wilson_running_included":False,
                "EFT1_LLSS_fixed_order_quartic_log_and_finite_joint_bubble_included":True,
                "native_mass_scalar_basis_map_certified":False},
            "sequence":"UV -> integrate out ALL heavy fermions at muF -> EFT1 with BOTH scalars -> integrate out BOTH scalars at muS -> SMEFT"
        }
        (staging/"conditional_scalar_matching.json").write_text(
            json.dumps(conditional,indent=2,allow_nan=False)+"\n",encoding="utf-8")
        (staging/"summary.json").write_text(
            json.dumps(summary,indent=2,allow_nan=False)+"\n",encoding="utf-8")
        (staging/"report.md").write_text(
            "# A_-4 Batch 4: one dense fitted benchmark, two thresholds (CONDITIONAL)\n\n"
            "All UV couplings are taken from the frozen fitted reference at 1e7 GeV. "
            "Canonical RGBeta evolves that same UV state to 1e5 GeV; the exact "
            "source PL tree/CG is projected with the resulting dense flavour matrices; "
            "EFT1 retains both scalars, then one joint scalar match uses their EFT1 "
            "mass and quartic state at 1e3 GeV.\n\n"
            f"- LLSS independent coordinates: {flav['new_canonical_coordinate_count']} "
            f"({flav['new_offdiagonal_flavour_coordinate_count']} off-diagonal flavour); "
            f"old diagonal example: {flav['old_canonical_coordinate_count']}.\n"
            "- Conditional neutrino masses (eV per |N_GammaCC|): " +
            ", ".join(f"{v:.9g}" for v in masses) + "\n"
            f"- Relative Weinberg projection residual: {max(response['weinberg_beta_projection_residual'],response['weinberg_bubble_projection_residual']):.3e}.\n\n"
            "The previous direct UV masses use the original UV state and are *not* "
            "a same-scale UV/EFT comparison to this evolved EFT result. "
            "GammaCC phase, quartic aliases, hard fermion terms and full Wilson "
            "mixing remain uncertified.\n",
            encoding="utf-8")
        require(not destination.exists(), "Output appeared during run; refusing overwrite")
        staging.rename(destination)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    return {"model": MODEL, "status": summary["status"],
            "output": str(destination),
            "offdiagonal_LLSS_rows":flav["new_offdiagonal_flavour_coordinate_count"],
            "conditional_neutrino_masses_eV_per_abs_N": masses.tolist(),
            "physical_C5_certified":False}


def preflight(root: Path, *, symbolic_run_id: str = RUN_ID) -> dict:
    root = root.resolve()
    frozen,tree,quartic,old = check_source(root,symbolic_run_id)
    from Numerical.orchestration.GeneralRunningBenchmark import load_canonical_general_exports
    uv,eft,info=load_canonical_general_exports(MODEL,symbolic_run_id)
    overrides = overrides_from_frozen(frozen,set(uv["report_betas"]))
    old_flavour=load(root/OLD_FLAVOUR)
    require(old_flavour["model"] == MODEL and len(old_flavour["operator_terms"]) == 1,
            "Missing one-source native provenance")
    require(old_flavour["operator_terms"][0]["source_term_inputform"] ==
            tree["TreeWilsonTerms"][1]["TermInputForm"], "Original source changed")
    return {"status":"PREFLIGHT_OK_NO_OUTPUT_WRITTEN","model":MODEL,
            "complete_frozen_UV_parameter_keys":sorted(overrides),
            "uv_symbolic_export_keys":sorted(uv["report_betas"]),
            "eft1_symbolic_export_keys":sorted(eft["report_betas"]),
            "native_tree_PL_sources":1,
            "both_scalars_retained_at_muF":True,
            "gammaCC_sign":None,
            "output":None}
