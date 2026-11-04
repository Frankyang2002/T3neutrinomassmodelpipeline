"""Qualitative physical heavy-fermion mass/threshold scan of the saved general T3 fits.

This is a LOCAL, non-fitting study.  UV (all fields) -> EFT1 (both scalars,
all fermions removed) -> SMEFT (all scalars removed).  The general Matchete C5
used here is a CONDITIONAL one-scale coefficient evaluated after those running
segments; no claim of a certified sequential hard Wilson matching is made.

Unlike the pre-existing fit evaluators, both UV MF and muF are varied.  Their
fixed-threshold evaluation functions are deliberately NOT called.

From repository root:
  python -X utf8 -m scripts.RunT3GeneralFermionThresholdScan --preflight
  python -X utf8 -m scripts.RunT3GeneralFermionThresholdScan --run-id pilot_v1
  python -X utf8 -m scripts.RunT3GeneralFermionThresholdScan --models all --run-id all16_v1
"""
from __future__ import annotations

import argparse
import csv
from copy import deepcopy
from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import traceback
from typing import Any, Mapping

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = (ROOT / "output/studies/general_scalar_extended/thesis_general_all16/"
            "continue_general16_v1/merged_all16.json")
OUTPUT = ROOT / "output/studies/thesis_qualitative/fermion_threshold"
PINNED_HEAD = "1b4258b19c8da91afe7f352c2df1b1665de61f26"
FACTORS = (0.5, 1.0, 2.0)
PILOT = ("A_-4", "B_-1")
REPS = {"A": (1, 3, 2), "B": (2, 2, 1), "C": (2, 2, 3),
        "D": (3, 1, 2), "E": (3, 3, 2)}
ALPHAS = {"A": (-4, -2, 0), "B": (-3, -1, 1),
          "C": (-3, -1, 1), "D": (-2, 0, 2), "E": (-4, -2, 0, 2)}
MODELS = tuple(f"{family}_{alpha:+d}" for family in ALPHAS for alpha in ALPHAS[family])
MIXED = frozenset({"B_-1", "C_-1"})
OBS = ("delta_m21_sq_ev2", "delta_m3l_sq_ev2", "sin2_theta12",
       "sin2_theta13", "sin2_theta23")
STUDY_KIND = "t3-general-physical-fermion-mass-and-threshold-scan-conditional-v1"


def check(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def read(path: Path) -> dict:
    result = json.loads(path.read_text(encoding="utf-8-sig"))
    check(isinstance(result, dict), f"Expected JSON object: {path}")
    return result


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path: Path, obj: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + ".tmp")
    temp.write_text(json.dumps(obj, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    temp.replace(path)


def tag(model: str) -> str:
    return model.replace("-", "m").replace("+", "p")


def repo_path(raw: str) -> Path:
    """Manifest paths were serialized on Windows; reject path traversal."""
    check(isinstance(raw, str) and bool(raw.strip()), "Invalid manifest result path")
    check(not re.match(r"^[a-zA-Z]:", raw) and not raw.startswith(("/", "\\")),
          "Manifest source must be project-relative")
    parts = raw.replace("\\", "/").split("/")
    check(all(part not in ("", ".", "..") for part in parts),
          "Manifest source contains empty or unsafe path segment")
    return ROOT.joinpath(*parts)


def matrix(raw: Any) -> np.ndarray:
    """Only fitted-config {real,imag} matrix encoding, never guesses about basis."""
    check(isinstance(raw, dict) and set(raw) == {"real", "imag"},
          "Expected fitted {real,imag} matrix")
    z = np.asarray(raw["real"], dtype=float) + 1j * np.asarray(raw["imag"], dtype=float)
    check(z.shape == (3, 3) and np.isfinite(z).all(), "Nonfinite/non-3x3 mass or Yukawa")
    return z


def matrix_json(z: np.ndarray) -> dict:
    z = np.asarray(z, dtype=complex)
    return {"real": z.real.tolist(), "imag": z.imag.tolist()}


def scale_configuration(cfg: dict, factor: float) -> dict:
    """Change physical UV MF and muF only; preserve all UV Yukawa textures."""
    check(factor in FACTORS, "Only the declared 0.5,1,2 factors are supported")
    changed = deepcopy(cfg)
    initial = matrix(cfg["base_state"]["ordinary"]["MF"])
    mu_f = float(cfg["scales"]["mu_fermion_threshold_gev"])
    changed["base_state"]["ordinary"]["MF"] = matrix_json(factor * initial)
    changed["scales"]["mu_fermion_threshold_gev"] = float(factor * mu_f)
    check(np.array_equal(matrix(cfg["base_state"]["ordinary"]["y1"]),
                         matrix(changed["base_state"]["ordinary"]["y1"])),
          "y1 changed during MF rescaling")
    check(np.array_equal(matrix(cfg["base_state"]["ordinary"]["y2"]),
                         matrix(changed["base_state"]["ordinary"]["y2"])),
          "y2 changed during MF rescaling")
    return changed


def scales(cfg: dict) -> tuple[float, float, float, float]:
    x = cfg["scales"]
    mu = (float(cfg["base_state"]["mu_gev"]),
          float(x["mu_fermion_threshold_gev"]),
          float(x["mu_scalar_threshold_gev"]), float(x["mu_low_gev"]))
    check(all(map(math.isfinite, mu)) and mu[0] > mu[1] > mu[2] > mu[3] > 0,
          "UV > fermion > scalar > low hierarchy is violated")
    # Existing B/C low-energy fit implements precisely these last two scales.
    check(mu[2] == 1e3 and mu[3] == 1e2,
          "Existing scalar matching/SMEFT evaluator requires muS=1000, low=100 GeV")
    return mu


def input_mass_spectrum(cfg: dict) -> list[float]:
    values = np.linalg.svd(matrix(cfg["base_state"]["ordinary"]["MF"]),
                           compute_uv=False)
    check(np.all(values > 0), "UV MF has non-positive singular value")
    return sorted(map(float, values))


def checkout_head() -> str | None:
    if not (ROOT / ".git").exists():
        return None
    proc = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                          text=True, capture_output=True, check=False)
    return proc.stdout.strip() if proc.returncode == 0 else None


def select_rows(manifest: dict, models: tuple[str, ...]) -> dict[str, dict]:
    check(manifest.get("accepted_count") == 16 and manifest.get("blocked_count") == 0,
          "Saved manifest does not declare all 16 accepted")
    rows = manifest.get("rows")
    check(isinstance(rows, list) and len(rows) == 16, "Manifest must list exactly 16 models")
    by = {r["model"]: r for r in rows}
    check(len(by) == 16 and set(by) == set(MODELS), "Manifest model set invalid")
    for model in models:
        row = by[model]
        check(row["status"] == "ExploratoryAcceptedFit" and
              all(k in row["prediction"] for k in OBS),
              f"{model}: not an accepted five-observable general fit")
    return by


@dataclass(frozen=True)
class Source:
    model: str
    row: dict
    directory: Path
    config: dict
    result: dict
    paths: dict[str, Path]


def load_source(model: str, row: dict) -> Source:
    path = repo_path(row["result_path"])
    check(path.name == "result.json" and path.parent.name == tag(model),
          f"{model}: manifest-selected fitted result path mismatch")
    source_result = read(path)
    check(source_result.get("model") == model and
          math.isfinite(float(source_result["chi2_final"])) and
          float(source_result["chi2_final"]) <= 5.86 and
          abs(float(source_result["chi2_final"]) - float(row["chi2_final"])) < 1e-8,
          f"{model}: manifest/source fitted chi2 mismatch")
    check(all(math.isclose(float(source_result["prediction"][k]),
                           float(row["prediction"][k]), rel_tol=1e-11, abs_tol=1e-13)
              for k in OBS), f"{model}: manifest/source prediction mismatch")
    cfg_path = path.with_name("fitted_config.json")
    cfg = read(cfg_path)
    rep = cfg["representation"]
    family, a = model.split("_")
    check((rep["d_s1"], rep["d_s2"], rep["d_f"]) == REPS[family]
          and rep["alpha"] == int(a) and rep["shared_scalar"] is False,
          f"{model}: fitted configuration representation mismatch")
    scales(cfg)
    parameters = source_result.get("best_parameters", {})
    for name in ("y1", "y2"):
        y = matrix(cfg["base_state"]["ordinary"][name])
        check(all(f"{name}_{i}{j}" in parameters and
                  np.isclose(y[i,j], float(parameters[f"{name}_{i}{j}"]),
                             rtol=0, atol=1e-13)
                  for i in range(3) for j in range(3)),
              f"{model}: fitted Yukawa differs from selected source best_parameters")
    check(source_result.get("full_general_sequential_matching_performed", 
                            source_result.get("full_sequential_general_scalar_matching_performed")) is False,
          f"{model}: unexpected matching certification/status change")
    paths = {"result": path, "fitted_config": cfg_path}
    if model not in MIXED:
        uv_fixed = path.with_name("general_uv_overrides.json")
        check(uv_fixed.is_file(), f"{model}: missing actual fitted UV overrides")
        paths["general_uv_overrides"] = uv_fixed
    return Source(model, row, path.parent, cfg, source_result, paths)


def ensure_mass_hierarchy(mu_f: float, mu_s: float,
                          heavy_masses: np.ndarray, scalar_masses: np.ndarray) -> None:
    h = np.asarray(heavy_masses, dtype=float)
    s = np.asarray(scalar_masses, dtype=float)
    check(h.shape == (3,) and np.isfinite(h).all() and (h > 0).all(),
          "Invalid evolved heavy-fermion masses")
    check(s.size >= 2 and np.isfinite(s).all() and (s > 0).all(),
          "Invalid scalar masses at scalar matching")
    check(mu_f > mu_s and min(h) > max(s),
          "Physical fermion/scalar mass hierarchy inverted: fermion-first invalid")


def fixed_point_result(fit: Any, low: Any, *, vev: float, ordering: str) -> dict:
    from Numerical.fitting.ParameterScan import (neutrino_mass_from_c5,
                                                charged_lepton_mass_basis_matrix)
    from physics.NeutrinoObservables import calculate_neutrino_observables
    mass = neutrino_mass_from_c5(low.K, vev_gev=vev)
    flavour = charged_lepton_mass_basis_matrix(mass, low.ye)
    neutrinos = calculate_neutrino_observables(flavour, ordering=ordering)
    obs = dict(zip(fit.observable_names, map(float, fit.prediction)))
    check(set(obs) == set(OBS), "Unsupported oscillation observable inventory")
    values = np.asarray(neutrinos.masses_ev, float)
    check(values.shape == (3,) and np.isfinite(values).all(), "Nonfinite neutrino eigenmasses")
    return {"masses_ev": values.tolist(), "sum_masses_ev": float(sum(values)),
            "observables": obs, "chi2_against_fixed_target": float(fit.chi2),
            "takagi_residual": float(neutrinos.takagi_residual)}


def baseline_agreement(actual: Mapping[str, float],
                       expected: Mapping[str, float], *, rtol: float) -> dict:
    """Saved general fit must reproduce before any varied-mass result is trusted."""
    rel = {k: abs(float(actual[k])-float(expected[k])) /
              max(abs(float(expected[k])), 1e-30) for k in OBS}
    return {"pass": all(v <= rtol for v in rel.values()),
            "relative_errors": rel, "relative_tolerance": rtol}


def _physical_input_check(cfg: dict, initial_state: Any) -> None:
    o = cfg["base_state"]["ordinary"]
    for key in ("MF", "y1", "y2"):
        check(np.allclose(np.asarray(getattr(initial_state, key), complex), matrix(o[key]),
                          rtol=0.0, atol=1e-9 if key == "MF" else 1e-14),
              f"UV {key} does not equal modified physical input")


def _nonmixed_heavy_mass_basis(uv_state: Any) -> Any:
    """Apply the final-C5 adapter's explicit default mass-basis tolerances.

    The underlying private routine has four required keyword-only arguments;
    the study must supply the same defaults as extract_final_c5_inputs.
    """
    from Numerical.running.FinalC5TrajectoryAdapter import _heavy_mass_basis_data

    return _heavy_mass_basis_data(
        uv_state,
        offdiagonal_atol=1.0e-10,
        imaginary_atol=1.0e-10,
        takagi_residual_rtol=1.0e-10,
        takagi_residual_atol=1.0e-10,
    )


def _mixed_fitted_scalar_fractions(result: Mapping[str, Any]) -> tuple[float, float]:
    """Retrieve unchanged B/C scalar boundaries from original or continued fits.

    Original fits store fractions at top level; continuation fits carry them
    under matching_and_rge_provenance. Do not silently alter scalar inputs.
    """
    provenance = result.get("matching_and_rge_provenance", {})
    check(isinstance(provenance, dict), "Malformed mixed-fit provenance")
    check(provenance.get("modified_scalar_boundary") is not True,
          "Continuation changed the scalar boundary")

    pairs = (
        ("extra_quartic_fraction", "quartic_fraction"),
        ("extra_mixed_scalar_bilinear_fraction", "mu2_fraction"),
    )
    values = []
    for direct, carried in pairs:
        present_direct = direct in result
        present_carried = carried in provenance
        check(present_direct or present_carried,
              f"Missing mixed-scalar boundary fraction {direct}/{carried}")
        value = float(result[direct] if present_direct else provenance[carried])
        check(math.isfinite(value), f"Nonfinite mixed-scalar fraction {direct}")
        if present_direct and present_carried:
            other = float(provenance[carried])
            check(math.isfinite(other) and
                  math.isclose(value, other, rel_tol=1e-12, abs_tol=1e-15),
                  f"Conflicting mixed-scalar boundary fraction {direct}/{carried}")
        values.append(value)

    quartic_fraction, mu2_fraction = values
    check(0.0 < quartic_fraction <= 1.0 and 0.0 < mu2_fraction <= 0.2,
          "Mixed-scalar fitted candidate boundary fractions out of range")
    return quartic_fraction, mu2_fraction


def _nonmixed(source: Source, cfg: dict, *, canonical_run_id: str,
              wolframscript: str, rtol: float, atol: float) -> dict:
    from Numerical.orchestration.GeneralRunningBenchmark import (
        initial_state, load_canonical_general_exports, _state_values)
    from Numerical.running.UVRunner import run_uv_segment
    from Numerical.running.IntermediateScalarRunner import run_intermediate_scalar_segment
    from Numerical.running.FermionThresholdBoundary import build_intermediate_scalar_boundary
    from Numerical.running.IntermediateScalarStateVector import pack_intermediate_scalar_state
    from Numerical.fitting.OscillationFit import OscillationFitTarget
    from Numerical.orchestration.RunSmallLambdaRefit import TARGET
    from scripts import RunT3GeneralAll16FromReference as general
    from scripts import RunT3GeneralScalarBCExploratoryFit as bc

    model = source.model
    uv_betas, eft_betas, provenance = load_canonical_general_exports(model, canonical_run_id)
    coefficient = general.ensure_general_matching(model, generate=False,
                                                    wolframscript=wolframscript)
    from Numerical.matching.GeneralT3FlavourC5 import GeneralT3FlavourC5
    kernel = GeneralT3FlavourC5(model, coefficient, wolframscript=wolframscript)
    fixed = read(source.paths["general_uv_overrides"])
    check(set(fixed) == set(uv_betas["report_betas"]),
          "Saved UV override keys differ from full-general beta keys")
    # The selected fitted source supplies the actual full-general UV couplings.
    # Only MF changes relative to that source; NO restricted-kernel substitution.
    fixed["MF"] = general.entry_matrix(matrix(cfg["base_state"]["ordinary"]["MF"]))
    mu_uv, mu_f, mu_s, _ = scales(cfg)
    which, start = initial_state(uv_betas, mu_uv=mu_uv, mu_f=mu_f,
                                 mu_s=mu_s, overrides=fixed)
    check(which == model, "UV beta model mismatch")
    _physical_input_check(cfg, start)
    uv = run_uv_segment(start, uv_betas, mu_f, rtol=rtol, atol=atol)
    check(uv.solver_success and math.isclose(uv.final_state.mu_gev, mu_f, rel_tol=1e-10),
          "UV RGE did not reach new muF")
    heavy = _nonmixed_heavy_mass_basis(uv.final_state)
    boundary = build_intermediate_scalar_boundary(uv.final_state)
    check(math.isclose(boundary.mu_gev, mu_f, rel_tol=1e-10),
          "Fermion boundary is not at new muF")
    _, layout = pack_intermediate_scalar_state(boundary)
    check(set(layout.names) == set(eft_betas["report_betas"]),
          "EFT1 beta/state inventory mismatch")
    eft = run_intermediate_scalar_segment(boundary, eft_betas, mu_s,
                                          rtol=rtol, atol=atol)
    check(eft.solver_success and math.isclose(eft.final_state.mu_gev, mu_s, rel_tol=1e-10),
          "EFT1 RGE did not reach scalar threshold")
    scalar_masses = np.sqrt(np.array([eft.final_state.mS1Sq, eft.final_state.mS2Sq]))
    ensure_mass_hierarchy(mu_f, mu_s, np.asarray(heavy.masses_gev), scalar_masses)
    c5 = np.asarray(kernel.evaluate(uv.final_state, eft.final_state), complex)
    check(c5.shape == (3,3) and np.isfinite(c5).all(), "Invalid general C5")
    c5 = 0.5*(c5+c5.T)
    target = OscillationFitTarget.from_json(TARGET)
    fit, low = bc.low_energy_fit(c5, eft.final_state, target,
                                 vev=float(cfg["numerical"]["vev_gev"]))
    values = fixed_point_result(fit, low, vev=float(cfg["numerical"]["vev_gev"]),
                                ordering=cfg["ordering"])
    values.update({
        "uv_endpoint_gev": float(uv.final_state.mu_gev),
        "eft1_initial_gev": float(boundary.mu_gev),
        "eft1_endpoint_gev": float(eft.final_state.mu_gev),
        "sm_endpoint_gev": float(low.mu_gev) if hasattr(low, "mu_gev") else 100.,
        "heavy_masses_at_muF_gev": list(map(float, heavy.masses_gev)),
        "scalar_masses_at_muS_gev": scalar_masses.tolist(),
        "uv_initial_mf_gev": matrix_json(np.asarray(start.MF)),
        "uv_endpoint_mf_gev": matrix_json(np.asarray(uv.final_state.MF)),
        "uv_nfev": int(uv.nfev), "eft1_nfev": int(eft.nfev),
        "c5_at_scalar_threshold_gev_inv": matrix_json(c5),
        "c5_at_low_scale_gev_inv": matrix_json(low.K),
        "matching_diagnostics": kernel.last_diagnostics,
        "general_symbolic_provenance": provenance,
        "matching_coefficient_sha256": digest(coefficient),
        "matching_scope": "full-general Matchete channel one-scale conditional; NOT certified sequential EFT1 Wilson matching",
    })
    return values


def _mixed(source: Source, cfg: dict, *, rtol: float, atol: float) -> dict:
    """Variable-muF native B/C quadratic RGE integration; no fixed-scale evaluator call."""
    from Numerical.experimental.GeneralT3ScalarRGState import (
        unpack_general_scalar_state, project_general_after_fermion, pack_general_scalar_state)
    from Numerical.fitting.OscillationFit import OscillationFitTarget
    from Numerical.orchestration.RunSmallLambdaRefit import TARGET
    from scripts import RunT3BCGeneralFromReference as bc_source
    from scripts import RunT3GeneralScalarBCExploratoryFit as bc
    from scripts import RunT3GeneralScalarBCPilot as pilot
    from scripts import RunT3ZeroExtraBoundaryBenchmark as native
    from scripts.RunGeneralT3NumericalC5 import DEFAULT_SOURCE
    from scripts.AnalyzeT3BCMatchingConventions import endpoint
    from Numerical.experimental.GeneralT3MixedScalarMassBasis import MixedDoubletMassBasis

    model = source.model
    # The reference constructor validates registry, representation and actual base
    # inputs; replace only with the manifest-selected FITTED config, not its ref Yukawas.
    record = bc_source.build_reference_source(model)
    record["cfg"] = deepcopy(cfg)
    record["lambdaT3"] = complex(float(cfg["base_state"]["ordinary"]["lambdaT3"]["real"]),
                                  float(cfg["base_state"]["ordinary"]["lambdaT3"]["imag"]))
    record["uv_on"] = True
    registry = record["registry"]
    model_quartic_fraction, model_mu2_fraction = _mixed_fitted_scalar_fractions(
        source.result)
    # The actual used native candidate RGEs, not canonical non-mixed kernels.
    bundle = bc.beta_bundle(model, record, registry)
    mu_uv, mu_f, mu_s, _ = scales(cfg)
    y0, ly0, extra, start_mu2 = pilot.prepare_initial(
        record, registry, "general_quartics_and_mixing",
        model_quartic_fraction, model_mu2_fraction)
    uv_initial = unpack_general_scalar_state(y0[:ly0.size], ly0, registry,
                                             uv=True, mu_gev=mu_uv).baseline
    _physical_input_check(cfg, uv_initial)
    uv_sol = native.integrate_segment(
        y0, lambda t,y: native.rhs_at(t,y,layout=ly0,registry=registry,
                                     betas=bundle["UV"],stage="UV"),
        mu_uv, mu_f, n_points=3)
    at_f = np.asarray(uv_sol.y[:,-1], float)
    uv_full = unpack_general_scalar_state(at_f[:ly0.size], ly0, registry,
                                          uv=True, mu_gev=mu_f)
    uv_end = uv_full.baseline
    heavy = bc.heavy_flavour_at_muF(uv_end)
    projected = project_general_after_fermion(uv_full, registry)
    ye, lye = pack_general_scalar_state(projected, registry)
    ye = np.concatenate((ye, at_f[ly0.size:]))
    check(np.array_equal(ye[lye.size:], at_f[ly0.size:]),
          "Fermion projection changed mixed quadratic boundary")
    new, reused = record["new"], record["reused"]
    eft_betas = native.check_payload(
        bundle["payloads"]["EFT1"],model=model,stage="EFT1",
        rep=lye.representation,base_names=tuple(lye.names[:-len(new)]),
        new_names=new,reused_names=reused)
    eft_sol = native.integrate_segment(
        ye, lambda t,y: native.rhs_at(t,y,layout=lye,registry=registry,
                                     betas=eft_betas,stage="EFT1"),
        mu_f, mu_s, n_points=3)
    y_s = np.asarray(eft_sol.y[:,-1], float)
    snap = native.snapshot(y_s, lye, registry,stage="EFT1",mu_gev=mu_s,
                           new=new,reused=reused)
    eft_state = unpack_general_scalar_state(y_s[:lye.size], lye, registry,
                                            uv=False, mu_gev=mu_s).baseline
    s = endpoint({"eft1_end": snap})
    mixing = MixedDoubletMassBasis.build(model, s["m1"], s["m2"], s["mu"])
    check(not mixing.tachyonic and not mixing.degenerate,
          "Invalid/degenerate B/C scalar mixing eigenbasis")
    scalar_masses = np.sqrt(np.asarray(mixing.eigenvalues_sq, dtype=float))
    ensure_mass_hierarchy(mu_f, mu_s, heavy[0], scalar_masses)
    coefficient = DEFAULT_SOURCE / (tag(model) +
        "_final_one_loop_probe_weinberg_ordered_coefficient.inputform.txt")
    compiled = bc.compile_source(coefficient)
    c5 = bc.c5_at_endpoint(model, snap, heavy, compiled)
    target = OscillationFitTarget.from_json(TARGET)
    fit, low = bc.low_energy_fit(c5, eft_state, target,
                                 vev=float(cfg["numerical"]["vev_gev"]))
    values = fixed_point_result(fit, low, vev=float(cfg["numerical"]["vev_gev"]),
                                ordering=cfg["ordering"])
    values.update({
        "uv_endpoint_gev": mu_f, "eft1_initial_gev": mu_f,
        "eft1_endpoint_gev": mu_s,
        "sm_endpoint_gev": 100.0,
        "heavy_masses_at_muF_gev": heavy[0].tolist(),
        "scalar_masses_at_muS_gev": scalar_masses.tolist(),
        "uv_initial_mf_gev": matrix_json(uv_initial.MF),
        "uv_endpoint_mf_gev": matrix_json(uv_end.MF),
        "heavy_mass_basis_residual": heavy[3],
        "uv_nfev": int(uv_sol.nfev), "eft1_nfev": int(eft_sol.nfev),
        "mixed_quadratic_initial_gev2": {"real": start_mu2.real, "imag": start_mu2.imag},
        "mixed_quadratic_at_muS_gev2": snap["quadratic_mu2_gev2"],
        "scalar_mixing_U": matrix_json(np.asarray(mixing.U)),
        "c5_at_scalar_threshold_gev_inv": matrix_json(c5),
        "c5_at_low_scale_gev_inv": matrix_json(low.K),
        "matching_coefficient_sha256": digest(coefficient),
        "matching_scope": "native mixed-scalar conditional signed-direct-quartic one-scale; NOT certified sequential EFT1 Wilson matching",
        "conditional_sign_map": "lambda_Matchete = -lambda_RGBeta for two direct quartics; not certified",
    })
    return values


def preflight(source: Source, *, canonical_run_id: str, wolframscript: str) -> dict:
    problems: list[str] = []
    checked = {key: {"project_path": str(p.relative_to(ROOT)), "sha256": digest(p)}
               for key, p in source.paths.items()}
    if source.model in MIXED:
        try:
            from scripts import RunT3BCGeneralFromReference as s
            from scripts import RunT3GeneralScalarBCExploratoryFit as mixed
            from scripts import RunT3ZeroExtraBoundaryBenchmark as b
            from scripts.RunGeneralT3NumericalC5 import DEFAULT_SOURCE
            p = DEFAULT_SOURCE / (tag(source.model)+
                "_final_one_loop_probe_weinberg_ordered_coefficient.inputform.txt")
            checked["mixed_coefficient"] = {"project_path": str(p.relative_to(ROOT)),
                                             "sha256": digest(p)}
            mixed.compile_source(p)
            for stage in ("UV", "EFT1"):
                bp = b.beta_path(source.model, stage)
                checked[f"native_{stage.lower()}_beta"] = {
                    "project_path": str(bp.relative_to(ROOT)), "sha256": digest(bp)}
            s.build_reference_source(source.model)
        except Exception as exc:
            problems.append(f"B/C native candidate dependencies: {type(exc).__name__}: {exc}")
    else:
        try:
            from Numerical.orchestration.GeneralRunningBenchmark import load_canonical_general_exports
            from scripts.RunT3GeneralAll16FromReference import ensure_general_matching
            uv, eft, info = load_canonical_general_exports(source.model, canonical_run_id)
            p = ensure_general_matching(source.model, generate=False,
                                        wolframscript=wolframscript)
            checked["general_matchete_c5"] = {"project_path": str(p.relative_to(ROOT)),
                                               "sha256": digest(p)}
            checked["canonical_uv_export_sha256"] = info["uv_sha256"]
            checked["canonical_eft1_export_sha256"] = info["eft1_sha256"]
            fixed = read(source.paths["general_uv_overrides"])
            check(set(fixed) == set(uv["report_betas"]),
                  "Saved UV couplings differ from full-general UV beta keys")
        except Exception as exc:
            problems.append(f"nonmixed full-general dependencies: {type(exc).__name__}: {exc}")
        if shutil.which(wolframscript) is None:
            problems.append(f"Wolfram executable not found: {wolframscript}")
    return {"model": source.model, "status": "Ready" if not problems else "Blocked",
            "problems": problems, "inputs": checked}


def trend(reference: dict, varying: dict) -> dict:
    result = {}
    for key in ("m1_ev", "m2_ev", "m3_ev", "sum_masses_ev", *OBS):
        a = (reference["masses_ev"][int(key[1])-1] if key.startswith("m") and
             key.endswith("_ev") and len(key)==5 else
             reference["observables"].get(key, reference.get(key)))
        b = (varying["masses_ev"][int(key[1])-1] if key.startswith("m") and
             key.endswith("_ev") and len(key)==5 else
             varying["observables"].get(key, varying.get(key)))
        a, b = float(a), float(b)
        rel = (b/a-1) if a != 0 else None
        result[key] = {"baseline": a, "value": b,
                       "absolute_change": b-a, "relative_change": rel,
                       "direction": "increase" if b>a else "decrease" if b<a else "unchanged"}
    return result



def representation_report(paths: list[dict], dest: Path) -> None:
    """Per-model directions, and family medians only across actually computed models."""
    rows = []
    family_data: dict[str, dict[str, list[float]]] = {}
    for item in paths:
        if item.get("status") != "ComputedConditional" or "relative" not in item:
            continue
        model = item["model"]
        factor = float(item["factor"])
        family = model.split("_")[0]
        for name, values in item["relative"].items():
            relative = values["relative_change"]
            rows.append({"model": model, "family": family,
                         "dS1_dS2_dF": "-".join(map(str, REPS[family])),
                         "factor": factor, "quantity": name,
                         "baseline": values["baseline"], "value": values["value"],
                         "absolute_change": values["absolute_change"],
                         "relative_change": relative,
                         "direction": values["direction"]})
            if relative is not None:
                family_data.setdefault(family, {}).setdefault(
                    f"MF_x{factor:g}:{name}", []).append(relative)
    with (dest/"trends.csv").open("w", newline="", encoding="utf8") as handle:
        fields = ("model", "family", "dS1_dS2_dF", "factor", "quantity",
                  "baseline", "value", "absolute_change", "relative_change", "direction")
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    agg = {family: {name: {"median_relative_change":float(np.median(values)),
                           "n_models":len(values)} for name,values in groups.items()}
           for family,groups in family_data.items()}
    save(dest/"representation_dependence.json", {
        "family_representation_dS1_dS2_dF": REPS,
        "conditional_matching_only": True,
        "family_medians_only_from_computed_models": agg,
        "rows_recorded": len(rows),
        "note": "Representation differences are conditional on distinct fitted UV textures, "
                "not isolated pure group-theory factors."})


def summary(paths: list[dict], dest: Path, *, head: str | None) -> None:
    save(dest/"summary.json", {"format": STUDY_KIND, "main_head": PINNED_HEAD,
                               "local_checkout_head": head, "rows": paths,
                               "fully_certified_sequential_matching": False})
    representation_report(paths, dest)
    with (dest/"summary.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(("model", "factor", "status", "muF_GeV", "reason", "m1_eV",
                    "m2_eV", "m3_eV", *OBS, "m3_relative_change"))
        for row in paths:
            data = row.get("prediction", {})
            ms = data.get("masses_ev", [None]*3)
            w.writerow((row["model"], row["factor"], row["status"], row.get("muF_gev"),
                        row.get("reason", ""), *ms,
                        *(data.get("observables",{}).get(k) for k in OBS),
                        row.get("relative",{}).get("m3_ev",{}).get("relative_change")))


def run(models: tuple[str, ...], *, run_id: str, preflight_only: bool,
        baseline_rtol: float, canonical_run_id: str, wolframscript: str,
        rtol: float, atol: float, allow_different_head: bool) -> int:
    check(re.fullmatch(r"[A-Za-z0-9_-]+", run_id) is not None,
          "--run-id supports only alphanumerics/underscore/hyphen")
    head = checkout_head()
    if head and head != PINNED_HEAD and not allow_different_head:
        raise RuntimeError(f"Local checkout {head} is not verified main {PINNED_HEAD}")
    by = select_rows(read(MANIFEST), models)
    dest = OUTPUT / run_id
    check(not dest.exists(), f"Output already exists; choose another --run-id: {dest}")
    dest.mkdir(parents=True, exist_ok=False)
    rows: list[dict] = []
    save(dest/"settings.json", {"format": STUDY_KIND, "models": models, "factors": FACTORS,
         "manifest_path": str(MANIFEST.relative_to(ROOT)),
         "manifest_sha256": digest(MANIFEST), "main_head": PINNED_HEAD,
         "checkout_head": head, "canonical_run_id": canonical_run_id,
         "rtol": rtol, "atol": atol, "baseline_rtol": baseline_rtol,
         "preflight_only": preflight_only,
         "matching_certification": "one-scale conditional, not fully sequential"})
    failed = False
    for model in models:
        try:
            source = load_source(model, by[model])
            before = input_mass_spectrum(source.config)
            pf = preflight(source, canonical_run_id=canonical_run_id,
                           wolframscript=wolframscript)
            save(dest/tag(model)/"preflight.json", pf)
            if pf["status"] != "Ready":
                raise RuntimeError("; ".join(pf["problems"]))
            if preflight_only:
                rows.append({"model": model, "factor": 1.0, "status": "Ready",
                             "muF_gev": scales(source.config)[1],
                             "source_result": str(source.paths["result"].relative_to(ROOT))})
                print(f"{model}: Ready", flush=True)
                continue
            baseline: dict | None = None
            for factor in (1., 0.5, 2.):  # Gate varied runs on a reproduced baseline.
                changed = scale_configuration(source.config, factor)
                mu_uv, mu_f, mu_s, mu_low = scales(changed)
                updated = input_mass_spectrum(changed)
                check(np.allclose(updated, factor*np.asarray(before), rtol=1e-12, atol=0),
                      "UV physical heavy mass spectrum did not rescale")
                if model in MIXED:
                    prediction = _mixed(source, changed, rtol=rtol, atol=atol)
                else:
                    prediction = _nonmixed(source, changed, canonical_run_id=canonical_run_id,
                                           wolframscript=wolframscript, rtol=rtol, atol=atol)
                check(math.isclose(prediction["uv_endpoint_gev"], mu_f, rel_tol=1e-10)
                      and math.isclose(prediction["eft1_initial_gev"], mu_f, rel_tol=1e-10)
                      and math.isclose(prediction["eft1_endpoint_gev"], mu_s, rel_tol=1e-10),
                      "Numerical UV/EFT1 endpoints do not track updated thresholds")
                record = {"model": model, "factor": factor, "status": "ComputedConditional",
                          "muF_gev": mu_f, "muS_gev": mu_s, "muUV_gev": mu_uv,
                          "muLow_gev": mu_low,
                          "uv_physical_MF_singular_values_gev": updated,
                          "source_result": str(source.paths["result"].relative_to(ROOT)),
                          "source_input_sha256": {k:digest(p) for k,p in source.paths.items()},
                          "prediction": prediction,
                          "representation": changed["representation"],
                          "eft_order": "UV(all fields)->EFT1(both scalars)->SMEFT"}
                if factor == 1.:
                    agreement = baseline_agreement(prediction["observables"],
                                                    source.row["prediction"],
                                                    rtol=baseline_rtol)
                    record["saved_baseline_validation"] = agreement
                    if not agreement["pass"]:
                        record["status"] = "BaselineMismatch"
                        save(dest/tag(model)/"factor_1.json", record)
                        rows.append(record)
                        raise RuntimeError("Saved general-fit baseline not reproduced; "
                                           "mass variation suppressed")
                    baseline = prediction
                else:
                    check(baseline is not None, "Baseline was not computed")
                    record["relative"] = trend(baseline, prediction)
                save(dest/tag(model)/f"factor_{str(factor).replace('.', 'p')}.json", record)
                rows.append(record)
                print(f"{model} MF factor={factor:g}: "
                      f"m_eV={np.asarray(prediction['masses_ev']).round(8).tolist()} "
                      f"chi2={prediction['chi2_against_fixed_target']:.6g}", flush=True)
        except Exception as exc:
            failed = True
            err = {"model": model, "factor": None, "status": "Blocked",
                   "reason": f"{type(exc).__name__}: {exc}",
                   "traceback": traceback.format_exc(limit=8)}
            rows.append(err)
            save(dest/tag(model)/"blocked.json", err)
            print(f"{model}: BLOCKED {err['reason']}", flush=True)
        summary(rows, dest, head=head)
    return 2 if failed else 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--models", default=",".join(PILOT),
                        help="Comma-separated 16-model identifiers or 'all'")
    parser.add_argument("--run-id", default="pilot_v1")
    parser.add_argument("--preflight", action="store_true")
    parser.add_argument("--canonical-run-id", default="full16_general_v1")
    parser.add_argument("--wolframscript", default="wolframscript")
    parser.add_argument("--baseline-rtol", type=float, default=0.005)
    parser.add_argument("--rtol", type=float, default=1e-7)
    parser.add_argument("--atol", type=float, default=1e-10)
    parser.add_argument("--allow-different-head", action="store_true")
    args = parser.parse_args(argv)
    requested = MODELS if args.models.lower() == "all" else tuple(
        x.strip() for x in args.models.split(",") if x.strip())
    if not requested or len(requested) != len(set(requested)) or set(requested)-set(MODELS):
        parser.error("--models must be nonrepeated A_-4,... or 'all'")
    if not (0 < args.baseline_rtol < 1 and 0 < args.rtol < 1 and args.atol > 0):
        parser.error("All tolerances must be positive and baseline-rtol < 1")
    return run(tuple(m for m in MODELS if m in requested), run_id=args.run_id,
               preflight_only=args.preflight, baseline_rtol=args.baseline_rtol,
               canonical_run_id=args.canonical_run_id,
               wolframscript=args.wolframscript, rtol=args.rtol, atol=args.atol,
               allow_different_head=args.allow_different_head)


if __name__ == "__main__":
    raise SystemExit(main())
