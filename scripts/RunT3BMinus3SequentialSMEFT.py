"""Numerical B(-3) handoff: saved fermion-first sequential C5 -> SMEFT -> observables.

This additive driver uses the *existing* successful B(-3) symbolic sequential
Matchete result and the *existing* general-fit UV/EFT1 numerical threshold
snapshots. It does not evaluate the unrelated one-scale C5, refit Yukawas,
or claim to re-run the threshold matching in this process.

Run from the repository root (PowerShell):
    python -X utf8 -m scripts.RunT3BMinus3SequentialSMEFT --check-only
    python -X utf8 -m scripts.RunT3BMinus3SequentialSMEFT --run-id first_v1

Output: output/general_sequential_pilot/B_m3_sequential_smeft_numerical/<run-id>/

This is a conditional first benchmark. It has not yet linked the fitted UV
trajectory *into* the symbolic Matchete threshold run on each evaluation, and
it does not constitute a general all-16 sequential matching implementation.
"""

from __future__ import annotations

import argparse
from dataclasses import fields
from hashlib import sha256
import json
import math
from pathlib import Path
import re
from typing import Any, Mapping

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
MODEL = "B_-3"
DEFAULT_C5 = (
    ROOT / "output/general_sequential_pilot/B_m3_matching/data/"
    "general_eft1_pilot/general_sequential_final_c5.json"
)
DEFAULT_FIT = (
    ROOT / "output/studies/general_scalar_extended/thesis_general_all16/"
    "thesis16_general_c5_full_v1/B_m3"
)
DEFAULT_OUTPUT = ROOT / "output/general_sequential_pilot/B_m3_sequential_smeft_numerical"
SCHEMA = "t3-b-minus-3-saved-sequential-to-smeft-numerical-v1"
SM_FIELDS = frozenset(("gY", "g2", "g3", "lambdaH", "yu", "yd", "ye"))


def read_object(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(f"Required input is absent: {path}")
    value = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def decode(value: Any) -> Any:
    """Decode canonical {real,imag} scalar/matrix entries without changing bases."""
    if isinstance(value, Mapping) and set(value) == {"real", "imag"}:
        a = np.asarray(value["real"], dtype=float)
        b = np.asarray(value["imag"], dtype=float)
        if a.shape != b.shape:
            raise ValueError("Real/imaginary shapes differ")
        if not np.isfinite(a).all() or not np.isfinite(b).all():
            raise ValueError("Non-finite encoded complex value")
        if a.ndim == 0:
            return complex(float(a), float(b))
        return a + 1j * b
    if isinstance(value, list):
        values = [decode(item) for item in value]
        if all(isinstance(item, (int, float, complex, np.number)) for item in values):
            return np.asarray(values, dtype=complex)
        if all(isinstance(item, np.ndarray) for item in values):
            return np.stack(values)
        return values
    if isinstance(value, Mapping):
        return {str(key): decode(item) for key, item in value.items()}
    if isinstance(value, (str, bool)) or value is None:
        return value
    if isinstance(value, (int, float)) and math.isfinite(float(value)):
        return value
    raise ValueError(f"Unsupported/non-finite numerical value: {value!r}")


def finite_real(value: Any, label: str) -> float:
    z = complex(decode(value))
    if not math.isfinite(z.real) or abs(z.imag) > 1e-10 * max(1.0, abs(z)):
        raise ValueError(f"{label}: expected finite real value, got {z!r}")
    return float(z.real)


def finite_matrix(value: Any, label: str) -> np.ndarray:
    result = np.asarray(decode(value), dtype=complex)
    if result.shape != (3, 3) or not (
        np.isfinite(result.real).all() and np.isfinite(result.imag).all()
    ):
        raise ValueError(f"{label}: expected a finite complex 3x3 matrix")
    return result


def encode_matrix(value: np.ndarray) -> dict[str, list]:
    a = np.asarray(value, dtype=complex)
    if a.shape != (3, 3) or not (
        np.isfinite(a.real).all() and np.isfinite(a.imag).all()
    ):
        raise ValueError("Cannot encode an invalid 3x3 matrix")
    return {"real": a.real.tolist(), "imag": a.imag.tolist()}


def source_preflight(c5_path: Path, fit_result_path: Path,
                     fit_config_path: Path) -> tuple[dict, dict, dict, dict]:
    c5 = read_object(c5_path)
    pilot = read_object(c5_path.with_name("pilot_summary.json"))
    if pilot.get("status") != "GeneralBMinus3SequentialEFT1Executed" or (
        pilot.get("resumed_scalar_threshold", {}).get("status") != "Success"
    ) or (
        pilot.get("resumed_scalar_threshold", {}).get("running_insertion_loaded") is not True
    ) or (
        pilot.get("rge", {}).get("weinberg_subspace_validation", {}).get(
            "matches_weinberg_subspace") is not True
    ):
        raise ValueError("B(-3) sequential Wilson running and scalar threshold were not completed")
    result = read_object(fit_result_path)
    config = read_object(fit_config_path)
    if c5.get("status") != "Success" or (
        c5.get("scheme") != "fixed_order_one_loop_MSbar_hierarchical"
    ):
        raise ValueError("B(-3) final C5 is not the saved sequential MSbar result")
    if c5.get("renormalized_hard_threshold", {}).get("status") != "Success":
        raise ValueError("Sequential C5 is missing its renormalised hard term")
    if not c5.get("combined", {}).get("ready_for_physical_majorana_numerics"):
        raise ValueError("Sequential final C5 is not ready for numerical substitution")
    if c5.get("checks", {}).get("equal_scale_running_vanishes") is not True:
        raise ValueError("Expected sequential leading-log matching checks absent")
    if result.get("model") != MODEL or not result.get("general_rges_used"):
        raise ValueError("Expected the saved B(-3) general UV/EFT1 numerical result")
    if config.get("representation") != {
        "d_s1": 2, "d_s2": 2, "d_f": 1, "alpha": -3,
        "shared_scalar": False,
    }:
        raise ValueError("Fitted configuration is not ordinary B(-3)")
    details = result.get("numerical_details", {})
    if (int(details.get("uv_nfev", 0)) <= 0 or
            int(details.get("eft1_nfev", 0)) <= 0):
        raise ValueError("Saved benchmark has no numerical UV/EFT1 integration evidence")
    if not isinstance(details.get("uv_at_fermion_threshold"), dict) or not (
        isinstance(details.get("eft1_at_scalar_threshold"), dict)
    ):
        raise ValueError("Saved fitted UV and EFT1 numerical endpoints are missing")
    scales = config.get("scales", {})
    if not (0 < float(scales["mu_scalar_threshold_gev"]) <
            float(scales["mu_fermion_threshold_gev"])):
        raise ValueError("Expected mu_F > mu_S > 0")
    if not (0 < float(scales["mu_low_gev"]) <
            float(scales["mu_scalar_threshold_gev"])):
        raise ValueError("Expected mu_S > mu_low > 0")
    return c5, result, config, pilot


def fitted_fermion_state(uv: Mapping[str, Any], config: Mapping[str, Any],
                         mu_f: float):
    """Recreate the saved mu_F state for the established heavy SVD convention."""
    from Numerical.core.State import SMNumericalState, T3Representation, T3UVState

    rep = T3Representation(**config["representation"]).validated()
    sm = SMNumericalState(
        **{name: finite_matrix(uv[name], name) if name in ("yu", "yd", "ye")
           else finite_real(uv[name], name) for name in SM_FIELDS}
    ).validated()
    known = {item.name for item in fields(T3UVState)} | SM_FIELDS
    base_names = ("y1", "y2", "MF")
    real_names = ("mS1Sq", "mS2Sq", "lambdaS1", "lambdaS2",
                  "lambdaH1", "lambdaH2", "lambda12")
    kwargs: dict[str, Any] = {
        "mu_gev": mu_f, "representation": rep, "sm": sm,
        **{k: finite_matrix(uv[k], k) for k in base_names},
        **{k: finite_real(uv[k], k) for k in real_names},
        "lambdaT3": complex(decode(uv["lambdaT3"])),
    }
    optional = ("lambdaH1Adj", "lambdaH2Adj", "lambdaS1Adj",
                "lambdaS2Adj", "lambda12Adj", "lambda12Cross")
    for key in optional:
        if key in uv:
            kwargs[key] = finite_real(uv[key], key)
    complex_optional = (
        "lambdaHHdagS2S2", "lambdaHHdagS1barS1bar", "lambdaS1bar2S2bar2",
        "lambdaS1barS2S2bar2", "lambdaS1S1bar2S2bar",
        "lambdaHHdagS1barS2barCross",
    )
    for key in complex_optional:
        if key in uv:
            kwargs[key] = complex(decode(uv[key]))
    kwargs["general_couplings"] = {
        key: decode(value) for key, value in uv.items() if key not in known
    }
    return T3UVState(**kwargs).validated()


def run_numerics(c5_path: Path, fit_result: dict, fit_config: dict,
                 target_path: Path) -> dict[str, Any]:
    from Numerical.running.FinalC5TrajectoryAdapter import (
        FinalC5EvaluationInputs, _heavy_mass_basis_data,
    )
    from Numerical.running.SMWeinbergStage import evaluate_final_weinberg_json
    from Numerical.running.SMWeinbergEvolution import (
        SMWeinbergInitialConditions, evolve_sm_weinberg,
    )
    from Numerical.fitting.OscillationFit import (
        OscillationFitTarget, evaluate_mass_matrix_fit,
    )
    from physics.NeutrinoTrajectory import (
        charged_lepton_mass_basis_matrix, neutrino_mass_from_c5,
    )
    from physics.NeutrinoObservables import calculate_neutrino_observables

    scales = fit_config["scales"]
    mu_f = float(scales["mu_fermion_threshold_gev"])
    mu_s = float(scales["mu_scalar_threshold_gev"])
    mu_low = float(scales["mu_low_gev"])
    vev = float(fit_config.get("numerical", {}).get("vev_gev", 246.22))
    detail = fit_result["numerical_details"]
    uv = detail["uv_at_fermion_threshold"]
    eft = detail["eft1_at_scalar_threshold"]

    heavy_state = fitted_fermion_state(uv, fit_config, mu_f)
    heavy = _heavy_mass_basis_data(
        heavy_state,
        offdiagonal_atol=1e-7, imaginary_atol=1e-9,
        takagi_residual_rtol=1e-8, takagi_residual_atol=1e-7,
    )
    s1sq = finite_real(eft["mS1Sq"], "mS1Sq(mu_S)")
    s2sq = finite_real(eft["mS2Sq"], "mS2Sq(mu_S)")
    if s1sq <= 0 or s2sq <= 0:
        raise ValueError("Scalar mass-squared at mu_S must be positive")
    inputs = FinalC5EvaluationInputs(
        heavy_masses_gev=np.asarray(heavy.masses_gev, dtype=float),
        y1=np.asarray(heavy.y1, dtype=complex),
        y2=np.asarray(heavy.y2, dtype=complex),
        scalar_mass1_gev=math.sqrt(s1sq),
        scalar_mass2_gev=math.sqrt(s2sq),
        grouped_scalar_matching_scale_gev=mu_s,
        lambda_t3=complex(decode(eft["lambdaT3"])),
    )
    # Existing symbolic adapter evaluates the saved hard MSbar + LLSS-induced
    # direct leading-log expression; no one-scale kernel is substituted.
    c5_at_s = np.asarray(
        evaluate_final_weinberg_json(c5_path, inputs.as_config()), dtype=complex
    )
    c5_at_s = 0.5 * (c5_at_s + c5_at_s.T)
    if c5_at_s.shape != (3, 3) or not np.isfinite(c5_at_s).all():
        raise ValueError("Sequential matched C5 is invalid")

    sm_initial = SMWeinbergInitialConditions(
        gY=finite_real(eft["gY"], "gY"),
        g2=finite_real(eft["g2"], "g2"),
        g3=finite_real(eft["g3"], "g3"),
        lambdaH=finite_real(eft["lambdaH"], "lambdaH"),
        yu=finite_matrix(eft["yu"], "yu"),
        yd=finite_matrix(eft["yd"], "yd"),
        ye=finite_matrix(eft["ye"], "ye"),
        K=c5_at_s,
    ).validated()
    numerical = fit_config.get("numerical", {})
    evolved = evolve_sm_weinberg(
        sm_initial, mu_s, mu_low,
        rtol=float(numerical.get("rtol", 1e-8)),
        atol=float(numerical.get("atol", 1e-11)),
    )
    mnu = neutrino_mass_from_c5(evolved.K, vev_gev=vev)
    in_lepton_basis = charged_lepton_mass_basis_matrix(mnu, evolved.ye)
    target = OscillationFitTarget.from_json(target_path)
    fit = evaluate_mass_matrix_fit(in_lepton_basis, target)
    obs = calculate_neutrino_observables(in_lepton_basis, ordering=target.ordering)
    if not all(math.isfinite(float(x)) for x in obs.masses_ev):
        raise ValueError("Non-finite neutrino masses")

    prediction = {name: float(value) for name, value in
                  zip(fit.observable_names, fit.prediction)}
    return {
        "model": MODEL,
        "status": "CompleteConditionalSequentialC5ToSMEFT",
        "saved_renormalisable_running": {
            "UV_nfev": int(detail["uv_nfev"]),
            "EFT1_nfev": int(detail["eft1_nfev"]),
            "rerun_by_this_driver": False,
        },
        "scales_gev": {"mu_F": mu_f, "mu_S": mu_s, "mu_low": mu_low},
        "sequential_matching": {
            "source_kind": "saved_B_minus_3_Matchete_F_then_joint_scalar_pilot",
            "hard_one_loop_MSbar_included": True,
            "direct_LLSS_to_Weinberg_leading_log_included": True,
            "numerical_C5_evaluated_from_saved_hard_plus_running_expression": True,
            "sequential_matching_reexecuted_in_this_driver": False,
            "full_seven_term_Wilson_running_used_in_final_C5": False,
            "heavy_mass_basis_residual": float(max(
                heavy.takagi_residual, heavy.basis_residual)),
            "fermion_masses_gev": [float(x) for x in heavy.masses_gev],
            "scalar_masses_gev": [inputs.scalar_mass1_gev, inputs.scalar_mass2_gev],
            "lambdaT3_at_muS": {
                "real": float(inputs.lambda_t3.real),
                "imag": float(inputs.lambda_t3.imag),
            },
        },
        "smeft": {
            "RGEs_evolved": True,
            "solver_success": bool(evolved.solver_success),
            "solver_nfev": int(evolved.nfev),
            "C5_at_muS_gev_inverse": encode_matrix(c5_at_s),
            "C5_at_muLow_gev_inverse": encode_matrix(evolved.K),
        },
        "neutrinos": {
            "vev_gev": vev,
            "ordering": obs.ordering,
            "masses_ev_m1_m2_m3": [float(x) for x in obs.masses_ev],
            "mass_matrix_gev": encode_matrix(mnu),
            "mass_matrix_charged_lepton_basis_gev": encode_matrix(in_lepton_basis),
            "PMNS_absolute": np.asarray(obs.pmns_abs, dtype=float).tolist(),
            "takagi_residual": float(obs.takagi_residual),
            "oscillation_observables": prediction,
            "oscillation_chi2": float(fit.chi2),
            "target_ordering": str(target.ordering),
        },
    }


def write_outputs(destination: Path, result: dict[str, Any]) -> None:
    destination.mkdir(parents=True, exist_ok=False)
    (destination / "result.json").write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    obs = result["neutrinos"]
    lines = ["quantity,value,unit"]
    for k, v in obs["oscillation_observables"].items():
        lines.append(f"{k},{v:.17g},dimensionless" if k.startswith("sin2") else
                     f"{k},{v:.17g},eV^2")
    for i, mass in enumerate(obs["masses_ev_m1_m2_m3"], 1):
        lines.append(f"m{i},{mass:.17g},eV")
    lines.append(f"chi2,{obs['oscillation_chi2']:.17g},dimensionless")
    (destination / "observables.csv").write_text("\n".join(lines) + "\n", encoding="utf-8")


def parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-id", default="first_v1")
    ap.add_argument("--check-only", action="store_true")
    ap.add_argument("--matching-json", type=Path, default=DEFAULT_C5)
    ap.add_argument("--fit-result", type=Path, default=DEFAULT_FIT / "result.json")
    ap.add_argument("--fit-config", type=Path, default=DEFAULT_FIT / "fitted_config.json")
    ap.add_argument("--target", type=Path, default=None,
                    help="Override the project's canonical oscillation target JSON")
    ap.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT)
    return ap


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]{0,63}", args.run_id):
        raise ValueError("Run ID must be one simple folder name")
    c5, fitted, config, pilot = source_preflight(
        args.matching_json, args.fit_result, args.fit_config
    )
    inputs = (args.matching_json, args.matching_json.with_name("pilot_summary.json"),
              args.fit_result, args.fit_config)
    print(f"B(-3) sequential source and fitted threshold snapshots: OK", flush=True)
    if args.check_only:
        print(f"Hard MSbar matching: {c5['renormalized_hard_threshold']['status']}")
        print("Direct EFT1 Weinberg contribution: " +
              str(c5.get("direct_running_contribution", {}).get("flavor_status")))
        print("Ready to evaluate SMEFT with the fitted B(-3) endpoint")
        return 0
    if args.target is None:
        from Numerical.orchestration.RunSmallLambdaRefit import TARGET
        target = Path(TARGET)
    else:
        target = args.target
    if not target.is_file():
        raise FileNotFoundError("Missing oscillation target: " + str(target))
    target = target.resolve()
    destination = args.output_root / args.run_id
    if destination.exists():
        raise FileExistsError(f"Output already exists; choose another --run-id: {destination}")
    result = run_numerics(args.matching_json, fitted, config, target)
    result["format"] = SCHEMA
    result["symbolic_sequential_pilot"] = {
        "status": pilot["status"],
        "source_Wilson_rge_nonzero_beta_components": pilot.get("rge", {}).get(
            "nonzero_beta_component_count"),
        "scalar_threshold_status": pilot["resumed_scalar_threshold"]["status"],
        "hard_pole_direct_log_relation_checked": pilot.get(
            "pole_rge_consistency", {}).get("direct_log_equals_twice_pole_residue"),
    }
    result["input_files"] = {
        "sequential_C5": str(args.matching_json.resolve()),
        "symbolic_pilot_summary": str(
            args.matching_json.with_name("pilot_summary.json").resolve()),
        "fitted_general_running": str(args.fit_result.resolve()),
        "fitted_UV_config": str(args.fit_config.resolve()),
        "target": str(target),
    }
    result["input_sha256"] = {
        str(path.resolve()): digest(path) for path in (*inputs, target)
    }
    result["scope"] = (
        "B(-3) existing symbolic F-first hard+direct-Wilson one-loop C5, "
        "numerical substitution from saved fitted UV/EFT1 endpoints, "
        "SMEFT one-loop RGE, neutrino masses and observables; not a new all-16 "
        "sequential threshold matching calculation"
    )
    write_outputs(destination, result)
    print(f"Wrote: {destination / 'result.json'}", flush=True)
    print(f"chi2 = {result['neutrinos']['oscillation_chi2']:.8g}", flush=True)
    print("masses (eV) = " + ", ".join(
        f"{m:.9g}" for m in result["neutrinos"]["masses_ev_m1_m2_m3"]
    ), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
