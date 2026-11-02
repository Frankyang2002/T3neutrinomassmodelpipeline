"""Record and plot an existing general-T3 *conditional one-scale* C5 RGE result.

Project destination: scripts/PlotGeneralT3C5Running.py

This post-processor does NOT perform UV or EFT1 Wilson matching. It reads the
numerical output of Numerical.matching.GeneralT3MultiModelNeutrinos and
reintegrates the existing SM+Weinberg RHS at sampled scales. It can explicitly
retry a previously failed public low-energy run with an absolute tolerance
appropriate to small C5 coefficients. A retry is recorded, never silently
reported as a success of the original run. On failure, diagnostics are saved.
The source's one-scale matching limitations are preserved in the summary.

Example:
 python -X utf8 -m scripts.PlotGeneralT3C5Running \
   --running-json output/general_observable/pilot_B_p1.json \
   --c5-json output/general_observable/B_p1_c5_conditional.json \
   --output-dir output/general_observable/B_p1_running
"""
from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

from Numerical.running.SMWeinbergEvolution import (
    SMWeinbergInitialConditions,
    _beta as sm_weinberg_rhs,
    _pack as sm_weinberg_pack,
    _unpack as sm_weinberg_unpack,
    evolve_sm_weinberg,
)
from physics.NeutrinoMass import numerical_neutrino_mass_matrix


LABELS = ((0, 0, "ee"), (0, 1, "emu"), (0, 2, "etau"),
          (1, 1, "mumu"), (1, 2, "mutau"), (2, 2, "tautau"))


def load_json(path: Path) -> dict:
    obj = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(obj, dict):
        raise ValueError(f"{path}: expected JSON object")
    return obj


def as_complex(value: object) -> complex:
    if isinstance(value, dict) and set(value) == {"real", "imag"}:
        return complex(float(value["real"]), float(value["imag"]))
    if isinstance(value, (int, float, complex)):
        return complex(value)
    raise ValueError(f"Not a complex numeric value: {type(value).__name__}")


def matrix3(value: object, label: str) -> np.ndarray:
    try:
        array = np.asarray([[as_complex(x) for x in row] for row in value],
                           dtype=np.complex128)
    except (ValueError, TypeError) as exc:
        raise ValueError(f"{label}: invalid complex 3x3 matrix") from exc
    if array.shape != (3, 3) or not np.all(np.isfinite(array.real)) or not np.all(np.isfinite(array.imag)):
        raise ValueError(f"{label}: must be a finite 3x3 complex matrix")
    return array


def encoded_matrix(value: np.ndarray) -> dict:
    return {"real": np.asarray(value).real.tolist(),
            "imag": np.asarray(value).imag.tolist()}


def positive(value: object, label: str) -> float:
    result = float(value)
    if not np.isfinite(result) or result <= 0:
        raise ValueError(f"{label} must be finite and positive")
    return result


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def save_plot(path: Path, scales: np.ndarray, curves: dict[str, np.ndarray],
              ylabel: str, title: str, *, logarithmic_y: bool = True) -> None:
    fig, ax = plt.subplots(figsize=(8, 4.8), constrained_layout=True)
    any_positive = False
    for label, curve in curves.items():
        values = np.asarray(curve, dtype=float)
        if values.shape != scales.shape or not np.all(np.isfinite(values)):
            raise ValueError(f"Nonfinite or wrong-size curve: {label}")
        any_positive |= bool(np.any(values > 0))
        ax.plot(scales, values, label=label, linewidth=1.7)
    ax.set_xscale("log")
    if logarithmic_y and any_positive:
        ax.set_yscale("log")
    ax.set_xlim(float(scales[0]), float(scales[-1]))
    ax.set_xlabel(r"Renormalisation scale $\mu$ (GeV)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(fontsize=8, ncol=2)
    ax.grid(alpha=0.25)
    fig.savefig(path, dpi=170)
    plt.close(fig)


def run(args: argparse.Namespace) -> dict:
    if args.points < 3:
        raise ValueError("--points must be at least 3")
    running = load_json(args.running_json)
    matched = load_json(args.c5_json)
    model = str(running.get("model", ""))
    if not model or matched.get("model") != model:
        raise ValueError("Running and C5 results do not have the same model")
    if running.get("status") != "NumericalRenormalisableRunningOnly":
        raise ValueError("Input is not the expected general renormalisable trajectory")
    if matched.get("status") != "ProvisionalDirectMatching_NotFullEFT1":
        raise ValueError("Not the expected conditional one-scale C5 output")
    previous_status = str(matched.get("low_energy_rge_status", "Unknown"))
    previous_errors = matched.get("low_energy_rge_errors", [])
    if not isinstance(previous_errors, list):
        previous_errors = [str(previous_errors)]
    scalar_scale = positive(running["scales_gev"]["S"], "scalar matching scale")
    fermion_scale = positive(running["scales_gev"]["F"], "fermion matching scale")
    # The older evaluator stores the *scalar* scale as observables_scale_gev
    # if its low-energy run failed. It does NOT store the requested mu_low.
    if args.mu_low is None and previous_status == "Success":
        low_scale = positive(matched["observables_scale_gev"], "saved low-energy scale")
    else:
        low_scale = positive(91.1876 if args.mu_low is None else args.mu_low,
                             "requested low-energy scale")
    if not low_scale < scalar_scale:
        raise ValueError("Need 0 < mu_low < mu_scalar")
    eft = running["eft1_at_scalar_threshold"]
    K0 = matrix3(matched["C5_at_scalar_threshold"], "C5 at scalar threshold")
    K_saved = (matrix3(matched["C5_at_low_scale"], "saved low-energy C5")
               if previous_status == "Success" and matched.get("C5_at_low_scale") is not None
               else None)
    if not np.allclose(K0, K0.T, rtol=1e-9, atol=1e-25):
        raise ValueError("Matched C5 is not complex symmetric")
    threshold_diagnostics = {
        "model": model,
        "matching_status": matched.get("status"),
        "original_low_energy_status": previous_status,
        "original_low_energy_errors": previous_errors,
        "mu_scalar_gev": scalar_scale,
        "mu_low_requested_gev": low_scale,
        "matching_is_sequential_eft": False,
        "matching_is_physically_certified": False,
        "c5_norm_at_scalar_gev_inv": float(np.linalg.norm(K0)),
        "SM_couplings_at_scalar": {name: float(eft[name]) for name in
                                    ("gY", "g2", "g3", "lambdaH")},
        "SM_yukawa_singular_max_at_scalar": {
            name: float(np.linalg.svd(matrix3(eft[name], name),
                                      compute_uv=False).max())
            for name in ("ye", "yu", "yd")},
        "solver_rtol": args.rtol,
        "solver_atol": args.atol,
        "solver": "existing SMWeinbergEvolution._beta / scipy DOP853",
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)

    def diagnostic_failure(stage: str, reason: str) -> None:
        payload = {**threshold_diagnostics, "status": "LowEnergyRGEFailed",
                   "failure_stage": stage, "failure_reason": str(reason)}
        (args.output_dir / "diagnostics.json").write_text(
            json.dumps(payload, indent=2, allow_nan=False) + "\n", encoding="utf-8")

    try:
        initial = SMWeinbergInitialConditions(
            gY=float(eft["gY"]), g2=float(eft["g2"]),
            g3=float(eft["g3"]), lambdaH=float(eft["lambdaH"]),
            ye=matrix3(eft["ye"], "ye"),
            yu=matrix3(eft["yu"], "yu"),
            yd=matrix3(eft["yd"], "yd"),
            K=K0,
        ).validated()
    except (RuntimeError, ValueError, FloatingPointError) as exc:
        diagnostic_failure("initial_conditions", str(exc))
        raise RuntimeError("SM+Weinberg initial state invalid; see diagnostics.json: "
                           + str(exc)) from exc

    # Reuse the SAME RHS and state-vector encoding as the existing public RGE
    # solver; the internal helpers are used only to request sampled output.
    t_samples = np.linspace(np.log(scalar_scale), np.log(low_scale), args.points)
    try:
        ode = solve_ivp(sm_weinberg_rhs, (t_samples[0], t_samples[-1]),
                        sm_weinberg_pack(initial), method="DOP853",
                        t_eval=t_samples, rtol=args.rtol, atol=args.atol)
        if not ode.success or ode.y.shape != (76, args.points):
            raise RuntimeError("Sampled SM+Weinberg integration failed: " + str(ode.message))
        # Independent public endpoint run checks the same RHS at matched tolerances.
        public = evolve_sm_weinberg(initial, scalar_scale, low_scale,
                                    rtol=args.rtol, atol=args.atol)
    except (RuntimeError, ValueError, FloatingPointError, OverflowError) as exc:
        diagnostic_failure("ode_integration", str(exc))
        raise RuntimeError("SM+Weinberg running still failed; see diagnostics.json: "
                           + str(exc)) from exc
    mu = np.geomspace(scalar_scale, low_scale, args.points)
    curves: dict[str, list[float]] = {name: [] for _, _, name in LABELS}
    masses: list[list[float]] = []
    rows: list[dict] = []
    for idx, scale in enumerate(mu):
        *_, K = sm_weinberg_unpack(ode.y[:, idx])
        K = 0.5 * (K + K.T)
        if not np.all(np.isfinite(K.real)) or not np.all(np.isfinite(K.imag)):
            raise RuntimeError(f"Nonfinite Weinberg coefficient at index {idx}")
        mass_matrix = numerical_neutrino_mass_matrix(K, vev_gev=args.vev)
        singular_masses_ev = sorted(float(x) * 1e9 for x in np.linalg.svd(mass_matrix, compute_uv=False))
        masses.append(singular_masses_ev)
        row = {"mu_gev": float(scale), "c5_frobenius_gev_inv": float(np.linalg.norm(K)),
               "m_lightest_ev": singular_masses_ev[0], "m_middle_ev": singular_masses_ev[1],
               "m_heaviest_ev": singular_masses_ev[2]}
        for i, j, name in LABELS:
            z = complex(K[i, j])
            row[f"c5_{name}_real_gev_inv"] = float(z.real)
            row[f"c5_{name}_imag_gev_inv"] = float(z.imag)
            row[f"c5_{name}_abs_gev_inv"] = float(abs(z))
            curves[name].append(float(abs(z)))
        rows.append(row)

    sample_end = sm_weinberg_unpack(ode.y[:, -1])[-1]
    # Compare to the public API using both relative and physically scaled norms.
    denom = max(float(np.linalg.norm(public.K)), 1e-40)
    internal_relative = float(np.linalg.norm(sample_end - public.K) / denom)
    saved_relative = (float(np.linalg.norm(sample_end - K_saved) / denom)
                      if K_saved is not None else None)
    if internal_relative > 1e-7:
        raise RuntimeError("Sampled and public Weinberg solvers disagree: "
                           + repr(internal_relative))

    mass_matrix_low = numerical_neutrino_mass_matrix(
        0.5 * (sample_end + sample_end.T), vev_gev=args.vev)
    min_fermion = float(np.linalg.svd(matrix3(running["uv_at_fermion_threshold"]["MF"], "MF"),
                                      compute_uv=False).min())
    scalar_squared = [float(eft["mS1Sq"]), float(eft["mS2Sq"])]
    scalar_masses = [float(np.sqrt(x)) if x > 0 else None for x in scalar_squared]
    issues = []
    if min_fermion <= fermion_scale:
        issues.append("Grouped-F matching scale lies at/above at least one running heavy-fermion mass; hierarchy check failed")
    if any(x is None or x <= scalar_scale for x in scalar_masses):
        issues.append("Scalar threshold lies at/above at least one unbroken-phase scalar mass")
    if saved_relative is not None and saved_relative > 1e-3:
        issues.append("Reintegrated low-energy C5 differs from stored public-run value (>0.1%); original script uses a looser absolute ODE tolerance")
    if previous_status != "Success":
        issues.append("Previous low-energy solver failed; a NEW RGE solution was computed with stated tolerances")

    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    csv_file = out / "c5_running.csv"
    with csv_file.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    save_plot(out / "c5_magnitudes.png", mu,
              {f"|C5_{name}|": np.asarray(values) for name, values in curves.items()},
              r"$|C_{5,\alpha\beta}|$ (GeV$^{-1}$)",
              f"{model} — conditional one-scale C5, SM running")
    m = np.asarray(masses, dtype=float)
    save_plot(out / "neutrino_masses.png", mu,
              {f"mass {i+1} (ascending)": m[:, i] for i in range(3)},
              "Takagi masses (eV)",
              f"{model} — conditional neutrino masses", logarithmic_y=True)
    save_plot(out / "c5_norm.png", mu,
              {"Frobenius norm": np.asarray([r["c5_frobenius_gev_inv"] for r in rows])},
              r"$\|C_5\|_F$ (GeV$^{-1}$)", f"{model} — Weinberg coefficient norm")

    summary = {
        "model": model,
        "status": ("RecomputedAfterEarlierFailure_ConditionalOneScaleC5_NotSequentialEFT"
                   if previous_status != "Success"
                   else "ConditionalOneScaleC5SampledSMWeinbergRGE_NotSequentialEFT"),
        "original_low_energy_rge_status": previous_status,
        "original_low_energy_rge_errors": previous_errors,
        "low_energy_rge_recomputed": True,
        "recomputed_solver_rtol": args.rtol,
        "recomputed_solver_atol": args.atol,
        "input_running_sha256": digest(args.running_json),
        "input_c5_sha256": digest(args.c5_json),
        "matching_method": "FullGeneralUnsplitMatcheteOneScaleCandidate",
        "source_c5_matching": str(matched.get("coefficient", "")),
        "mu_scalar_gev": scalar_scale, "mu_low_gev": low_scale,
        "samples": args.points, "rge_nfev": int(ode.nfev),
        "c5_at_scalar_gev_inv": encoded_matrix(K0),
        "c5_at_low_reintegrated_gev_inv": encoded_matrix(0.5 * (sample_end + sample_end.T)),
        "neutrino_mass_at_low_gev": encoded_matrix(mass_matrix_low),
        "neutrino_takagi_masses_low_ev_ascending": masses[-1],
        "source_neutrino_masses_ev": matched.get("neutrino_masses_ev"),
        "sampled_vs_public_endpoint_rel": internal_relative,
        "sampled_vs_original_saved_endpoint_rel": saved_relative,
        "fermion_min_running_mass_at_F_gev": min_fermion,
        "scalar_unbroken_mass_parameters_at_S_gev": scalar_masses,
        "diagnostic_issues": issues,
        "matchete_rgbeta_coupling_map_certified": False,
        "fermion_threshold_wilson_matching_performed": False,
        "eft1_wilson_running_performed": False,
        "scalar_threshold_sequential_matching_performed": False,
        "physical_general_neutrino_prediction_certified": False,
        "c5_running_csv": "c5_running.csv",
        "plots": ["c5_magnitudes.png", "neutrino_masses.png", "c5_norm.png"],
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n",
                                      encoding="utf-8")
    print(f"{model}: {args.points} C5 samples, low mu={low_scale:g} GeV")
    if previous_status != "Success":
        print("Previous low-energy run FAILED:", previous_errors)
        print("NEW low-energy run SUCCEEDED at atol=", args.atol,
              "; this is separately recorded in summary.json")
    print("Sampled/public endpoint relative difference:", internal_relative)
    print("Sampled/previous stored endpoint relative difference:", saved_relative)
    for issue in issues:
        print("WARNING:", issue)
    print("Recorded:", csv_file)
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--running-json", type=Path, required=True)
    parser.add_argument("--c5-json", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--points", type=int, default=41)
    parser.add_argument("--mu-low", type=float,
                        help="Requested low scale in GeV; defaults to previous successful scale or 91.1876")
    parser.add_argument("--vev", type=float, default=246.22)
    parser.add_argument("--rtol", type=float, default=1e-9)
    parser.add_argument("--atol", type=float, default=1e-21)
    args = parser.parse_args()
    run(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
