"""One-command fixed-order pipeline with optional exploratory T3-B coupled LL.

This is an additive front end to the official
Numerical.orchestration.PipelineNumericalResults.run_pipeline_numerical_results.
Without --enable-coupled it invokes ONLY the existing production pipeline.
With the flag it additionally writes data/t3b_coupled_ll_exploratory.json.

All three heavy fermions decouple together at mu_F. Their separate masses are
retained in the common-scale matching-log boundary. The new output is a
leading-log-resummed two-operator approximation, NOT a validated NLO result.
No existing C5, figures, or running_diagnostics.json are replaced by the
exploratory result.

Requirement: previously installed T3BExploratoryCoupledBenchmark,
T3BPerturbativeAccuracy, T3BHeavyLogAlignment, T3BC12FlavorRGE modules.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np


MODEL = (2, 2, 1, -1, False)
_TINY = 1e-35


def _relative(a, b):
    a, b = np.asarray(a, dtype=complex), np.asarray(b, dtype=complex)
    return float(np.linalg.norm(a - b) / max(np.linalg.norm(b), _TINY))


def _matrix_json(a):
    a = np.asarray(a, dtype=complex)
    if a.shape != (3, 3) or not np.all(np.isfinite(a)):
        raise ValueError("Expected a finite 3x3 complex matrix")
    return {"real": a.real.tolist(), "imag": a.imag.tolist(),
            "frobenius_norm": float(np.linalg.norm(a))}


def _matrix_read(obj):
    a = np.asarray(obj["real"], dtype=float) + 1j * np.asarray(obj["imag"], dtype=float)
    if a.shape != (3, 3) or not np.all(np.isfinite(a)):
        raise ValueError("Invalid saved C5 matrix")
    return a


def _model_check(payload):
    rep = payload["representation"]
    actual = (int(rep["d_s1"]), int(rep["d_s2"]),
              int(rep["d_f"]), int(rep["alpha"]),
              bool(rep.get("shared_scalar", False)))
    if actual != MODEL:
        raise ValueError("Exploratory coupled LL is only implemented for ordinary "
                         f"T3-B alpha=-1; received {actual}")
    if payload.get("benchmark_search", {}).get("enabled", False):
        raise ValueError("Opt-in coupled LL requires a fixed (non-search) benchmark "
                         "to ensure identical UV inputs in both paths")
    return actual


def _make_record(payload, output_dir):
    from common.RunRecords import RunRecord
    dims = _model_from_payload(payload)
    return RunRecord(name=Path(output_dir).name,
                     alpha=dims[3], d_s1=dims[0], d_s2=dims[1],
                     d_f=dims[2], shared_scalar=dims[4],
                     return_code=0, summary={}, output_dir=Path(output_dir))


def _model_from_payload(payload):
    rep = payload["representation"]
    return (int(rep["d_s1"]), int(rep["d_s2"]), int(rep["d_f"]),
            int(rep["alpha"]), bool(rep.get("shared_scalar", False)))


def _scalar_complex_json(value):
    z = complex(value)
    return [float(z.real), float(z.imag)]


def _intermediate_samples(renormalisable):
    """Use the already defined SM+S1+S2 background, exclusively mu_F -> mu_S."""
    result = []
    for i in range(renormalisable.eft1.n_points):
        s = renormalisable.eft1.state_at_index(i)
        result.append({
            "mu_gev": float(s.mu_gev),
            "g2": float(s.sm.g2),
            "lambdaH": float(s.sm.lambdaH),
            "lambda12": float(s.lambda12),
            "lambda12Adj": float(s.lambda12Adj),
            "lambdaT3": _scalar_complex_json(s.lambdaT3),
            "ye": _matrix_json(s.sm.ye),
            "yu": _matrix_json(s.sm.yu),
            "yd": _matrix_json(s.sm.yd),
        })
    if len(result) < 2:
        raise ValueError("Scalar-only background requires at least two points")
    if not all(result[i]["mu_gev"] > result[i + 1]["mu_gev"]
               for i in range(len(result) - 1)):
        raise ValueError("Intermediate scales must be strictly descending")
    return result


def _neutrino_summary(k, ye, *, vev_gev, ordering):
    from physics.NeutrinoObservables import calculate_neutrino_observables
    from physics.NeutrinoTrajectory import (
        charged_lepton_mass_basis_matrix, neutrino_mass_from_c5,
    )
    from Numerical.diagnostics.RunningDiagnostics import mixing_angles_from_pmns_abs

    m = neutrino_mass_from_c5(k, vev_gev=vev_gev)
    m = charged_lepton_mass_basis_matrix(m, ye)
    o = calculate_neutrino_observables(m, ordering=ordering)
    if not np.isfinite(o.takagi_residual) or o.takagi_residual > 1e-7:
        raise RuntimeError("Bad Takagi diagonalisation")
    s12, s13, s23 = mixing_angles_from_pmns_abs(o.pmns_abs)
    angles = {"sin2_theta12": float(s12), "sin2_theta13": float(s13),
              "sin2_theta23": float(s23),
              "theta12_deg": float(np.degrees(np.arcsin(np.sqrt(s12)))),
              "theta13_deg": float(np.degrees(np.arcsin(np.sqrt(s13)))),
              "theta23_deg": float(np.degrees(np.arcsin(np.sqrt(s23))))}
    return {"ordering": o.ordering, "masses_ev": np.asarray(o.masses_ev).tolist(),
            "delta_m21_sq_ev2": float(o.delta_m21_sq_ev2),
            "delta_m31_sq_ev2": float(o.delta_m31_sq_ev2),
            "delta_m32_sq_ev2": float(o.delta_m32_sq_ev2),
            "pmns_abs": np.asarray(o.pmns_abs).tolist(),
            "mixing_angles": angles, "takagi_residual": float(o.takagi_residual)}


def _low_energy_comparison(renormalisable, fixed_ms, coupled_ms, *, mu_low,
                           vev_gev, ordering, production_low):
    from Numerical.running.ScalarThresholdBoundary import build_sm_weinberg_initial_conditions
    from Numerical.running.SMWeinbergEvolution import evolve_sm_weinberg

    mu_s = float(renormalisable.mu_scalar_threshold_gev)
    if not 0 < mu_low < mu_s:
        raise ValueError("Low scale must be between zero and MS")
    boundary = renormalisable.final_sm_boundary
    results = {}
    low_matrices = {}
    for label, c5 in (("fixed_order", fixed_ms), ("exploratory_coupled", coupled_ms)):
        initial = build_sm_weinberg_initial_conditions(boundary, c5)
        end = evolve_sm_weinberg(initial, mu_s, mu_low)
        results[label] = _neutrino_summary(end.K, end.ye,
                                          vev_gev=vev_gev, ordering=ordering)
        low_matrices[label] = end.K

    official_low = _matrix_read(production_low["c5"])
    baseline_residual = _relative(low_matrices["fixed_order"], official_low)
    if baseline_residual > 2e-5:
        raise ValueError("Rebuilt fixed-order low-scale C5 disagrees with production: "
                         f"{baseline_residual:.4g}")
    f = results["fixed_order"]
    c = results["exploratory_coupled"]
    fm, cm = np.asarray(f["masses_ev"]), np.asarray(c["masses_ev"])
    return {
        "scale_gev": float(mu_low),
        "fixed_order": f,
        "exploratory_coupled": c,
        "relative_c5_difference_at_low_scale": _relative(
            low_matrices["exploratory_coupled"], low_matrices["fixed_order"]),
        "production_low_c5_reconstruction_residual": baseline_residual,
        "differences": {
            "masses_ev": (cm - fm).tolist(),
            "relative_masses": [None if a == 0 else float((b - a) / a)
                                for a, b in zip(fm, cm)],
            "delta_m21_sq_ev2": c["delta_m21_sq_ev2"] - f["delta_m21_sq_ev2"],
            "delta_m31_sq_ev2": c["delta_m31_sq_ev2"] - f["delta_m31_sq_ev2"],
            "theta12_deg": c["mixing_angles"]["theta12_deg"] - f["mixing_angles"]["theta12_deg"],
            "theta13_deg": c["mixing_angles"]["theta13_deg"] - f["mixing_angles"]["theta13_deg"],
            "theta23_deg": c["mixing_angles"]["theta23_deg"] - f["mixing_angles"]["theta23_deg"],
        },
    }


def compute_opt_in_coupled(*, config: dict, uv_rgbeta: Path,
                           eft1_rgbeta: Path, final_c5: Path,
                           production: dict, max_order: int = 4,
                           n_samples: int = 120) -> dict:
    """Attach an independent, clearly non-production coupled result.

    Production's already calculated C5 (hard + direct) is the fixed-order
    anchor. No normalisation fit is used. All heavy generations are matched
    together at mu_F; individual heavy masses enter matching logarithms.
    """
    _model_check(config)
    if not 2 <= n_samples <= 10000:
        raise ValueError("n_samples must be between 2 and 10000")
    if not 1 <= max_order <= 8:
        raise ValueError("max_order must be in [1,8]")
    if production.get("status") != "Success":
        raise ValueError("Official fixed-order numerical stage has not succeeded")

    from Numerical.core.RGBetaEvaluator import load_rgbeta_payload
    from Numerical.fitting.ScanCLI import build_uv_state_from_config
    from Numerical.running.T3Trajectory import run_renormalisable_t3_trajectory
    from Numerical.running.FinalC5TrajectoryAdapter import extract_final_c5_inputs
    from Numerical.diagnostics.FinalC5ContributionDiagnostics import evaluate_final_c5_contributions
    from Numerical.diagnostics.T3BHeavyLogAlignment import matrix_from_json
    from Numerical.diagnostics.T3BExploratoryCoupledBenchmark import (
        calculate, fixed_order_log_from_common_boundary,
    )
    from Numerical.diagnostics.T3BPerturbativeAccuracy import perturbative_accuracy

    scales = config["scales"]
    mu_f, mu_s, mu_low = (float(scales[k]) for k in (
        "mu_fermion_threshold_gev", "mu_scalar_threshold_gev", "mu_low_gev"))
    nconf = config.get("numerical", {})
    rtol = float(nconf.get("rtol", 1e-8))
    atol = float(nconf.get("atol", 1e-11))
    vev = float(nconf.get("vev_gev", 246.22))
    ordering = str(config.get("ordering", "NO")).upper()
    if ordering not in ("NO", "IO"):
        raise ValueError("Coupled opt-in requires explicit NO or IO mass ordering")

    initial = build_uv_state_from_config(
        SimpleNamespace(raw={**config, "scan": config.get("scan", {"bindings": {}})}),
        {},
    ).validated()
    rep = initial.representation
    if (rep.d_s1, rep.d_s2, rep.d_f, rep.alpha, rep.shared_scalar) != MODEL:
        raise ValueError("Loaded UV state is not ordinary T3-B alpha=-1")
    ren = run_renormalisable_t3_trajectory(
        initial, load_rgbeta_payload(Path(uv_rgbeta)),
        load_rgbeta_payload(Path(eft1_rgbeta)),
        mu_fermion_threshold_gev=mu_f, mu_scalar_threshold_gev=mu_s,
        eft1_save_scales_gev=np.geomspace(mu_f, mu_s, n_samples),
        rtol=rtol, atol=atol,
    )
    samples = _intermediate_samples(ren)
    heavy = extract_final_c5_inputs(ren)
    breakdown = evaluate_final_c5_contributions(Path(final_c5), ren)
    fixed = np.asarray(breakdown.combined, dtype=complex)
    direct = np.asarray(breakdown.direct_running, dtype=complex)
    if np.linalg.norm(direct) <= _TINY or np.linalg.norm(fixed) <= _TINY:
        raise ValueError("Expected nonzero C5 direct and fixed-order matrices")
    official = _matrix_read(production["scalar_threshold"]["c5"])
    fixed_residual = _relative(fixed, official)
    if fixed_residual > 2e-6:
        raise ValueError("Matched C5 differs from the official numerical stage: "
                         f"{fixed_residual:.4g}")
    predicted = fixed_order_log_from_common_boundary(
        samples, heavy.y1, heavy.y2, heavy.heavy_masses_gev)
    normalisation_error = _relative(predicted, direct)
    if normalisation_error > 1e-8:
        raise ValueError("Unfitted common-scale direct term does not match production: "
                         f"{normalisation_error:.4g}")

    control = {
        "status": "FixedOrderControlComplete",
        "mu_f_gev": mu_f,
        "mu_s_gev": mu_s,
        "intermediate_running": {"samples": samples},
        "fixed_order_c5_at_ms": _matrix_json(fixed),
    }
    normalization = {
        "status": "UnfittedNormalizationPassed",
        "production_direct_c5": _matrix_json(direct),
        "predicted_direct_c5": _matrix_json(predicted),
    }
    coupled = calculate(control, normalization, heavy.y1, heavy.y2,
                        heavy.heavy_masses_gev,
                        threshold_scheme="common_scale_log_boundary",
                        compare_legacy_threshold=False)
    coupled_ms = matrix_from_json(coupled["exploratory_total_c5_at_ms"])
    low = _low_energy_comparison(
        ren, fixed, coupled_ms,
        mu_low=mu_low, vev_gev=vev, ordering=ordering,
        production_low=production["low_energy"],
    )
    audit = perturbative_accuracy(
        control, coupled, heavy.y1, heavy.y2,
        heavy.heavy_masses_gev, max_order=max_order,
    )
    return {
        "status": "OptionalExploratoryCoupledLLComplete",
        "model": "T3-B alpha=-1",
        "matching_scheme": "simultaneous heavy-fermion decoupling at mu_F",
        "threshold_scheme": "common_scale_log_boundary",
        "mu_f_gev": mu_f, "mu_s_gev": mu_s, "mu_low_gev": mu_low,
        "heavy_masses_gev": np.asarray(heavy.heavy_masses_gev).tolist(),
        "production_c5_overwritten": False,
        "validated_physical_prediction_available": False,
        "perturbative_accuracy": {
            "scope": "one-loop Wilson RGE leading-log resummation; not complete NLO",
            "bookkeeping": audit,
        },
        "matching_checks": {
            "production_fixed_c5_residual": fixed_residual,
            "unfitted_direct_leading_log_residual": normalisation_error,
        },
        "coupled_at_ms": coupled,
        "low_energy": low,
        "notes": [
            "The standard fixed-order production C5 and figures remain authoritative.",
            "Complete higher-order matching and full three-generation operator-space closure are not established.",
            "No sequential heavy-fermion EFTs are used or required by this common-scale scheme.",
        ],
    }


def _load_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def run(*, config: Path, output_dir: Path, enable_coupled: bool = False,
        max_order: int = 4, samples: int = 120) -> dict:
    """Run the OFFICIAL numerical stage, then optionally the coupled sidecar."""
    from Numerical.orchestration.PipelineNumericalResults import run_pipeline_numerical_results

    config, output_dir = Path(config), Path(output_dir)
    payload = _load_json(config)
    record = _make_record(payload, output_dir)
    ok = run_pipeline_numerical_results(record, config)
    if not ok:
        raise RuntimeError("Standard pipeline numerical stage failed: "
                           + str(record.summary.get("PipelineNumericalResultsError", "unknown error")))
    official_path = output_dir / "data" / "running_diagnostics.json"
    report = {
        "status": "IntegratedFixedOrderComplete",
        "production_status": "Success",
        "production_diagnostics_file": str(official_path),
        "coupled_ll_enabled": bool(enable_coupled),
        "production_default_unchanged": True,
        "validated_physical_prediction_available": False if enable_coupled else None,
    }
    if enable_coupled:
        ext = compute_opt_in_coupled(
            config=payload,
            uv_rgbeta=output_dir / "data" / "uv_rgbeta_rge.json",
            eft1_rgbeta=output_dir / "data" / "eft1_rgbeta_rge.json",
            final_c5=output_dir / "data" / "final_weinberg_coefficient.json",
            production=_load_json(official_path),
            max_order=max_order,
            n_samples=samples,
        )
        target = output_dir / "data" / "t3b_coupled_ll_exploratory.json"
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(ext, indent=2) + "\n", encoding="utf-8")
        report["status"] = "IntegratedFixedAndExploratoryCoupledComplete"
        report["exploratory_coupled_file"] = str(target)
        report["relative_c5_change_at_ms"] = ext["coupled_at_ms"]["relative_difference_from_fixed_order"]
        report["relative_c5_change_at_low_scale"] = ext["low_energy"]["relative_c5_difference_at_low_scale"]
        report["reverse_feedback_relative"] = ext["coupled_at_ms"]["relative_reverse_feedback"]
        terms = ext["perturbative_accuracy"]["bookkeeping"]["terms"]
        report["order2_relative_to_order1"] = terms.get("order_2", {}).get("norm_relative_to_order1")
        report["order3_reverse_relative_to_order1"] = terms.get("order_3", {}).get("reverse_contribution_relative_to_order1")

        # Presentation only: refresh HTML companions of already-built analytic
        # GroupFactors reports. No Matchete, RGBeta or group-factor derivation.
        # This opt-in branch does not touch the authoritative production C5.
        from Reports.GroupFactorCoupledHTML import write_group_factor_html_companions
        project_root = Path(__file__).resolve().parents[2]
        html_files = []
        for reports_root in (
            project_root / "Reports" / "output" / "analytical",
            project_root / "Reports" / "output",
        ):
            tex_files = sorted((reports_root / "GroupFactors").glob("GF_*.tex"))
            if tex_files:
                html_files.extend(write_group_factor_html_companions(
                    [], tex_files, report_root=reports_root, explicit_t3b=True,
                    representation_data_dir=output_dir / "data",
                    preserve_existing_html=True,
                ))
        report["group_factor_html_files"] = [str(p) for p in html_files]
        if not html_files:
            report["group_factor_html_note"] = (
                "No saved GroupFactors TeX found; generate analytic GroupFactors "
                "reports first or use the HTML writer on their report root."
            )
    return report


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--config", type=Path, required=True,
                   help="Existing numerical benchmark configuration JSON")
    p.add_argument("--output-dir", type=Path, required=True,
                   help="Existing matching/RGE output directory (contains data/*.json)")
    p.add_argument("--enable-coupled", action="store_true",
                   help="Opt into exploratory coupled C12-C5 LL alongside standard fixed-order results")
    p.add_argument("--max-order", type=int, default=4)
    p.add_argument("--samples", type=int, default=120)
    args = p.parse_args(argv)
    report = run(config=args.config, output_dir=args.output_dir,
                 enable_coupled=args.enable_coupled,
                 max_order=args.max_order, samples=args.samples)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
