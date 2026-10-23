"""Opt-in T3-B coupled leading-log sidecar for an EXISTING production run.

Consumes the independent production fixed-order control and previous coupled
benchmark. Leaves the standard pipeline, its C5 and its fitted outputs untouched.
The sidecar is NOT promoted to the physical C5 output or a standard prediction.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from Numerical.diagnostics.T3BPerturbativeAccuracy import perturbative_accuracy


def _load(path: Path) -> dict:
    value = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object in {path}")
    return value


def check_low_energy(coupled: dict, low: dict) -> dict:
    """Check that an existing low-energy comparison belongs to this coupled run."""
    if low.get("status") != "ExploratoryLowEnergyComparisonComplete":
        raise ValueError("Low energy data has an unexpected status")
    if low.get("threshold_scheme") != "common_scale_log_boundary":
        raise ValueError("Low energy result uses the wrong threshold scheme")
    recorded = float(low["relative_c5_difference_at_ms"])
    actual = float(coupled["relative_difference_from_fixed_order"])
    if not np.isclose(recorded, actual, rtol=2e-9, atol=1e-11):
        raise ValueError("Low energy data is inconsistent with coupled input")
    return {
        "mu_low_gev": float(low["mu_low_gev"]),
        "relative_c5_difference_at_low_scale": float(low["relative_c5_difference_at_low_scale"]),
        "low_energy_status": low["status"],
        "low_energy_is_exploratory": True,
    }


def run(*, config: Path, uv_rgbeta: Path, control: Path,
        normalization: Path, coupled: Path,
        low_energy: Path | None = None, max_order: int = 4) -> dict:
    """Create an optional sidecar with explicit loop-order accounting."""
    from Numerical.core.RGBetaEvaluator import load_rgbeta_payload
    from Numerical.fitting.ScanCLI import build_uv_state_from_config
    from Numerical.running.UVRunner import run_uv_segment
    from Numerical.running.FinalC5TrajectoryAdapter import _majorana_takagi_mass_basis

    cfg, ctl, nrm, cpl = (_load(p) for p in (config, control, normalization, coupled))
    if nrm.get("status") != "UnfittedNormalizationPassed" or float(
            nrm.get("relative_frobenius_residual", 1.0)) > 1e-8:
        raise ValueError("Absolute first-order production matching must pass")
    if cpl.get("validated_physical_prediction_available") is not False:
        raise ValueError("Coupled data must remain explicitly exploratory")
    raw = dict(cfg)
    raw.setdefault("scan", {"bindings": {}})
    initial = build_uv_state_from_config(SimpleNamespace(raw=raw), {}).validated()
    rep = initial.representation
    if (rep.d_s1, rep.d_s2, rep.d_f, rep.alpha, rep.shared_scalar) != (2, 2, 1, -1, False):
        raise ValueError("This optional stage supports ordinary T3-B alpha=-1 only")
    mu_f = float(ctl["mu_f_gev"])
    if not np.isclose(mu_f, float(cfg["scales"]["mu_fermion_threshold_gev"])):
        raise ValueError("Input common fermion matching scales differ")
    uv = run_uv_segment(initial, load_rgbeta_payload(uv_rgbeta), mu_f)
    basis = _majorana_takagi_mass_basis(
        uv.final_state, symmetry_atol=1e-10,
        residual_rtol=1e-9, residual_atol=1e-10)
    if "heavy_masses_gev" in cpl and not np.allclose(
            basis.masses_gev, cpl["heavy_masses_gev"], rtol=1e-9, atol=1e-6):
        raise ValueError("Heavy masses disagree with the saved coupled run")
    audit = perturbative_accuracy(ctl, cpl, basis.y1, basis.y2,
                                  basis.masses_gev, max_order=max_order)
    output = {
        "status": "OptionalExploratoryCoupledLLSidecarComplete",
        "t3_model": "T3-B alpha=-1",
        "simultaneous_fermion_decoupling": True,
        "production_pipeline_modified": False,
        "production_c5_overwritten": False,
        "validated_physical_prediction_available": False,
        "relative_c5_difference_at_ms": cpl["relative_difference_from_fixed_order"],
        "relative_reverse_feedback": cpl["relative_reverse_feedback"],
        "perturbative_bookkeeping": audit,
        "inputs": {"config": str(config), "control": str(control),
                   "normalization": str(normalization), "coupled": str(coupled)},
    }
    if low_energy is not None:
        output["low_energy_comparison"] = check_low_energy(cpl, _load(low_energy))
        output["inputs"]["low_energy"] = str(low_energy)
    return output


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("config", "uv-rgbeta", "control", "normalization", "coupled", "output"):
        parser.add_argument("--" + name, required=True, type=Path)
    parser.add_argument("--low-energy", type=Path, default=None)
    parser.add_argument("--max-order", type=int, default=4)
    args = parser.parse_args(argv)
    result = run(config=args.config, uv_rgbeta=args.uv_rgbeta,
                 control=args.control, normalization=args.normalization,
                 coupled=args.coupled, low_energy=args.low_energy,
                 max_order=args.max_order)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    summary = {"status": result["status"], "relative_c5_difference_at_ms": result["relative_c5_difference_at_ms"],
               "order2_relative_to_order1": result["perturbative_bookkeeping"]["terms"].get("order_2", {}).get("norm_relative_to_order1"),
               "order3_reverse_relative_to_order1": result["perturbative_bookkeeping"]["terms"].get("order_3", {}).get("reverse_contribution_relative_to_order1"),
               "resummed_reconstruction_residual": result["perturbative_bookkeeping"]["resummed_c5_reconstruction_residual"],
               "physical_prediction": False, "output": str(args.output)}
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
