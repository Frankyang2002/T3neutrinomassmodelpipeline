"""Conditional B_-1/C_-1 hybrid benchmark: archived fixed-scale C5 + new RG background.

Runs existing canonical general UV/EFT1 numerical ODEs using the same improved
UV overrides as the 14-model study. Treats a separately archived fixed-scale C5
as the SM+Weinberg boundary at mu_s. This is an *assumption*, not an EFT matching
calculation, and the archived C5 parameters need not match the evolved benchmark.
No existing production code or saved 14-model outputs are modified.

PowerShell:
  python -X utf8 -m scripts.RunArchivedBCWeinberg --canonical-run-id full16_general_v1
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np

from Numerical.orchestration.GeneralRunningBenchmark import (
    load_canonical_general_exports, run as run_general,
)
from Numerical.running.SMWeinbergEvolution import (
    SMWeinbergInitialConditions, evolve_sm_weinberg,
)
from physics.NeutrinoTrajectory import (
    neutrino_mass_from_c5, charged_lepton_mass_basis_matrix,
)
from physics.NeutrinoObservables import calculate_neutrino_observables
from scripts.RunGeneralT3ImprovedRunning import improved_overrides

MODELS = {"B_-1": "B_m1", "C_-1": "C_m1"}
ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ARCHIVE = ROOT / "Reports/validation/bc_fixed_order_archived_reference.json"
DEFAULT_OUT = ROOT / "output/general_numerical_BC_archived_hybrid"


def decode(value):
    if isinstance(value, dict) and set(value) == {"real", "imag"}:
        return complex(float(value["real"]), float(value["imag"]))
    if isinstance(value, list):
        return [decode(x) for x in value]
    return value


def encode(value):
    if isinstance(value, np.ndarray):
        return encode(value.tolist())
    if isinstance(value, (list, tuple)):
        return [encode(x) for x in value]
    if isinstance(value, complex):
        return {"real": float(value.real), "imag": float(value.imag)}
    if isinstance(value, np.generic):
        return encode(value.item())
    if isinstance(value, dict):
        return {str(k): encode(v) for k, v in value.items()}
    return value


def archived_matrix(value):
    arr = np.asarray(value, dtype=float)
    if arr.shape != (3, 3, 2) or not np.isfinite(arr).all():
        raise ValueError("Archived C5 must be a finite 3x3 matrix of [real,imag] pairs")
    result = arr[:, :, 0] + 1j * arr[:, :, 1]
    rel = np.linalg.norm(result - result.T) / max(np.linalg.norm(result), 1e-90)
    if rel > 1e-8:
        raise ValueError(f"Archived C5 is not symmetric: relative residual {rel}")
    return (result + result.T) / 2.0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--model", choices=tuple(MODELS), action="append")
    ap.add_argument("--canonical-run-id", default="full16_general_v1")
    ap.add_argument("--archive", type=Path, default=DEFAULT_ARCHIVE)
    ap.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--mu-uv", type=float, default=1e7)
    ap.add_argument("--mu-f", type=float, default=1e5)
    ap.add_argument("--mu-s", type=float, default=1e3)
    ap.add_argument("--mu-low", type=float, default=91.1876)
    ap.add_argument("--vev", type=float, default=246.22)
    ap.add_argument("--ordering", choices=("AUTO", "NO", "IO"), default="AUTO")
    opts = ap.parse_args(argv)
    if not opts.mu_uv > opts.mu_f > opts.mu_s > opts.mu_low > 0:
        ap.error("Require mu_uv > mu_f > mu_s > mu_low > 0")
    archive = json.loads(opts.archive.read_text(encoding="utf-8-sig"))
    if archive.get("format") != "t3-bc-native-direct-one-loop-fixed-scale-mass-v1":
        raise ValueError("Unexpected archived reference schema")
    if archive.get("full_rg_evolved_prediction_certified") is not False:
        raise ValueError("Unexpected archived reference certification")
    opts.output_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for model in (tuple(dict.fromkeys(opts.model)) if opts.model else tuple(MODELS)):
        row = {"model": model}
        try:
            ref = archive["models"][MODELS[model]]
            if ref.get("model") != model or ref.get("status") != "FIXED_SCALE_ONE_LOOP_COMPUTED_NOT_RG_CERTIFIED":
                raise ValueError("Wrong or changed archived matching status")
            c5 = archived_matrix(ref["C5_GeV_inverse_complex_pairs"])
            uv, eft1, provenance = load_canonical_general_exports(model, opts.canonical_run_id)
            running = run_general(uv, eft1, mu_uv=opts.mu_uv, mu_f=opts.mu_f,
                                  mu_s=opts.mu_s, overrides=improved_overrides(model))
            if running.get("model") != model or running.get("status") != "NumericalRenormalisableRunningOnly":
                raise ValueError("Unexpected general running status/model")
            if not (running.get("uv_solver", {}).get("success") and running.get("eft1_solver", {}).get("success")):
                raise ValueError("General UV/EFT1 solver incomplete")
            sm = running["eft1_at_scalar_threshold"]
            init = SMWeinbergInitialConditions(
                gY=float(sm["gY"]), g2=float(sm["g2"]), g3=float(sm["g3"]),
                lambdaH=float(sm["lambdaH"]),
                ye=np.asarray(decode(sm["ye"]), dtype=complex),
                yu=np.asarray(decode(sm["yu"]), dtype=complex),
                yd=np.asarray(decode(sm["yd"]), dtype=complex), K=c5,
            ).validated()
            evolved = evolve_sm_weinberg(init, opts.mu_s, opts.mu_low)
            mass = charged_lepton_mass_basis_matrix(
                neutrino_mass_from_c5(evolved.K, vev_gev=opts.vev), evolved.ye)
            obs = calculate_neutrino_observables(mass, ordering=opts.ordering)
            result = {
                "model": model,
                "status": "ConditionalArchivedC5HybridToMZ",
                "physical_neutrino_prediction_certified": False,
                "mixed_scalar_matching_performed": False,
                "archived_C5_assumed_at_mu_s": True,
                "archived_C5_parameter_set_matches_running_background": False,
                "archived_C5_parameter_consistency_verified": False,
                "mixed_scalar_bilinear_rge_included": False,
                "higher_dimensional_EFT1_running_performed": False,
                "note": "Fixed-scale archived C5 is imposed at mu_s as an external boundary; no scale/parameter or matching consistency is claimed.",
                "archive_source": str(opts.archive),
                "archived_fixed_scale_masses_ev": ref["masses_eV_ascending"],
                "running_symbolic_provenance": provenance,
                "scales_gev": {"uv":opts.mu_uv,"f":opts.mu_f,"s_assumed":opts.mu_s,"low":opts.mu_low},
                "C5_assumed_at_mu_s_gev_inverse": encode(c5),
                "C5_at_MZ_gev_inverse": encode(evolved.K),
                "mass_matrix_at_MZ_gev": encode(mass),
                "neutrino_masses_ev": encode(obs.masses_ev),
                "delta_m21_sq_ev2": float(obs.delta_m21_sq_ev2),
                "delta_m31_sq_ev2": float(obs.delta_m31_sq_ev2),
                "delta_m32_sq_ev2": float(obs.delta_m32_sq_ev2),
                "pmns_abs": encode(obs.pmns_abs),
                "ordering": obs.ordering,
                "takagi_residual": float(obs.takagi_residual),
                "sm_weinberg_solver_nfev": evolved.nfev,
                "uv_solver": running["uv_solver"],
                "eft1_solver": running["eft1_solver"],
            }
            dest = opts.output_dir / f"{MODELS[model]}_conditional_neutrinos.json"
            dest.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
            row.update(status="SuccessConditional", masses_ev=result["neutrino_masses_ev"],
                       output_json=str(dest))
            print(f"{model}: conditional SM+Weinberg to MZ complete; masses [eV]: {row['masses_ev']}", flush=True)
        except Exception as exc:
            row.update(status="Failed", error=f"{type(exc).__name__}: {exc}")
            print(f"{model}: FAILED: {row['error']}", file=sys.stderr, flush=True)
        rows.append(row)
        (opts.output_dir / "summary.json").write_text(json.dumps({
            "scope": "conditional_hybrid_archived_fixed_scale_C5_plus_improved_SM_background",
            "physical_predictions_certified": False,
            "models": rows,
        }, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    good = sum(x["status"] == "SuccessConditional" for x in rows)
    print(f"Conditional B/C low-energy evolution: {good}/{len(rows)}; summary: {opts.output_dir / 'summary.json'}")
    return 0 if good == len(rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())
