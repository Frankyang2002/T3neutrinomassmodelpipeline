"""Standalone T3-B end-to-end CONTROL run using existing full-flavour pipeline.

This runs the renormalisable UV/EFT1 trajectory and existing fixed-order C5
matching without changing production outputs. It does NOT implement the
missing full-flavour C12<->C5 evolution or replace threshold matching.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from Numerical.core.RGBetaEvaluator import load_rgbeta_payload
from Numerical.fitting.ScanCLI import build_uv_state_from_config
from Numerical.running.T3Trajectory import run_renormalisable_t3_trajectory
from Numerical.running.FinalC5TrajectoryAdapter import FinalC5TrajectoryEvaluator
from Numerical.diagnostics.FinalC5ContributionDiagnostics import evaluate_final_c5_contributions


def _matrix(value):
    a = np.asarray(value, dtype=complex)
    return {"real": a.real.tolist(), "imag": a.imag.tolist(),
            "shape": list(a.shape), "frobenius_norm": float(np.linalg.norm(a))}


def run(*, config: Path, uv_rgbeta: Path, eft1_rgbeta: Path,
        final_c5: Path, n_samples: int = 40):
    """Return an actual full-flavour fixed-order control and EFT1 trajectory."""
    payload = json.loads(Path(config).read_text(encoding="utf-8"))
    raw = dict(payload)
    raw.setdefault("scan", {"bindings": {}})
    initial = build_uv_state_from_config(SimpleNamespace(raw=raw), {}).validated()
    rep = initial.representation
    expected = (2, 2, 1, -1, False)
    actual = (rep.d_s1, rep.d_s2, rep.d_f, rep.alpha, rep.shared_scalar)
    if actual != expected:
        raise ValueError(f"Only T3-B alpha=-1 supported: got {actual}")
    if not (2 <= n_samples <= 10000):
        raise ValueError("n_samples must lie between 2 and 10000")
    scales = payload["scales"]
    mu_f = float(scales["mu_fermion_threshold_gev"])
    mu_s = float(scales["mu_scalar_threshold_gev"])
    numerical = payload.get("numerical", {})
    rtol = float(numerical.get("rtol", 1e-8))
    atol = float(numerical.get("atol", 1e-11))
    ren = run_renormalisable_t3_trajectory(
        initial, load_rgbeta_payload(uv_rgbeta), load_rgbeta_payload(eft1_rgbeta),
        mu_fermion_threshold_gev=mu_f, mu_scalar_threshold_gev=mu_s,
        eft1_save_scales_gev=np.geomspace(mu_f, mu_s, n_samples),
        rtol=rtol, atol=atol)
    rows = []
    for i in range(ren.eft1.n_points):
        state = ren.eft1.state_at_index(i)
        rows.append({"mu_gev": float(state.mu_gev),
                     "g2": float(state.sm.g2),
                     "lambdaH": float(state.sm.lambdaH),
                     "lambda12": float(state.lambda12),
                     "lambda12Adj": float(state.lambda12Adj),
                     "lambdaT3": [float(complex(state.lambdaT3).real),
                                  float(complex(state.lambdaT3).imag)],
                     "yu": _matrix(state.sm.yu),
                     "yd": _matrix(state.sm.yd),
                     "ye": _matrix(state.sm.ye)})
    c5 = np.asarray(FinalC5TrajectoryEvaluator(final_c5)({}, ren), dtype=complex)
    if c5.shape != (3, 3):
        raise ValueError(f"Expected full flavour C5 (3,3), got {c5.shape}")
    return {"status": "FixedOrderControlComplete",
            "model": "T3-B alpha=-1",
            "scope": "renormalisable UV/EFT1 and existing full-flavour fixed-order C5",
            "mu_f_gev": mu_f, "mu_s_gev": mu_s,
            "intermediate_running": {"n_points": ren.eft1.n_points,
                                     "solver_success": ren.eft1.solver_success,
                                     "samples": rows},
            "fixed_order_c5_at_ms": _matrix(c5),
            "full_flavour_c12_c5_coupled_rge": "NOT_IMPLEMENTED",
            "scalar_threshold_resummed_matching": "NOT_IMPLEMENTED",
            "neutrino_mass_comparison": "NOT_AVAILABLE_UNTIL_COUPLED_C5_VERIFIED",
            "note": "Do not feed the one-generation coupled solver with 3x3 matrices."
            }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--config", type=Path, required=True)
    p.add_argument("--uv-rgbeta", type=Path, required=True)
    p.add_argument("--eft1-rgbeta", type=Path, required=True)
    p.add_argument("--final-c5", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--samples", type=int, default=40)
    args=p.parse_args()
    result=run(config=args.config, uv_rgbeta=args.uv_rgbeta,
               eft1_rgbeta=args.eft1_rgbeta, final_c5=args.final_c5,
               n_samples=args.samples)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2), encoding="utf-8")
    print(json.dumps({"status":result["status"],
        "intermediate_points":result["intermediate_running"]["n_points"],
        "fixed_order_c5_frobenius_norm": result["fixed_order_c5_at_ms"]["frobenius_norm"],
        "output":str(args.output)},indent=2))

if __name__=="__main__":main()
