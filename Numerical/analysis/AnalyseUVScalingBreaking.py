"""Measure T1/T2/T3 breaking in exported UV RGBeta beta functions.

Evaluates the existing symbolic report_betas at the saved fitted UV boundary.
No RGE integration, refit, or symbolic regeneration. The three sector
contributions are *differences under rescaling*, not classifications into
gauge/Yukawa/scalar monomials.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from Numerical.core.RGBetaEvaluator import (
    evaluate_inputform_expression, load_rgbeta_payload, state_environment,
)
from Numerical.fitting.BenchmarkSensitivity import scaled_config_payload
from Numerical.fitting.ScanCLI import build_uv_state_from_config

ROOT = Path(__file__).resolve().parents[2]
MODELS = ("optimised", "shared_a_m2")
DIRECTIONS = ("T1", "T2", "T3")
BETAS = ("lambdaT3", "y1", "y2")


def state_from(raw):
    return build_uv_state_from_config(
        SimpleNamespace(raw={**raw, "scan": {"bindings": {}}}), {}
    )


def terms(state, payload):
    env = state_environment(state)
    missing = set(BETAS) - set(payload.get("report_betas", {}))
    if missing:
        raise ValueError("Missing report_betas: " + repr(sorted(missing)))
    beta = {name: evaluate_inputform_expression(payload["report_betas"][name], env)
            for name in BETAS}
    a, b, lam = state.y1, state.y2, state.lambdaT3
    q = a @ b.T + b @ a.T
    dl = beta["lambdaT3"] * q
    dy1 = lam * (beta["y1"] @ b.T + b @ beta["y1"].T)
    dy2 = lam * (a @ beta["y2"].T + beta["y2"] @ a.T)
    return {"lambda": dl, "y1": dy1, "y2": dy2,
            "total": dl + dy1 + dy2, "P": lam * q}


def norm(x):
    return float(np.linalg.norm(x, ord="fro"))


def run(gammas, modes, output):
    rows, errors = [], []
    base = ROOT / "output/studies/lambda_5e-9_refit/comparison"
    uvroot = ROOT / "output/full/comparison/smallY_smallL"
    for mode in modes:
        for config_path in sorted((base / mode / "models").glob("*/shared_config.json")):
            model = config_path.parent.name
            try:
                raw = json.loads(config_path.read_text(encoding="utf-8-sig"))
                matches = list((uvroot / model).glob("*/data/uv_rgbeta_rge.json"))
                if len(matches) != 1:
                    raise FileNotFoundError(f"Expected one RGBeta export, found {len(matches)}")
                payload = load_rgbeta_payload(matches[0])
                original = state_from(raw)
                original_terms = terms(original, payload)
                pnorm = norm(original_terms["P"])
                if pnorm == 0:
                    raise ValueError("Zero baseline P; relative diagnostics undefined")
                # The UV initial boundary is the only place where this
                # expression-level comparison is performed.
                for direction in DIRECTIONS:
                    for gamma in gammas:
                        varied = state_from(scaled_config_payload(raw, direction, gamma))
                        modified = terms(varied, payload)
                        per = {key: modified[key] - original_terms[key]
                               for key in ("lambda", "y1", "y2", "total")}
                        recombined = per["lambda"] + per["y1"] + per["y2"]
                        if norm(recombined - per["total"]) > 1e-10 * max(norm(per["total"]), 1e-25):
                            raise AssertionError("Product-rule sector sum failed")
                        pchange = norm(modified["P"] - original_terms["P"]) / pnorm
                        if pchange > 1e-9:
                            raise AssertionError(f"P is not invariant at input: {pchange:.3g}")
                        # beta convention: 16 pi^2 dP/dln(mu) = beta_P
                        row = {"mode": mode, "model": model, "direction": direction,
                               "gamma": gamma, "uv_boundary_gev": float(original.mu_gev),
                               "rgbeta_source": str(matches[0].relative_to(ROOT)),
                               "relative_input_P_error": pchange,
                               "baseline_P_norm": pnorm,
                               "baseline_betaP_norm": norm(original_terms["total"]),
                               "delta_betaP_norm": norm(per["total"]),
                               "delta_betaP_over_P": norm(per["total"]) / pnorm,
                               "delta_betaP_over_16pi2_P": norm(per["total"]) / (16 * np.pi**2 * pnorm),
                               "delta_beta_lambda_sector_norm": norm(per["lambda"]),
                               "delta_beta_y1_sector_norm": norm(per["y1"]),
                               "delta_beta_y2_sector_norm": norm(per["y2"]),
                               "sector_cancellation_ratio": norm(per["total"]) /
                               max(sum(norm(per[k]) for k in ("lambda", "y1", "y2")), 1e-300)}
                        rows.append(row)
                print(mode, model, "OK", flush=True)
            except Exception as exc:
                errors.append({"mode": mode, "model": model, "error": f"{type(exc).__name__}: {exc}"})
                print(mode, model, "FAILED", str(exc), flush=True)
    output.mkdir(parents=True, exist_ok=True)
    if rows:
        with (output / "uv_scaling_breaking.csv").open("w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    summary = []
    for mode in modes:
        for direction in DIRECTIONS:
            subset = [r for r in rows if r["mode"] == mode and r["direction"] == direction]
            if subset:
                biggest = max(subset, key=lambda r: r["delta_betaP_over_16pi2_P"])
                summary.append({"mode": mode, "direction": direction,
                                "n_model_gamma_points": len(subset),
                                "largest_normalised_derivative": biggest["delta_betaP_over_16pi2_P"],
                                "model": biggest["model"], "gamma": biggest["gamma"]})
    if summary:
        with (output / "summary.csv").open("w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(summary[0]))
            writer.writeheader()
            writer.writerows(summary)
    (output / "summary.json").write_text(json.dumps({
        "definition": "P=lambdaT3*(y1@y2.T+y2@y1.T)",
        "beta_convention": "16*pi^2*dP/dln(mu)=beta_P",
        "interpretation": "Pointwise symbolic-beta evaluation at UV initial conditions; not integrated sensitivity or gauge/yukawa category decomposition",
        "gammas": list(gammas), "rows": len(rows), "summaries": summary, "errors": errors,
    }, indent=2) + "\n", encoding="utf-8")
    print("Saved:", output, "rows:", len(rows), "errors:", len(errors))
    return 1 if errors else 0


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gammas", nargs="+", type=float, default=[0.9, 1.1])
    ap.add_argument("--mode", choices=["all", *MODELS], default="all")
    ap.add_argument("--output", type=Path, default=ROOT / "output/studies/lambda_5e-9_refit/uv_scaling_breaking")
    args = ap.parse_args()
    if any(not np.isfinite(g) or g <= 0 for g in args.gammas):
        ap.error("Gamma must be positive and finite")
    return run(args.gammas, MODELS if args.mode == "all" else (args.mode,), args.output)


if __name__ == "__main__":
    raise SystemExit(main())
