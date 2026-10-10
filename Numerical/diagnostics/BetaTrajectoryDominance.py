"""Read-only scale-dependent UV beta-sector diagnostic.

Reintegrates only the canonical UV ODE using saved RGBeta JSON and the fixed
benchmark configuration. No Matchete, Mathematica, matching, or full pipeline.

Usage:
 python -m Numerical.diagnostics.BetaTrajectoryDominance --strict
 python -m Numerical.diagnostics.BetaTrajectoryDominance --models all --strict
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from types import SimpleNamespace

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from Numerical.core.RGBetaEvaluator import (
    evaluate_rgbeta_payload, load_rgbeta_payload, state_environment,
)
from Numerical.diagnostics.BetaTermDominance import (
    SCENARIOS, DEFAULT_CONFIG, DEFAULT_RAW, DEFAULT_REPORT,
    analyse_beta, find_inputs,
)
from Numerical.fitting.ScanCLI import build_uv_state_from_config
from Numerical.running.UVRunner import run_uv_segment

REPRESENTATIVE = {
    "A": "T3_dS1_1_dS2_3_dF_2_alpha_m2",
    "B": "T3_dS1_2_dS2_2_dF_1_alpha_m1",
    "C": "T3_dS1_2_dS2_2_dF_3_alpha_m1",
    "D": "T3_dS1_3_dS2_1_dF_2_alpha_m2",
    "E": "T3_dS1_3_dS2_3_dF_2_alpha_m2",
}
CATEGORIES = ("gauge", "yukawa", "scalar", "mixed", "mass", "other")
COUPLINGS = ("y1", "y2", "lambdaT3")
COLOURS = {"gauge": "tab:blue", "yukawa": "tab:orange",
           "scalar": "tab:green", "mixed": "tab:purple",
           "mass": "tab:brown", "other": "tab:gray"}


def analyse_states(result, payload: dict) -> list[dict]:
    """Evaluate full and component betas for every saved *complete* UV state."""
    values = []
    for index in range(result.n_points):
        state = result.state_at_index(index)
        expected = evaluate_rgbeta_payload(payload, state)
        environment = state_environment(state)
        beta_results = {
            coupling: analyse_beta(
                str(payload["report_betas"][coupling]),
                environment, expected[coupling], target=coupling,
            )
            for coupling in COUPLINGS
        }
        values.append({
            "mu_gev": float(result.mu_gev[index]),
            "betas": beta_results,
        })
    return values


def analyse_case(config: Path, rgbeta: Path, *, points: int = 33,
                 rtol: float = 1e-8, atol: float = 1e-11) -> dict:
    raw = json.loads(config.read_text(encoding="utf-8-sig"))
    raw.setdefault("scan", {"bindings": {}})
    state = build_uv_state_from_config(SimpleNamespace(raw=raw), {}).validated()
    scales = raw.get("scales", {})
    threshold = float(scales["mu_fermion_threshold_gev"])
    if not 0 < threshold < state.mu_gev:
        raise ValueError("Invalid UV/fermion-threshold scales")
    if points < 2:
        raise ValueError("points must be at least 2")
    payload = load_rgbeta_payload(rgbeta)
    grid = np.geomspace(state.mu_gev, threshold, points)
    result = run_uv_segment(state, payload, threshold, save_scales_gev=grid,
                            rtol=rtol, atol=atol)
    samples = analyse_states(result, payload)
    return {
        "config": str(config), "rgbeta": str(rgbeta),
        "uv_initial_gev": float(state.mu_gev),
        "fermion_threshold_gev": threshold,
        "solver_nfev": result.nfev, "samples": samples,
    }


def _write_figures(cases: list[dict], figures: Path) -> None:
    figures.mkdir(parents=True, exist_ok=True)
    by_model: dict[str, list[dict]] = {}
    for case in cases:
        by_model.setdefault(case["model"], []).append(case)
    for model, model_cases in by_model.items():
        for coupling in COUPLINGS:
            fig, ax = plt.subplots(figsize=(8.7, 5.3))
            for category in CATEGORIES:
                for case in model_cases:
                    scenario = case["scenario"]
                    x = [p["mu_gev"] for p in case["samples"]]
                    y = [p["betas"][coupling]["categories"][category]["dominance"]
                         for p in case["samples"]]
                    if not any(v > 1e-9 for v in y):
                        continue
                    ax.plot(x, y, color=COLOURS[category],
                            ls="--" if "largeY" in scenario else "-",
                            alpha=0.55 if "largeL" in scenario else 1.0,
                            lw=1.7, label=f"{category} · {scenario}")
            ax.set_xscale("log")
            ax.invert_xaxis() if not ax.xaxis_inverted() else None
            ax.set_ylim(-0.025, 1.025)
            ax.set_xlabel(r"Renormalisation scale $\mu$ [GeV]")
            ax.set_ylabel(r"$D_{x,a}=\|\beta_{x,a}\|/\sum_b\|\beta_{x,b}\|$")
            ax.set_title(f"{model} — {coupling} UV beta-sector fractions")
            ax.grid(True, alpha=.23)
            handles, labels = ax.get_legend_handles_labels()
            if handles:
                ax.legend(handles, labels, fontsize=7, ncol=2,
                          loc="upper center", bbox_to_anchor=(.5, -.15))
            fig.tight_layout()
            fig.savefig(figures / f"{model}_{coupling}_fractions.png",
                        dpi=170, bbox_inches="tight")
            plt.close(fig)

            fig, ax = plt.subplots(figsize=(8.7, 4.5))
            for case in model_cases:
                x = [p["mu_gev"] for p in case["samples"]]
                y = [p["betas"][coupling]["cancellation_ratio"]
                     for p in case["samples"]]
                ax.plot(x, y, lw=1.9, label=case["scenario"])
            ax.set_xscale("log")
            ax.invert_xaxis() if not ax.xaxis_inverted() else None
            ax.set_xlabel(r"Renormalisation scale $\mu$ [GeV]")
            ax.set_ylabel(r"$R_x=\|\beta_x\|/\sum_b\|\beta_{x,b}\|$")
            ax.set_ylim(-0.025, 1.025)
            ax.set_title(f"{model} — {coupling} beta cancellation")
            ax.grid(True, alpha=.23)
            ax.legend(fontsize=8)
            fig.tight_layout()
            fig.savefig(figures / f"{model}_{coupling}_cancellation.png",
                        dpi=170, bbox_inches="tight")
            plt.close(fig)


def run_analysis(config_root: Path = DEFAULT_CONFIG, raw_root: Path = DEFAULT_RAW,
                 output_root: Path = DEFAULT_REPORT / "uv_trajectory",
                 *, models: tuple[str, ...] = tuple(REPRESENTATIVE.values()),
                 scenarios: tuple[str, ...] = SCENARIOS,
                 points: int = 33, strict: bool = False) -> dict:
    cases: list[dict] = []
    errors: list[str] = []
    for model in models:
        for scenario in scenarios:
            try:
                config, rgbeta = find_inputs(config_root, raw_root, scenario, model)
                data = analyse_case(config, rgbeta, points=points)
                cases.append({"model": model, "scenario": scenario, **data})
                print(f"  {model} / {scenario}: {len(data['samples'])} UV points")
            except (OSError, ValueError, TypeError, KeyError, RuntimeError) as exc:
                message = f"{model} / {scenario}: {exc}"
                errors.append(message)
                if strict:
                    raise RuntimeError(message) from exc
    if not cases:
        raise RuntimeError("No valid UV dominance trajectories")
    output_root.mkdir(parents=True, exist_ok=True)
    (output_root / "uv_trajectory_beta_terms.json").write_text(
        json.dumps({"method": "full UV ODE reintegration + grouped RGBeta terms",
                    "beta_convention": "16*pi^2*dX/dln(mu) = beta1",
                    "samples_per_case": points,
                    "models": list(models), "scenarios": list(scenarios),
                    "cases": cases, "errors": errors}, indent=2),
        encoding="utf-8")
    rows = []
    for case in cases:
        for sample in case["samples"]:
            for coupling, beta in sample["betas"].items():
                for category, block in beta["categories"].items():
                    rows.append({
                        "model": case["model"], "scenario": case["scenario"],
                        "mu_gev": sample["mu_gev"], "coupling": coupling,
                        "category": category,
                        "norm_16pi2_beta": block["norm_16pi2_beta"],
                        "dominance": block["dominance"],
                        "full_norm_16pi2_beta": beta["full_norm_16pi2_beta"],
                        "cancellation_ratio": beta["cancellation_ratio"],
                    })
    with (output_root / "uv_trajectory_beta_dominance.csv").open(
        "w", newline="", encoding="utf-8"
    ) as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    _write_figures(cases, output_root / "figures")
    return {"cases": len(cases), "errors": errors, "output": str(output_root)}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--configs", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--raw", type=Path, default=DEFAULT_RAW)
    parser.add_argument("--output", type=Path, default=DEFAULT_REPORT / "uv_trajectory")
    parser.add_argument("--models", choices=("representative", "all"),
                        default="representative")
    parser.add_argument("--model", action="append", help="Explicit full model key; overrides --models")
    parser.add_argument("--scenario", action="append", choices=SCENARIOS)
    parser.add_argument("--points", type=int, default=33)
    parser.add_argument("--strict", action="store_true")
    args = parser.parse_args(argv)
    if args.model:
        models = tuple(args.model)
    elif args.models == "all":
        models = tuple(sorted(p.stem for p in
                 (args.configs / SCENARIOS[0] / "input").glob("T3_dS1_*.json")))
    else:
        models = tuple(REPRESENTATIVE.values())
    result = run_analysis(
        args.configs, args.raw, args.output, models=models,
        scenarios=tuple(args.scenario or SCENARIOS),
        points=args.points, strict=args.strict,
    )
    print(f"Cases: {result['cases']}, errors: {len(result['errors'])}, "
          f"output: {result['output']}")
    return 1 if result["errors"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
