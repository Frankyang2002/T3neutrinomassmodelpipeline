"""UV beta-term dominance for the 5e-9 optimised and shared A(-2) studies.

Runs only the canonical UV RGE (10^7 -> 10^5 GeV) using each study's
existing model-specific RGBeta JSON. No fitting, matching, or symbolic export.

Outputs per-mode CSV, detailed JSON, figures, and a model summary.
The gauge/yukawa/scalar/mixed attribution reuses BetaTermDominance's
coupling-dependency heuristic: it is not a diagram-by-diagram separation.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from Numerical.core.RGBetaEvaluator import (
    evaluate_inputform_expression, evaluate_rgbeta_payload,
    load_rgbeta_payload, state_environment,
)
from Numerical.diagnostics.BetaTermDominance import (
    CATEGORIES, additive_terms, category,
)
from Numerical.fitting.ScanCLI import build_uv_state_from_config
from Numerical.orchestration import RunLambda5e9Comparisons as comparisons
from Numerical.orchestration import RunLambda5e9Refit as study
from Numerical.running.UVRunner import run_uv_segment

COUPLINGS = ("y1", "y2", "lambdaT3")
PALETTE = {"gauge":"tab:blue", "yukawa":"tab:orange",
           "scalar":"tab:green", "mixed":"tab:purple",
           "mass":"tab:brown", "other":"tab:gray"}


def save_json(path: Path, obj: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(obj, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    temporary.replace(path)


def decomposition(expression: str, environment: dict, expected: object,
                  coupling: str, physical_value: object) -> dict:
    """Sum top-level terms by family and verify exact beta reconstruction."""
    total = np.asarray(expected, dtype=complex)
    accum = {k: np.zeros_like(total, dtype=complex) for k in CATEGORIES}
    per_term = []
    for expression_term in additive_terms(expression):
        family = category(expression_term, target=coupling)
        term = np.asarray(evaluate_inputform_expression(expression_term, environment), dtype=complex)
        if term.shape != total.shape:
            raise ValueError(f"{coupling}: term shape {term.shape} != {total.shape}")
        accum[family] += term
        per_term.append({"expression": expression_term, "category": family,
                         "norm_16pi2_beta": float(np.linalg.norm(term))})
    reconstruction = sum(accum.values(), np.zeros_like(total, dtype=complex))
    if not np.allclose(reconstruction, total, rtol=1e-9, atol=1e-13):
        raise ValueError(f"{coupling}: beta terms fail canonical reconstruction")
    norms = {k: float(np.linalg.norm(v)) for k, v in accum.items()}
    sum_norms = sum(norms.values())
    total_norm = float(np.linalg.norm(total))
    result = {"full_norm_16pi2_beta": total_norm,
              "cancellation_ratio": total_norm / sum_norms if sum_norms else None,
              "categories": {}, "terms": per_term}
    x = np.asarray(physical_value, dtype=complex)
    if coupling in ("y1", "y2"):
        if x.shape != total.shape:
            raise ValueError(f"{coupling}: physical Yukawa shape mismatch")
        norm_sq = float(np.vdot(x, x).real)
        # Contribution to d ||Y||_F^2 / d ln(mu) (NOT per log10 scale).
        result["yukawa_norm_sq"] = norm_sq
        result["net_norm_sq_slope"] = float(2 * np.vdot(x, total).real / (16*np.pi**2))
        result["net_norm_slope"] = (result["net_norm_sq_slope"] / (2*np.sqrt(norm_sq))) if norm_sq else None
    for family, value in accum.items():
        entry = {"norm_16pi2_beta": norms[family],
                 "dominance": norms[family]/sum_norms if sum_norms else 0.0}
        if coupling in ("y1", "y2"):
            slope = float(2 * np.vdot(x, value).real / (16*np.pi**2))
            entry["norm_sq_slope"] = slope
            entry["norm_slope"] = slope / (2*np.sqrt(norm_sq)) if norm_sq else None
        result["categories"][family] = entry
    if coupling in ("y1", "y2"):
        reconstructed_slope = sum(result["categories"][k]["norm_sq_slope"] for k in CATEGORIES)
        if not np.isclose(reconstructed_slope, result["net_norm_sq_slope"], atol=1e-12, rtol=1e-9):
            raise ValueError(f"{coupling}: projected slopes fail reconstruction")
    return result


def make_case(model: str, mode: str, config: dict, rgbeta: Path, points: int) -> dict:
    raw = json.loads(json.dumps(config))
    raw["scan"] = {"bindings": {}}
    state = build_uv_state_from_config(SimpleNamespace(raw=raw), {}).validated()
    threshold = float(raw["scales"]["mu_fermion_threshold_gev"])
    if not 0 < threshold < float(state.mu_gev):
        raise ValueError("Invalid UV and fermion-threshold scales")
    payload = load_rgbeta_payload(rgbeta)
    report_betas = payload["report_betas"]
    for coupling in COUPLINGS:
        if coupling not in report_betas:
            raise KeyError(f"{model}: missing {coupling} symbolic beta")
    grid = np.geomspace(float(state.mu_gev), threshold, points)
    result = run_uv_segment(state, payload, threshold, save_scales_gev=grid,
                            rtol=float(raw.get("numerical",{}).get("rtol",1e-8)),
                            atol=float(raw.get("numerical",{}).get("atol",1e-11)))
    samples = []
    for i in range(result.n_points):
        sample_state = result.state_at_index(i)
        expected = evaluate_rgbeta_payload(payload, sample_state)
        env = state_environment(sample_state)
        betas = {coupling: decomposition(str(report_betas[coupling]), env,
                                        expected[coupling], coupling,
                                        getattr(sample_state, coupling))
                 for coupling in COUPLINGS}
        samples.append({"mu_gev": float(result.mu_gev[i]),
                        "y1_norm": float(np.linalg.norm(sample_state.y1)),
                        "y2_norm": float(np.linalg.norm(sample_state.y2)),
                        "lambdaT3_abs": float(abs(sample_state.lambdaT3)),
                        "betas": betas})
    return {"model": model, "mode": mode,
            "uv_rgbeta": str(rgbeta), "uv_rgbeta_sha256": hashlib.sha256(rgbeta.read_bytes()).hexdigest(),
            "uv_initial_gev": float(state.mu_gev), "fermion_threshold_gev": threshold,
            "solver_nfev": int(result.nfev), "samples": samples}


def draw(case: dict, figures: Path) -> None:
    figures.mkdir(parents=True, exist_ok=True)
    samples = case["samples"]
    x = np.array([p["mu_gev"] for p in samples])
    model = case["model"]
    for coupling in COUPLINGS:
        for kind in ("fractions", "cancellation"):
            fig, ax = plt.subplots(figsize=(8.8, 4.8))
            if kind == "fractions":
                for family in CATEGORIES:
                    vals = [p["betas"][coupling]["categories"][family]["dominance"] for p in samples]
                    if max(vals) > 1e-9:
                        ax.plot(x, vals, label=family, color=PALETTE[family])
                ax.set_ylabel(r"$\|\beta_k\| / \sum_j \|\beta_j\|$")
                ax.set_ylim(-0.025, 1.025)
                ax.legend(fontsize=8, ncol=3)
            else:
                vals = [p["betas"][coupling]["cancellation_ratio"] for p in samples]
                ax.plot(x, vals, label="net / sum of magnitudes")
                ax.set_ylabel("Cancellation ratio")
                ax.set_ylim(-0.025, 1.025)
            ax.set_xscale("log");ax.invert_xaxis(); ax.grid(alpha=0.2)
            ax.set_xlabel(r"Renormalisation scale $\mu$ [GeV]")
            ax.set_title(f"{model} | {coupling} | {case['mode']} | {kind}")
            fig.tight_layout()
            fig.savefig(figures / f"{model}_{coupling}_{kind}.png", dpi=165)
            plt.close(fig)
        if coupling not in ("y1", "y2"):
            continue
        fig, ax = plt.subplots(figsize=(8.8, 4.8))
        for family in CATEGORIES:
            vals = [p["betas"][coupling]["categories"][family]["norm_sq_slope"] for p in samples]
            if max(abs(v) for v in vals) > 1e-16:
                ax.plot(x, vals, label=family, color=PALETTE[family])
        net = [p["betas"][coupling]["net_norm_sq_slope"] for p in samples]
        ax.plot(x, net, color="black", ls="--", lw=1.8, label="net slope")
        ax.axhline(0, color="gray", lw=0.8)
        ax.set_xscale("log");ax.invert_xaxis(); ax.grid(alpha=0.2)
        ax.set_xlabel(r"Renormalisation scale $\mu$ [GeV]")
        ax.set_ylabel(r"$d\|Y\|_F^2 / d\ln\mu$")
        ax.set_title(f"{model} | {coupling} norm-squared slope | {case['mode']}")
        ax.legend(fontsize=8, ncol=3)
        fig.tight_layout()
        fig.savefig(figures / f"{model}_{coupling}_norm_slope.png", dpi=165)
        plt.close(fig)
        # Compare the full UV norm trajectory against integrating the
        # reconstructed beta-projected slope; this is a diagnostic identity,
        # not a separate leading-term approximation.
        t = np.log(x)
        y = np.array([p[coupling + "_norm"]**2 for p in samples])
        reconstructed = np.array(net)
        approximated = np.empty_like(y); approximated[0] = y[0]
        for j in range(1, len(y)):
            approximated[j] = approximated[j-1] + 0.5 * (reconstructed[j-1] + reconstructed[j]) * (t[j]-t[j-1])
        fig, ax = plt.subplots(figsize=(8.8,4.8))
        ax.plot(x, np.sqrt(np.maximum(0,y)), label="full UV ODE", lw=2)
        ax.plot(x, np.sqrt(np.maximum(0,approximated)), ls="--", label="integrated symbolic beta slope")
        ax.set_xscale("log");ax.invert_xaxis();ax.grid(alpha=.2)
        ax.set_xlabel(r"Renormalisation scale $\mu$ [GeV]")
        ax.set_ylabel(r"$\|Y\|_F$")
        ax.set_title(f"{model} | {coupling} slope reconstruction | {case['mode']}")
        ax.legend()
        fig.tight_layout()
        fig.savefig(figures / f"{model}_{coupling}_slope_check.png", dpi=165)
        plt.close(fig)


def run_mode(mode: str, models: list, points: int, out: Path, resume: bool, figures: bool) -> None:
    out.mkdir(parents=True, exist_ok=True)
    all_cases = []
    rows = []
    summaries = []
    for n, (row, config, artifacts) in enumerate(models, 1):
        name = row["model"]
        cache = out / "models" / name / "beta_dominance.json"
        try:
            payload = {"mode":mode,"model":name,"config":config,"points":points,
                       "uv_rgbeta_sha256":hashlib.sha256(Path(artifacts[0]).read_bytes()).hexdigest()}
            digest = hashlib.sha256(json.dumps(payload,sort_keys=True).encode()).hexdigest()
            if resume and cache.is_file():
                case = json.loads(cache.read_text(encoding="utf-8"))
                if case.get("input_sha256") == digest and len(case.get("samples", [])) == points:
                    print(f"[{n}/{len(models)}] REUSED {mode} {name}",flush=True)
                else:
                    case = None
            else:
                case = None
            if case is None:
                print(f"[{n}/{len(models)}] RUN {mode} {name}",flush=True)
                case = make_case(name, mode, config, Path(artifacts[0]), points)
                case["input_sha256"] = digest
                save_json(cache,case)
            if figures:
                draw(case,out / "figures")
            all_cases.append({"model":name,"path":str(cache)})
            for sample in case["samples"]:
                for coupling in COUPLINGS:
                    beta = sample["betas"][coupling]
                    for family in CATEGORIES:
                        c = beta["categories"][family]
                        rows.append({"mode":mode,"model":name,"mu_gev":sample["mu_gev"],
                                     "coupling":coupling,"category":family,
                                     "norm_16pi2_beta":c["norm_16pi2_beta"],
                                     "dominance":c["dominance"],
                                     "cancellation_ratio":beta["cancellation_ratio"],
                                     "norm_sq_slope":c.get("norm_sq_slope", ""),
                                     "net_norm_sq_slope":beta.get("net_norm_sq_slope","")})
            summary = {"mode":mode,"model":name}
            for coupling in COUPLINGS:
                at_uv = case["samples"][0]["betas"][coupling]
                family = max(CATEGORIES,key=lambda k:at_uv["categories"][k]["dominance"])
                summary[coupling+"_dominant_at_uv"] = family
                summary[coupling+"_dominance_fraction_at_uv"] = at_uv["categories"][family]["dominance"]
                summary[coupling+"_cancellation_at_uv"] = at_uv["cancellation_ratio"]
            summaries.append(summary)
        except Exception as exc:
            print(f"[{n}/{len(models)}] FAILED {mode} {name}: {type(exc).__name__}: {exc}",flush=True)
            raise
    if rows:
        with (out/"beta_dominance_trajectory.csv").open("w",newline="",encoding="utf-8") as f:
            writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    if summaries:
        with (out/"beta_dominance_summary.csv").open("w",newline="",encoding="utf-8") as f:
            writer=csv.DictWriter(f,fieldnames=list(summaries[0]));writer.writeheader();writer.writerows(summaries)
    save_json(out/"manifest.json",{"mode":mode,"points":points,"cases":all_cases,
                                  "convention":"16*pi^2*dX/dln(mu)=beta1",
                                  "caveat":"Dependency-based beta term groups; only UV region; no matching included"})
    print(f"{mode}: complete {len(all_cases)} models -> {out}",flush=True)


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("optimised","shared","both"),default="both")
    parser.add_argument("--points",type=int,default=33)
    parser.add_argument("--resume",action="store_true")
    parser.add_argument("--no-figures",action="store_true")
    parser.add_argument("--artifact-map",type=Path)
    parser.add_argument("--output",type=Path,default=study.OUT/"beta_dominance")
    args=parser.parse_args(argv)
    if args.points<2:parser.error("--points must be >= 2")
    for mode in (("optimised","shared") if args.mode=="both" else (args.mode,)):
        _, models=comparisons.load_inputs(args.artifact_map,mode)
        out=args.output/("shared_a_m2" if mode=="shared" else "optimised")
        run_mode(mode,models,args.points,out,args.resume,not args.no_figures)


if __name__=="__main__":
    main()
