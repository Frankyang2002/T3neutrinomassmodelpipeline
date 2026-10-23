"""Term-level UV scaling-breaking decomposition for the saved lambda=5e-9 refit.

Evaluates the repository's RGBeta report_betas without re-running the RGEs.
Categories use symbols occurring in each additive beta-function term. The
individual contributions are complex matrices and are added before norms.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path

import numpy as np

from Numerical.analysis.AnalyseUVScalingBreaking import state_from, norm
from Numerical.core.RGBetaEvaluator import (
    evaluate_inputform_expression, load_rgbeta_payload, state_environment,
)
from Numerical.fitting.BenchmarkSensitivity import scaled_config_payload

ROOT = Path(__file__).resolve().parents[2]
MODES = ("optimised", "shared_a_m2")
DIRECTIONS = ("T1", "T2", "T3")
COUPLINGS = ("lambdaT3", "y1", "y2")
CATEGORIES = ("gauge", "yukawa", "scalar", "mixed", "other")
SYMBOLS = re.compile(r"[A-Za-z_$][A-Za-z0-9_$]*")
YUKAWAS = {"y1", "y2", "yu", "yd", "ye", "h"}
GAUGES = {"gY", "g2", "g3"}


def split_additive(expression: str) -> list[str]:
    """Split only at top-level +/-, preserving nested Matrix/Tr arguments."""
    s = expression.strip()
    pieces: list[str] = []
    start = 0
    depth = 0
    for i, char in enumerate(s):
        if char in "([":
            depth += 1
        elif char in ")]":
            depth -= 1
            if depth < 0:
                raise ValueError("Unbalanced RGBeta brackets")
        elif char in "+-" and i > start and depth == 0:
            # Exclude exponent notation (e.g. 1e-9).
            if s[i - 1] in "eE" and i >= 2 and (s[i - 2].isdigit() or s[i - 2] == "."):
                continue
            pieces.append(s[start:i].strip())
            start = i
    if depth:
        raise ValueError("Unbalanced RGBeta brackets")
    pieces.append(s[start:].strip())
    if not all(pieces):
        raise ValueError("Empty RGBeta additive term")
    return pieces


def classify(term: str) -> str:
    symbols = set(SYMBOLS.findall(term))
    contains_gauge = bool(symbols & GAUGES)
    contains_yuk = bool(symbols & YUKAWAS)
    contains_scalar = any(name.startswith("lambda") for name in symbols)
    # A scalar beta always contains lambdaT3; classify the *additional*
    # factor(s) so lambdaT3*g2^2 is gauge rather than mixed.
    symbols.discard("lambdaT3")
    contains_scalar = any(name.startswith("lambda") for name in symbols)
    tags = [contains_gauge, contains_yuk, contains_scalar]
    if sum(tags) > 1:
        return "mixed"
    if contains_gauge:
        return "gauge"
    if contains_yuk:
        return "yukawa"
    if contains_scalar:
        return "scalar"
    return "other"


def beta_by_category(state, payload):
    env = state_environment(state)
    result = {}
    for coupling in COUPLINGS:
        raw = payload["report_betas"][coupling]
        shape = np.asarray(getattr(state, coupling)).shape
        zero = np.zeros(shape, dtype=complex)
        groups = {category: zero.copy() for category in CATEGORIES}
        terms = []
        for expression in split_additive(raw):
            category = classify(expression)
            value = np.asarray(evaluate_inputform_expression(expression, env), dtype=complex)
            if value.shape != shape:
                raise ValueError(f"{coupling}: term has shape {value.shape}; expected {shape}: {expression}")
            groups[category] += value
            terms.append((expression, category, value))
        whole = np.asarray(evaluate_inputform_expression(raw, env), dtype=complex)
        reconstructed = sum(groups.values())
        scale = max(float(np.linalg.norm(whole)), 1e-25)
        if float(np.linalg.norm(reconstructed - whole)) > 5e-11 * scale:
            raise ValueError(f"{coupling}: additive-term reconstruction mismatch")
        result[coupling] = (groups, terms)
    return result


def product_contribution(state, betas):
    a, b, lam = state.y1, state.y2, state.lambdaT3
    dl = complex(np.asarray(betas["lambdaT3"]).item())
    d1 = np.asarray(betas["y1"], dtype=complex)
    d2 = np.asarray(betas["y2"], dtype=complex)
    return (dl * (a @ b.T + b @ a.T)
            + lam * (d1 @ b.T + b @ d1.T + a @ d2.T + d2 @ a.T))


def write_csv(path: Path, rows: list[dict]):
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def run(gammas, modes, destination):
    root = ROOT / "output/studies/lambda_5e-9_refit/comparison"
    uvroot = ROOT / "output/full/comparison/smallY_smallL"
    detail, terms_output, errors = [], [], []
    for mode in modes:
        configs = sorted((root / mode / "models").glob("*/shared_config.json"))
        if len(configs) != 16:
            raise FileNotFoundError(f"Expected 16 {mode} configurations, found {len(configs)}")
        for config in configs:
            model = config.parent.name
            try:
                raw = json.loads(config.read_text(encoding="utf-8-sig"))
                matches = list((uvroot / model).glob("*/data/uv_rgbeta_rge.json"))
                if len(matches) != 1:
                    raise FileNotFoundError(f"{model}: expected one UV RGBeta export, found {len(matches)}")
                payload = load_rgbeta_payload(matches[0])
                base_state = state_from(raw)
                base_groups = beta_by_category(base_state, payload)
                base_p = base_state.lambdaT3 * (base_state.y1 @ base_state.y2.T + base_state.y2 @ base_state.y1.T)
                p_norm = norm(base_p)
                if p_norm <= 0:
                    raise ValueError("Zero P norm")
                for direction in DIRECTIONS:
                    for gamma in gammas:
                        varied_state = state_from(scaled_config_payload(raw, direction, gamma))
                        varied_groups = beta_by_category(varied_state, payload)
                        p_varied = varied_state.lambdaT3 * (varied_state.y1 @ varied_state.y2.T + varied_state.y2 @ varied_state.y1.T)
                        p_error = norm(p_varied - base_p) / p_norm
                        if p_error > 1e-9:
                            raise ValueError(f"Input P not invariant: {p_error}")
                        by_category = {}
                        by_coupling = {}
                        for cat in CATEGORIES:
                            original_beta = {name: base_groups[name][0][cat] for name in COUPLINGS}
                            modified_beta = {name: varied_groups[name][0][cat] for name in COUPLINGS}
                            by_category[cat] = (product_contribution(varied_state, modified_beta)
                                                - product_contribution(base_state, original_beta))
                        for name in COUPLINGS:
                            original_beta = {k: np.zeros_like(base_groups[k][0]["gauge"]) for k in COUPLINGS}
                            modified_beta = {k: np.zeros_like(varied_groups[k][0]["gauge"]) for k in COUPLINGS}
                            original_beta[name] = sum(base_groups[name][0].values())
                            modified_beta[name] = sum(varied_groups[name][0].values())
                            by_coupling[name] = (product_contribution(varied_state, modified_beta)
                                                 - product_contribution(base_state, original_beta))
                        total = sum(by_category.values())
                        full_before = {k: sum(base_groups[k][0].values()) for k in COUPLINGS}
                        full_after = {k: sum(varied_groups[k][0].values()) for k in COUPLINGS}
                        full_delta = product_contribution(varied_state, full_after) - product_contribution(base_state, full_before)
                        if norm(total - full_delta) > 5e-10 * max(norm(full_delta), 1e-25):
                            raise AssertionError("Category contributions do not reconstruct full betaP difference")
                        if norm(sum(by_coupling.values()) - full_delta) > 5e-10 * max(norm(full_delta), 1e-25):
                            raise AssertionError("Coupling sectors do not reconstruct full betaP difference")
                        denominator = 16 * math.pi**2 * p_norm
                        row = dict(mode=mode, model=model, direction=direction, gamma=float(gamma),
                                   uv_boundary_gev=float(base_state.mu_gev), relative_input_P_error=p_error,
                                   delta_betaP_over_16pi2_P=norm(total)/denominator,
                                   sum_category_norms_over_16pi2_P=sum(norm(v) for v in by_category.values())/denominator,
                                   category_alignment_ratio=norm(total)/max(sum(norm(v) for v in by_category.values()), 1e-300))
                        row.update({f"{cat}_over_16pi2_P": norm(by_category[cat])/denominator for cat in CATEGORIES})
                        row.update({f"{name}_sector_over_16pi2_P": norm(by_coupling[name])/denominator for name in COUPLINGS})
                        detail.append(row)
                        # Terms are source-verifiable strings, grouped by their own symbols.
                        for name in COUPLINGS:
                            for expression, category, baseline_value in base_groups[name][1]:
                                term_group = [x for x in varied_groups[name][1] if x[0] == expression]
                                if len(term_group) != 1:
                                    raise ValueError("Cannot uniquely match additive term across states")
                                zero_original = {k: np.zeros_like(base_groups[k][0]["gauge"]) for k in COUPLINGS}
                                zero_modified = {k: np.zeros_like(varied_groups[k][0]["gauge"]) for k in COUPLINGS}
                                zero_original[name] = baseline_value
                                zero_modified[name] = term_group[0][2]
                                contribution = (product_contribution(varied_state, zero_modified)
                                                - product_contribution(base_state, zero_original))
                                terms_output.append(dict(mode=mode,model=model,direction=direction,gamma=float(gamma),
                                                         beta_coupling=name,category=category,expression=expression,
                                                         term_breaking_norm_over_16pi2_P=norm(contribution)/denominator))
                print(mode, model, "OK", flush=True)
            except Exception as exc:
                errors.append(dict(mode=mode, model=model, error=f"{type(exc).__name__}: {exc}"))
                print(mode, model, "FAILED:", exc, flush=True)
    destination.mkdir(parents=True, exist_ok=True)
    write_csv(destination / "category_breaking.csv", detail)
    write_csv(destination / "symbolic_term_breaking.csv", terms_output)
    summary = []
    for mode in modes:
        for direction in DIRECTIONS:
            subset = [r for r in detail if r["mode"] == mode and r["direction"] == direction]
            if not subset:
                continue
            strongest = max(subset, key=lambda r: r["delta_betaP_over_16pi2_P"])
            summary.append(dict(mode=mode, direction=direction, count=len(subset),
                                largest_breaking=strongest["delta_betaP_over_16pi2_P"],
                                largest_model=strongest["model"],gamma=strongest["gamma"],
                                **{f"{c}_at_largest":strongest[f"{c}_over_16pi2_P"] for c in CATEGORIES}))
    write_csv(destination / "summary.csv", summary)
    (destination / "summary.json").write_text(json.dumps(dict(
        description="Instantaneous UV RGBeta term-level scaling breaking; separate term norms are not additive",
        gamma_values=list(gammas), rows=len(detail), symbolic_term_rows=len(terms_output),
        errors=errors, summary=summary, category_definition={
            "gauge":"term contains gauge coupling and no other nontrivial coupling family",
            "yukawa":"term contains Yukawa coupling and no other nontrivial coupling family",
            "scalar":"term contains scalar quartic(s) besides lambdaT3 only",
            "mixed":"term contains multiple coupling families",
            "other":"no nontrivial gauge, Yukawa or extra scalar factor"},
        caution="The lambdaT3 prefactor is excluded when categorising beta_lambdaT3 terms; norms of matrices do not add linearly."
    ), indent=2) + "\n",encoding="utf-8")
    print("Saved:",destination,"rows:",len(detail),"terms:",len(terms_output),"errors:",len(errors))
    return 1 if errors else 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gammas",type=float,nargs="+",default=[0.9,1.1])
    parser.add_argument("--mode",choices=[*MODES,"all"],default="all")
    args = parser.parse_args()
    if any(not math.isfinite(g) or g <= 0 or g == 1 for g in args.gammas):
        parser.error("Gamma must be finite, positive and not 1")
    modes = MODES if args.mode == "all" else (args.mode,)
    dest = ROOT / "output/studies/lambda_5e-9_refit/uv_scaling_term_decomposition"
    raise SystemExit(run(args.gammas,modes,dest))


if __name__ == "__main__":
    main()
