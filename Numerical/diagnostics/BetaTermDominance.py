"""Read-only UV RGBeta one-loop beta-term diagnostics at fixed full-study UV benchmarks.

Run: python -m Numerical.diagnostics.BetaTermDominance
No RG integration, Mathematica, Matchete, or matching is performed.

Terms are classified by which couplings appear in a top-level additive term.
This is a transparent coupling-dependency decomposition, not a diagrammatic
separation. In particular a gauge*Yukawa term is recorded as 'mixed', not
arbitrarily assigned wholly to gauge or Yukawa.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from Numerical.core.RGBetaEvaluator import (
    evaluate_inputform_expression, evaluate_rgbeta_payload, load_rgbeta_payload,
    state_environment,
)
from Numerical.fitting.ScanCLI import build_uv_state_from_config

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG = PROJECT_ROOT / 'configs' / 'generated_model_sets' / 'comparison'
DEFAULT_RAW = PROJECT_ROOT / 'output' / 'full' / 'comparison'
DEFAULT_REPORT = PROJECT_ROOT / 'Reports' / 'output' / 'full' / 'beta_dominance'
SCENARIOS = ('smallY_smallL', 'smallY_largeL', 'largeY_smallL', 'largeY_largeL')
GAUGE = {'gY', 'g2', 'g3'}
YUKAWA = {'yu', 'yd', 'ye', 'y1', 'y2', 'h'}
SCALAR = {'lambdaH', 'lambdaT3', 'lambdaS', 'lambda3', 'lambda4', 'lambda5',
          'lambdaS1', 'lambdaS2', 'lambdaH1', 'lambdaH2', 'lambda12',
          'lambdaH1Adj', 'lambdaH2Adj', 'lambdaS1Adj', 'lambdaS2Adj',
          'lambda12Adj', 'lambda12Cross', 'lambdaHHdagS2S2',
          'lambdaHHdagS1barS1bar', 'lambdaS1bar2S2bar2',
          'lambdaS1barS2S2bar2', 'lambdaS1S1bar2S2bar',
          'lambdaHHdagS1barS2barCross'}
MASS = {'MF', 'mS1Sq', 'mS2Sq', 'mSSq'}
SYMBOL_RE = re.compile(r'[A-Za-z_$][A-Za-z0-9_$]*')
CATEGORIES = ('gauge', 'yukawa', 'scalar', 'mixed', 'mass', 'other')


def additive_terms(expression: str) -> list[str]:
    """Split at top-level +/- only, keeping unary minus and powers intact."""
    expression = expression.strip()
    # Remove only parentheses enclosing the entire expression.
    while expression.startswith('(') and expression.endswith(')'):
        depth = 0
        outer = True
        for i, ch in enumerate(expression):
            if ch == '(': depth += 1
            elif ch == ')': depth -= 1
            if depth == 0 and i != len(expression) - 1:
                outer = False
                break
        if not outer: break
        expression = expression[1:-1].strip()
    terms: list[str] = []
    paren = bracket = 0
    start = 0
    for i, ch in enumerate(expression):
        if ch == '(' : paren += 1
        elif ch == ')': paren -= 1
        elif ch == '[': bracket += 1
        elif ch == ']': bracket -= 1
        elif ch in '+-' and i > start and paren == 0 and bracket == 0:
            previous = expression[i - 1]
            if previous in 'eE' and i >= 2 and expression[i-2].isdigit():
                continue
            if previous in '+-*/^.,([':
                continue
            terms.append(expression[start:i].strip())
            start = i
    if paren != 0 or bracket != 0:
        raise ValueError('Unbalanced RGBeta brackets')
    terms.append(expression[start:].strip())
    return [term for term in terms if term]


def category(term: str, target: str | None = None) -> str:
    names = SYMBOL_RE.findall(term)
    symbols = set(names)
    # A single external coupling factor (e.g. g2^2*y1 in beta_y1)
    # does not make the gauge correction a gauge-Yukawa mixed loop.
    # This is a dependency heuristic; terms not factorisable this way remain mixed.
    if target and names.count(target) == 1:
        symbols.discard(target)
    present = [name for name, family in (
        ('gauge', GAUGE), ('yukawa', YUKAWA), ('scalar', SCALAR),
        ('mass', MASS),
    ) if symbols & family]
    if len(present) == 1: return present[0]
    if len(present) > 1: return 'mixed'
    return 'other'


def norm(value: object) -> float:
    return float(np.linalg.norm(np.asarray(value, dtype=complex).ravel()))


def analyse_beta(expression: str, environment: dict, expected: object,
                 *, target: str | None = None) -> dict:
    buckets: dict[str, object] = {}
    lines = []
    for term in additive_terms(expression):
        group = category(term, target=target)
        value = evaluate_inputform_expression(term, environment)
        buckets[group] = value if group not in buckets else buckets[group] + value
        lines.append({'expression': term, 'category': group,
                      'norm_16pi2_beta': norm(value)})
    total = evaluate_inputform_expression(expression, environment)
    summed = sum((np.asarray(x, dtype=complex) for x in buckets.values()),
                 np.zeros_like(np.asarray(total, dtype=complex)))
    if not np.allclose(summed, total, rtol=1e-9, atol=1e-13):
        raise ValueError('Additive-term reconstruction does not reproduce exported beta')
    if not np.allclose(total, expected, rtol=1e-9, atol=1e-13):
        raise ValueError('Term evaluator disagrees with canonical RGBeta evaluator')
    norms = {k: norm(v) for k, v in buckets.items()}
    denominator = sum(norms.values())
    return {
        'full_norm_16pi2_beta': norm(total),
        'sum_group_norms': denominator,
        'cancellation_ratio': norm(total) / denominator if denominator else None,
        'categories': {k: {'norm_16pi2_beta': norms.get(k, 0.0),
                           'dominance': norms.get(k, 0.0) / denominator if denominator else 0.0}
                       for k in CATEGORIES},
        'terms': lines,
    }


def find_inputs(config_root: Path, raw_root: Path, scenario: str, model_key: str):
    config = config_root / scenario / 'input' / f'{model_key}.json'
    if not config.is_file():
        raise FileNotFoundError(f'Missing full-study benchmark configuration: {config}')
    files = sorted((raw_root / scenario / model_key).glob('T3_*/data/uv_rgbeta_rge.json'))
    if len(files) != 1:
        raise FileNotFoundError(f'Expected exactly one UV RGBeta JSON for {scenario}/{model_key}; found {len(files)}')
    return config, files[0]


def run_analysis(config_root: Path, raw_root: Path, report_root: Path,
                 *, selected_scenarios: tuple[str, ...] = SCENARIOS,
                 selected_models: tuple[str, ...] = (), strict: bool = False) -> dict:
    rows: list[dict] = []
    reports = []
    errors: list[str] = []
    for scenario in selected_scenarios:
        config_dir = config_root / scenario / 'input'
        configs = sorted(config_dir.glob('T3_dS1_*.json'))
        if selected_models:
            configs = [p for p in configs if p.stem in selected_models]
        for config in configs:
            model_key = config.stem
            try:
                config_path, uv_path = find_inputs(config_root, raw_root, scenario, model_key)
                raw = json.loads(config_path.read_text(encoding='utf-8-sig'))
                raw.setdefault('scan', {'bindings': {}})
                state = build_uv_state_from_config(SimpleNamespace(raw=raw), {}).validated()
                uv = load_rgbeta_payload(uv_path)
                evaluated = evaluate_rgbeta_payload(uv, state)
                env = state_environment(state)
                betas = uv['report_betas']
                result = {}
                for name in ('y1', 'y2', 'lambdaT3'):
                    if name not in betas:
                        raise KeyError(f'Missing expected UV beta {name!r}; keys={list(betas)}')
                    result[name] = analyse_beta(str(betas[name]), env, evaluated[name], target=name)
                    for family, entry in result[name]['categories'].items():
                        rows.append({'scenario': scenario, 'model': model_key,
                                     'mu_gev': state.mu_gev, 'beta': name,
                                     'category': family, **entry,
                                     'cancellation_ratio': result[name]['cancellation_ratio']})
                reports.append({'scenario': scenario, 'model': model_key,
                                'config': str(config_path), 'rgbeta': str(uv_path),
                                'mu_gev': state.mu_gev, 'betas': result})
            except (OSError, KeyError, ValueError, TypeError) as exc:
                message = f'{scenario}/{model_key}: {exc}'
                errors.append(message)
                if strict: raise RuntimeError(message) from exc
    if not reports:
        raise RuntimeError('No beta diagnostics generated. ' + '; '.join(errors[:4]))
    report_root.mkdir(parents=True, exist_ok=True)
    (report_root / 'uv_initial_beta_terms.json').write_text(
        json.dumps({'method': 'top-level InputForm additive terms; mixed terms are not assigned to a single sector',
                    'loop_convention': '16*pi^2*dX/dln(mu) = beta_1',
                    'cases': reports, 'errors': errors}, indent=2), encoding='utf-8')
    with (report_root / 'uv_initial_beta_dominance.csv').open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return {'cases': len(reports), 'errors': errors, 'path': str(report_root)}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--configs', type=Path, default=DEFAULT_CONFIG)
    parser.add_argument('--raw', type=Path, default=DEFAULT_RAW)
    parser.add_argument('--output', type=Path, default=DEFAULT_REPORT)
    parser.add_argument('--scenario', action='append', choices=SCENARIOS)
    parser.add_argument('--model', action='append', help='full T3_dS1_... model key')
    parser.add_argument('--strict', action='store_true')
    args = parser.parse_args(argv)
    result = run_analysis(args.configs, args.raw, args.output,
        selected_scenarios=tuple(args.scenario or SCENARIOS),
        selected_models=tuple(args.model or ()), strict=args.strict)
    print(f"UV beta diagnostics: {result['cases']} cases; {len(result['errors'])} errors; {result['path']}")
    for error in result['errors']: print('  WARNING:', error)
    return 1 if result['errors'] else 0


if __name__ == '__main__':
    raise SystemExit(main())
