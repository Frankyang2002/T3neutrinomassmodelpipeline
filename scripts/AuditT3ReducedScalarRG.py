"""Batch 5: symbolic one-loop RG generation orders in the reduced T3 scalar sector.

Reads the saved 16-model candidate RGBeta UV/EFT1 beta exports and Stage29
operator registry. All 50 additional quartic slots and 8 allowed quadratic
slots are set zero at the boundary. The full beta expression is reduced as a
formal polynomial, including cancellations of identical monomials.

Computes whether an absent coupling is sourced immediately at one RG loop,
then determines the minimum number of *iterated one-loop beta insertions*
needed to generate it. Does not do an RG integration, calculate C5, or certify
Matchete<->RGBeta normalization. Saved data are never overwritten.

    python -m scripts.AuditT3ReducedScalarRG
"""
from __future__ import annotations

import argparse
import ast
import csv
from fractions import Fraction
from hashlib import sha256
import json
import math
from pathlib import Path
from typing import Any

from scripts.AuditT3QuadraticRGClosure import (
    MODELS, EXPECTED, MASS_SYMBOLS, parse_beta, tag,
)

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'output' / 'general_t3_quadratic_eft1'
REGISTRY = ROOT / 'RGE' / 'general' / 'generated' / 'general_rgbeta_quartic_registry.json'
MANIFEST = ROOT / 'docs' / 'general_scalar_trial_operator_manifest.json'
MATCHING = ROOT / 'output' / 'general_matching_symbolic_report.json'
OUT = ROOT / 'output' / 'scalar_relevance'
ONE_LOOP = 16 * math.pi**2

# Monomials are sorted factor tuples: ('lambdaT3', 'lambdaT3', '__TR_0').
# The Tr[...] token represents a formal, otherwise arbitrary flavour trace.
Poly = dict[tuple[str, ...], Fraction]


def fail(msg: str) -> None:
    raise ValueError('Batch 5 fail-closed: ' + msg)


def read(path: Path) -> dict[str, Any]:
    if not path.is_file():
        fail('required input missing: ' + str(path))
    value = json.loads(path.read_text(encoding='utf-8-sig'))
    if not isinstance(value, dict):
        fail('invalid JSON object: ' + str(path))
    return value


def fingerprint(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _clean(p: Poly) -> Poly:
    return {m: c for m, c in p.items() if c}


def _sum(a: Poly, b: Poly, sign: int = 1) -> Poly:
    out = dict(a)
    for m, c in b.items():
        out[m] = out.get(m, Fraction(0)) + sign*c
    return _clean(out)


def _times(a: Poly, b: Poly) -> Poly:
    if len(a)*len(b) > 100000:
        fail('polynomial expansion exceeded complexity cap')
    out: Poly = {}
    for k, x in a.items():
        for l, y in b.items():
            m = tuple(sorted(k+l))
            out[m] = out.get(m, Fraction(0)) + x*y
    return _clean(out)


def polynomial(node: ast.AST) -> Poly:
    """Exact rational polynomial from the restricted RGBeta InputForm AST."""
    if isinstance(node, ast.Constant) and type(node.value) in (int, float):
        return _clean({(): Fraction(str(node.value))})
    if isinstance(node, ast.Name):
        return {(node.id,): Fraction(1)}
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.USub, ast.UAdd)):
        p = polynomial(node.operand)
        return {m: -c for m, c in p.items()} if isinstance(node.op, ast.USub) else p
    if not isinstance(node, ast.BinOp):
        fail('non-polynomial beta expression: ' + type(node).__name__)
    left, right = polynomial(node.left), polynomial(node.right)
    if isinstance(node.op, ast.Add):
        return _sum(left, right)
    if isinstance(node.op, ast.Sub):
        return _sum(left, right, -1)
    if isinstance(node.op, ast.Mult):
        return _times(left, right)
    if isinstance(node.op, ast.Div):
        if set(right) != {()} or right[()] == 0:
            fail('beta expression divides by a coupling or zero')
        return _clean({m: c/right[()] for m, c in left.items()})
    if isinstance(node.op, ast.Pow):
        if set(right) != {()} or right[()].denominator != 1:
            fail('unsupported non-integer coupling exponent')
        n = int(right[()])
        if not 0 <= n <= 5:
            fail('unsupported coupling exponent')
        result: Poly = {(): Fraction(1)}
        for _ in range(n):
            result = _times(result, left)
        return result
    fail('unsupported beta expression operation')


def expand_beta(expression: str) -> Poly:
    if not isinstance(expression, str) or not expression.strip() or any(
            x in expression for x in ('BetaTerm[', '$Failed', '$Aborted', 'Missing[')):
        fail('unresolved or absent RGBeta beta expression')
    return polynomial(parse_beta(expression))


def bare(symbol: str) -> str:
    return symbol[4:] if symbol.startswith('BAR_') else symbol


def source_terms(p: Poly, omitted: set[str]) -> Poly:
    return {m: c for m, c in p.items() if not any(bare(x) in omitted for x in m)}


def minimum_orders(betas: dict[str, Poly], omitted: set[str]) -> dict[str, int | None]:
    """Formal iterated one-loop insertion count with zero omitted boundaries.

    Retained symbols/traces have zero insertion cost. Each omitted symbol has
    initially infinite cost. A beta monomial inherits the sum of costs of its
    omitted factors, plus one for the RG derivative. Repeated relaxation
    handles chains (quartic -> quartic -> quadratic) and nonlinear terms.
    """
    if set(betas) != omitted:
        fail('beta inventory does not match omitted-coupling inventory')
    unreachable = 10**8
    orders: dict[str, int] = {name: unreachable for name in omitted}
    for _ in range(len(omitted) + 1):
        updated = False
        for name, p in betas.items():
            choices = []
            for monomial in p:
                cost = 1
                for factor in monomial:
                    f = bare(factor)
                    if f in omitted:
                        cost += orders[f]
                if cost < unreachable:
                    choices.append(cost)
            if choices and min(choices) < orders[name]:
                orders[name] = min(choices)
                updated = True
        if not updated:
            break
    return {name: (None if cost >= unreachable else cost) for name, cost in orders.items()}


def _expected_direct_symbols(model: str, registry_model: dict) -> set[str]:
    """Map the four Matchete Wick-eligible slots to RGBeta via counts.

    This maps only field structures, NOT CG normalization or physical C5.
    """
    wanted = {(0,2,0,0,2,0), (0,2,0,2,0,0)} if model in ('B_-1','C_-1') else set()
    pairs = []
    for item in registry_model['new_rgbeta_symbols']:
        pairs.append((item['name'], tuple(item['counts'])))
    for item in registry_model['preexisting_rgbeta_symbols']:
        pairs.append((item['existing_name'], tuple(item['counts'])))
    got = {n for n,c in pairs if c in wanted}
    if len(got) != len(wanted):
        fail(model + ': four directly matched B/C quartic slots not represented in RGBeta registry')
    return got


def _validate_saved_matching(model: str, report: dict) -> None:
    tag_name = tag(model)
    matches = [r for r in report.get('results', []) if r.get('model') == tag_name]
    if len(matches) != 1:
        fail('saved symbolic matching missing/duplicate model: '+model)
    detected = set(matches[0].get('variants', {}).get('extended', {}).get('extra_couplings_observed', []))
    exp = ({'lambdaT3ExtraQ020020I1','lambdaT3ExtraQ020200I1'}
           if model in ('B_-1','C_-1') else set())
    if detected != exp:
        fail(model + ': saved C5 coupling dependency differs from 4-slot audit')
    if matches[0].get('extended_log_operator_count_matches_manifest') is not True:
        fail(model + ': incomplete saved extended matching interaction inventory')


def check_export(model: str, stage: str, raw: dict, spec: dict,
                 has_quadratic: bool) -> tuple[dict[str, Poly], dict[str, Any]]:
    kind = EXPECTED.get(model, 'none')
    meta = raw.get('metadata', {})
    if (raw.get('format') != 'general-t3-quadratic-rgbeta-experimental-v1'
            or raw.get('status') != 'SuccessCandidateMassBeta'
            or not isinstance(meta, dict)
            or raw.get('physical_mapping_certified') is not False):
        fail(model+'/'+stage+': saved RGBeta beta certificate invalid')
    expected = {
        'GeneralScalarModel': model, 'GeneralScalarStage': stage,
        'GeneralScalarBasisStatus': 'CandidateNativeRGBetaUnmappedToMatchete',
        'QuadraticKind': kind, 'QuadraticMassSymbol': MASS_SYMBOLS[kind],
        'PhysicalMatcheteRGBetaMassBasisCertified': False,
        'ReportBetaConvention': '16*pi^2*dX/dln(mu); gauge BetaTerm divided by 2*g',
    }
    if any(meta.get(k) != v for k,v in expected.items()):
        fail(model+'/'+stage+': saved beta metadata inconsistent')
    extra = {x['name'] for x in spec['new_rgbeta_symbols']}
    reused = {x['existing_name'] for x in spec['preexisting_rgbeta_symbols']}
    quartics = extra | reused
    if (set(meta.get('GeneralScalarNewBetaKeys', [])) != extra
            or meta.get('GeneralScalarAdditionalQuartics') != spec['total_additional_quartics']
            or len(quartics) != spec['total_additional_quartics']
            or int(has_quadratic) != int(kind != 'none')):
        fail(model+'/'+stage+': model-specific quartic or quadratic inventory mismatch')
    omitted = quartics | ({MASS_SYMBOLS[kind]} if kind != 'none' else set())
    report = raw.get('report_betas')
    if not isinstance(report, dict) or not omitted.issubset(report):
        # The no-additional-operator models have an empty omitted set. We still
        # require the complete nonempty beta inventory from the saved solver.
        if not isinstance(report, dict) or not report:
            fail(model+'/'+stage+': no beta functions')
        if not omitted.issubset(report):
            fail(model+'/'+stage+': missing omitted beta functions')
    source = {name: expand_beta(report[name]) for name in omitted}
    for name, poly in source.items():
        for mon in poly:
            unknown = {bare(f) for f in mon if (bare(f).startswith('lambdaT3GeneralQ')
                        or bare(f).startswith('muT3General')) and bare(f) not in omitted}
            if unknown:
                fail(model+'/'+stage+'/'+name+': unexpected absent extra symbol '+repr(sorted(unknown)))
    mass_name = MASS_SYMBOLS[kind]
    if kind == 'none':
        if any(k.startswith('muT3General') for k in report):
            fail(model+'/'+stage+': unexpected mass mixing beta')
    else:
        if source_terms(source[mass_name], omitted):
            fail(model+'/'+stage+': quadratic beta directly sourced with ALL extras zero')
    return source, {'quartics':sorted(quartics), 'omitted':sorted(omitted),
                    'direct': sorted(_expected_direct_symbols(model, spec)),
                    'mass': mass_name}


def build(registry_path: Path = REGISTRY, manifest_path: Path = MANIFEST,
          matching_path: Path = MATCHING, data_dir: Path = DATA) -> dict:
    registry, manifest, matching = read(registry_path), read(manifest_path), read(matching_path)
    if (registry.get('totals', {}).get('matchete_extra_quartics') != 50
            or registry.get('totals', {}).get('quadratic_couplings_blocked') != 8
            or manifest.get('model_count') != 16
            or set(registry.get('models', {})) != set(MODELS)
            or set(manifest.get('models', {})) != set(MODELS)):
        fail('16-model general scalar registry or manifest incomplete')
    hashes = {str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p): fingerprint(p)
              for p in (registry_path, manifest_path, matching_path)}
    models = {}
    total_quartics, total_masses, uv_quartic_one, uv_mass_order2 = 0,0,0,0
    for model in MODELS:
        spec = registry['models'][model]
        manifest_row = manifest['models'][model]
        _validate_saved_matching(model, matching)
        if (spec['total_additional_quartics'] != manifest_row['extra_quartic_complex_couplings']
                or spec['quadratics_not_implemented'] != manifest_row['quadratic_complex_couplings_not_enabled']):
            fail(model + ': operator manifest and registry disagree')
        total_quartics += spec['total_additional_quartics']
        has_mass = spec['quadratics_not_implemented'] == 1
        total_masses += int(has_mass)
        stages = {}
        for stage in ('UV', 'EFT1'):
            path = data_dir / f'{tag(model)}_{stage.lower()}_quadratic_live_rgbeta.json'
            raw = read(path)
            hashes[str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)] = fingerprint(path)
            betas, fields = check_export(model,stage,raw,spec,has_mass)
            orders = minimum_orders(betas,set(fields['omitted']))
            quartics = fields['quartics']
            immediate = sorted(q for q in quartics if source_terms(betas[q],set(fields['omitted'])))
            if stage=='UV':
                uv_quartic_one += len(immediate)
                if has_mass and orders[fields['mass']] == 2:
                    uv_mass_order2 += 1
            stages[stage.lower()] = {
                'immediately_generated_quartics': immediate,
                'not_immediately_generated_quartics':sorted(set(quartics)-set(immediate)),
                'one_loop_RGE_insertion_orders':orders,
                'direct_C5_quartics': fields['direct'],
                'direct_C5_quartic_orders':{q:orders[q] for q in fields['direct']},
                'mass_symbol':fields['mass'],
                'quadratic_first_generation_order':orders.get(fields['mass']) if has_mass else None,
                'zero_extra_boundary_quadratic_beta_exactly_zero':True if has_mass else None,
                'symbolic_source_terms':{
                    q:[{'coefficient':str(c),'factors':list(m)} for m,c in sorted(source_terms(betas[q],set(fields['omitted'])).items())]
                    for q in immediate},
            }
        models[model] = {'extra_quartics':len(stages['uv']['one_loop_RGE_insertion_orders'])-int(has_mass),
                         'extra_quadratic':has_mass,'stages':stages}
    if (total_quartics,total_masses)!=(50,8):
        fail('incorrect global 50 quartics / 8 quadratic slots')
    return {
        'status':'CandidateReducedScalarOneLoopRGGenerationOrderAudit',
        'scope':'symbolic UV/EFT1 saved RGBeta beta functions at zero additional quartic and quadratic boundary',
        'models_checked':16, 'quartic_slots_checked':50, 'quadratic_slots_checked':8,
        'uv_immediately_sourced_quartic_slots':uv_quartic_one,
        'uv_quadratic_first_generated_at_order_2_count':uv_mass_order2,
        'model_results':models,
        'input_sha256':hashes,
        'physics_limits':[
            'Counts use an RGBeta candidate interaction basis, not a certified Matchete/RGBeta physical coupling map.',
            'Radiative generation order counts iterated one-loop RG insertions, not one-loop matching diagrams.',
            'A zero one-loop beta at the boundary does not imply exact RG closure or a protecting symmetry.',
            'Separate UV and EFT1 chains do not implement threshold matching, field redefinitions, or a two-stage integration.',
            'No physical neutrino mass, C5 correction, or quantitative relative error bar has been computed.',
        ],
    }


def write(result: dict, destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    (destination/'reduced_scalar_rg_order_audit.json').write_text(
        json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    with (destination/'reduced_scalar_rg_order_audit.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f)
        w.writerow(['model','stage','coupling','kind','generation_order','direct_C5_vertex','immediate_source'])
        for model, item in result['model_results'].items():
            for stage, block in item['stages'].items():
                for name, order in sorted(block['one_loop_RGE_insertion_orders'].items()):
                    w.writerow([model,stage,name,'quadratic' if name==block['mass_symbol'] else 'quartic',
                                'not generated at this iterated one-loop order' if order is None else order,
                                name in block['direct_C5_quartics'],name in block['immediately_generated_quartics']])
    lines = [
        '# T3 reduced scalar-sector RG generation orders — Batch 5','',
        f"- Models: **{result['models_checked']}**; additional quartic slots: **{result['quartic_slots_checked']}**; bilinear slots: **{result['quadratic_slots_checked']}**.",
        f"- UV quartic slots sourced immediately with all additional couplings zero: **{result['uv_immediately_sourced_quartic_slots']}**.",
        f"- UV additional bilinears first reachable via two iterated one-loop beta insertions: **{result['uv_quadratic_first_generated_at_order_2_count']}**.",
        '- These counts refer to candidate RGBeta beta functions, not complete C5 matching or exact RG closure.','',
        '| Model | UV immediate quartics | UV direct C5 quartic orders | UV bilinear order | EFT1 immediate quartics | EFT1 bilinear order |',
        '|---|---:|---|---:|---:|---:|',
    ]
    for model, item in result['model_results'].items():
        uv,eft = item['stages']['uv'],item['stages']['eft1']
        d=', '.join(str(v) if v is not None else 'not sourced' for v in uv['direct_C5_quartic_orders'].values()) or '—'
        m=lambda v: '—' if v is None else str(v)
        lines.append(f"| {model} | {len(uv['immediately_generated_quartics'])} | {d} | {m(uv['quadratic_first_generation_order'])} | {len(eft['immediately_generated_quartics'])} | {m(eft['quadratic_first_generation_order'])} |")
    lines += ['', '## Interpretation','',
        'For a UV boundary where all additional couplings vanish, a first-order induced quartic has leading-log scaling',
        r'`delta_lambda ~ beta_lambda^(1) log(mu/muUV)/(16*pi^2)`.',
        'Insertion into a one-loop neutrino-mass diagram is formally of second loop order if the boundary value is zero.',
        'A bilinear generated through a sourced quartic appears after a second beta-function insertion and has a double-log scaling,',
        r'`delta_mu2 ~ mS^2 * (log(mu/muUV)/(16*pi^2))^2` times source coefficients.',
        'These are formal orders, not numerical error estimates. Accidental cancellations or large logarithms require separate analysis.','',
        f"For each of the UV and EFT1 segments, the loop-log factor for a 100:1 scale ratio is `{math.log(100)/ONE_LOOP:.6f}`.",
        '', '## Limitations','']
    lines += ['- '+s for s in result['physics_limits']]
    (destination/'reduced_scalar_rg_order_summary.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')


def main() -> None:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,default=OUT)
    args=parser.parse_args()
    result=build()
    destination=args.output_dir if args.output_dir.is_absolute() else ROOT/args.output_dir
    write(result,destination)
    print('Batch 5 symbolic RG audit:', result['models_checked'], 'models;',
          result['uv_immediately_sourced_quartic_slots'], 'UV quartic sources;',
          result['uv_quadratic_first_generated_at_order_2_count'], 'UV two-beta quadratic sectors')
    print('Input files:',len(result['input_sha256']),'SHA-256 fingerprinted')
    print('Result:',destination/'reduced_scalar_rg_order_summary.md')


if __name__ == '__main__':
    main()
