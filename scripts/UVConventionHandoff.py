"""Stage 29: convention-labelled numerical UV boundary values for general quartics.

Converts the exact 46 Stage23 *conditional* signed matrices into 16 small
numerical UV quartic dictionaries. Checks the recorded Stage26--28 evidence and
the current production ODE symbol coverage. Does NOT perform beta evolution,
modify production, or certify a physical Matchete/RGBeta component map.

    python -m scripts.UVConventionHandoff --extra-quartic 0.01
"""
from __future__ import annotations

import argparse
import ast
from collections import Counter
from decimal import Decimal, InvalidOperation
import hashlib
import json
import math
from pathlib import Path

import sympy as sp

from scripts.QuarticPrefactorAudit import exact_matrix, run as stage23_run
from scripts.RGBetaTripletNativeStructureAudit import run as stage28_run

ROOT = Path(__file__).resolve().parents[1]
DATA = Path('output/general_t3_stage22')
REGISTRY = ROOT / 'RGE/general/generated/general_rgbeta_quartic_registry.json'
PRODUCTION = {
    'state': ROOT / 'Numerical/core/State.py',
    'vector': ROOT / 'Numerical/core/StateVector.py',
    'evaluator': ROOT / 'Numerical/core/RGBetaEvaluator.py',
}
EXACT_MODELS = frozenset(
    ['A_+0', 'A_-2', 'A_-4', 'B_+1', 'B_-1', 'B_-3',
     'C_+1', 'C_-1', 'C_-3', 'D_+0', 'D_+2', 'D_-2',
     'E_+0', 'E_+2', 'E_-2', 'E_-4']
)
REQUIRED_COUNTS = {'slots': 50, 'blocks': 46}


def fail(msg: str) -> None:
    raise ValueError('Stage29 fail-closed: ' + msg)


def exact_input(value: str) -> sp.Rational:
    """Positive decimal (not an arbitrary parsed algebraic Python expression)."""
    try:
        d = Decimal(value)
    except (InvalidOperation, ValueError, TypeError):
        fail('Expected a positive decimal --extra-quartic input')
    if not d.is_finite() or d <= 0 or d > 1:
        fail('Illustrative extra quartic must satisfy 0 < lambda <= 1')
    return sp.Rational(str(d))


def text_value(expr: sp.Expr) -> dict:
    expr = sp.simplify(expr)
    if expr.free_symbols or sp.simplify(sp.im(expr)) != 0:
        fail('Expected a real, numerically evaluable coupling')
    numerical = float(sp.N(expr, 17))
    if not math.isfinite(numerical):
        fail('Nonfinite numerical coupling')
    return {'exact': str(expr), 'float': numerical}


def _read(path: Path) -> dict:
    if not path.is_file():
        fail('Missing ' + str(path))
    result = json.loads(path.read_text(encoding='utf-8-sig'))
    if not isinstance(result, dict):
        fail('JSON root is not a dictionary: ' + str(path))
    return result


def source_symbols(paths: dict[str, Path]) -> dict[str, set[str]]:
    """Inspect exact local Python source, without importing the numerical solver."""
    for name in ('state', 'vector', 'evaluator'):
        if name not in paths or not paths[name].is_file():
            fail('Missing production source ' + name)
    trees = {}
    for name, path in paths.items():
        try:
            trees[name] = ast.parse(path.read_text(encoding='utf-8-sig'), filename=str(path))
        except (OSError, UnicodeError, SyntaxError) as exc:
            fail('Cannot parse production source ' + name + ': ' + str(exc))
    classes = [x for x in trees['state'].body
               if isinstance(x, ast.ClassDef) and x.name == 'T3UVState']
    if len(classes) != 1:
        fail('Cannot identify unique production T3UVState dataclass')
    annotated = {x.target.id for x in classes[0].body
                 if isinstance(x, ast.AnnAssign) and isinstance(x.target, ast.Name)}
    if not {'lambdaS1', 'lambdaT3', 'lambdaH1', 'lambda12'} <= annotated:
        fail('Unexpected production T3UVState contract')
    result = {'state': annotated}
    for name in ('vector', 'evaluator'):
        result[name] = {x.value for x in ast.walk(trees[name])
                        if isinstance(x, ast.Constant) and isinstance(x.value, str)}
    return result


def _registry_inventory(registry: dict) -> tuple[dict[str, dict[str, set[str]]], set[tuple[str, str]], set[tuple[str, str]]]:
    totals = registry.get('totals', {})
    if (totals.get('matchete_extra_quartics') != 50 or
        totals.get('rgbeta_new_vertices') != 38 or
        totals.get('rgbeta_preexisting_structural') != 12):
        fail('Unexpected general quartic registry totals')
    models = registry.get('models')
    if not isinstance(models, dict) or set(models) != EXACT_MODELS:
        fail('Registry model inventory differs from reviewed 16 models')
    by_model, new, reused = {}, set(), set()
    for model, entry in models.items():
        nn = [x['name'] for x in entry.get('new_rgbeta_symbols', [])]
        rr = [x['existing_name'] for x in entry.get('preexisting_rgbeta_symbols', [])]
        if len(set(nn + rr)) != len(nn + rr):
            fail('Duplicate registry symbol: ' + model)
        if len(nn + rr) != entry.get('total_additional_quartics'):
            fail('Registry model total mismatch: ' + model)
        if any(not s.startswith('lambdaT3GeneralQ') for s in nn):
            fail('New quartic inventory contains an unexpected symbol')
        by_model[model] = {'new': set(nn), 'reused': set(rr)}
        new.update((model, x) for x in nn)
        reused.update((model, x) for x in rr)
    if len(new) != 38 or len(reused) != 12 or new & reused:
        fail('Invalid 38/12 registry partition')
    return by_model, new, reused


def validate_evidence(s23: dict, s26: dict, s27: dict, s28: dict) -> None:
    if (s23.get('totals', {}).get('slots') != 50 or
        s23.get('totals', {}).get('blocks') != 46 or
        s23.get('signed_physical_coupling_map_certified') is not False or
        len(s23.get('blocks', [])) != 46):
        fail('Stage23 conditional matrices missing or falsely certified')
    if (s26.get('models_checked') != 16 or s26.get('slots_checked') != 50 or
        s26.get('blocks_checked') != 46 or
        s26.get('signed_polynomial_identity_conditional_verified') is not True or
        s26.get('physical_coupling_mapping_certified') is not False or
        len(s26.get('blocks', [])) != 46):
        fail('Stage26 signed action evidence missing')
    if (s27.get('source_native_matchete_tensor_count') != 6 or
        s27.get('su2_commutators_casimir_hermiticity_intertwiner_exact') is not True or
        s27.get('physical_coupling_mapping_certified') is not False):
        fail('Stage27 triplet witness missing')
    if (s28.get('models_checked') != 6 or
        s28.get('rgbeta_delS2_orthogonality_verified') is not True or
        s28.get('matchete_triplet_invariant_metric_and_generators_exact') is not True or
        s28.get('native_rgbeta_numerical_triplet_generators_exported') is not False or
        s28.get('physical_coupling_mapping_certified') is not False):
        fail('Stage28 native S2 evidence missing or falsely certified')
    actual = Counter((x['model'], x['block']) for x in s23['blocks'])
    checked = Counter((x['model'], x['block']) for x in s26['blocks'])
    if actual != checked or any(v != 1 for v in actual.values()):
        fail('Stage23/26 block inventories are inconsistent')
    if any(x.get('exact_signed_action_identity_in_provisional_basis') is not True
           for x in s26['blocks']):
        fail('Stage26 block lacks its exact conditional action identity')


def build(s23: dict, s26: dict, s27: dict, s28: dict, registry: dict,
          source: dict[str, set[str]], quartic: sp.Rational) -> dict:
    validate_evidence(s23, s26, s27, s28)
    inventory, new, reused = _registry_inventory(registry)
    if not isinstance(quartic, sp.Rational) or not 0 < quartic <= 1:
        fail('Benchmark quartic value is not an exact positive rational <=1')
    if not {'state', 'vector', 'evaluator'} <= source.keys():
        fail('Numerical source coverage absent')
    models = {name: {'matchete_quartics': {}, 'rgbeta_quartics': {}, 'blocks': []}
              for name in sorted(EXACT_MODELS)}
    visited_mat, visited_rgb = set(), set()
    status_counts = Counter()
    block_shapes = Counter()
    max_abs = 0.0
    for block in s23['blocks']:
        name = block['model']
        if name not in models:
            fail('Unexpected Stage23 model ' + str(name))
        mn, rn = block['matchete_couplings'], block['rgbeta_couplings']
        if not isinstance(mn, list) or not isinstance(rn, list) or len(mn) != len(rn) or len(mn) not in (1, 2):
            fail('Invalid Stage23 block coupling size')
        if len(set(mn)) != len(mn) or len(set(rn)) != len(rn):
            fail('Repeated coupling within block')
        A = exact_matrix(block['r_in_m'])
        B = exact_matrix(block['lambda_rgbeta_from_lambda_matchete_if_assumptions_hold'])
        if A.shape != (len(mn), len(mn)) or B.shape != A.shape:
            fail('Block matrix dimension mismatch')
        if any(sp.simplify(z) != 0 for z in A*B + sp.eye(len(mn))):
            fail('Conditional signed coupling matrix does not invert geometry')
        if block.get('physical_mapping_certified') is not False or \
           block.get('assumed_relative_potential_sign') != '-1':
            fail('Physical certification or sign metadata unexpectedly changed')
        vals = sp.Matrix([quartic] * len(mn))
        converted = B * vals
        if any(sp.simplify(z) != 0 for z in A*converted + vals):
            fail('Numerical signed coupling reconstruction failed')
        for symbol, value in zip(mn, vals):
            key = (name, symbol)
            if key in visited_mat:
                fail('Duplicate Matchete UV boundary symbol')
            visited_mat.add(key)
            models[name]['matchete_quartics'][symbol] = text_value(value)
        for symbol, value in zip(rn, converted):
            key = (name, symbol)
            if key in visited_rgb or symbol not in inventory[name]['new'] | inventory[name]['reused']:
                fail('Duplicate or unregistered RGBeta UV boundary symbol')
            visited_rgb.add(key)
            obj = text_value(value)
            obj['registry_class'] = ('new' if symbol in inventory[name]['new'] else 'reused')
            models[name]['rgbeta_quartics'][symbol] = obj
            max_abs = max(max_abs, abs(obj['float']))
        models[name]['blocks'].append({
            'block': block['block'],
            'geometry_class': block['geometry_class'],
            'matchete_symbols': mn, 'rgbeta_symbols': rn,
            'signed_map_exact': block['lambda_rgbeta_from_lambda_matchete_if_assumptions_hold'],
            'reconstruction_residual_exact': '0',
        })
        status_counts[block['geometry_class']] += len(mn)
        block_shapes[len(mn)] += 1
    expected = new | reused
    if visited_rgb != expected or len(visited_mat) != 50 or len(visited_rgb) != 50:
        fail('UV handoff does not exactly exhaust 50 quartic slots')
    # Record the *actual* implementation gap; never imply imported values have
    # entered the production real ODE state or its beta evaluator.
    coverage = {}
    for layer in ('state', 'vector', 'evaluator'):
        fields = source[layer]
        coverage[layer] = {
            'new_quartic_slots_recognised': sum(sym in fields for _, sym in new),
            'reused_quartic_slots_recognised': sum(sym in fields for _, sym in reused),
            'missing_new_symbols': sorted({sym for _, sym in new if sym not in fields}),
        }
    return {
        'format': 'stage29-conditional-uv-quartic-handoff-v1',
        'status': 'NumericalBoundaryCandidatesOnlyNoRGEIntegration',
        'models_checked': 16, 'slots_checked': 50, 'blocks_checked': 46,
        'new_quartic_slots': 38, 'reused_quartic_slots': 12,
        'benchmark_lambda_matchete_per_additional_slot_exact': str(quartic),
        'benchmark_scale_gev': 1.0e7,
        'benchmark_scale_source': 'docs/PHYSICS.md chosen UV scale (not an RGE result)',
        'lambda_t3_topology_coupling_modified': False,
        'class_E_topology_lambda_t3_factor_applied_to_additional_quartics': False,
        'triplet_component_convention': 'Stage27 single Matchete-CG reference convention, provisionally identified with abstract RGBeta S2',
        'triplet_barred_field_physical_intertwiner_certified': False,
        'source_signed_action_conditional_verified': True,
        'conditional_block_reconstruction_verified': True,
        'physical_coupling_mapping_certified': False,
        'one_loop_beta_functions_for_new_quartics_evaluated': False,
        'rge_evolution_performed': False,
        'production_state_modified': False,
        'production_ode_symbol_coverage': coverage,
        'geometry_class_slot_counts': dict(sorted(status_counts.items())),
        'block_dimension_counts': {str(k): v for k, v in sorted(block_shapes.items())},
        'max_abs_conditional_rgbeta_quartic_uv': max_abs,
        'models': models,
        'next_step': 'Export 38 new one-loop quartic beta functions and extend numerical state, RGBeta evaluator and vector layout before real RGE evolution.',
    }


def run(folder: Path = DATA, output: Path | None = None, *,
        root: Path = ROOT, extra_quartic: str = '0.01', check_upstream: bool = True) -> dict:
    stored23 = _read(folder / 'stage23_quartic_prefactor_candidates.json')
    stored26 = _read(folder / 'stage26_conditional_quartic_action.json')
    stored27 = _read(folder / 'stage27_triplet_generator_witness.json')
    stored28 = _read(folder / 'stage28_triplet_native_structure_audit.json')
    if check_upstream:
        # Re-run the source-linked upstream checks to reject stale JSONs.
        if stage23_run(folder, None) != stored23:
            fail('Stage23 stored signed map is stale')
        if stage28_run(folder, None, check_stage27=True) != stored28:
            fail('Stage28 stored native structure is stale')
        # stage28_run rechecks Stage27 -> Stage26 -> Stage25 and their sources.
    paths = {name: root / path.relative_to(ROOT) for name, path in PRODUCTION.items()}
    inputs = [folder / s for s in (
        'stage23_quartic_prefactor_candidates.json',
        'stage26_conditional_quartic_action.json',
        'stage27_triplet_generator_witness.json',
        'stage28_triplet_native_structure_audit.json',
    )] + [root / 'RGE/general/generated/general_rgbeta_quartic_registry.json'] + list(paths.values())
    report = build(stored23, stored26, stored27, stored28,
                   _read(root / 'RGE/general/generated/general_rgbeta_quartic_registry.json'),
                   source_symbols(paths), exact_input(extra_quartic))
    report['input_sha256'] = {str(p.relative_to(root)) if p.is_relative_to(root) else p.name:
                              hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}
    report['upstream_full_freshness_checked'] = bool(check_upstream)
    if output is not None:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, default=DATA)
    parser.add_argument('--output', type=Path,
                        default=DATA / 'stage29_convention_uv_handoff.json')
    parser.add_argument('--extra-quartic', default='0.01',
                        help='Illustrative Matchete coefficient for each of 50 slots (default 0.01)')
    ns = parser.parse_args()
    try:
        r = run(ns.input_dir, ns.output, extra_quartic=ns.extra_quartic)
    except (KeyError, ValueError, TypeError, OSError, json.JSONDecodeError, sp.SympifyError) as exc:
        parser.exit(2, str(exc) + '\n')
    print('Stage29 convention-labelled numerical UV values: 16/16 models, 50/50 slots, 46/46 blocks')
    print('New slots needing production one-loop state/beta support: 38')
    print('Largest |lambda_RGBeta| at chosen UV point:', format(r['max_abs_conditional_rgbeta_quartic_uv'], '.8g'))
    print('Actual RGE integration performed: NO')
    print('Physical Matchete/RGBeta coupling mapping: UNCERTIFIED')
    print('Report:', ns.output)


if __name__ == '__main__':
    main()
