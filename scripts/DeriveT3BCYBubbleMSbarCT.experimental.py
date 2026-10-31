#!/usr/bin/env python3
"""Workstream C: independent UV/MSbar subtraction for the *scalar-bubble sector*.

No Matchete matching is called. Reconstruct K_YY and K_XY from the saved EFT1
LLSS/HHSS tree vertices and exported native CG tensors (v4), then derive the
UV singularity of a two-scalar loop from its momentum integral, NOT by
negating the native Matchete pole. Compare the resulting -K/eps CT to the
independent native Y/LLHH pole as an audit.

SCOPE IS LIMITED: an ordered neutral-component, native GammaCC/Proj basis,
with common hbar=1/(16*pi^2), 1/MF factored out. Not the full EFT1/EFT2
counterterm difference; not UV/IR certification for other EFT diagrams;
not a normalized, symmetric physical Weinberg C5.

Run at the repository root:
    python scripts/DeriveT3BCYBubbleMSbarCT.experimental.py --self-test
    python scripts/DeriveT3BCYBubbleMSbarCT.experimental.py --output output/audits/t3_bc_y_scalar_bubble_msbar_v1/summary.json
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from hashlib import sha256
import importlib.util
import json
import math
from pathlib import Path
import sys

import sympy as sp

STEMS = {'B_m1': 'B_-1', 'C_m1': 'C_-1'}
BASE = Path('output/general_t3_mass_eigen_matchete')
AUDITS = Path('output/audits')
EXPECTED_KEYS = {'YY': 18, 'XY': 32}


class BubbleCTError(ValueError):
    pass


def require(condition, message):
    if not condition:
        raise BubbleCTError(message)


def sibling(script: str, modname: str):
    path = Path(__file__).resolve().with_name(script)
    require(path.is_file(), 'Missing prerequisite complete script: ' + str(path))
    spec = importlib.util.spec_from_file_location(modname, path)
    require(spec is not None and spec.loader is not None, 'Cannot load ' + str(path))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def equals(a: dict, b: dict):
    """Exact comparison after combining SU2 components and symbolic factors."""
    keys = set(a) | set(b)
    differences = []
    for key in sorted(keys, key=str):
        residual = sp.simplify(a.get(key, 0) - b.get(key, 0))
        if residual != 0:
            differences.append({'key': repr(key), 'difference': str(residual),
                                'lhs': str(a.get(key, 0)),
                                'rhs': str(b.get(key, 0))})
    return differences


def sum_maps(a: dict, b: dict) -> dict:
    return {k: sp.simplify(a.get(k, 0) + b.get(k, 0))
            for k in set(a) | set(b)
            if sp.simplify(a.get(k, 0) + b.get(k, 0)) != 0}


def invert_map(a: dict) -> dict:
    return {key: -value for key, value in a.items()}


def save_coefficients(a: dict) -> list:
    return [{'coupling_factors': list(key[0]), 'chirality': key[1],
             'H_bar_count': key[2], 'coefficient_without_hbar_and_inverse_MF': str(val)}
            for key, val in sorted(a.items(), key=lambda item: str(item[0]))]


def scalar_bubble_cutoff_uv_slope(cutoff_squared, mass1_squared, mass2_squared):
    """d/d ln(Lambda^2) of the dimensionless Euclidean bubble cutoff integral.

    J = int_0^(Lambda^2) dt t/((t+m1^2)(t+m2^2)).
    Its UV slope is independent of both masses and equals one.
    """
    L = cutoff_squared
    return sp.factor(L**2 / ((L + mass1_squared) * (L + mass2_squared)))


def bubble_finite_msbar(r: float) -> float:
    """F_XY at mu=mY, full massive bubble, for reference only."""
    require(math.isfinite(r) and 0.0 <= r <= 1.0, 'Require 0 <= r <= 1')
    if r == 0.0:
        return 1.0
    if 1.0-r < 1e-5:
        d = 1.0-r
        return sum(d**n / (n*(n+1)) for n in range(1, 8))
    return 1.0 + r * math.log(r) / (1.0-r)


def feynman_parameter_quadrature(r, n=40000):
    """Independent direct integral -int_0^1 dx ln(r+(1-r)*x)."""
    require(0 < r <= 1, 'Quadrature requires r>0')
    return -sum(math.log(r + (1-r)*(i+.5)/n) for i in range(n))/n


def audit_one(stem, v3, v2, v4, root, native_dir, ledger_dir, old_v4):
    model = STEMS[stem]
    registry = native_dir / f'{stem}_native_CG.json'
    kernels, counts, hashes = v4.get_pair_kernels(v3, v2, root, stem, registry)
    for channel, expected in EXPECTED_KEYS.items():
        require(len(kernels[channel]) == expected,
                f'{model}: changed nonzero {channel} kernel inventory; '
                f'expected {expected}, got {len(kernels[channel])}')
    require(not (set(kernels['YY']) & set(kernels['XY'])),
            f'{model}: YY and XY coefficient supports overlap unexpectedly')
    # These W/Q kernels come entirely from EFT1 tree vertices and native CGs,
    # without using the matched Matchete pole to assign the CT coefficients.
    independent_uv = sum_maps(kernels['YY'], kernels['XY'])
    independent_msbar_ct = invert_map(independent_uv)

    pole_file = ledger_dir / f'{stem}_y_llhh_raw_pole_residue.inputform.txt'
    finite_file = (ledger_dir /
                   f'{stem}_y_llhh_formal_finite_at_muY.NOT_PHYSICAL_C5.inputform.txt')
    native = v3.read_registry(registry, model)
    pole, pole_sha = v4.native_sector(v3, v2, pole_file, native, 50)
    finite, finite_sha = v4.native_sector(v3, v2, finite_file, native, 32)
    require(len(pole) == 50 and len(finite) == 32,
            f'{model}: Native sector inventory changed')
    pole_diff = equals(independent_uv, pole)
    finite_diff = equals(kernels['XY'], finite)
    # Observe two distinct identities, no assignment from a native residue:
    # (i) bare bubble UV K independently computed, (ii) CT = -K from MSbar.
    ct_cancel = equals(sum_maps(independent_uv, independent_msbar_ct), {})
    native_cancel = equals(sum_maps(pole, independent_msbar_ct), {})

    prior = old_v4.get('models', {}).get(model)
    require(isinstance(prior, dict) and
            prior.get('format') == 't3bc-y-scalar-bubble-finitemass-diagnostic-v1',
            f'{model}: Missing v4 audit; run prior complete script first')
    expected_hashes = {**hashes, 'pole': pole_sha, 'formal_finite': finite_sha}
    require(prior.get('source_sha256') == expected_hashes,
            f'{model}: Original source / CG / ledger hashes differ from v4 audit')
    require(prior.get('tree_UV_kernel_equals_native_pole') is True and
            prior.get('tree_XY_kernel_equals_formal_finite_at_muY') is True,
            f'{model}: Earlier native kernel tests not passed')

    masses = root / f'{stem}_final_parameters.json'
    p = json.loads(masses.read_text(encoding='utf-8'))
    require(p.get('model') == model, f'{model}: Parameter model mismatch')
    m2 = p['heavy_mass_squared_gev2']
    require(len(m2) == 2 and 0 < m2[0] < m2[1],
            f'{model}: Bad scalar mass hierarchy')
    r = float(m2[0]/m2[1])
    full_value = bubble_finite_msbar(r)
    require(abs(full_value - feynman_parameter_quadrature(r)) < 1e-7,
            f'{model}: Independent Feynman parameter integral fails')
    require(abs(full_value - prior['scalar_bubble_msbar_finite_at_muY']['XY_full_massive']) < 1e-12,
            f'{model}: Unexpected v4 finite massive-bubble result')

    # For the strict mX=0, p=0 bubble the IR endpoint is integrable:
    #  dt t/((t+M^2)t) = dt/(t+M^2) near t=0.
    # For YY, t/(t+M^2)^2 is also regular near t=0.
    # The Y-related subset is UV-divergent but has no IR pole.
    result = {
        'format': 't3bc-y-native-scalar-bubble-msbar-subtraction-v1',
        'model': model,
        'source_sha256': expected_hashes,
        'counts': counts,
        'nonzero_component_monomials': {
            'YY': len(kernels['YY']), 'XY': len(kernels['XY']),
            'native_pole': len(pole), 'native_formal_finite': len(finite),
        },
        'mass_ratio_r': r,
        'bubble_finite_reference_at_muY': {
            'YY': 0, 'XY_massless_X': 1,
            'XY_both_masses_kept': full_value,
            'strict_dimension5_vs_full_mass_difference': 1-full_value,
        },
        'uv_proof': {
            'euclidean_integral':
                'J(L;a,b)=integral_0^L [t dt / ((t+a)(t+b))], L=Lambda^2',
            'uv_log_slope': 'lim_(L->infinity) dJ/dlnL = 1 for a>0,b>=0',
            'massless_X_IR_integrand':
                't/((t+MY^2)*t) -> 1/MY^2 at t->0, no IR pole',
            'YY_IR_integrand':
                't/(t+MY^2)^2 -> 0 at t->0, no IR pole',
            'scheme': 'MSbar; 1/epsilon_bar subtracted, +hbar loop convention',
            'MSbar_pole_coefficient_per_Wick_contracted_pair': '1',
            'independent_of_native_Matchete_pole': True,
        },
        'tree_WQ_kernel_equals_native_Y_pole': not pole_diff,
        'independent_bubble_MSbar_CT_cancels_tree_WQ_UV': not ct_cancel,
        'independent_bubble_MSbar_CT_cancels_native_Y_pole': not native_cancel,
        'Y_formal_finite_equals_XY_kernel_at_muY': not finite_diff,
        'UV_comparison_mismatches': pole_diff,
        'CT_cancel_mismatches': ct_cancel,
        'native_CT_cancel_mismatches': native_cancel,
        'finite_comparison_mismatches': finite_diff,
        'independent_bubble_sector_CT_without_hbar_and_inverse_MF':
            save_coefficients(independent_msbar_ct),
        'independent_bubble_sector_MSbar_hard_finite_at_muY_without_hbar_and_inverse_MF':
            save_coefficients(kernels['XY']),
        'operator_convention':
            'native ordered LLHH neutral component, GammaCC/Proj[+/-1], hbar/MF factored out',
        'derived_CT_operator_definition':
            'delta L_LLHH = (hbar/MF) * (-K_YY-K_XY) / epsilon_bar, '
            'in the native ordered neutral component at order quartic*LLSS',
        'derived_hard_finite_operator_definition':
            'Delta L_LLHH^bubble,MSbar(muY) = (hbar/MF)*K_XY, '
            'strict massless X at dimension-5; YY term vanishes at muY',
        'scalar_bubble_sector_CT_derived_independently': True,
        'UV_IR_separated_for_YY_XY_scalar_bubbles': True,
        'EFT1_all_operator_counterterms_derived': False,
        'EFT1_EFT2_full_counterterm_difference_derived': False,
        'EFT2_light_loop_soft_overlap_checked': False,
        'higher_dim_mass_corrections_included': False,
        'Majorana_GammaCC_to_symmetric_Weinberg_normalization_certified': False,
        'physical_C5_Y': None,
        'physical_matching_certified': False,
    }
    return result


def self_test():
    L, a, b = sp.symbols('L a b', positive=True)
    uv_slope = scalar_bubble_cutoff_uv_slope(L,a,b)
    require(sp.limit(uv_slope,L,sp.oo) == 1,
            'Independent UV cutoff integral does not have unit residue')
    # Test one massless line and two equal-mass lines.
    require(sp.limit(scalar_bubble_cutoff_uv_slope(L,a,sp.Integer(0)),L,sp.oo)==1,
            'Mixed massless-X/massive-Y UV residue is not one')
    require(sp.limit(scalar_bubble_cutoff_uv_slope(L,a,a),L,sp.oo)==1,
            'Equal-mass YY UV residue is not one')
    t=sp.symbols('t',nonnegative=True)
    require(sp.limit(t/((t+a)*t),t,0,dir='+') == 1/a,
            'Mixed bubble has an unexpected IR endpoint')
    require(sp.limit(t/((t+a)**2),t,0,dir='+') == 0,
            'YY bubble has an unexpected IR endpoint')
    for r in (0.01,0.5,0.76,0.99999):
        require(abs(bubble_finite_msbar(r)-feynman_parameter_quadrature(r,4000)) < 1e-6,
                'Feynman parameter quadrature disagrees with closed form')
    require(sp.simplify(sp.Integer(1)+sp.Rational(1,2)*sp.log(sp.Rational(1,2)) /
                        (1-sp.Rational(1,2)) - (1-sp.log(2)))==0,
            'Feynman parameter exact test fails')
    require(bubble_finite_msbar(0)==1 and bubble_finite_msbar(1)==0,
            'Endpoint finite convention changed')
    example={('A','PR',0):sp.Rational(1,2),('B','PL',2):sp.Integer(1)}
    cancel= sum_maps(example,invert_map(example))
    require(cancel=={},'UV CT subtraction failed')
    require(len(equals(example,{'B':sp.Integer(1)}))>0,
            'Tampered coefficient was not rejected')
    print('PASS: independent cutoff-integral UV slope and both IR endpoints')
    print('PASS: MSbar massive-bubble finite Feynman-parameter integral')
    print('PASS: independent source-kernel CT algebra and tamper rejection')
    print('PASS: no full EFT1/EFT2 CT difference or physical C5 certification')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--self-test', action='store_true')
    ap.add_argument('--output', type=Path)
    ap.add_argument('--root', type=Path, default=BASE)
    ap.add_argument('--native-dir', type=Path,
                    default=AUDITS/'t3_bc_y_native_cg_components_v1')
    ap.add_argument('--ledger-dir', type=Path,
                    default=AUDITS/'t3_bc_y_renormalization_ledger_v1')
    ap.add_argument('--v4-summary', type=Path,
                    default=AUDITS/'t3_bc_y_finite_scalar_bubble_v1'/'summary.json')
    args = ap.parse_args()
    try:
        if args.self_test:
            self_test()
        if args.output:
            require(args.v4_summary.is_file(), 'Missing successful v4 summary: '+str(args.v4_summary))
            prior = json.loads(args.v4_summary.read_text(encoding='utf-8'))
            require(prior.get('format') == 't3bc-y-bubble-mass-finite-audit-v1',
                    'Unexpected v4 summary format: '+str(prior.get('format')))
            v3 = sibling('CompareT3BCYNativeSU2Coefficients.experimental.py','t3bc_v3_uv')
            v2 = sibling('AuditT3BCYMixedScalarPoleSupport.experimental.py','t3bc_v2_uv')
            v4 = sibling('CheckT3BCYFiniteScalarBubble.experimental.py','t3bc_v4_uv')
            summary={
                'format':'t3bc-y-scalar-bubble-msbar-counterterm-audit-v1',
                'model_results':{},
                'diagnostic_scope':
                    'independent source-vertex scalar bubble CT only; NOT full EFT matching',
                'physical_C5_Y':None,
                'full_EFT1_EFT2_counterterm_difference_certified':False,
            }
            for stem in STEMS:
                result=audit_one(stem,v3,v2,v4,args.root,args.native_dir,
                                 args.ledger_dir,prior)
                summary['model_results'][result['model']]=result
                good=(result['tree_WQ_kernel_equals_native_Y_pole'] and
                      result['independent_bubble_MSbar_CT_cancels_native_Y_pole'] and
                      result['Y_formal_finite_equals_XY_kernel_at_muY'])
                print(f"{result['model']}: "
                      f"{result['nonzero_component_monomials']['YY']} YY + "
                      f"{result['nonzero_component_monomials']['XY']} XY; "
                      f"source-derived bubble CT cancels native pole={good}; "
                      "full EFT counterterms=NOT CERTIFIED")
                require(good, result['model']+': Source-UV/native UV/finite mismatch')
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(json.dumps(summary,indent=2,ensure_ascii=False)+'\n',
                                   encoding='utf-8')
            print('PASS: independent scalar-bubble CT diagnostic -> '+str(args.output))
        if not (args.self_test or args.output):
            ap.error('Provide --self-test and/or --output')
        return 0
    except (BubbleCTError, ValueError, TypeError, KeyError, IndexError,
            RuntimeError, AssertionError, FileNotFoundError, OSError) as e:
        print('FAIL (no matching certification): '+str(e),file=sys.stderr)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
