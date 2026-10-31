#!/usr/bin/env python3
"""Workstream C: explicit hard/soft regions for the zero-momentum XY scalar bubble.

This is a dimension-five power-counting and finite-mass accuracy audit, NOT a
full EFT2 method-of-regions amplitude, operator counterterm or physical C5.
The inputs are the actual B(-1)/C(-1) parameter JSONs and the previously
successful v4 and v6 local reports. No synthetic physical couplings are used.

Run from the repository root:
  python scripts/AnalyzeT3BCYXYRegions.experimental.py --self-test
  python scripts/AnalyzeT3BCYXYRegions.experimental.py --output output/audits/t3_bc_y_xy_regions_v1/summary.json
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
import sys

ROOT = Path('output/general_t3_mass_eigen_matchete')
AUDITS = Path('output/audits')
STEMS = {'B_m1': 'B_-1', 'C_m1': 'C_-1'}
V4 = AUDITS / 't3_bc_y_finite_scalar_bubble_v1/summary.json'
V6 = AUDITS / 't3_bc_y_eft2_light_xx_overlap_v1/summary.json'

class RegionsError(ValueError):
    pass


def require(condition: bool, message: str):
    if not condition:
        raise RegionsError(message)


def load(path: Path) -> tuple[dict, str]:
    require(path.is_file(), 'Missing local prerequisite: ' + str(path))
    raw = path.read_bytes()
    return json.loads(raw), sha256(raw).hexdigest()


def finite_bubble(r: float, mu_over_M: float = 1.0) -> float:
    """Zero-momentum MSbar B0^fin(M^2,m^2;mu), excluding 1/(16pi^2)."""
    require(math.isfinite(r) and 0 < r < 1, 'Require nondegenerate 0<r<1')
    require(math.isfinite(mu_over_M) and mu_over_M > 0, 'Invalid mu/M')
    # With t=1-r, avoid loss of significance when nearly degenerate.
    t = 1.0-r
    quotient = (r * math.log1p(-t))/t
    return 1.0 + quotient + 2.0*math.log(mu_over_M)


def hard_leading(mu_over_M: float = 1.0) -> float:
    """Finite piece of hard, order-r^0 integral; UV pole is separately +1."""
    return 1.0 + 2.0*math.log(mu_over_M)


def region_at_order(r: float, n: int, mu_over_M: float = 1.0) -> dict:
    """For n>=1, hard IR and soft UV poles are separate until summed.

    Exact MSbar finite terms at that order:
      H_n = r^n [ +1/eps_IR + 1 + log(mu^2/M^2) ]
      S_n = r^n [ -1/eps_UV - 1 - log(mu^2/m^2) ]
      sum finite = r^n log(r).
    The common 1/(16pi^2) is suppressed. No soft n=0 term.
    """
    require(isinstance(n, int) and n >= 1, 'n>=1 required')
    require(0 < r < 1 and mu_over_M > 0, 'Invalid mass/scale ratio')
    a = r**n
    hard = a*(1.0 + 2.0*math.log(mu_over_M))
    soft = -a*(1.0 + 2.0*math.log(mu_over_M) - math.log(r))
    return {
        'power_n': n,
        'power_r_n': a,
        'hard_IR_pole_residue': a,
        'soft_UV_pole_residue': -a,
        'pole_sum_if_epsilon_IR_identified_with_epsilon_UV': 0.0,
        'hard_finite': hard,
        'soft_finite': soft,
        'combined_finite': hard+soft,
        'combined_finite_expected': a*math.log(r),
    }


def truncated(r: float, n: int, mu_over_M: float = 1.0) -> tuple[float, float]:
    """Sum hard leading and n high-power hard+soft combinations; exact tail."""
    require(isinstance(n, int) and n >= 0, 'n must be nonnegative')
    leading = hard_leading(mu_over_M)
    # Avoid loss of significance from the explicit geometric sum.
    partial = leading + math.log(r)*r*(-math.expm1(n*math.log(r)))/(1-r)
    tail = math.log(r)*r**(n+1)/(1-r)
    return partial, tail


def model_result(stem: str, root: Path, v4: dict, v6: dict,
                 v4sha: str, v6sha: str, orders: int) -> dict:
    name = STEMS[stem]
    a4 = v4.get('models', {}).get(name)
    a6 = v6.get('models', {}).get(name)
    require(isinstance(a4, dict) and isinstance(a6, dict), name + ': missing v4/v6 benchmark record')
    require(a4.get('tree_UV_kernel_equals_native_pole') is True and
            a4.get('tree_XY_kernel_equals_formal_finite_at_muY') is True,
            name + ': v4 failed source/finite comparison')
    checks6 = (
        'exact_LLXX_component_equality', 'exact_HHXX_component_equality',
        'exact_XX_loop_UV_component_equality',
        'EFT1_minus_EFT2_scalar_bubble_UV_equals_Y_vertices',
        'EFT1_minus_EFT2_scalar_bubble_CT_equals_minus_Y',
        'EFT1_EFT2_scalar_bubble_CT_difference_cancels_native_Y_pole',
        'XX_finite_MSbar_log_cancelled_between_EFTs',
        'independent_Y_kernel_equals_native_pole',
        'previous_v5_CT_matches_v6_CT_difference',
    )
    require(all(a6.get(flag) is True for flag in checks6), name + ': v6 prerequisite failed')
    require(a6.get('complete_EFT1_EFT2_operator_counterterms_certified') is False and
            a6.get('physical_C5_Y') is None,
            name + ': unexpected full renormalization certification in v6')
    param_file = root / f'{stem}_final_parameters.json'
    par, phash = load(param_file)
    require(par.get('model') == name, name + ': wrong parameter JSON')
    require(par.get('threshold_plan', {}).get('decreasing_scale_order') ==
            ['NewScalarY', 'NewScalarX'], name + ': wrong threshold order')
    m2 = par['heavy_mass_squared_gev2']
    m = par['heavy_masses_gev']
    require(len(m2)==2 and len(m)==2 and all(float(x)>0 for x in m2+m),
            name + ': invalid mass array')
    r = float(m2[0])/float(m2[1])
    require(0<r<1, name+': heavy mass ordering wrong')
    require(abs(float(m[0])**2/float(m[1])**2-r) < 5e-12,
            name + ': squared masses do not match masses')
    require(abs(r-float(a4['squared_mass_ratio_r']))<1e-12,
            name + ': v4 finite-bubble mass ratio changed')
    f = finite_bubble(r)
    prev = a4['scalar_bubble_msbar_finite_at_muY']['XY_full_massive']
    require(abs(f-float(prev))<2e-13, name + ': v4 full-mass bubble disagrees')
    require(a6['independent_Y_kernel_equals_native_pole'] is True,
            name + ': native source coefficient check failed')
    rows = [region_at_order(r, n) for n in range(1,orders+1)]
    require(all(math.isclose(t['combined_finite'],t['combined_finite_expected'],
                             rel_tol=3e-15,abs_tol=2e-16) for t in rows),
            name + ': hard/soft expansion mismatch')
    cumulative=[]
    for n in [0,1,2,3,5,10,20,30,orders]:
        if n>orders or n in [t['order'] for t in cumulative]:continue
        value,tail=truncated(r,n)
        require(abs((value+tail)-f)<2e-13, name+': expansion tail identity failed')
        cumulative.append({'order':n,'finite_up_to_order':value,
                           'tail_to_full_bubble':tail,
                           'abs_tail':abs(tail),
                           'tail_relative_to_full_bubble':abs(tail)/abs(f)})
    return {
        'format':'t3bc-y-xy-method-of-regions-scalar-integral-v1',
        'model':name,
        'source_sha256':{'parameters':phash,'v4_summary':v4sha,'v6_summary':v6sha},
        'mass_gev':{'X':float(m[0]),'Y':float(m[1])},
        'mass_squared_ratio_r':r,
        'log_r':math.log(r),
        'hard_leading_dimension5_finite_muY':1.0,
        'full_mass_bubble_finite_muY':f,
        'full_minus_dimension5':f-1.0,
        'first_power_correction_r_log_r':r*math.log(r),
        'full_mass_over_leading_dim5_factor':f,
        'first_soft_region_occurs_at_order_r':True,
        'first_local_operator_required_in_EFT2':
            'schematic dim-7 LLHHXX/MY^2 from tree Y exchange, plus h.c.; not exported/matched',
        'region_per_order':rows,
        'cumulative_power_expansion':cumulative,
        'v4_source_kernel_and_finite_agreement_verified_from_prior_report':True,
        'v6_XX_operator_UV_and_CT_difference_verified_from_prior_report':True,
        'full_mass_factor_must_not_be_used_as_physical_C5_multiplier':True,
        'gammaCC_weinberg_normalization_certified':False,
        'dim7_EFT2_bare_operator_matching_performed':False,
        'dim7_EFT2_renormalization_performed':False,
        'complete_EFT1_EFT2_operator_counterterms_certified':False,
        'physical_C5_Y':None,
    }


def self_test():
    # Algebraic first UV-IR-cancelled term and exact geometric tail.
    for r in (0.01,0.25,0.759189881435,0.999):
        exact = finite_bubble(r)
        for n in (0,1,3,8,30):
            v,t = truncated(r,n)
            assert abs(v+t-exact) < 1e-12, (r,n,v,t,exact)
        for n in (1,2,4,10):
            row = region_at_order(r,n)
            assert abs(row['hard_IR_pole_residue']+row['soft_UV_pole_residue'])<1e-15
            assert math.isclose(row['combined_finite'],r**n*math.log(r),rel_tol=1e-10,abs_tol=1e-15)
    # Exact benchmark limit: r=1/2 -> 1-ln(2).
    assert abs(finite_bubble(0.5) - (1.0-math.log(2.0)))<1e-14
    assert abs(finite_bubble(0.5,2.0)-(1.0-math.log(2.0)+2*math.log(2)))<1e-14
    # Nonzero dimension-seven correction: r ln r, not 0 or full remainder.
    assert abs(truncated(0.5,1)[0]-(1.0-0.5*math.log(2))) < 1e-14
    assert abs(truncated(0.5,1)[1]+0.5*math.log(2)) < 1e-14
    try:
        finite_bubble(1.2)
    except RegionsError:pass
    else:raise AssertionError('Bad heavy-scalar ordering accepted')
    print('PASS: hard n=0 UV pole, n>=1 hard-IR/soft-UV opposite residues')
    print('PASS: exact full-mass bubble, power-series tail and scale dependence')
    print('PASS: dimension-7 soft onset; physical C5 and dim-7 matching uncertified')


def main() -> int:
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--self-test',action='store_true')
    ap.add_argument('--output',type=Path)
    ap.add_argument('--root',type=Path,default=ROOT)
    ap.add_argument('--v4-summary',type=Path,default=V4)
    ap.add_argument('--v6-summary',type=Path,default=V6)
    ap.add_argument('--max-order',type=int,default=30)
    args=ap.parse_args()
    try:
        require(args.max_order>=10 and args.max_order<=100,
                'Use a max order between 10 and 100')
        if args.self_test:self_test()
        if args.output:
            a,sha4=load(args.v4_summary)
            b,sha6=load(args.v6_summary)
            require(a.get('format')=='t3bc-y-bubble-mass-finite-audit-v1',
                    'Wrong v4 summary format')
            require(b.get('format')=='t3bc-y-e1-e2-light-xx-overlap-report-v1',
                    'Wrong v6 summary format')
            result={'format':'t3bc-y-xy-regions-report-v1','models':{},
                    'scope':'Zero-momentum scalar bubble only, at mu=MY',
                    'region_pole_assignment':'Hard n>=1 poles are IR; soft n>=1 poles are UV',
                    'full_eft2_dim7_matching_performed':False,
                    'full_EFT1_EFT2_renormalization_certified':False,
                    'physical_C5_Y':None}
            for stem,name in STEMS.items():
                one=model_result(stem,args.root,a,b,sha4,sha6,args.max_order)
                result['models'][name]=one
                print(f"{name}: r={one['mass_squared_ratio_r']:.9f}; "
                      f"hard(dim5)={one['hard_leading_dimension5_finite_muY']:.9f}; "
                      f"full massive={one['full_mass_bubble_finite_muY']:.9f}; "
                      f"combined power remainder={one['full_minus_dimension5']:.9f}")
            args.output.parent.mkdir(parents=True,exist_ok=True)
            args.output.write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
            print('PASS: scalar XY hard/soft region diagnostic -> '+str(args.output))
            print('Physical C5: NOT CERTIFIED')
        else:
            require(args.self_test,'Specify --self-test and/or --output')
        return 0
    except (RegionsError,ValueError,KeyError,TypeError,AssertionError,IOError,
            OSError,OverflowError,ZeroDivisionError) as exc:
        print('FAIL (physical C5 not certified): '+str(exc),file=sys.stderr)
        return 1

if __name__=='__main__':
    raise SystemExit(main())
