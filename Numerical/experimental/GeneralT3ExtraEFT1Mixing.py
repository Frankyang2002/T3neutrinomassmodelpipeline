"""Experimental B/C(-1) extra-quartic contribution to EFT1 LLSS -> Weinberg RGE.

The real basis is (Re H1,Im H1,Re H2,Im H2,
                   Re S1_1,Im S1_1,Re S1_2,Im S1_2,
                   Re S2_1,Im S2_1,Re S2_2,Im S2_2).
The quartic potential convention is V4 = lambda_abcd phi_a phi_b phi_c phi_d/4!.
Only the scalar-pair term (not a complete EFT1 Wilson beta or matching) is
calculated, with the existing tree-level Matchete PL Wilson seed. No beta or
C5 is inserted into the production numerical trajectory.
"""
from __future__ import annotations

from collections import Counter
from itertools import combinations_with_replacement
from pathlib import Path
from typing import Any, Mapping
import json

import sympy as sp

ACTIVE_RGBETA = ('lambdaHHdagS2S2', 'lambdaHHdagS1barS1bar')
ACTIVE_MAT = ('lambdaT3ExtraQ020020I1', 'lambdaT3ExtraQ020200I1')
HIGGS_REAL_INDICES = (1, 2, 3, 4)


def _real_fields():
    return sp.symbols('r1:13', real=True)


def extra_potential_polynomial(*, include_s2: bool = True,
                               include_s1: bool = True) -> sp.Expr:
    """Use certified T=(delta delta + exchange)/2 or eps-eps tensor.

    λ here is RGBeta's potential coefficient, not Matchete's Lagrangian
    parameter.  The Stage15 source check identified λ_M = -λ_RGBeta.
    """
    r = _real_fields()
    z = tuple((r[2*k] + sp.I*r[2*k+1])/sp.sqrt(2) for k in range(6))
    h, s1, s2 = z[:2], z[2:4], z[4:6]
    lam_s2, lam_s1 = (sp.Symbol(x) for x in ACTIVE_RGBETA)
    monomials = []
    if include_s2:
        # Hbar Hbar S2 S2: P_abcd Hbar_a Hbar_b S2_c S2_d
        mono = (sp.conjugate(h[0])*s2[0] + sp.conjugate(h[1])*s2[1])**2
        monomials.append(lam_s2*mono + sp.conjugate(lam_s2*mono))
    if include_s1:
        # Hbar Hbar S1bar S1bar: eps_12=+1, project symmetrically.
        mono = (sp.conjugate(h[0])*sp.conjugate(s1[1]) -
                sp.conjugate(h[1])*sp.conjugate(s1[0]))**2
        monomials.append(lam_s1*mono + sp.conjugate(lam_s1*mono))
    return sp.expand(sum(monomials, sp.S.Zero))


def real_quartic_entries(*, include_s2: bool = True,
                         include_s1: bool = True) -> dict[tuple[int,int,int,int], sp.Expr]:
    """Exact symmetric fourth derivatives of V4; indices are one-based."""
    r = _real_fields()
    expression = extra_potential_polynomial(include_s2=include_s2,
                                            include_s1=include_s1)
    if expression == 0:
        return {}
    polynomial = sp.Poly(expression, *r)
    result = {}
    for exponents, coefficient in polynomial.terms():
        index = tuple(k+1 for k, power in enumerate(exponents) for _ in range(power))
        if len(index) != 4:
            raise ValueError('Unexpected degree in quartic potential')
        factor = sp.prod(sp.factorial(power) for power in exponents)
        result[index] = sp.simplify(coefficient*factor)
    return {key:val for key,val in result.items() if val != 0}


def extra_scalar_pair_component(quartic: Mapping[tuple[int,int,int,int],sp.Expr],
                                wilson, i: int, j: int, a: int, b: int) -> sp.Expr:
    """Exact extra scalar-pair term from the implemented master tensor RGE.

    Master equation: -1/2 Σ_{σ(a,b),c,d}(2 gauge - lambda_abcd) C_ijcd.
    Only its delta-lambda piece is evaluated here; no gauge/yukawa terms change.
    """
    total = sp.S.Zero
    for aa,bb in ((a,b),(b,a)):
        for c in range(1,13):
            for d in range(1,13):
                key = tuple(sorted((aa,bb,c,d)))
                tensor_value = quartic.get(key, sp.S.Zero)
                if tensor_value != 0:
                    total += tensor_value * wilson(i,j,c,d)/2
    return sp.factor(sp.expand(total))


def _expr_to_input(expr: sp.Expr) -> str:
    return sp.sstr(sp.factor(expr))


def _assert_linear_in_quartics(expr: sp.Expr) -> bool:
    """Treat barred and unbarred λ as independent formal variables."""
    atoms = []
    for name in ACTIVE_RGBETA:
        symbol = sp.Symbol(name)
        atoms.extend([symbol, sp.conjugate(symbol)])
    independent = sp.symbols('p0:4')
    p = sp.expand(expr.xreplace(dict(zip(atoms, independent))))
    poly = sp.Poly(p, *independent)
    return poly.total_degree() <= 1 and poly.coeff_monomial(1) == 0


def _classify_flavor(expr: sp.Expr) -> list[str]:
    names = []
    for symbol in sp.symbols('y1 y2'):
        if expr.has(symbol) or expr.has(sp.conjugate(symbol)):
            names.append(str(symbol))
    return names


def calculate_extra_eft1_mixing(model: str, seed: Mapping[str, Any],
                                eft1_rgbeta: Mapping[str, Any],
                                source_convention: Mapping[str, Any]) -> dict[str, Any]:
    """Export exact extra-quartic scalar-pair beta for Weinberg real components."""
    if model not in ('B_m1','C_m1'):
        raise ValueError('Only B/C alpha=-1 supported')
    if source_convention.get('status') != 'SourceLagrangianConventionDerived' or not source_convention.get('source_checks_passed'):
        raise ValueError('Missing certified Stage15 quartic convention')
    certified = {(x['model'], x['rgbeta_coupling']):x
                 for x in source_convention.get('results',[])}
    for coupling in ACTIVE_RGBETA:
        c=certified.get((model,coupling))
        if c is None or c.get('relative_sign') != -1 or c.get('requires_complex_conjugation') is not False:
            raise ValueError(f'Uncertified field/sign convention: {coupling}')
    if seed.get('Level') != 1 or seed.get('IntegratedFields') != ['F'] or seed.get('ActiveHeavyFields') != ['S1','S2']:
        raise ValueError('EFT1 tree Wilson seed is not F -> SM+S1+S2')
    meta = eft1_rgbeta.get('metadata',{})
    if eft1_rgbeta.get('status') != 'Success' or any(meta.get(k)!=v for k,v in {
        'dS1':2,'dS2':2,'dF':1 if model=='B_m1' else 3,'alpha':-1,
        'IntegratedField':'F'}.items()):
        raise ValueError('EFT1 RGBeta metadata disagreement')
    if 'report_betas' not in eft1_rgbeta or not all(x in eft1_rgbeta['report_betas'] for x in ACTIVE_RGBETA):
        raise ValueError('Missing scalar-EFT RG beta functions for active quartics')
    from RGE.running.intermediate.ScalarOnlyTensorAdapters import build_eft1_wilson_tensor
    wilson = build_eft1_wilson_tensor(seed,d_s1=2,d_s2=2,chirality='PL',
                                      lepton_indices=(1,2),project_operator_symmetry=True)
    if not wilson.entries:
        raise ValueError('No projected PL EFT1 Wilson tensor')
    quartic = real_quartic_entries()
    rows=[]
    for a,b in combinations_with_replacement(HIGGS_REAL_INDICES, 2):
        expr = extra_scalar_pair_component(quartic,wilson,1,2,a,b)
        if expr == 0:
            continue
        if not _assert_linear_in_quartics(expr):
            raise ValueError('Nonlinear quartic dependence in scalar-pair beta')
        rows.append({'higgs_real_indices':[a,b],
                     'beta_16pi2_expression':_expr_to_input(expr),
                     'yukawa_families':_classify_flavor(expr)})
    if not rows:
        raise ValueError('No nonzero Weinberg scalar-pair mixing found')
    all_yuk = set(sum((x['yukawa_families'] for x in rows),[]))
    if all_yuk != {'y1','y2'}:
        raise ValueError(f'Expected both y1^2 and y2^2 EFT1 Wilson sectors; got {all_yuk}')
    return {
        'status':'ExtraEFT1ScalarPairMixingDerived',
        'model':model,
        'source':'Matchete EFT1 PL tree seed plus RGBeta potential invariants and project master psi2phi2 tensor RGE',
        'operator':'PL two-lepton/two-real-scalar tensor C_ijab; i,j represent the two SU(2) components of one lepton generation',
        'quartic_potential':'lambda_abcd real scalar fourth derivatives of V4; V4=-L4',
        'source_convention':'lambda_Matchete=-lambda_RGBeta, same named scalar field basis',
        'master_rge_term':'delta_beta_ijab^(1) = (1/2) sum_swap(a,b),c,d delta_lambda_abcd C_ijcd',
        'quartic_real_nonzero_components':len(quartic),
        'pl_wilson_seed_terms':sum(t.get('Chirality')=='PL' for t in seed.get('TreeWilsonTerms',[])),
        'nonzero_weinberg_real_components':len(rows),
        'components':rows,
        'full_eft1_wilson_running_implemented':False,
        'scalar_threshold_extra_matching_implemented':False,
        'extra_low_energy_neutrino_mass':None,
        'complete_two_loop_result':False,
        'limitation':'Only extra quartic-dependent scalar-pair anomalous-dimension term. No full flavor Wilson ODE or scalar-threshold finite matching is included.',
    }


def run_file(model: str, seed: Path, eft: Path, convention: Path) -> dict:
    return calculate_extra_eft1_mixing(model,
        json.loads(Path(seed).read_text(encoding='utf-8-sig')),
        json.loads(Path(eft).read_text(encoding='utf-8-sig')),
        json.loads(Path(convention).read_text(encoding='utf-8-sig')))
