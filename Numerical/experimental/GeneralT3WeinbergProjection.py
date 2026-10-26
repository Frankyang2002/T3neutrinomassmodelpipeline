"""Stage 19: certified projection of the *selected* Stage 18 EFT1 LLHH source.

Three independent checks precede numerical transport:
 (1) exact SU(2) completion from the full Matchete PL seed agrees with
     the production Weinberg real-field tensor (including its /4 convention);
 (2) Stage17's mixed-lepton source agrees with the exact completed tensor;
 (3) the production direct-running C12->C5 coefficient agrees with its
     established Weinberg-tensor anomalous dimension.

Only the single-insertion extra-quartic source is transported. Homogeneous
EFT1 Wilson mixing, finite scalar matching, and complete two-loop physics are
not computed. This module does not alter any production matching artifact.
"""
from __future__ import annotations

from itertools import combinations_with_replacement
from typing import Any, Mapping

import numpy as np
import sympy as sp

from Numerical.experimental.GeneralT3ExtraEFT1Mixing import (
    extra_scalar_pair_component, real_quartic_entries,
)

HIGGS = (1, 2, 3, 4)
LEPTONS = (1, 2)
STAGE18_STATUS = 'ExtraEFT1LLHHTensorSourceIntegrated'


def _analytic_mixed_template() -> dict[tuple[int,int], complex]:
    """Production polynomial: kappa*(nu*H0-e*H+)**2, C=d^4/4.

    H+=(phi1+i phi2)/sqrt(2), H0=(phi3+i phi4)/sqrt(2).
    For lepton SU2 pair (nu,e), C_{12,ab} is -kappa/4 times
    (1,i,i,-1) on (13,14,23,24).
    """
    return {(1, 3): -.25+0j, (1, 4): -.25j,
            (2, 3): -.25j, (2, 4): .25+0j}


def _parse_scalar(expression: str) -> sp.Expr:
    names=('MF','y1','y2','lambdaT3','lambdaHHdagS2S2','lambdaHHdagS1barS1bar')
    symbols={n:sp.Symbol(n, nonzero=True) if n=='MF' else sp.Symbol(n) for n in names}
    return sp.sympify(expression, locals={**symbols,'conjugate':sp.conjugate,
                                         'sqrt':sp.sqrt,'log':sp.log, 'I':sp.I})


def verify_production_c5_normalization(wilson_rge: Mapping[str,Any],
                                       final_weinberg: Mapping[str,Any]) -> dict[str, Any]:
    """Cross-check established tensor beta against final direct-running C5.

    The production final adapter takes one-generation C12 = conjugate(y1)*
    conjugate(y2)/MF. The expression's log coefficient must match the
    Weinberg-subspace beta_kappa exactly (without its 1/16pi2).
    """
    info=wilson_rge.get('weinberg_subspace_validation',{})
    direct=final_weinberg.get('direct_running_contribution',{})
    if (wilson_rge.get('status')!='Success' or
        not info.get('matches_weinberg_subspace') or
        final_weinberg.get('status')!='Success' or
        not final_weinberg.get('combined',{}).get('ready_for_physical_majorana_numerics') or
        direct.get('manifestly_symmetric_under_pq') is not True):
        raise ValueError('Existing Weinberg tensor or final C5 certificate is unavailable')
    beta=info.get('beta_kappa_16pi2')
    raw=direct.get('raw_expression')
    if not isinstance(raw,str) or raw.count('C12[p,q]')!=1 or not isinstance(beta,str):
        raise ValueError('Unrecognised direct C12 running expression')
    if not ('log(MS/MF)' in raw):
        raise ValueError('Expected production matching log(MS/MF)')
    text=raw.replace('C12[p,q]','(conjugate(y1)*conjugate(y2)/MF)')
    diff=sp.cancel(sp.expand(_parse_scalar(text)).coeff(sp.log(sp.Symbol('MS')/sp.Symbol('MF',nonzero=True))) - _parse_scalar(beta))
    # The explicit log coefficient is independent of the matching log itself.
    if diff != 0:
        raise ValueError(f'Production direct C5 normalization mismatch: {diff}')
    return {'status':'ProductionDirectC5NormalizationVerified',
            'beta_kappa_16pi2':beta,
            'direct_running_log_coefficient_equals_beta_kappa':True}


def _production_template():
    """Read the actual one-generation Weinberg tensor from production source."""
    from RGE.matching.WeinbergTensorAdapter import build_sm_eft, build_weinberg_wilson_tensor
    scalars,fermions=build_sm_eft()
    if tuple(fermions.global_index('L',k) for k in LEPTONS)!=(1,2):
        raise ValueError('Unexpected lepton SU2 index ordering')
    if tuple(scalars.block('H').indices)!=HIGGS:
        raise ValueError('Unexpected Higgs real basis ordering')
    template=build_weinberg_wilson_tensor(scalars,fermions,sp.S.One)
    for pair,value in _analytic_mixed_template().items():
        if sp.simplify(template.get((1,2,*pair),0)-sp.sympify(value)) != 0:
            raise ValueError(f'Production Weinberg mixed tensor changed at {pair}')
    return template


def verify_full_seed_weinberg_projection(model:str, wilson_seed:Mapping[str,Any],
                                         stage17:Mapping[str,Any]) -> dict[str,Any]:
    """Check ALL 30 independent (L,L,H,H) components, not just the four mixed.

    This guards against extrapolating the four Stage18 norms into a physical
    Weinberg operator without checking its full SU2 multiplet.
    """
    if model not in ('B_m1','C_m1') or stage17.get('model')!=model or stage17.get('status')!='ExtraEFT1ScalarPairMixingDerived':
        raise ValueError('Missing successful matching Stage17 model certificate')
    from RGE.running.intermediate.ScalarOnlyTensorAdapters import build_eft1_wilson_tensor
    wilson=build_eft1_wilson_tensor(wilson_seed,d_s1=2,d_s2=2,chirality='PL',
                                    lepton_indices=LEPTONS,project_operator_symmetry=True)
    template=_production_template()
    quartic=real_quartic_entries()
    keys=[(i,j,a,b) for i,j in combinations_with_replacement(LEPTONS,2)
          for a,b in combinations_with_replacement(HIGGS,2)]
    data={k:extra_scalar_pair_component(quartic,wilson,*k) for k in keys}
    reference=(1,2,1,3)
    if sp.simplify(template.get(reference,0))==0 or sp.simplify(data[reference])==0:
        raise ValueError('Reference extra Weinberg component absent')
    beta=sp.factor(sp.cancel(data[reference]/template[reference]))
    mismatches={}
    for key in keys:
        residual=sp.cancel(sp.expand(data[key]-beta*template.get(key,0)))
        if sp.simplify(residual)!=0:
            mismatches[','.join(map(str,key))]=sp.sstr(residual)
    if mismatches:
        raise ValueError('Extra EFT1 scalar-pair beta NOT pure Weinberg: '+str(mismatches))
    actual_rows={tuple(row['higgs_real_indices']):_parse_scalar(row['beta_16pi2_expression'])
                 for row in stage17.get('components',[])}
    if len(actual_rows)!=4 or set(actual_rows)!=set(_analytic_mixed_template()):
        raise ValueError('Stage17 mixed-lepton components do not match Weinberg support')
    for pair,expr in actual_rows.items():
        if sp.simplify(expr-data[(1,2,*pair)])!=0:
            raise ValueError(f'Stage17 tensor component disagrees with full seed at {pair}')
    return {'status':'FullExtraWeinbergTensorProjectionVerified',
            'model':model,'independent_components_verified':len(keys),
            'mixed_components_matched_to_stage17':4,
            'beta_kappa_16pi2_expression':sp.sstr(beta),
            'all_component_residuals_zero':True,
            'template':'RGE.matching.WeinbergTensorAdapter.build_weinberg_wilson_tensor(kappa=1)',
            'normalization':'C_ijab=(1/4)*partial_i partial_j partial_a partial_b [kappa*(nu*H0-e*H+)**2]'}


def _decode_matrix(payload:Any) -> np.ndarray:
    if not isinstance(payload,Mapping):
        raise ValueError('Missing Wilson component matrix')
    m=np.asarray(payload['real'],dtype=float)+1j*np.asarray(payload['imag'],dtype=float)
    if m.shape!=(3,3) or not np.isfinite(m).all() or not np.allclose(m,m.T,rtol=1e-10,atol=1e-40):
        raise ValueError('Invalid complex-symmetric 3x3 Wilson matrix')
    return m


def project_stage18_integral(stage18:Mapping[str,Any], model:str,
                             certificate:Mapping[str,Any],
                             normalization:Mapping[str,Any]) -> tuple[np.ndarray,dict[str,Any]]:
    """Extract *selected-source* delta K by SU2-template least squares.

    Need all four complex components and an independent full 30-component
    symbolic certificate. The normalized projection is based on the actual
    production Weinberg tensor; the result is NOT a complete matching C5.
    """
    if (stage18.get('status')!=STAGE18_STATUS or stage18.get('model')!=model or
        stage18.get('extra_source_only') is not True or
        stage18.get('eft1_wilson_seed_frozen_at_mu_F') is not True or
        stage18.get('full_eft1_wilson_anomalous_dimension_implemented') is not False):
        raise ValueError('Unrecognised Stage18 single-insertion source')
    if (certificate.get('status')!='FullExtraWeinbergTensorProjectionVerified' or
        certificate.get('model')!=model or
        certificate.get('independent_components_verified')!=30 or
        certificate.get('all_component_residuals_zero') is not True or
        normalization.get('status')!='ProductionDirectC5NormalizationVerified'):
        raise ValueError('Full-tensor and physical normalization certificates required')
    rows=stage18.get('higgs_real_components',[])
    if len(rows)!=4:
        raise ValueError('Expected exactly four mixed Higgs components')
    observed={}
    for row in rows:
        key=tuple(row['higgs_real_indices'])
        if key in observed:
            raise ValueError('Repeated Higgs real component')
        observed[key]=_decode_matrix(row['integrated_delta_C_ijab_gev_inverse'])
    template=_analytic_mixed_template()
    if set(observed)!=set(template):
        raise ValueError('Stage18 Higgs support differs from certified Weinberg tensor')
    denom=sum(abs(t)**2 for t in template.values())
    delta=sum(np.conj(t)*observed[k] for k,t in template.items())/denom
    delta=(delta+delta.T)/2
    total=sum(np.linalg.norm(v,'fro')**2 for v in observed.values())**.5
    residual=sum(np.linalg.norm(observed[k]-t*delta,'fro')**2 for k,t in template.items())**.5
    relative=residual/max(total,1e-300)
    if relative>2e-8 or (total>0 and residual>1e-38+2e-8*total):
        raise ValueError(f'Integrated tensor failed complex Weinberg phase/sign test: relative={relative:.6g}')
    return delta, {'status':'IntegratedExtraWeinbergSourceProjected',
                   'model':model,'component_count':4,
                   'projection_relative_residual':float(relative),
                   'extra_delta_kappa_norm_gev_inverse':float(np.linalg.norm(delta,'fro')),
                   'loop_order_statement':'one-loop LLSS-source insertion; O(hbar^2) when extra quartics are radiatively generated from zero',
                   'scalar_threshold_projection':'tree-level LLHH survives removing S1/S2, finite scalar matching omitted',
                   'complete_eft1_wilson_running':False,
                   'complete_two_loop_prediction':False}


def encode_matrix(m: np.ndarray, *, unit: str) -> dict[str,Any]:
    z=np.asarray(m,dtype=complex)
    if z.shape!=(3,3) or not np.isfinite(z).all():
        raise ValueError('Nonfinite or incorrectly shaped output matrix')
    return {'real':z.real.tolist(),'imag':z.imag.tolist(),
            'frobenius_norm':float(np.linalg.norm(z,'fro')),'unit':unit}
