"""Batch 5: explicit one-loop massive UV triangle / EFT1-soft matching ledger.

This calculates a genuine subtraction identity for the *already implemented*
heavy-fermion--two-scalar triangle, at zero external momentum:

    C_UV_full = C_hard_dim5(mu) + C_EFT1_soft_MSbar(mu) + C_dim7plus.

The last term is power suppressed, NOT a higher-loop term. The hard and soft
terms individually depend on the arbitrary renormalization scale mu, but their
sum does not at one loop with couplings held fixed. All three fermions and the
complete scalar mass matrix are retained; no scalar is removed at MF.

Separately, reconstruct the raw native Matchete PL Laurent pole/constant/log
coefficient from the Batch-4 numerical results. A deliberately hypothetical
native-as-bare-UV subtraction is returned for comparison only. Matchete's
saved MatchedEFT output is not established as a bare UV amplitude, so this is
NOT a certified physical one-loop hard coefficient or an MSbar determination.

No arbitrary counterterm or missing amplitude is assigned a value.
"""
from __future__ import annotations

from math import isfinite, log, pi
import numpy as np

SCHEMA='t3-f-triangle-UV-minus-EFT1-soft-and-native-pole-ledger-batch5-v1'
LOOP=1/(16*pi*pi)
SCALE_MULTIPLIERS=(0.5,1.,2.)


def _rel(a,b):
    return float(np.linalg.norm(np.asarray(a)-np.asarray(b)) /
                 max(np.linalg.norm(np.asarray(b)),1e-100))


def _finite_symmetric(a,name):
    x=np.asarray(a,complex)
    if x.shape!=(3,3) or not np.isfinite(x).all() or not np.allclose(
            x,x.T,rtol=1e-10,atol=1e-24):
        raise ValueError(name+': require finite complex symmetric 3x3')
    return (x+x.T)/2


def native_laurent_from_evaluated(per_kernel_and_fermion:np.ndarray,
                                   masses:np.ndarray,mu:float)->dict:
    """Undo the 1+ln(mu^2/Mk^2) multiplier on each native PL source kernel.

    Batch 4's native per-kernel numbers already contain one hbar factor.
    The exact pole=constant=log kernel identities are verified in Batch 4;
    this is only their *numerical reconstruction* and does not drop poles.
    """
    p=np.asarray(per_kernel_and_fermion,complex)
    m=np.asarray(masses,float)
    if (p.ndim!=4 or p.shape[1:]!=(3,3,3) or p.shape[0]<1 or
            m.shape!=(3,) or not np.isfinite(p).all() or
            not np.isfinite(m).all() or np.min(m)<=0 or not isfinite(mu) or mu<=0):
        raise ValueError('Unexpected native kernel/flavour/fermion geometry')
    factor=1+2*np.log(mu/m)
    if np.min(abs(factor))<=0.05:
        raise ValueError('Cannot invert a nearly vanishing native finite multiplier')
    k=p/factor[np.newaxis,:,np.newaxis,np.newaxis]
    poles=k.sum(axis=(0,1))
    logs=np.einsum('kfij,f->ij',k,2*np.log(mu/m),optimize=True)
    finite=p.sum(axis=(0,1))
    if _rel(poles+logs,finite)>1e-12:
        raise ValueError('Native Laurent pole/constant/log decomposition failed')
    for z in (poles,logs,finite):_finite_symmetric(z,'native Laurent')
    return {'per_native_kernel_and_fermion_constant':k,
            'pole_residue':poles,'finite_constant':poles.copy(),
            'finite_scale_log':logs,'formal_pole_dropped':finite,
            'lambda_epsilon_pole_removed_physically':False,
            'native_finite_is_MSbar_hard':False}


def triangle_soft_ledger(c_by_fermion:np.ndarray,masses:np.ndarray,
                         quartic_hhss:np.ndarray,mass_squared:np.ndarray,
                         matching_scale:float=1e5,*,multipliers=SCALE_MULTIPLIERS,
                         projection_tolerance:float=1e-7)->dict:
    """Evaluate exact UV and EFT1 soft loop in the SAME tensor convention.

    Existing hard_matching.py defines w=C_LLSS lambda/(2*16pi2),
    kh=1+ln(mu^2/Mk^2), ks=-B0bar(m_a2,m_b2;mu),
    k_full=triangle_ratio(Mk^2,m_a2,m_b2).

    Exact UV minus renormalized EFT1 soft includes higher-dimensional terms.
    The mass-independent local hard dim5 coefficient is *not* replaced by them.
    Changing mu here changes only the SUBTRACTION scale in the kernel; it is
    not a second heavy-fermion threshold and it does not rerun couplings.
    """
    from Numerical.lo_certification.completeness.hard_matching import (
        evaluate_hard_scalar_channel,weinberg_tensor)
    from Numerical.parallel_t3.threshold.matching import master_weinberg_template

    if not np.allclose(weinberg_tensor(),2*master_weinberg_template(),atol=1e-14):
        raise ValueError('Hard-to-master Weinberg normalization changed')
    mf=float(matching_scale)
    if not isfinite(mf) or mf<=0:raise ValueError('Invalid physical fermion threshold')
    multipliers=tuple(float(x) for x in multipliers)
    if (not multipliers or len(set(multipliers))!=len(multipliers) or
            any(not isfinite(x) or x<=0 for x in multipliers) or 1. not in multipliers):
        raise ValueError('Need unique positive scales including common MF')
    from_c = np.asarray(c_by_fermion,complex)
    masses=np.asarray(masses,float)
    quartic_hhss=np.asarray(quartic_hhss,float)
    mass_squared=np.asarray(mass_squared,float)
    if from_c.ndim!=7 or from_c.shape[:5]!=(3,3,3,2,2) or masses.shape!=(3,):
        raise ValueError('Exactly three heavy mass eigenstates and full LLSS tensor required')
    n=from_c.shape[-1]
    if (from_c.shape[-2]!=n or quartic_hhss.shape!=(4,4,n,n) or
            mass_squared.shape!=(n,n) or not np.isfinite(from_c).all()):
        raise ValueError('Complete heavy scalar dimensions missing')
    result={}
    for factor in multipliers:
        scale=mf*factor
        # Equal mu_f,mu_s in this call selects a single common logarithmic
        # renormalization scale; the scalar masses remain the IR regulator.
        r=evaluate_hard_scalar_channel(
            from_c,masses,quartic_hhss,mass_squared,scale,scale,
            C_total_muF=from_c.sum(axis=0),projection_tolerance=projection_tolerance)
        if not r['projection_passed']:
            raise ValueError(f'Full triangle Weinberg projection failed at scale factor {factor}: {r["residuals"]}')
        mat={key:2*np.asarray(a,complex) for key,a in r['matrices'].items()}
        for key in mat:mat[key]=_finite_symmetric(mat[key],key)
        if _rel(mat['hard_muF']+mat['scalar_threshold_muS']+
                mat['higher_dimension_remainder'],mat['full_triangle'])>1e-8:
            raise ValueError('Analytic full-minus-soft exact triangle closure failed')
        # C_hard_full_minus_soft is a finite mass-dependent object and is not
        # matched as a local dimension-five coefficient.
        pole=np.zeros((3,3),complex)
        for row in r['per_pair_terms']:
            z=row['projected_C5_base_not_independently_gauge_invariant']
            pole += 2*(np.asarray(z['real'])+1j*np.asarray(z['imag']))
        pole=_finite_symmetric(pole,'EFT1 scalar bubble pole base')
        # EFT1 bubble = -w * B0, whose UV pole is -w/epsilon_bar.
        if abs(factor-1.)<1e-14:
            base_pole=pole
            terms=r['per_pair_terms']
        result[factor]={
            'scale_gev':scale,
            'hard_local_dim5':mat['hard_muF'],
            'EFT1_soft_finite':mat['scalar_threshold_muS'],
            'full_UV_triangle':mat['full_triangle'],
            'exact_UV_minus_soft':mat['full_triangle']-mat['scalar_threshold_muS'],
            'power_suppressed_dim7plus_remainder':mat['higher_dimension_remainder'],
            'hard_local_plus_soft':mat['hard_muF']+mat['scalar_threshold_muS'],
            'scalar_RGE_log_at_equal_scales':mat['eft1_scalar_log'],
            'eft_soft_1over_epsilon_pole':-pole,
            'projection_residuals':r['residuals'],
        }
    ref=result[1.]
    max_scale_variation=max(_rel(s['hard_local_plus_soft'],ref['hard_local_plus_soft'])
                             for s in result.values())
    max_full_variation=max(_rel(s['full_UV_triangle'],ref['full_UV_triangle'])
                            for s in result.values())
    max_exact_variation=max(_rel(s['exact_UV_minus_soft']+s['EFT1_soft_finite'],
                             ref['full_UV_triangle']) for s in result.values())
    if max(max_scale_variation,max_full_variation,max_exact_variation)>2e-8:
        raise ValueError('UV/soft matching scale cancellation failed')
    # Ratio of dimension-7+ corrections is useful but *never* used to change C5.
    ratio=_rel(ref['power_suppressed_dim7plus_remainder'],ref['full_UV_triangle'])
    return {'schema':SCHEMA,'scales':result,
            'matching_scale_cancellation_relative':max_scale_variation,
            'full_UV_scale_invariance_relative':max_full_variation,
            'UV_minus_EFT1_soft_exact_closure_relative':max_exact_variation,
            'dimension7plus_relative_to_full_triangle':ratio,
            'eft_soft_pole_at_MF':-base_pole,
            'fermion_count':len(masses),
            'heavy_scalar_real_count':n,
            'scalar_pairs_evaluated':len(terms),
            'mass_eigenvalues_squared_gev2':np.linalg.eigvalsh(mass_squared),
            'source_one_loop_formula_verified':'hard_matching.py original full massive triangle + B0bar',
            'massive_remainder_inserted_in_dim5':False,
            'full_one_loop_matching_certified':False}


def comparison_with_native(triangle:dict,native:dict)->dict:
    """Numerical Laurent-pole comparison; never assert native UV status."""
    r=triangle['scales'][1.]
    n=_finite_symmetric(native['formal_pole_dropped'],'native formal coefficient')
    pole=_finite_symmetric(native['pole_residue'],'native raw pole')
    eft=_finite_symmetric(triangle['eft_soft_pole_at_MF'],'EFT1 soft pole')
    tri=_finite_symmetric(r['hard_local_dim5'],'analytic hard local')
    soft=_finite_symmetric(r['EFT1_soft_finite'],'EFT1 soft finite')
    hypo=_finite_symmetric(n-soft,'hypothetical native UV minus EFT soft')
    residual=pole+eft
    # This pole test is informative, NOT a mandatory cancellation: the source
    # is Matchete MatchedEFT and its counterterms are not reconstructed.
    return {'hypothetical_hard_MF':hypo,
            'native_finite_formal_MF':n,
            'EFT1_soft_finite_at_MF':soft,
            'native_raw_pole_residue':pole,
            'analytic_EFT1_soft_pole_residue':eft,
            'unaccounted_pole_or_counterterm_residual':residual,
            'native_vs_analytic_hard_relative':_rel(n,tri),
            'pole_cancellation_relative':float(np.linalg.norm(residual)/max(
                np.linalg.norm(pole),np.linalg.norm(eft),1e-100)),
            'hypothesis_not_certified_reason':(
                'Matchete F-only MatchedEFT raw Laurent output is not established '
                'as an unsubtracted UV amplitude; EFT1 soft and UV counterterms '
                'must be independently identified before this subtraction is physical'),
            'hypothesis_is_MSbar_hard':False,
            'counterterms_applied':False,
            'missing_hard_operator_topologies_completed':False}
