"""Source-resolved one-loop perturbative bookkeeping for all 16 T3 models.

A genuinely first-order EFT1 endpoint uses the *order-zero* MF LLSS tensor
and the seven-channel candidate anomalous dimension evaluated at MF:

  C_strict(MS) = C_tree(MF) + C_F_hard^1(MF)
                 + log(MS/MF)/(16 pi^2) beta[C_tree(MF); g(MF)].

The hard boundary is NOT evolved at this order. C5 generation is one part of
full C evolution and does not receive an independently inserted leading log.
Existing Workstream C uses the order-zero LLSS tensor for its scalar bubble,
not the RGE-evolved LLSS, in strict_one_loop mode.

This is a first-order calculation for the INCLUDED one-loop contributions:
other first-threshold hard topologies, scalar matching topologies, full UV
parameter one-loop threshold corrections and operator completion are missing.
"""
from __future__ import annotations

import math
from typing import Mapping

import numpy as np

LOOP = 16.0 * math.pi**2
DEFAULT_ADM = 'fonseca_eq22_six_term'
SCOPE = 'strict_one_loop_included_topologies_not_complete_physical_matching'


def first_order_endpoint(tree: np.ndarray, hard: np.ndarray,
                         beta_channels: Mapping[str,np.ndarray], log_ratio: float
                         ) -> tuple[np.ndarray,dict[str,np.ndarray]]:
    """Pure full-tensor O(hbar) composition; never evolve O(hbar) hard once more.

    beta_channels are 16*pi**2*dC/dln(mu) evaluated ONLY at tree level.
    All arrays have the same tensor shape; no operator is projected/discarded.
    """
    tree=np.asarray(tree,complex)
    hard=np.asarray(hard,complex)
    if tree.ndim != 4 or tree.shape!=hard.shape or tree.shape[0]!=tree.shape[1] or tree.shape[2]!=tree.shape[3]:
        raise ValueError('Expected matching full psi2phi2 tree and hard tensors')
    if not np.isfinite(log_ratio) or not log_ratio < 0:
        raise ValueError('Need finite descending EFT1 logarithm ln(MS/MF)<0')
    if not np.isfinite(tree).all() or not np.isfinite(hard).all():
        raise ValueError('Nonfinite tree or hard Wilson')
    if not beta_channels:
        raise ValueError('Missing the complete candidate ADM channel inventory')
    changes={}
    for name, arr in beta_channels.items():
        if not isinstance(name,str) or not name or name in changes:
            raise ValueError('Invalid/duplicate beta channel name')
        beta=np.asarray(arr,complex)
        if beta.shape!=tree.shape or not np.isfinite(beta).all():
            raise ValueError(f'Incorrect ADM channel {name}')
        changes[name]=(log_ratio/LOOP)*beta
    total=sum(changes.values(), np.zeros_like(tree))
    strict=tree+hard+total
    if not np.isfinite(strict).all():
        raise FloatingPointError('Nonfinite strict one-loop Wilson')
    return strict,changes


def _norm(x:np.ndarray)->float:
    return float(np.linalg.norm(x))


def construct_strict_eft1(case, hard_c5:np.ndarray, *,
                          adm_variant:str=DEFAULT_ADM) -> tuple[dict,dict]:
    """Construct the one-loop-truncated full EFT1 Wilson endpoint.

    case is the *actual* Batch-1 MF input with its sourced hard triangle.
    The tree-level coefficient is reconstructed by removing ONLY its known
    Weinberg seed, leaving every LLSS component unchanged.
    """
    from Numerical.parallel_t3.wilson.adm import beta, SYMMETRIC_TERMS
    from Numerical.parallel_t3.wilson.coordinates import (
        embed_c5, project_c5, symmetry_residual)
    from Numerical.parallel_t3.wilson.ward_singlet_v7 import _gauge_norm

    m_f=float(case.mf)
    m_s=float(case.ms)
    if not 0<m_s<m_f:
        raise ValueError('Incorrect EFT order MF>MS')
    hard_c5=np.asarray(hard_c5,complex)
    if hard_c5.shape!=(3,3) or not np.isfinite(hard_c5).all():
        raise ValueError('Hard MF contribution must be finite 3x3')
    if not np.allclose(hard_c5,hard_c5.T,rtol=1e-10,atol=1e-25):
        raise ValueError('Hard Weinberg coefficient must be complex symmetric')
    c_f=np.asarray(case.initial,complex)
    h=embed_c5(hard_c5,case.basis)
    c0=c_f-h
    c5_before,_=project_c5(c_f,case.basis)
    c5_tree,_=project_c5(c0,case.basis)
    if _norm(c5_before-hard_c5)>max(1e-23,1e-9*_norm(hard_c5)):
        raise ValueError('Batch-1 input MF hard boundary does not match stored hard C5')
    if _norm(c5_tree)>max(1e-23,1e-9*_norm(hard_c5)):
        raise ValueError('Reconstructed MF tree has a separate unrecognized Weinberg seed')
    if symmetry_residual(c0)>max(1e-24,1e-10*np.max(abs(c0))):
        raise ValueError('Reconstructed tree violates Wilson exchange symmetry')
    args=case.background.at(math.log(m_f))
    parts,backend=beta(c0,*args,variant=adm_variant)
    if set(parts)!=set(SYMMETRIC_TERMS):
        raise ValueError('ADM channel inventory changed')
    strict,changes=first_order_endpoint(c0,h,parts,math.log(m_s/m_f))
    gauges=_gauge_norm(strict,case.basis)
    if max(gauges.values())>1e-7:
        raise ValueError('Strict first-order endpoint violates gauge Ward identity '+str(gauges))
    c5_run=np.zeros((3,3),complex)
    channels={}
    for key in SYMMETRIC_TERMS:
        value,_res=project_c5(changes[key],case.basis)
        channels[key]=value
        c5_run+=value
    c5_strict,projection_residual=project_c5(strict,case.basis)
    c5_pred=c5_tree+hard_c5+c5_run
    if _norm(c5_strict-c5_pred)>max(1e-23,1e-9*_norm(c5_strict)):
        raise ValueError('Strict one-loop projection differs from hard + seven-channel running')
    data={'C_psi2phi2_tree0_mf':c0,
          'C_psi2phi2_hard_one_loop_mf':h,
          'C_psi2phi2_ms_strict_1loop':strict,
          'C5_MF_hard_fermion_triangle':hard_c5,
          'C5_MF_tree0':c5_tree,
          'C5_EFT1_strict_running_delta':c5_run,
          'C5_EFT1_MS_strict':c5_strict}
    for key,val in channels.items():
        data['C5_EFT1_strict_running_'+key]=val
    manifest={
        'mode':'strict_one_loop', 'scope':SCOPE,
        'mf_gev':m_f,'ms_gev':m_s,'log_ms_over_mf':math.log(m_s/m_f),
        'beta_model':adm_variant,'beta_backend':backend,
        'formula':'C0_MF + C_hard1_MF + log(MS/MF)/(16pi2)*beta(C0_MF)',
        'C0_is_MF_LLSS_from_B_without_one_loop_hard':True,
        'F_hard_C5_included_once':True,
        'F_hard_not_evolved_as_order_zero':True,
        'EFT1_renormalisable_running_background_order':'MF_frozen_for_strict_first_order_Wilson',
        'SMEFT_running_after_scalar_threshold':'existing_RGE_improved_SM_Weinberg',
        'C5_strict_vs_explicit_terms_relative':_norm(c5_strict-c5_pred)/max(_norm(c5_strict),1e-100),
        'C5_projection_residual':float(projection_residual),
        'ward_relative':gauges,
        'C5_norms':{k:_norm(data[k]) for k in (
            'C5_MF_hard_fermion_triangle','C5_EFT1_strict_running_delta','C5_EFT1_MS_strict')},
        'incomplete_finite_F_matching':True,'incomplete_finite_S_matching':True,
        'unverified_full_ADM':True,'physical_C5_certified':False}
    return data,manifest
