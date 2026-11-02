"""Strict same-input, same-scale UV/EFT numerical comparison; no EFT imports.

No comparison is made against a different scale, input point or convention.
The EFT expansion control is m_scalar^2/m_fermion^2, not exact equality.
"""
from __future__ import annotations
from typing import Mapping,Any
import numpy as np
from .core import _complex,takagi,encoded


def compare_same_point(direct: Mapping[str,Any], eft: Mapping[str,Any],
                       *, scale_rtol:float=2e-12) -> dict:
    if eft.get('physical_dimension')!='mass_matrix_gev':
        raise ValueError('EFT reference must specify physical_dimension=mass_matrix_gev')
    if eft.get('scheme') is None:
        raise ValueError('EFT reference missing scheme/convention')
    source=direct.get('source_provenance') or {}
    if not source.get('benchmark_sha256') or source['benchmark_sha256']!=eft.get('benchmark_sha256'):
        raise ValueError('Different or undocumented UV benchmark inputs')
    mu=float(direct['renormalization_scale_gev'])
    nu=float(eft['renormalization_scale_gev'])
    if abs(mu-nu)>scale_rtol*max(mu,nu):
        return {'status':'not_compared_different_scales',
                'uv_direct_scale_gev':mu,'eft_scale_gev':nu,
                'note':'Run both at the same scale with corresponding parameter values; do not require exact agreement'}
    if eft.get('same_renormalized_uv_parameters') is not True:
        return {'status':'not_compared_different_or_unknown_input_point'}
    uv=_complex(direct['total']['complex_symmetric_matrix_gev'],'direct matrix',(3,3))
    em=_complex(eft['matrix_gev'],'eft matrix',(3,3))
    _,um,_=takagi(uv)
    _,ee,_=takagi(em)
    denom=max(float(np.linalg.norm(uv)),float(np.linalg.norm(em)),1e-90)
    err=float(np.linalg.norm(uv-em)/denom)
    ratios=np.full(3,np.nan)
    for i in range(3):
        if abs(ee[i])>1e-20*max(float(max(ee)),1e-100):
            ratios[i]=float(um[i]/ee[i])
    return {'status':'same_scale_comparison_completed',
            'relative_mass_matrix_frobenius_difference':err,
            'direct_masses_ev':(um*1e9).tolist(),
            'eft_masses_ev':(ee*1e9).tolist(),
            'direct_over_eft_mass_ratios':[
                None if not np.isfinite(x) else float(x) for x in ratios],
            'eft_scheme':eft['scheme'],
            'max_ms2_over_mf2':direct['mass_hierarchy']['max_scalar_squared_to_fermion_squared_ratio'],
            'interpretation':'Discrepancies need not vanish when m_scalar^2/m_fermion^2 is not small; matching/order/scheme also matter.',
            'difference_matrix_gev':encoded(uv-em)}
