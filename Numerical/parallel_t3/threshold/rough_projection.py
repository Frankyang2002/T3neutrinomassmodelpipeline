"""Source-preserving Weinberg-purity and SU(2)xU(1) Ward diagnostics.

This is *not* a repaired EFT1 beta function or a new matching prescription.
Only read the full evolved Wilson coefficient; do not replace its LLHH entries.
Every candidate reindexing is diagnostic ONLY and may not be applied without an
independently established field map/provenance.
"""
from __future__ import annotations

import json
from pathlib import Path
import numpy as np


FORMAT = 't3-rough-llhh-source-projection-diagnostic-v1'
HIGGS = ('H[1,R]', 'H[1,I]', 'H[2,R]', 'H[2,I]')


def _complex(x: np.ndarray | complex) -> dict:
    a = np.asarray(x, dtype=np.complex128)
    return {'real': a.real.tolist(), 'imag': a.imag.tolist()}


def _weyl_doublet_generators() -> dict[str, np.ndarray]:
    return {
        'SU2_T1': np.array([[0., .5], [.5, 0.]], complex),
        'SU2_T2': np.array([[0., -.5j], [.5j, 0.]], complex),
        'SU2_T3': np.diag([.5, -.5]).astype(complex),
        'U1_Y': -.5*np.eye(2, dtype=complex),
    }


def _real_higgs_generator(generator: np.ndarray) -> np.ndarray:
    """Same interleaved-real convention as wilson.basis.real_scalar_generator."""
    g = np.zeros((4, 4), dtype=np.complex128)
    for i in range(2):
        for j in range(2):
            a, b = generator[i,j].real, generator[i,j].imag
            g[2*i:2*i+2, 2*j:2*j+2] = 1j*np.array([[b,a],[-a,b]])
    return g


def _ward(a: np.ndarray, f: np.ndarray, s: np.ndarray) -> np.ndarray:
    """Coefficient variation, with C_ijAB psi_i psi_j H_A H_B (all sums)."""
    return (np.einsum('ki,pqkjAB->pqijAB', f, a, optimize=True)
          + np.einsum('kj,pqikAB->pqijAB', f, a, optimize=True)
          + np.einsum('CA,pqijCB->pqijAB', s, a, optimize=True)
          + np.einsum('CB,pqijAC->pqijAB', s, a, optimize=True))


def _project(a: np.ndarray, template: np.ndarray) -> tuple[np.ndarray, float, np.ndarray]:
    c5 = np.einsum('pqijAB,ijAB->pq', a, template.conj(), optimize=True) / np.vdot(template,template)
    c5 = (c5 + c5.T) / 2
    reconstructed = np.einsum('pq,ijAB->pqijAB',c5,template,optimize=True)
    error = a - reconstructed
    return c5, float(np.linalg.norm(error)/max(float(np.linalg.norm(a)),1e-100)), error


def diagnose_llhh(a: np.ndarray, template: np.ndarray, *,
                  declared_c5: np.ndarray | None=None,
                  top_components: int=12) -> dict:
    """Quantify the exact C Weinberg-template residual before any projection.

    a[p,q,weak_i,weak_j,H_A,H_B] has shape (3,3,2,2,4,4).
    The caller must extract a from the full *actual* A tensor using its
    explicit lepton/Higgs indices (not a newly assumed position list).
    """
    a = np.asarray(a, dtype=np.complex128)
    template = np.asarray(template, dtype=np.complex128)
    if a.shape != (3,3,2,2,4,4) or template.shape != (2,2,4,4):
        raise ValueError('Expected source LLHH (3,3,2,2,4,4) and C template (2,2,4,4)')
    if not np.isfinite(a).all() or not np.isfinite(template).all() or not np.linalg.norm(template):
        raise ValueError('Nonfinite/zero source or template')
    if not np.allclose(a,a.transpose(1,0,3,2,4,5),rtol=1e-9,atol=1e-100):
        # The existing C preflight owns hard rejection of inconsistent tensors.
        raise ValueError('Full LLHH Weyl-pair symmetry violated')
    if not np.allclose(a,a.swapaxes(-1,-2),rtol=1e-9,atol=1e-100):
        raise ValueError('Full LLHH real-scalar-pair symmetry violated')
    c5,residual,error = _project(a,template)
    c5_component = 4*a[:,:,0,0,2,2]
    norm = float(np.linalg.norm(a))
    strength = float(np.linalg.norm(c5))
    template_norm = float(np.linalg.norm(template))
    wards = {}
    for name, gl in _weyl_doublet_generators().items():
        gs = _real_higgs_generator(.5*np.eye(2) if name=='U1_Y' else gl)
        ward = _ward(a,gl,gs)
        wards[name] = {'relative_to_LLHH':float(np.linalg.norm(ward)/max(norm,1e-100)),
                       'absolute_gev_inverse':float(np.linalg.norm(ward))}
    # These are NOT permitted field transformations: they are ranked solely
    # to identify a systematic index/sign convention discrepancy for review.
    weak_flip = [1,0]
    h_flip = [2,3,0,1]
    imag_flip = np.array([1,-1,1,-1],float)
    t_weak = np.take(np.take(template,weak_flip,axis=0),weak_flip,axis=1)
    t_higgs = np.take(np.take(template,h_flip,axis=2),h_flip,axis=3)
    alternatives = {
        'declared_C_Weinberg_template':template,
        'H_imaginary_components_conjugated':(template*imag_flip[None,None,:,None]*imag_flip[None,None,None,:]),
        'both_lepton_weak_components_swapped':t_weak,
        'H_doublet_components_swapped':t_higgs,
        'lepton_and_H_doublet_components_swapped':np.take(np.take(t_weak,h_flip,axis=2),h_flip,axis=3),
    }
    trials = sorted(({'hypothesis_only':name, 'relative_residual':_project(a,t)[1]}
                     for name,t in alternatives.items()),key=lambda d:d['relative_residual'])
    most=[]
    for flat in np.argsort(np.abs(error).ravel())[::-1][:max(0,int(top_components))]:
        p,q,i,j,A,B = (int(v) for v in np.unravel_index(int(flat),a.shape))
        v = a[p,q,i,j,A,B]
        diff = error[p,q,i,j,A,B]
        most.append({'flavours_zero_based':[p,q], 'weak_one_based':[i+1,j+1],
                     'higgs_real_labels':[HIGGS[A],HIGGS[B]],
                     'source_coefficient_gev_inverse':_complex(v),
                     'difference_from_Weinberg_gev_inverse':_complex(diff),
                     'absolute_difference_gev_inverse':float(abs(diff))})
    ans = {
        'format':FORMAT,
        'status':'SOURCE_DIAGNOSTIC_NOT_MATCHING',
        'source_LLHH_norm_gev_inverse':norm,
        'source_LLHH_maxabs_gev_inverse':float(np.max(abs(a))),
        'projected_C5_norm_gev_inverse':strength,
        'fitted_Weinberg_component_norm_gev_inverse':strength*template_norm,
        'residual_absolute_norm_gev_inverse':float(np.linalg.norm(error)),
        'relative_Weinberg_purity_residual':residual,
        'pure_Weinberg_at_tolerance_1e-8':bool(residual <= 1e-8),
        'C5_least_squares_projection_gev_inverse':_complex(c5),
        'C5_neutrino_neutral_component_gev_inverse':_complex(c5_component),
        'C5_component_vs_LS_relative_to_larger':float(np.linalg.norm(c5_component-c5)/max(
            strength,float(np.linalg.norm(c5_component)),1e-100)),
        'ward_SU2xU1':wards,
        'alternative_coordinate_hypotheses_DO_NOT_APPLY':trials,
        'largest_nonWeinberg_components':most,
        'interpretation':(
            'A large residual can indicate field-index conventions or an EFT1 tensor '
            'that is not gauge-invariant in the stated basis. Alternative hypotheses '
            'are diagnostics only. Do not drop non-Weinberg components silently.'),
    }
    if declared_c5 is not None:
        supplied = np.asarray(declared_c5,dtype=np.complex128)
        if supplied.shape != (3,3) or not np.isfinite(supplied).all():
            raise ValueError('Supplied C5_EFT1_ms is not a finite 3x3 matrix')
        ans['stream2_separately_declared_C5'] = {
            'matrix_gev_inverse':_complex(supplied),
            'relative_to_C_projector':float(np.linalg.norm(supplied-c5)/max(
                float(np.linalg.norm(supplied)),strength,1e-100)),
            'relative_to_neutral_component':float(np.linalg.norm(supplied-c5_component)/max(
                float(np.linalg.norm(supplied)),float(np.linalg.norm(c5_component)),1e-100))}
    return ans


def diagnose_pair(wilson_manifest:Path, *, output:Path, model:str) -> dict:
    """Inspect adapted genuine A endpoint without changing a single coefficient."""
    from .matching import extract_lepton_sectors,master_weinberg_template
    from .models import n_real, normalize_model
    from .pipeline import _load_pair,W_FORMAT

    model = normalize_model(model)
    w,wa,_,_ = _load_pair(Path(wilson_manifest),W_FORMAT)
    if normalize_model(w.get('model','')) != model:
        raise ValueError('Mismatch between diagnostic model and A source model')
    hh,_ = extract_lepton_sectors(wa['C_psi2phi2_ms'],
        np.asarray(w['lepton_indices'],dtype=int),
        np.asarray(w['higgs_source_indices'],dtype=int),
        np.asarray(w['heavy_source_indices'],dtype=int),
        number_of_real_scalars=n_real(model)+4)
    result = diagnose_llhh(hh, master_weinberg_template(),
                           declared_c5=wa.get('C5_EFT1_ms'))
    # Source-available order-zero MF comparison: this is NOT a new LLSS log
    # and is not a substitute for the full one-loop beta history.
    mf_full = wa.get('C_psi2phi2_tree_mf')
    if mf_full is None:
        result['source_MF_orderzero_witness'] = {
            'availability':'not_in_Stream2_endpoint',
            'interpretation':'Cannot localise when Weinberg purity was lost from MS data alone'}
    else:
        mf_hh, _ = extract_lepton_sectors(mf_full,
            np.asarray(w['lepton_indices'],dtype=int),
            np.asarray(w['higgs_source_indices'],dtype=int),
            np.asarray(w['heavy_source_indices'],dtype=int),
            number_of_real_scalars=n_real(model)+4)
        mf = diagnose_llhh(mf_hh, master_weinberg_template())
        mf_norm = mf['source_LLHH_norm_gev_inverse']
        result['source_MF_orderzero_witness'] = {
            'availability':'C_psi2phi2_tree_mf',
            'source_LLHH_norm_gev_inverse':mf_norm,
            'zero_source_no_relative_purity_information':mf_norm<1e-100,
            'relative_Weinberg_purity_residual':mf['relative_Weinberg_purity_residual'] if mf_norm>=1e-100 else None,
            'ward_SU2xU1':mf['ward_SU2xU1'] if mf_norm>=1e-100 else None,
            'interpretation':'Order-zero boundary witness only; cannot certify the one-loop ADM'}
    wards = result['ward_SU2xU1']
    maximum_ward = max(x['relative_to_LLHH'] for x in wards.values())
    best = result['alternative_coordinate_hypotheses_DO_NOT_APPLY'][0]
    mf = result['source_MF_orderzero_witness']
    if result['relative_Weinberg_purity_residual'] < 1e-8 and maximum_ward<1e-8:
        hint = 'pure_Weinberg_in_declared_basis'
    elif best['relative_residual']<1e-7 and best['hypothesis_only']!='declared_C_Weinberg_template':
        hint = 'possible_component_frame_discrepancy_verify_actual_source_map'
    elif maximum_ward>1e-7:
        hint = 'source_LLHH_breaks_gauge_Ward_in_declared_basis_or_basis_map_wrong'
    else:
        hint = 'invariant_but_not_projected_check_operator_conventions'
    result['diagnostic_summary'] = {
        'weinberg_purity_relative':result['relative_Weinberg_purity_residual'],
        'SU2_T1_relative':wards['SU2_T1']['relative_to_LLHH'],
        'SU2_T2_relative':wards['SU2_T2']['relative_to_LLHH'],
        'SU2_T3_relative':wards['SU2_T3']['relative_to_LLHH'],
        'U1_Y_relative':wards['U1_Y']['relative_to_LLHH'],
        'max_gauge_ward_relative':maximum_ward,
        'best_diagnostic_template':best['hypothesis_only'],
        'best_diagnostic_template_residual':best['relative_residual'],
        'source_MF_witness_status':mf['availability'],
        'source_MF_relative_residual':mf.get('relative_Weinberg_purity_residual'),
        'first_investigation_hint_NOT_A_REPAIR':hint}
    result.update(model=model, source_manifest=str(wilson_manifest),
                  source_C5_status='evolved_EFT1_at_MS_not_recomputed',
                  scalar_threshold_matching_performed=False)
    path=Path(output)
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    return result
