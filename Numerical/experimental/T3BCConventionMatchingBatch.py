"""Read-only B/C restricted direct-triangle convention comparison.

No physical matching is certified. Does not change the source evaluator.
"""
from __future__ import annotations
import numpy as np
from Numerical.experimental.T3BCRunningOneScaleTriangle import assemble_ordered_c5, packed, require


def relative(a,b):
    return float(np.linalg.norm(a-b)/max(np.linalg.norm(b),1e-250))


def spectrum(c5,vev):
    return np.sort(np.linalg.svd(-0.5*vev**2*c5*1e9,compute_uv=False))[::-1]


def evaluate(boundary,evaluator,vev,archived):
    """Recompute both rotations with exactly the same coefficients/input state."""
    basis=boundary['scalar_basis']
    mx,my=map(float,np.sqrt(basis.eigenvalues_sq))
    V=boundary['heavy_basis']; y1_old=boundary['y1']; y2=boundary['y2']
    y1_same=boundary['y1_same_V_alternative']
    require(np.linalg.norm(V.conj().T@V-np.eye(3))<1e-10,'nonunitary heavy rotation')
    # Real (or effectively real) Takagi rotations make the two prescriptions
    # indistinguishable. This is a legitimate benchmark outcome, not an error.
    yukawa_rotation_relative=relative(y1_same,y1_old)
    distinguishable=yukawa_rotation_relative>1e-12
    args=(boundary['mf'], None,y2,mx,my,basis.U,boundary['quartics'],evaluator)
    old,old_terms=assemble_ordered_c5(args[0],y1_old,*args[2:])
    same,same_terms=assemble_ordered_c5(args[0],y1_same,*args[2:])
    stored=np.asarray(archived['candidate_ordered_source_C5_3x3_gev_inverse_NOT_PHYSICAL']['real'],float)+1j*np.asarray(archived['candidate_ordered_source_C5_3x3_gev_inverse_NOT_PHYSICAL']['imag'],float)
    parity=relative(old,stored)
    require(parity<2e-8,'original production convention does not reproduce archived result')
    require(relative(old,old.T)<1e-10 and relative(same,same.T)<1e-10,'non-symmetric C5')
    old_m=spectrum(old,vev); same_m=spectrum(same,vev)
    return {
       'status':'TwoConditionalYukawaRotationsRecomputedAgainstSameNativeRestrictedSource',
       'model':boundary['model'],'matching_scale_gev':boundary['mu_gev'],
       'mf_gev':list(map(float,boundary['mf'])),'mx_gev':mx,'my_gev':my,
       'vev_reference_gev':float(vev),
       'production_baseline':{'rotation':'y1@conj(V), y2@V',
         'C5_GeV_inverse_NOT_PHYSICAL':packed(old),
         'neutrino_masses_descending_ev_NOT_PHYSICAL':old_m.tolist(),
         'per_heavy_terms':old_terms},
       'same_U_weyl_candidate':{'rotation':'y1@V, y2@V',
         'C5_GeV_inverse_NOT_PHYSICAL':packed(same),
         'neutrino_masses_descending_ev_NOT_PHYSICAL':same_m.tolist(),
         'per_heavy_terms':same_terms},
       'heavy_yukawa_rotation_relative_difference':yukawa_rotation_relative,
       'rotation_prescriptions_numerically_distinguishable':bool(distinguishable),
       'rotation_comparison_interpretation':(
           'Distinct rotations numerically compared in this benchmark'
           if distinguishable else
           'Current running benchmark does not distinguish V from conj(V) in y1; '
           'zero C5 sensitivity here does not certify the cross-package convention'),
       'relative_C5_same_U_minus_production':relative(same,old),
       'relative_neutrino_spectrum_same_U_minus_production':relative(same_m,old_m),
       'recomputed_production_vs_archived_C5_relative':parity,
       'both_cases_use_same_running_endpoint_and_same_coefficient_function':True,
       'native_matchete_rgbeta_coupling_dictionary_certified':False,
       'physical_C5_certified':False,'physical_neutrino_masses_certified':False,
       'complete_sequential_EFT_threshold_matching':False,
    }


def check_gate(source,archived,native,action):
    require(archived.get('physical_C5_certified') is False,'archived result wrongly physical')
    require(archived.get('full_fermion_X_Y_sequential_RG_transport') is False,'source scope changed')
    require(source.get('physical_matchete_rgbeta_mapping_certified') is False,'source scope changed')
    require(native.get('status')=='NativeAbstractChiralChecksPassedPhysicalDictionaryPending',
            'native Matchete result missing or changed')
    for model in ('B_m1','C_m1'):
        row=native['models'][model]
        require(row.get('passed') is True and all(row.get('checks',{}).values()),
                'native abstract identities failed '+model)
    require(native.get('rgbeta_cross_package_dictionary_certified') is False,
            'native result scope changed')
    require(action.get('physical_C5_certified') is False,
            'action report is wrongly labelled physical')
    require(action.get('status')=='MajoranaWeylActionCovariancePassedNativeDictionaryPending',
            'action result missing or changed')
    return {'native_abstract_matchete_check':'PassedBothModels',
            'weyl_action_check':'Passed',
            'cross_package_absolute_dictionary':'Uncertified'}
