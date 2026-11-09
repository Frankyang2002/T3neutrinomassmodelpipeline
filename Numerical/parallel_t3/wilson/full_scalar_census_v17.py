"""v17 real scalar 150-slot census (read-only, conditional source frame).

A complete quartic *census* is not a physical coupling dictionary. Cross-check
v9 portal 58 + v16 nonportal 76 + 16 SM Higgs norms against the independently
executed RGBeta v14 150 real projector diagonals. Complex slots stay separate.
"""
from __future__ import annotations
from hashlib import sha256
from pathlib import Path
import numpy as np
from .nonportal_native_dictionary_v16 import (
    read, sha, audit_from_paths as audit_v16, require,
    EXPECTED_MODELS, DIMENSIONS, FORMAT as V16_FORMAT, STATUS as V16_STATUS,
)
from .portal_real_action_v9 import STATUS as V9_STATUS, FORMAT as V9_FORMAT
from .native_context_quartic_v14 import (
    audit_from_paths as audit_v14,
    STATUS as V14_STATUS,
)

FORMAT='t3-parallel-a-full-real-scalar-census-conditional-v17'
STATUS='RealQuartic150of150SourceFrameCensus_NativeProjectorsVerified_PhysicalDictionaryConditional'

def audit_from_paths(root:Path,v16_native:Path,v5_portal:Path,native14:Path,
                     v16_saved:Path,v9_saved:Path,v14_saved:Path)->dict:
    root=Path(root).resolve()
    v16_new=audit_v16(v16_native,v5_portal,native14,root,samples=5)
    v16_old=read(v16_saved)
    require(v16_old.get('format')==v16_new['format']==V16_FORMAT and
            v16_old.get('status')==v16_new['status']==V16_STATUS and
            v16_old.get('input_sha256')==v16_new.get('input_sha256'),
            'v16 saved report missing or changed source provenance')
    require(v16_old.get('nonportal_coupling_slots')==76 and
            v16_old.get('nonportal_sectors')==48 and
            v16_old.get('certification')==v16_new.get('certification') and
            float(v16_old.get('max_relative_residual',1.))<1e-10 and
            set(v16_old.get('models',{}))==EXPECTED_MODELS,
            'v16 report does not cover 76/48 sectors')
    for tag in EXPECTED_MODELS:
        for fam in ('self1','self2','cross'):
            left=v16_old['models'][tag][fam]
            right=v16_new['models'][tag][fam]
            require(left.get('matchete_couplings')==right['matchete_couplings'] and
                    left.get('rgbeta_couplings')==right['rgbeta_couplings'] and
                    np.allclose(left.get('conditional_rgbeta_from_matchete_real_matrix',[]),
                                right['conditional_rgbeta_from_matchete_real_matrix'],atol=1e-10),
                    f'{tag}/{fam}: saved v16 map altered')
    v14_current=audit_v14(root,native14,v9_saved)
    v14_old=read(v14_saved)
    require(v14_old.get('status')==v14_current.get('status')==V14_STATUS and
            v14_old.get('input_sha256')==v14_current.get('input_sha256') and
            v14_old.get('source_census')==v14_current.get('source_census') and
            v14_old.get('real_projector_conditionally_verified') is True and
            v14_old.get('physical_Wilson_evolution_authorized') is False and
            v14_old.get('full_scalar_potential_dictionary_certified') is False and
            v14_old.get('native_numeric_S2_frame_verified') is False,
            'saved native RGBeta v14 summary is stale, incomplete or blocked')
    census=v14_current['source_census']
    require(census.get('real_registered')==150 and
            census.get('real_diagonal_unit')==150 and
            census.get('portal_response_numeric')==110,
            'native RGBeta full real projectors not 150/110')
    v9=read(v9_saved)
    require(v9.get('format')==V9_FORMAT and v9.get('status')==V9_STATUS and
            v9.get('native_Matchete_portal_coupling_slots')==58 and
            v9.get('models_checked')==16 and
            set(v9.get('models',{}))==EXPECTED_MODELS and
            v9.get('certification',{}).get('physical_Wilson_evolution_authorized') is False,
            'saved v9 portal action 58 evidence incomplete')
    models={};total=0
    for tag in sorted(EXPECTED_MODELS):
        nr=sum(len(v16_new['models'][tag][fam]['rgbeta_couplings'])
               for fam in ('self1','self2','cross'))
        nh=sum(len(s['rgbeta_registered_couplings']) for s in v9['models'][tag])
        total+=1+nh+nr
        models[tag]={'SM_Higgs_real_couplings':['lambdaH'],
                     'Higgs_portal_real_coupling_slots':nh,
                     'nonportal_scalar_real_coupling_slots':nr,
                     'total_real_quartic_slots':1+nh+nr,
                     'native_RGBeta_real_projector_diagonals_verified':True,
                     'physical_coupling_dictionary_certified':False}
    require(total==150,'real quartic 16+58+76 != 150')
    return {'format':FORMAT,'status':STATUS,
            'native_RGBeta_real_quartic_slots':150,
            'native_RGBeta_real_projector_unit_diagonals':150,
            'native_RGBeta_portal_response_matrix_entries':110,
            'Matchete_Higgs_quartic_slots':16,
            'Matchete_Higgs_portal_slots':58,
            'Matchete_nonportal_scalar_slots':76,
            'Matchete_nonportal_sectors':48,
            'native_RGBeta_complex_quartics_exploratory_only':28,
            'models':models,
            'input_sha256':{name:sha(p) for name,p in (
                ('Matchete_nonportal_v16',v16_native),('Matchete_portals_v5',v5_portal),
                ('native_RGBeta_v14',native14),('nonportal_dictionary_v16',v16_saved),
                ('portal_real_action_v9',v9_saved),('native_quartic_v14_summary',v14_saved))},
            'certification':{
                'real_quartic_registration_census_complete':True,
                'native_RGBeta_real_projector_unit_responses_verified':True,
                'Matchete_nonportal_real_CG_action_in_source_aligned_SU2_span':True,
                'independent_Matchete_RGBeta_scalar_component_frame_verified':False,
                'all_complex_matchete_rgbeta_actions_certified':False,
                'BC_minus1_mixed_scalar_mass_diagonalization_certified':False,
                'signed_GammaCCWeyl_certified':False,
                'physical_Wilson_evolution_authorized':False,
                'EFT_order':'UV -> integrate out ALL fermions -> integrate out ALL scalars -> SMEFT'
            }}
