"""v20: immutable 178-slot real+complex source-conditioned scalar inventory.

Only the registration census is complete. No unique RGBeta S2 numerical axes,
physical Matchete/RGBeta coupling dictionary, or mixed-mass basis is claimed.
"""
from __future__ import annotations
from pathlib import Path
from .complex_scalar_bridge_v19 import (
    FORMAT as V19_FORMAT, STATUS as V19_STATUS, audit as audit_v19,
    load, sha, require,
)
from .full_scalar_census_v17 import FORMAT as V17_FORMAT, STATUS as V17_STATUS
from .native_context_quartic_v14 import audit_witness, read as read_v14

FORMAT='t3-parallel-a-178-real-complex-scalar-registry-census-v20'
STATUS='Real150Complex28NativeRegistrationAndConditionalSourceActionCensusComplete_PhysicsBlocked'

def audit(root:Path, v17:Path, v19:Path, native14:Path,
          native18_topology:Path, native18_bc:Path, v15:Path,
          mat50:Path,rgb50:Path,v3:Path)->dict:
    root=root.resolve()
    new=audit_v19(root,native18_topology,native18_bc,native14,v15,mat50,rgb50,v3,samples=3)
    saved=load(v19);real=load(v17)
    require(saved.get('format')==new['format']==V19_FORMAT and
            saved.get('status')==new['status']==V19_STATUS and
            saved.get('input_sha256')==new['input_sha256'] and
            saved.get('certification')==new['certification'] and
            saved.get('native_Matchete_topology_CGs')==16 and
            saved.get('native_Matchete_BC_extra_CGs_rechecked')==14 and
            saved.get('baseline_EFT1_complex_slots')==28 and
            saved.get('max_native_topology_6_over_24_real_action_relative_residual',1.)<1e-9,
            'saved v19 missing, stale or inconsistently certified')
    for tag in new['models']:
        require(saved['models'][tag]['baseline_complex_quartic_slots']==
                new['models'][tag]['baseline_complex_quartic_slots'],
                tag+': saved v19 complex census tampered')
    # v17 was produced by a source-guarded read-only audit. Recheck every
    # exact witness-file hash rather than relying on its status string alone.
    v17_inputs={
        'Matchete_nonportal_v16':'matchete_nonportal_native_v16.json',
        'Matchete_portals_v5':'matchete_portals_native_v5.json',
        'native_RGBeta_v14':'rgbeta_quartic_projectors_native_v14.json',
        'nonportal_dictionary_v16':'nonportal_dictionary_conditional_v16.json',
        'portal_real_action_v9':'portal_real_action_conditional_v9.json',
        'native_quartic_v14_summary':'native_quartic_projectors_conditional_v14.json',
    }
    recorded=real.get('input_sha256',{})
    require(isinstance(recorded,dict) and set(recorded)==set(v17_inputs),
            'v17 recorded input-hash inventory incomplete')
    for label,name in v17_inputs.items():
        witness=root/'output/parallel_t3/wilson'/name
        require(witness.is_file() and recorded[label]==sha(witness),
                'v17 source-input witness changed: '+name)
    require(real.get('format')==V17_FORMAT and real.get('status')==V17_STATUS and
            real.get('native_RGBeta_real_quartic_slots')==150 and
            real.get('native_RGBeta_real_projector_unit_diagonals')==150 and
            real.get('native_RGBeta_portal_response_matrix_entries')==110 and
            real.get('Matchete_nonportal_scalar_slots')==76 and
            real.get('Matchete_Higgs_portal_slots')==58 and
            real.get('Matchete_Higgs_quartic_slots')==16 and
            real.get('certification',{}).get('physical_Wilson_evolution_authorized') is False,
            'saved v17 real-quartic census invalid')
    native=read_v14(native14)
    check=audit_witness(native,root)
    require(check['source_census'].get('quartic_registered')==178 and
            check['source_census'].get('real_registered')==150 and
            check['source_census'].get('complex_registered')==28 and
            check['source_census'].get('real_diagonal_unit')==150 and
            check['complex_probe_census'].get('complex_direct_unit_conj_zero')==28,
            'actual installed RGBeta native quartic census changed')
    models={}
    for tag in sorted(new['models']):
        rr=real['models'][tag]['total_real_quartic_slots']
        cc=new['models'][tag]['baseline_complex_quartic_slots']
        models[tag]={'native_RGBeta_real_slots':rr,'native_RGBeta_complex_slots':cc,
                     'quartic_slots_total':rr+cc,
                     'physical_Matchete_RGBeta_dictionary_certified':False}
    require(sum(m['quartic_slots_total'] for m in models.values())==178 and
            sum(m['native_RGBeta_real_slots'] for m in models.values())==150 and
            sum(m['native_RGBeta_complex_slots'] for m in models.values())==28,
            'full scalar registry tally is not 178=150+28')
    return {'format':FORMAT,'status':STATUS,'models_checked':16,
        'total_EFT1_scalar_quartic_registration_slots':178,
        'native_RGBeta_real_quartics':150,'native_RGBeta_complex_quartics':28,
        'Matchete_native_Higgs_portal_CG_slots':58,
        'Matchete_native_nonportal_real_CG_slots':76,
        'Matchete_native_topology_complex_CG_slots':16,
        'Matchete_BC_current_extra_complex_CG_slots':14,
        'real_native_RGBeta_projector_diagonal_checks':150,
        'real_native_RGBeta_portal_response_checks':110,
        'complex_native_RGBeta_direct_conjugate_checks':28,
        'models':models,
        'input_sha256':{x:sha(p) for x,p in (
            ('real_census_v17',v17),('complex_bridge_v19',v19),
            ('native_rgbeta_v14',native14))},
        'certification':{
            'native_EFT1_quartic_registration_inventory_complete':True,
            'real_scalar_source_conditional_dictionary_checked':True,
            'complex_quartic_baseline_subspace_source_conditioned_checked':True,
            'native_RGBeta_triplet_numerical_frame_certified':False,
            'complete_physical_scalar_coupling_dictionary_certified':False,
            'BC_minus1_mixed_scalar_mass_basis_certified':False,
            'signed_GammaCCWeyl_certified':False,
            'one_loop_Wilson_ADM_physical_certified':False,
            'physical_Wilson_evolution_authorized':False,
            'EFT_order':'UV -> integrate out ALL fermions -> integrate out ALL scalars -> SMEFT'}}
