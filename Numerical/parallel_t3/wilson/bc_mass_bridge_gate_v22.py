"""A v22: immutable B/C native complex CG mass-transport + 178-slot readiness gate.

The report is a CHECKPOINT: it explicitly does not authorize a physical EFT run.
"""
from __future__ import annotations
from pathlib import Path
from .bc_mass_eigen_native_bridge_v21 import (
    FORMAT as V21_FORMAT, STATUS as V21_STATUS,
    BC, BCBridgeBlocked, read_json, digest, require,
)
from .full_scalar_census_v20 import (
    FORMAT as V20_FORMAT, STATUS as V20_STATUS,
)

FORMAT='t3-parallel-a-bc-source-conditional-178-slot-readiness-gate-v22'
STATUS='Real150Complex28AndBCMassEigenNativeCGTransportChecked_PhysicalMatchingAndADMBlocked'

BLOCKERS={
    'native_RGBeta_S2_component_intertwiner':'not certified: S2 generators remain an abstract representation',
    'complete_physical_Matchete_RGBeta_quartic_dictionary':'not certified: conditional source component frames',
    'matchete_mass_eigenstate_redefined_Lagrangian':'not certified: existing scalar rotation is experimental',
    'B_C_mixed_scalar_HardMatching':'not certified: no native matched heavy-scalar eigenfield action',
    'signed_GammaCC_to_Weyl':'not certified by v20/v21 scalar witnesses',
    'physical_one_loop_dim5_ADM':'unresolved loop counterterm/normalization',
    'physical_fermion_hard_boundary':'not provided by current scalar audit',
    'scalar_matching_and_SMEFT_running':'not performed',
}

def audit(root:Path,v20:Path,v21:Path)->dict:
    r20=read_json(v20);r21=read_json(v21)
    require(r20.get('format')==V20_FORMAT and r20.get('status')==V20_STATUS and
            r20.get('total_EFT1_scalar_quartic_registration_slots')==178 and
            r20.get('native_RGBeta_real_quartics')==150 and
            r20.get('native_RGBeta_complex_quartics')==28 and
            r20.get('certification',{}).get('physical_Wilson_evolution_authorized') is False,
            'invalid v20 178-slot census')
    require(r21.get('format')==V21_FORMAT and r21.get('status')==V21_STATUS and
            r21.get('models_checked')==2 and
            r21.get('native_Matchete_complex_CG_sources_rotated')==16 and
            r21.get('native_topology_CGs_rotated')==2 and
            r21.get('native_BC_optional_extra_CGs_rotated')==14 and
            r21.get('baseline_RGBeta_special_complex_quartic_slots_conditionally_mapped')==12 and
            r21.get('max_complex_quartic_action_transport_relative_residual',1.)<2e-11,
            'invalid v21 native Matchete mass-rotation audit')
    flags=r21.get('certification',{})
    for k in ('matchete_redefined_mass_eigenstate_lagrangian_certified',
              'complete_physical_Matchete_RGBeta_coupling_map_certified',
              'signed_GammaCCWeyl_certified',
              'one_loop_Wilson_ADM_physical_certified',
              'physical_Wilson_evolution_authorized'):
        require(flags.get(k) is False,'v21 unsafe certification flag '+k)
    model=r21.get('models',{})
    require(isinstance(model,dict) and set(model)==set(BC),'v21 B/C model set wrong')
    for tag in BC:
        row=model[tag]
        require(row.get('source_native_complex_CG_count')==8 and
                row.get('baseline_RGBeta_special_complex_slots')==6 and
                row.get('baseline_six_from_seven_rank')==6 and
                row.get('optional_seventh_source_direction_action_nonzero') is True and
                row.get('optional_seventh_source_direction_worst_holdout_residual',1.)<2e-11 and
                row.get('max_source_native_action_relative_residual',1.)<2e-11 and
                row.get('physical_component_coupling_map_certified') is False and
                row.get('physical_scalar_matching_performed') is False,
                tag+': B/C source-conditional rank/action guard not met')
    recorded=r21.get('source_sha256',{})
    source_names={
        'v18_topology':'matchete_topology_native_v18.json',
        'v18_bc_extra':'matchete_bc_extra_native_v18.json',
        'v19_complex_bridge':'complex_scalar_bridge_conditional_v19.json',
        'v20_quartic_census':'full_scalar_census_conditional_v20.json',
    }
    for k,name in source_names.items():
        witness=root/'output/parallel_t3/wilson'/name
        require(witness.is_file() and recorded.get(k)==digest(witness),
                'v21 source evidence stale '+name)
    mass_folder=root/'output/general_t3_mass_basis'
    for tag in BC:
        name=tag[0]+'_m1_final_mass_basis.json'
        witness=mass_folder/name
        require(witness.is_file() and recorded.get(tag+'_mass_basis')==digest(witness),
                'v21 source mass eigenstate evidence stale '+name)
    require(digest(v20)==recorded['v20_quartic_census'],
            'v21 v20 hash changed')
    return {'format':FORMAT,'status':STATUS,'models_checked':2,
        'native_scalar_real_quartic_census':150,
        'native_scalar_complex_quartic_census':28,
        'full_scalar_quartic_slots':178,
        'native_BC_complex_Matchete_CGs_mass_eigen_transported':16,
        'source_conditioned_BC_six_from_seven_special_quartics':12,
        'source_derived_from_saved_BC_mass_rotations':True,
        'unresolved_physical_blockers':BLOCKERS,
        'source_sha256':{'v20':digest(v20),'v21':digest(v21)},
        'certification':{
            'source_conditioned_scalar_quartic_inventory_complete':True,
            'native_BC_source_CG_mass_eigen_transport_checked':True,
            'Matchete_physical_mass_eigenstate_matching_certified':False,
            'physical_Matchete_RGBeta_coupling_map_certified':False,
            'signed_GammaCCWeyl_certified':False,
            'one_loop_Wilson_ADM_physical_certified':False,
            'physical_Wilson_evolution_authorized':False,
            'EFT_order':'UV -> integrate out ALL fermions -> integrate out ALL scalars -> SMEFT'}}
