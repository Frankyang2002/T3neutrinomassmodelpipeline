"""v27: explicit B/C all-fermion-first, all-scalar-second threshold itinerary.

No RG evolution, no scalar threshold coefficient, and no SMEFT Wilson C5.
Existing 1 TeV EFT1 endpoints lie below the two scalar eigenmasses; this is
recorded as a source-scale incompatibility for direct matching at either mass.
"""
from __future__ import annotations
from pathlib import Path
import math
from .bc_native_mass_eigen_reconcile_v25 import BC,GateError,load,digest,require
from .bc_one_scale_probe_audit_v26 import FORMAT as F26,STATUS as S26
from .bc_native_mass_eigen_reconcile_v25 import FORMAT as F25,STATUS as S25

FORMAT='t3-parallel-a-bc-fermion-first-scalar-threshold-itinerary-v27'
STATUS='BCMatcheteEigenbasisLinkedToSourceEFT1_ThresholdReevaluationRequired_PhysicsBlocked'
ORDER='UV -> integrate out ALL fermions -> integrate out ALL scalars -> SMEFT'

def audit(root:Path,v22_path:Path,v25_path:Path,v26_path:Path)->dict:
    root=Path(root).resolve()
    v22=load(v22_path);v25=load(v25_path);v26=load(v26_path)
    require(v25.get('format')==F25 and v25.get('status')==S25 and
            v25.get('models_checked')==2,'invalid v25 evidence')
    require(v26.get('format')==F26 and v26.get('status')==S26 and
            v26.get('source_sha256',{}).get('v25')==digest(v25_path),
            'v26 one-step evidence stale')
    require(v25.get('source_sha256',{}).get('v22')==digest(v22_path),
            'v25 source v22 hash stale')
    from .bc_mass_bridge_gate_v22 import FORMAT as F22,STATUS as S22
    require(v22.get('format')==F22 and v22.get('status')==S22 and
            v22.get('full_scalar_quartic_slots')==178 and
            v22.get('certification',{}).get('physical_Wilson_evolution_authorized') is False,
            'v22 scalar census invalid')
    out={}
    for model in BC:
        m=v25['models'][model];loop=v26['models'][model]
        require(loop.get('one_step_all_heavy_not_fermion_first') is True and
                loop.get('actual_EFT1_scalar_threshold_matching_performed') is False,
                model+': one-step probe mislabeled as EFT1')
        mx,my=m['mass_gev']; endpoint=m['source_eft1_endpoint_gev']
        require(all(math.isfinite(v) and v>0 for v in (mx,my,endpoint)) and
                endpoint<mx<my, model+': physical threshold ordering changed')
        require(m.get('two_scalar_thresholds_above_saved_eft1_endpoint') is True,
                model+': endpoint classification wrong')
        loggap=math.log(my/mx)
        log_from_endpoint=[math.log(mx/endpoint),math.log(my/endpoint)]
        require(loggap>0 and all(x>0 for x in log_from_endpoint),model+': threshold logs nonpositive')
        out[model]={
            'model':model,'scalar_mass_eigenstates':{'X_gev':mx,'Y_gev':my},
            'saved_candidate_EFT1_scale_gev':endpoint,
            'source_candidate_EFT1_state_is_at_scalar_matching_scale':False,
            'log_my_over_mx':loggap,
            'log_mx_over_saved_eft1_scale':log_from_endpoint[0],
            'log_my_over_saved_eft1_scale':log_from_endpoint[1],
            'source_UV_fermion_first_matching_verified':False,
            'joint_scalar_matching_completed':False,
            'sequential_eigenstate_scalar_matching_completed':False,
            'direct_one_scale_matchete_probe_not_valid_EFT1_scalar_matching':True,
            'required_next_numerical_inputs':[
                'all-heavy-fermion EFT1 Wilson boundary at fermion threshold',
                'source-conditioned EFT1 running couplings at chosen scalar matching scale',
                'signed LLSS-to-C5 scalar loop matching in the X/Y eigenfield basis',
                'physical one-loop EFT1 Wilson anomalous dimension'],
        }
    return {'format':FORMAT,'status':STATUS,'models_checked':2,'models':out,
        'full_scalar_quartic_registration_slots':178,
        'evidence':{'native_BC_mass_eigenfield_Matchete_checks':2,
                    'native_BC_complex_CG_mass_rotations':16,
                    'existing_direct_one_scale_exploratory_Weinberg_terms':472},
        'source_sha256':{'v22':digest(v22_path),'v25':digest(v25_path),'v26':digest(v26_path)},
        'certification':{'EFT_order':ORDER,
            'physical_fermion_hard_matching_boundary':False,
            'physical_Matchete_RGBeta_coupling_map':False,
            'actual_scalar_matching_at_threshold':False,
            'one_loop_dim5_ADM':False,
            'signed_GammaCC_to_Weyl':False,
            'SMEFT_running_completed':False,
            'physical_neutrino_mass_prediction':False,
            'physical_Wilson_evolution_authorized':False}}
