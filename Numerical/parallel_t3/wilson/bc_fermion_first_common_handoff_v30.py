"""v30: guard v27->v29 EFT-order handoff. No Wilson matching or SMEFT output."""
from __future__ import annotations
from pathlib import Path
from .bc_joint_scalar_scale_v28 import ORDER, sha, need, read
from .bc_joint_scalar_scale_audit_v29 import FORMAT as F29, STATUS as S29

FORMAT='t3-parallel-a-bc-all-fermion-first-all-scalar-joint-scale-handoff-v30'
STATUS='BCJointScalarRenormalizableScaleInputsReady_WilsonBoundaryMatchingAndADMBlocked'


def audit(v27_path:Path,v29_path:Path, *,root:Path)->dict:
    from .bc_threshold_order_ledger_v27 import FORMAT as F27, STATUS as S27
    v27=read(v27_path);v29=read(v29_path)
    need(v27.get('format')==F27 and v27.get('status')==S27 and
         v27.get('models_checked')==2 and
         v27.get('certification',{}).get('EFT_order')==ORDER,
         'invalid prior fermion-first scalar ledger')
    need(v29.get('format')==F29 and v29.get('status')==S29 and
         v29.get('models_checked')==2 and
         v29.get('certification',{}).get('EFT_order')==ORDER and
         v29.get('certification',{}).get('physical_Wilson_evolution_authorized') is False,
         'invalid new renormalisable common-scale audit')
    # The v27 ledger is not a new physics certificate. Revalidate its prior
    # SHA256 chain against the three immutable v22/v25/v26 source witnesses.
    upstream={'v22':'bc_quartic_mass_readiness_conditional_v22.json',
              'v25':'bc_native_matchete_mass_eigen_concordance_v25.json',
              'v26':'bc_existing_one_scale_matchete_probe_v26.json'}
    recorded=v27.get('source_sha256')
    need(isinstance(recorded,dict) and set(recorded)==set(upstream),
         'v27 source hash inventory incomplete')
    for key,name in upstream.items():
        need(recorded[key]==sha(v27_path.parent/name),
             'v27 upstream file changed: '+name)
    for model in ('B_-1','C_-1'):
        prior=v27['models'][model];fresh=v29['models'][model]
        need(prior.get('source_UV_fermion_first_matching_verified') is False and
             prior.get('joint_scalar_matching_completed') is False and
             prior.get('direct_one_scale_matchete_probe_not_valid_EFT1_scalar_matching') is True,
             model+': old matching claim changed')
        need(fresh.get('candidate_renormalizable_EFT1_state_revalidated') is True and
             fresh.get('physical_C5_boundary_missing') is True,
             model+': claimed physical C5 too early')
    return {'format':FORMAT,'status':STATUS,'models_checked':2,
      'source_sha256':{'v27':sha(v27_path),'v29':sha(v29_path)},
      'models':{name:{'candidate_common_muS_gev':v29['models'][name]['common_muS_gev'],
                      'beta_backend':v29['models'][name]['backend'],
                      'source_1_TeV_not_used_as_matching_boundary':True,
                      'all_scalars_integrated_out_at_one_common_scale':False,
                      'remaining_tasks':['signed all-heavy-fermion LLSS Wilson boundary',
                         'verified full one-loop EFT1 Wilson anomalous dimension',
                         'joint X/Y scalar matching into SMEFT at common muS']}
                for name in ('B_-1','C_-1')},
      'certification':{'EFT_order':ORDER,
           'renormalizable_candidate_scale_prepared':True,
           'physical_Wilson_boundary':False,'physical_one_loop_ADM':False,
           'signed_GammaCC_to_Weyl':False,
           'all_heavy_scalar_threshold_matching':False,
           'SMEFT_matching':False,
           'physical_neutrino_mass_prediction':False,
           'physical_Wilson_evolution_authorized':False}}
