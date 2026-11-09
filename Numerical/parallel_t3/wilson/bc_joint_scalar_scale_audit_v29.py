"""v29: independently recheck saved common scalar threshold state and provenance.

No mixing of v25 and new v28 numerical benchmark inputs. Only both scalar
states at ONE joint muS are accepted. One-step all-heavy C5 is never input.
"""
from __future__ import annotations

from pathlib import Path
import math
import re
import numpy as np
from .bc_joint_scalar_scale_v28 import (
    FORMAT as F28, STATUS as S28, ORDER, MODELS, need, sha, read, checked_complex,
)

FORMAT='t3-parallel-a-bc-joint-threshold-state-readiness-v29'
STATUS='TwoScalarCommonEFT1ScaleStatesReauthenticated_C5AndScalarMatchingBlocked'


def source_paths(root: Path, report: dict)->dict[str,Path]:
    root=Path(root).resolve()
    paths=report.get('source_files');hashes=report.get('source_sha256')
    need(isinstance(paths,dict) and set(paths)==set(hashes or {}) and len(paths)>=3,
         'missing source path/digest inventory')
    result={}
    for key,value in paths.items():
        need(isinstance(key,str) and isinstance(value,str)
             and len(value)<350 and '\\' not in value and ':' not in value,
             'invalid provenance source path')
        parts=value.split('/')
        need(all(part not in ('','.','..') for part in parts) and not value.startswith('/'),
             'source escapes repository')
        path=(root / Path(*parts)).resolve()
        need(path.is_relative_to(root), 'source path resolves outside repository')
        recorded=hashes.get(key)
        need(isinstance(recorded,str) and bool(re.fullmatch('[0-9a-f]{64}',recorded)) and
             sha(path)==recorded,'stale source hash: '+value)
        result[key]=path
    return result


def audit_model(root:Path,report:dict)->dict:
    model=report.get('model')
    need(report.get('format')==F28 and report.get('status')==S28 and
         model in MODELS and report.get('EFT_order')==ORDER,
         'unknown or incorrect EFT1 source report')
    certification=report.get('certification')
    for field in ('physical_Matchete_RGBeta_coupling_dictionary',
                  'fermion_threshold_Wilson_boundary_computed',
                  'signed_GammaCC_to_Weyl','full_one_loop_dim5_ADM',
                  'joint_two_scalar_threshold_Wilson_matching_performed',
                  'physical_SMEFT_C5_generated','physical_Wilson_evolution_authorized'):
        need(isinstance(certification,dict) and certification.get(field) is False,
             'unsupported physics authorization: '+field)
    need(certification.get('all_two_scalar_eigenstates_retained_to_common_matching_scale') is True and
         report.get('native_scalar_component_tensor_transport_repeated_here') is False and
         report.get('independent_real_hessian_spectrum_verified') is True and
         report.get('scalar_eigenprojectors_complete') is True,
         'missing full-scalar source-frame checks')
    need(report.get('backend') in ('five_mass','quadratic_only') and
         report.get('backend_full_higgs_mass_beta') is (report['backend']=='five_mass'),
         'incorrect beta backend class')
    if report['backend']=='quadratic_only':
        need(report.get('state_at_common_scale',{}).get('mHSq_gev2') is None,
             'quadratic-only audit pretends to evolve Higgs mass')
        need(report.get('saved_quadratic_endpoint_relative_residual') is not None and
             0<=report['saved_quadratic_endpoint_relative_residual']<2e-6,
             'missing saved quadratic endpoint replay check')
    else:
        need(math.isfinite(report.get('state_at_common_scale',{}).get('mHSq_gev2',float('nan'))),
             'full Higgs mass beta missing from five-mass state')
    muF=float(report['mu_fermion_first_gev']);mu0=float(report['mu_saved_candidate_endpoint_gev'])
    need(muF>mu0>0 and muF==1e5 and mu0==1e3,'candidate EFT1 scale hierarchy changed')
    data=report.get('common_scalar_matching_scale');other=report.get('other_joint_scale_scheme')
    need(isinstance(data,dict) and isinstance(other,dict) and
         {data.get('scheme'),other.get('scheme')}=={'geometric','upper'},
         'both common-scale prescriptions must be sampled')
    for item in (data,other):
        mu=float(item.get('matching_scale_gev',float('nan')))
        mx=float(item.get('mX_gev',float('nan')))
        my=float(item.get('mY_gev',float('nan')))
        need(all(math.isfinite(x) and x>0 for x in (mu,mx,my)) and
             mu0<mu<muF and mx<=my,'nonpositive/unordered scalar spectrum')
        need(abs(float(item['self_consistent_root_log_residual']))<1e-9 and
             abs(math.log(mu/mx)-item['log_mu_over_mX'])<1e-10 and
             abs(math.log(mu/my)-item['log_mu_over_mY'])<1e-10,
             'joint crossing logs inconsistent')
        if item['scheme']=='geometric':
            need(abs(math.log(mu)-.5*(math.log(mx)+math.log(my)))<1e-9,
                 'geometric common scale fails self-consistency')
        else:
            need(abs(math.log(mu/my))<1e-9,'upper common scale fails self-consistency')
        need(item.get('two_scalar_eigenstates_both_retained_in_EFT1_until_matching') is True and
             item.get('intermediate_single_scalar_EFT2_used') is False,
             'invalid EFT scalar decoupling order')
    actual=report.get('state_at_common_scale',{})
    need(abs(actual.get('mu_gev',-1)-data['matching_scale_gev'])<2e-8*data['matching_scale_gev'],
         'state at wrong RGE scale')
    basis=actual.get('mass_basis',{})
    eig=np.asarray(basis.get('masses_squared_gev2',[]),float)
    u=basis.get('unitary_flavour_rotation_U',{})
    need(eig.shape==(2,) and np.isfinite(eig).all() and np.all(eig>0),
         'mass basis eigenvalues invalid')
    need(np.allclose(eig,data['scalar_mass_squared_eigenvalues_gev2'],rtol=2e-11,atol=3e-7),
         'mass basis eigenvalues differ from self-consistent crossing')
    need(basis.get('matchete_field_redefinition_installed') is False and
         basis.get('threshold_matching_performed') is False,
         'new matching falsely reported')
    need(isinstance(u,dict) and set(u)=={'real','imag'},'invalid eigenbasis U representation')
    U=np.asarray(u['real'],float)+1j*np.asarray(u['imag'],float)
    need(U.shape==(2,2) and np.isfinite(U).all() and
         np.allclose(U.conj().T@U,np.eye(2),atol=3e-12),
         'joint scalar mass basis U not unitary')
    A=basis['mass_squared_flavour_matrix_gev2']
    amat=np.asarray(A['real'],float)+1j*np.asarray(A['imag'],float)
    need(amat.shape==(2,2) and np.isfinite(amat).all() and
         np.allclose(U.conj().T@amat@U,np.diag(eig),atol=5e-5,rtol=2e-11),
         'joint eigenbasis does not diagonalise saved mass matrix')
    need(float(report['epsilon_SU2_covariance_max_absolute_error'])<3e-12,
         'SU2 epsilon covariance failed')
    paths=source_paths(root,report)
    return {'model':model,'backend':report['backend'],
        'common_muS_gev':data['matching_scale_gev'],
        'scalar_masses_gev':[data['mX_gev'],data['mY_gev']],
        'max_log':data['max_absolute_scalar_matching_log'],
        'upper_muS_gev':other['matching_scale_gev'],
        'source_file_count':len(paths),
        'candidate_renormalizable_EFT1_state_revalidated':True,
        'physical_C5_boundary_missing':True}


def audit(root:Path,a:dict,b:dict)->dict:
    need({a.get('model'),b.get('model')}==set(MODELS),'need B_-1 and C_-1 jointly')
    results={d['model']:audit_model(root,d) for d in (a,b)}
    return {'format':FORMAT,'status':STATUS,'models':results,'models_checked':2,
            'single_joint_scalar_threshold_for_each_model':True,
            'intermediate_single_eigenstate_EFT2_not_used':True,
            'certification':{'EFT_order':ORDER,'physical_fermion_Wilson_boundary':False,
                'full_one_loop_Wilson_ADM':False,'signed_Weyl_normalization':False,
                'all_scalar_matching_done':False,'SMEFT_C5_generated':False,
                'physical_Wilson_evolution_authorized':False}}
