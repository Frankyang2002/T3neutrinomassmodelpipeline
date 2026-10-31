"""B/C(-1) restricted-source flavour, running, and SM-only proxy bridge.

This module NEVER identifies a candidate fixed-order Matchete coefficient with
renormalized physical C5.  The explicit SM-only continuation is a proxy: it
incorrectly treats X/Y/F as absent above their physical thresholds, and must
not be combined with the historical sequential Y/X Wilson trajectory.
"""
from __future__ import annotations

from hashlib import sha256
import json
import math
from pathlib import Path
from typing import Any
import numpy as np

MODELS = ('B_-1', 'C_-1')
FORMAT = 't3-bc-mega-running-source-to-sm-proxy-v1'
UV_SOURCE_FORMAT = 'workstream-b-bc-coupled-uv-through-eft1-five-mass-v1'
INPUT_FORMAT = 't3-bc-common-uv-scale-restricted-direct-matchete-triangle-v1'


def require(ok: bool, reason: str) -> None:
    if not ok:
        raise ValueError('T3 B/C mega bridge: '+reason)


def digest(path: Path) -> str:
    require(path.is_file(), 'missing '+str(path))
    return sha256(path.read_bytes()).hexdigest()


def load(path: Path) -> dict:
    require(path.is_file(), 'missing '+str(path))
    result = json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(result, dict), 'JSON root must be an object')
    return result


def unpack_matrix(obj: Any, shape=(3,3)) -> np.ndarray:
    require(isinstance(obj,dict) and set(obj)=={'real','imag'},'complex matrix schema')
    a=np.asarray(obj['real'],float);b=np.asarray(obj['imag'],float)
    require(a.shape==shape and b.shape==shape and np.isfinite(a).all()
            and np.isfinite(b).all(), 'nonfinite/incorrect matrix shape')
    return a+1j*b


def pack_matrix(a: np.ndarray) -> dict:
    v=np.asarray(a,complex)
    require(v.ndim==2 and np.isfinite(v).all(),'invalid result matrix')
    return {'real':v.real.tolist(),'imag':v.imag.tolist()}


def relative(a: np.ndarray,b: np.ndarray) -> float:
    return float(np.linalg.norm(a-b)/max(np.linalg.norm(b),1e-40))


def check_topology_orientation(path:Path)->dict:
    """Check the EXACT source field heads before using conditional phases.

    The identity CConj[V chi] = V* CConj[chi] is algebraic.  It does not
    independently establish 4-component Majorana source normalization or CGs.
    """
    content=path.read_text(encoding='utf-8')
    required=(
      'useCConj = which === 1;',
      'scalarBar = which === 2;',
      'If[useCConj, CConj[NewFermion[r]], NewFermion[r]]',
      'If[useCConj, CConj[NewFermion[j, r]], NewFermion[j, r]]',
      'y[p, r] NCM[Bar[l[i, p]], f] scalar CG[cg, labels]'
    )
    require(all(s in content for s in required),'source Yukawa orientation changed')
    return {'source_sha256':digest(path),
            'y1_source_heavy_spinor':'CConj[F]',
            'y2_source_heavy_spinor':'F',
            'field_redefinition_identity':'F=V chi => CConj[F]=conj(V) CConj[chi]',
            'coupling_transforms':'y1_mass=y1_UV @ conj(V); y2_mass=y2_UV @ V',
            'native_matchete_physically_certified':False}


def complex_spinor_action_test(*,seed:int=218, trials:int=32)->dict:
    """Exact algebraic redefinition of Yukawa interactions, not a Matchete test.

    Contract y1 with conjugated spinor and y2 with un-conjugated spinor in
    arbitrary complex flavour configurations. Include a deliberately WRONG
    y1@V alternative to check that the test actually distinguishes them.
    """
    rng=np.random.default_rng(seed)
    worst=0.; wrong=0.; mass_error=0.
    for _ in range(trials):
        z=lambda *n:rng.normal(size=n)+1j*rng.normal(size=n)
        Q,_=np.linalg.qr(z(3,3))
        y1,y2=z(3,3),z(3,3)
        chi=z(3);chi_c=z(3); ell=z(3)
        # The conjugate and unconjugate spinors are independent algebraic
        # symbols before imposing the charge-conjugation relation.
        old1=ell.T@y1@(Q.conj()@chi_c)
        old2=ell.T@y2@(Q@chi)
        new1=ell.T@(y1@Q.conj())@chi_c
        new2=ell.T@(y2@Q)@chi
        scale=max(1.,abs(old1),abs(old2))
        worst=max(worst,abs(old1-new1)/scale,abs(old2-new2)/scale)
        wrong=max(wrong,abs(old1-ell.T@(y1@Q)@chi_c)/scale)
        # A complex symmetric Majorana matrix in original weak basis.
        masses=np.array([2.,3.,7.]);M=Q.conj()@np.diag(masses)@Q.conj().T
        mass_error=max(mass_error,relative(Q.T@M@Q,np.diag(masses)))
    require(worst<1e-12 and mass_error<1e-12,'complex spinor action covariance failed')
    require(wrong>0.1,'negative control is not distinguishing conjugate rotation')
    return {'status':'ExactFieldRedefinitionAlgebraChecked',
            'trials':trials,'max_correct_action_relative_residual':float(worst),
            'max_majorana_mass_reconstruction_relative_residual':float(mass_error),
            'wrong_same_V_negative_control_max_relative':float(wrong),
            'source_field_redefinition_checked':True,
            'independent_live_matchete_4spinor_test_performed':False,
            'native_quartic_yukawa_CG_dictionary_certified':False}


def check_running_pair(model:str, source:dict, running:dict, *,
                       source_path:Path, running_path:Path)->tuple[dict,dict]:
    require(model in MODELS,'model unsupported')
    require(source.get('format')==UV_SOURCE_FORMAT and
            source.get('status')=='CandidateCoupledUVAndEFT1FiveRealMassesIntegrated' and
            source.get('model')==model and source.get('solver_success') is True,
            'unexpected upstream B/C five-mass integration')
    require(running.get('format')==INPUT_FORMAT and
            running.get('status')=='RunningCommonScaleRestrictedTriangleEvaluatedNotPhysicalC5' and
            running.get('model')==model,'unexpected running one-scale source')
    require(running.get('physical_C5_certified') is False and
            source.get('physical_matchete_rgbeta_mapping_certified') is False and
            running.get('full_fermion_X_Y_sequential_RG_transport') is False,
            'upstream certificate semantics changed; review required')
    require(running.get('running_bc_source_sha256')==digest(source_path),
            'running input provenance differs from saved B/C UV source')
    mu=float(running['matching_input_scale_gev'])
    require(mu==float(source['uv_threshold']['mu_gev'])==1e5,
            'not same running UV endpoint')
    parity=running['running_input_direct_vs_sequential_parity']
    require(parity['status']=='RunningInputRestrictedFullMassTriangleParityPassed' and
            float(parity['max_relative_residual'])<2e-8,
            'existing restricted triangle parity not passed')
    require(float(running['frozen_source_parity']['archived_frozen_C5_parity_relative'])<2e-8,
            'frozen native coefficient reference parity not passed')
    c5=unpack_matrix(running['candidate_ordered_source_C5_3x3_gev_inverse_NOT_PHYSICAL'])
    require(relative(c5,c5.T)<1e-11,'candidate C5 not symmetric')
    masses=np.asarray(running['MF_takagi_masses_gev'],float)
    require(masses.shape==(3,) and np.all(np.diff(masses)>0),
            'nondegenerate ordered Majorana masses required')
    info={'model':model,'same_uv_matching_scale_gev':mu,
          'five_mass_input_sha256':digest(source_path),
          'one_scale_input_sha256':digest(running_path),
          'running_parity_relative':parity['max_relative_residual'],
          'c5_source_symmetric_relative':relative(c5,c5.T),
          'full_sequential_transport':False,
          'physical_C5_certified':False}
    return info, source['uv_threshold']['state']['sm']


def take_sm_initial(sm:dict,K:np.ndarray):
    from Numerical.running.SMWeinbergEvolution import SMWeinbergInitialConditions
    return SMWeinbergInitialConditions(
        gY=float(sm['gY']),g2=float(sm['g2']),g3=float(sm['g3']),
        lambdaH=float(sm['lambdaH']),
        ye=unpack_matrix(sm['ye']),yu=unpack_matrix(sm['yu']),
        yd=unpack_matrix(sm['yd']),K=K).validated()


def mass_summary(K:np.ndarray,vev:float)->dict:
    M=-0.5*float(vev)**2*np.asarray(K,complex)*1e9
    require(np.isfinite(M).all() and relative(M,M.T)<1e-11,
            'invalid symmetric reference neutrino matrix')
    masses=np.sort(np.linalg.svd(M,compute_uv=False))
    return {'mnu_reference_matrix_ev':pack_matrix(M),
            'masses_ascending_ev':masses.tolist(),
            'masses_descending_ev':masses[::-1].tolist(),
            'delta_m21_sq_ev2':float(masses[1]**2-masses[0]**2),
            'delta_m31_sq_ev2':float(masses[2]**2-masses[0]**2)}


def run_sm_only_proxy(model:str,uv_source:dict,one_scale:dict,*,
                      low_scale:float=100.,rtol:float=1e-8,atol:float=1e-21)->dict:
    """Pure SM(+C5) numerical proxy evolved from 1e5 to 100 GeV.

    F, X, Y are actually active in parts of the interval; this is NOT the
    valid low-energy EFT trajectory. Use a consistent SAME-UV-scale initial
    SM state, not an unrelated historical frozen initial condition.
    """
    from Numerical.running.SMWeinbergEvolution import evolve_sm_weinberg
    from RGE.running.weinberg.WeinbergRGE import sm_weinberg_beta,LOOP_FACTOR
    require(1e5>low_scale>0 and 1e-12<rtol<1e-3,'invalid proxy controls')
    c5=unpack_matrix(one_scale['candidate_ordered_source_C5_3x3_gev_inverse_NOT_PHYSICAL'])
    state=take_sm_initial(uv_source['uv_threshold']['state']['sm'],c5)
    low=evolve_sm_weinberg(state,1e5,low_scale,rtol=rtol,atol=atol)
    require(low.solver_success and np.isfinite(low.K).all() and
            relative(low.K,low.K.T)<1e-12,'SM-only reference RGE integration invalid')
    *_,beta_k=sm_weinberg_beta(gY=state.gY,g2_value=state.g2,g3=state.g3,
                         lambdaH_value=state.lambdaH,ye=state.ye,
                         yu=state.yu,yd=state.yd,K=state.K)
    leading=c5+math.log(low_scale/1e5)*beta_k/LOOP_FACTOR
    vev=float(one_scale['reference_EW_vev_gev'])
    result={'format':FORMAT,'status':'SMOnlyOneScaleReferenceComputedNotSequentialEFT',
            'model':model,'initial_mu_gev':1e5,'low_mu_gev':float(low_scale),
            'vev_gev':vev,'solver_rtol':rtol,'solver_atol':atol,
            'solver_nfev':int(low.nfev),
            'UV_matching_scale_reference':mass_summary(c5,vev),
            'SM_only_100gev_reference':mass_summary(low.K,vev),
            'reference_C5_at_low_gev_inverse':pack_matrix(low.K),
            'one_step_frozen_beta_approximation_error':relative(leading,low.K),
            'relative_C5_change_SM_only_proxy':relative(low.K,c5),
            'source_SM_state_mu_is_same_as_C5_mu':True,
            'actual_active_F_X_Y_in_interval_ignored':True,
            'physical_neutrino_prediction':False,
            'complete_EFT1_EFT2_running':False,
            'finite_matching_included':False}
    return result


def blocked_physical_gates()->dict:
    return {
        'verified_within_scope':[
          'all-running UV parameters sampled at same common endpoint',
          'restricted direct versus sequential full-mass triangle',
          'exact mixed-scalar quadratic mass diagonalization and propagator',
          'native seven extra-quartic tensor transport',
          'SM+Weinberg numerical evolution as isolated reference'],
        'not_implemented_in_this_campaign':[
          'source action Matchete-RGBeta full invariant and finite normalization dictionary',
          'live Matchete 4-spinor complex heavy-flavour native equality',
          'new native one-loop Match with independent nonzero mixed quadratic',
          'physical F/Y/X finite field, parameter and EOM threshold matching',
          'Wilsonian EFT1/EFT2 running on the SAME new five-mass trajectory',
          'complete all-16 one-loop general Yukawa matching'],
        'physical_C5_certified':False,
        'same_trajectory_sequential_matching_performed':False,
    }
