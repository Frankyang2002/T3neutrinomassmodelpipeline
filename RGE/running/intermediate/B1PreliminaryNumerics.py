"""B_-1 source-backed lambdaT3-only pilot; explicitly conditional and incomplete.

Inputs are existing repository artifacts. No physical source GammaCC/Weyl factor,
Matchete/RGBeta quartic dictionary or full EFT1 running is inferred.
UV -> all heavy fermions at muF -> EFT1(S1,S2) -> both scalars at muS -> SMEFT.
"""
from __future__ import annotations

from hashlib import sha256
from math import isfinite, log, pi
from pathlib import Path
import json
import numpy as np

FORMAT = 't3-b1-lambdaT3-native-source-preliminary-v1'
LOOP = 1 / (16 * pi**2)
MODEL = 'B_-1'


class PilotError(ValueError):
    pass


def require(ok: bool, why: str):
    if not ok:
        raise PilotError(why)


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def digest(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def complex_array(obj, shape):
    a = np.asarray(obj['real'], dtype=float) + 1j*np.asarray(obj['imag'], dtype=float)
    require(a.shape == shape and np.isfinite(a).all(), f'invalid complex matrix shape {shape}')
    return a


def cjson(a):
    x = np.asarray(a, dtype=complex)
    require(np.isfinite(x).all(), 'nonfinite complex result')
    return {'real': x.real.tolist(), 'imag': x.imag.tolist()}


def real_lift_b1(U):
    """s_original = R s_XY, with S2=-eps (U21 X+U22 Y)*.

    Real order per doublet: (component1 R,I,component2 R,I). Both
    matrices use complex fields = (real+i imaginary)/sqrt(2).
    """
    U = np.asarray(U, complex)
    require(U.shape == (2, 2) and np.isfinite(U).all() and
            np.allclose(U.conj().T @ U, np.eye(2), atol=2e-11, rtol=2e-11),
            'mixed-doublet U must be unitary 2x2')
    eps = np.array([[0, 1], [-1, 0]], dtype=complex)
    R = np.zeros((8, 8), dtype=float)
    for col in range(8):
        phi = np.zeros((2, 2), dtype=complex)  # X,Y complex doublets
        phi[col // 4, (col % 4) // 2] = 1 if col % 2 == 0 else 1j
        s1 = U[0, 0]*phi[0] + U[0, 1]*phi[1]
        s2 = -eps @ (U[1, 0]*phi[0]+U[1, 1]*phi[1]).conj()
        for i, z in enumerate(np.r_[s1, s2]):
            R[2*i, col] = z.real
            R[2*i+1, col] = z.imag
    require(np.allclose(R.T @ R, np.eye(8), atol=2e-11, rtol=2e-11),
            'complex mixed-doublet map is not real-orthogonal')
    return R


def mass_original_from_eigenbasis(U_s, squared):
    """Use actual muS scalar eigenvalues, transporting both full doublets."""
    vals = np.asarray(squared, dtype=float)
    require(vals.shape == (2,) and np.isfinite(vals).all() and np.min(vals) > 0,
            'muS scalar mass eigenvalues must be two positive numbers')
    R_s = real_lift_b1(U_s)
    out = R_s @ np.diag(np.repeat(vals, 4)) @ R_s.T
    require(np.allclose(out, out.T, atol=1e-7, rtol=1e-12), 'real scalar mass not symmetric')
    require(np.allclose(np.sort(np.linalg.eigvalsh(out)),
                        np.sort(np.repeat(vals,4)), rtol=2e-11, atol=1e-7),
            'full eight-mode scalar eigenvalues mismatch')
    return out


def validate_basis_record(obj, model, mu):
    require(obj.get('format')=='t3-mixed-doublet-mass-basis-transport-v1'
            and obj.get('model')==model and abs(float(obj.get('scale_gev',-1))/mu-1)<1e-11,
            'mixed-doublet source record wrong model, scale, or format')
    require(obj.get('quadratic_hessian_agreement') is True and
            obj.get('max_su2_covariance_absolute_error',1) < 1e-9,
            'source mass basis lacks Hessian/SU2 consistency')
    b=obj['basis']
    require(b.get('basis')=='(S1, epsilon*conjugate(S2)) = U @ (X,Y)' and
            b.get('epsilon')==[[0.0,1.0],[-1.0,0.0]] and
            b.get('unitarity_verified') is True and b.get('diagonal_mass_verified') is True,
            'basis/epsilon convention differs from source')
    U=complex_array(b['unitary_flavour_rotation_U'], (2,2))
    A=complex_array(b['mass_squared_flavour_matrix_gev2'],(2,2))
    vals=np.asarray(b['masses_squared_gev2'],float)
    require(np.allclose(U.conj().T@A@U,np.diag(vals),atol=1e-6,rtol=1e-10),
            'source mass eigenvectors fail to diagonalize mass matrix')
    return U,vals


def q_in_fermion_basis(channels, R):
    """quartic original -> F-source X/Y real basis, with both legs rotated."""
    q=np.asarray(channels,dtype=float)
    require(q.ndim == 5 and q.shape[1:]==(4,4,8,8) and
            np.isfinite(q).all(), 'quartic channel tensor wrong shape')
    require(np.allclose(q,q.swapaxes(1,2),atol=1e-13,rtol=1e-12) and
            np.allclose(q,q.swapaxes(3,4),atol=1e-13,rtol=1e-12),
            'native quartic channels are not symmetric')
    return np.einsum('kABab,ai,bj->kABij',q,R,R,optimize=True)


def evaluate_response(C, quartic, mass2, mu_f, mu_s, *, v=246.22):
    """Source polynomial coefficient per unknown N_Gamma, lambdaT3-only.

    Delegates finite scalar bubble and Weinberg projection to the existing
    joint_scalar_matching/core.py. Never assumes full EFT1 operator closure.
    """
    from Numerical.joint_scalar_matching.core import (
        scalar_bubble_matching, project_weinberg, check_tensors,
    )
    C,q=check_tensors(C,quartic,8)
    require(isfinite(mu_f) and isfinite(mu_s) and 0 < mu_s < mu_f,
            'scales must satisfy 0 < muS < muF')
    beta=np.einsum('fgijab,ABab->fgijAB',C,q,optimize=True)
    beta5, wr=project_weinberg(beta,strict=False)
    bubbles=scalar_bubble_matching(C,q,mass2,mu_s,reference_scale=mu_s)
    bubble5, br=project_weinberg(bubbles['amplitude'],strict=False)
    leading_log=beta5*log(mu_s/mu_f)*LOOP
    result=leading_log+bubble5
    # Gauge-covariant result is necessary to interpret a pure Weinberg C5.
    ward_ok=wr<=1e-8 and br<=1e-8
    # An independent change in matching scale, holding inputs fixed:
    other_mu=mu_s*np.exp(.03125)
    other=scalar_bubble_matching(C,q,mass2,other_mu,reference_scale=mu_s)
    other5,_=project_weinberg(other['amplitude'],strict=False)
    new_sum=beta5*log(other_mu/mu_f)*LOOP+other5
    cancellation=float(np.linalg.norm(new_sum-result)/max(np.linalg.norm(result),
                                 np.linalg.norm(bubble5),np.linalg.norm(leading_log),1e-30))
    require(cancellation <= 2e-9, f'joint scalar matching-scale cancellation failed ({cancellation})')
    flavour_sym=float(np.linalg.norm(result-result.T)/max(np.linalg.norm(result),1e-30))
    require(flavour_sym<=1e-8,'source-induced projected flavour tensor is not symmetric')
    out={'status':'PreliminaryNativeQuarticOneChannel_PerUnknownGammaCC' if ward_ok
          else 'NotWeinbergInvariant_DoNotInterpretAsC5',
         'full_SU2_Weinberg_projection_pass':bool(ward_ok),
         'beta_SU2_residual_relative':float(wr),
         'joint_bubble_SU2_residual_relative':float(br),
         'matching_scale_cancellation_relative':cancellation,
         'beta_C5_master_16pi2_per_N_gev_inverse':cjson(2*beta5),
         'C5_master_leading_log_per_N_gev_inverse':cjson(2*leading_log),
         'C5_master_joint_scalar_bubble_per_N_gev_inverse':cjson(2*bubble5),
         'C5_master_total_projected_per_N_gev_inverse':cjson(2*result),
         'joint_to_master_Weinberg_coefficient_factor':2,
         'scalar_masses_gev':np.sqrt(np.linalg.eigvalsh(mass2)).tolist(),
         'physical_C5_certified':False}
    if ward_ok:
        # (1/2) C5 (L.H)^2 and H0=v/sqrt(2): m_nu=v^2*C5/2.
        # Singular values give Takagi masses for complex symmetric Majorana m.
        eig=np.linalg.svd(.5*v*v*(2*result)*1e9,compute_uv=False)
        out['neutrino_mass_singular_values_ev_per_abs_NGamma']=sorted(float(x) for x in eig)
    else:
        out['neutrino_mass_singular_values_ev_per_abs_NGamma']=None
    return out


def produce(root: Path) -> dict:
    """Use source-backed real B1 inputs, never substitutes UV benchmark at muS."""
    from Numerical.joint_scalar_matching.fermion_boundary_bridge import load_source
    from Numerical.joint_scalar_matching.source_quartic_v4 import (
        build_source_channels, source_seed_path, evaluate_channels,
    )
    root=Path(root)
    sr=root/'output/general_observable/lo_fermion_boundary/releases/source_021dcc1e8fe8ba0f/B_m1/boundary.json'
    sf=root/'output/general_observable/mixed_f_synchronised_v1/B_m1/mass_basis/B_m1_final_mass_basis.json'
    ss=root/'output/general_t3_mass_basis/B_m1_final_mass_basis.json'
    ordered=root/'output/general_observable/ordered_wilson_tensor_v4/B_m1/ordered_tensor.json'
    fparams=root/'output/general_observable/mixed_f_synchronised_v1/B_m1/matchete/B_m1_final_parameters.json'
    quartic=root/'output/general_observable/lo_joint_scalar_matching/source_quartic_v4/native_quartic_channels/B_m1.json'
    source=load_source(sr,expected_model=MODEL)
    require(source.mu_f_gev==100000.0 and not source.normalization_certified,
            'only the documented unnormalized muF boundary is accepted')
    fj=read(sf);sj=read(ss); oj=read(ordered); qj=read(quartic); fp=read(fparams)
    require(oj.get('format')=='t3-mixed-native-matchete-ordered-tensor-v1' and
            oj.get('model')==MODEL and oj.get('mu_F_gev')==source.mu_f_gev,
            'ordered source/MF basis record mismatch')
    # The original Matchete parameter handoff retains the source SHA256.
    # Check that handoff, not a newly serialized copy of the source JSON.
    require(fp.get('format')=='t3-experimental-mass-eigen-matchete-parameters-v1'
            and fp.get('model')==MODEL
            and fp.get('source_eft1_scale_gev')==source.mu_f_gev
            and fp.get('source_sha256')==oj.get('source_mass_basis_sha256')
            and fp.get('source_path')==str(sf.relative_to(root)).replace('\\','/'),
            'Matchete F-source parameters do not identify the recorded muF mass basis')
    require(fj.get('status')=='ExperimentalMassBasisTensorTransportNoPhysicalMatching' and
            sj.get('status')=='ExperimentalMassBasisTensorTransportNoPhysicalMatching',
            'mass records no longer have known provisional status')
    Uf,_=validate_basis_record(fj,MODEL,source.mu_f_gev)
    require(np.allclose(Uf,complex_array(fp['unitary_U'],(2,2)),rtol=0,atol=1e-13),
            'Matchete F-source U differs from the recorded muF eigenbasis')
    mu_s=float(sj['scale_gev']);Us,eval_s=validate_basis_record(sj,MODEL,mu_s)
    require(mu_s==1000.0,'unexpected common scalar threshold for this pilot')
    Rf=real_lift_b1(Uf)
    mass_original=mass_original_from_eigenbasis(Us,eval_s)
    mass_boundary=Rf.T@mass_original@Rf
    require(np.allclose(np.sort(np.linalg.eigvalsh(mass_boundary)),
                        np.sort(np.repeat(eval_s,4)),rtol=2e-11,atol=1e-7),
            'muF X/Y to muS X/Y mass transport changes eigenvalues')
    seed=source_seed_path(root/'output/full16_symbolic',MODEL)
    ch=build_source_channels(seed,MODEL)
    require(qj.get('format')=='t3-source-hhss-quartic-channel-export-v4' and
            qj.get('model')==MODEL and
            qj.get('seed_quartic_sha256')==ch.seed_sha256 and
            qj.get('source_boundary_sha256')==source.hash and
            sorted(ch.by_name)==qj.get('couplings_required_at_muS'),
            'native v4 quartic export does not match the scalar seed/source')
    require(tuple(source.labels)==('NewScalarX',)*4+('NewScalarY',)*4,
            'unexpected F-source scalar component ordering')
    candidate=sj['quartics']['lambdaT3']['source_rgbeta_coupling']
    lam=complex(float(candidate['real']),float(candidate['imag']))
    require(np.isfinite(lam) and abs(lam)>0, 'nonfinite/zero muS lambdaT3 candidate')
    # Deliberate truncation: other native source couplings are NOT claimed zero.
    hypotheses={name:({ 'real':lam.real,'imag':lam.imag} if name=='lambdaT3' else 0.)
                for name in sorted(ch.by_name)}
    raw,names=evaluate_channels(ch,hypotheses)
    qchannels=q_in_fermion_basis(raw,Rf)
    require(names==tuple(sorted(hypotheses)),'source native quartic channels reordered')
    index=names.index('lambdaT3')
    q=qchannels[index]
    resp=evaluate_response(source.raw,q,mass_boundary,source.mu_f_gev,mu_s)
    source_responses=[]
    summed=np.zeros((3,3),complex)
    # Decomposition is about the *projected amplitudes*, not a claim that
    # every individual PL source is separately gauge invariant.
    from Numerical.joint_scalar_matching.core import project_weinberg
    for k,idx in enumerate(source.source_indices):
        c=source.by_source[k]
        beta=np.einsum('fgijab,ABab->fgijAB',c,q,optimize=True)
        proj,_=project_weinberg(beta,strict=False)
        source_responses.append({'matchete_PL_source_index':idx,
                                 'beta_C5_master_16pi2_per_N_gev_inverse':cjson(2*proj)})
        summed+=2*proj
    target=(np.asarray(resp['beta_C5_master_16pi2_per_N_gev_inverse']['real'])+
            1j*np.asarray(resp['beta_C5_master_16pi2_per_N_gev_inverse']['imag']))
    require(np.linalg.norm(summed-target)<=2e-9*max(np.linalg.norm(target),1e-30),
            'per-source beta channels fail to reconstruct the all-source result')
    return {'format':FORMAT,'model':MODEL,
            'status':resp['status'],
            'input_source_sha256':source.hash,
            'mass_basis_muF_current_file_sha256':digest(sf),
            'mass_basis_muF_native_matchete_source_sha256':fp['source_sha256'],
            'native_F_source_parameter_file_sha256':digest(fparams),
            'mass_basis_muS_sha256':digest(ss),
            'native_quartic_seed_sha256':ch.seed_sha256,
            'native_quartic_export_sha256':digest(quartic),
            'mu_F_gev':source.mu_f_gev,'mu_S_gev':mu_s,
            'source_basis':'NewScalarX/NewScalarY at muF (source-matched)',
            'muS_mass_basis':'NewScalarX/NewScalarY at muS (rotated into muF source coordinates)',
            'muS_source_rgbeta_lambdaT3_as_matchete_hypothesis':cjson(lam),
            'native_quartic_channels':list(names),
            'only_included_quartic_channel':'lambdaT3',
            'other_native_quartic_channels_excluded_not_physically_zero':
                [name for name in names if name!='lambdaT3'],
            'unknown_GammaCC_to_Weyl_factor':None,
            'matchete_rgbeta_lambdaT3_dictionary_certified':False,
            'full_EFT1_seven_term_running_included':False,
            'hard_fermion_one_loop_matching_included':False,
            'other_dim5_operators_included':False,
            'joint_scalar_threshold_performed_once':True,
            'both_scalar_multiplets_retained_until_muS':True,
            'source_count':len(source.source_indices),
            'per_source_beta':source_responses,
            'response':resp,
            'physical_C5_computed':False,
            'physical_neutrino_mass_computed':False,
            'interpretation':('Preliminary lambdaT3-only contribution per unit unknown signed '
              'GammaCC-to-Weyl factor; uses actual recorded muF LLSS and '
              'muS scalar mass and RGBeta lambdaT3 candidate. Assigning the '
              'RGBeta lambdaT3 number to the Matchete scalar source is an '
              'explicit uncertified mapping hypothesis. All other quartic '
              'channels are omitted for this diagnostic, not physically zero. '
              'This is not a complete one-loop or physical neutrino mass.')}
