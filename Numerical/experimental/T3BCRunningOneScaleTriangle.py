"""B/C(-1) common-UV-scale one-loop triangle evaluator, NOT physical C5.

This module takes the *same running UV state* for all y1, y2, MF, scalar
masses, mixing, and the three source-selected quartics. It reuses the archived
Matchete direct ordered A5 expression via the caller-provided coefficient
function. Heavy Majorana fields are Takagi-rotated before summation.

Limitations: only archived restricted O(y^2 lambda) one-loop topology, no
independent finite EFT field/EOM/counterterm matching, no EFT1/EFT2 transport,
no global Matchete<->RGBeta action certificate, and no automatic validity for
arbitrary extra mixed scalar quadratic. Never set physical_C5_certified true.
"""
from __future__ import annotations

import math
from typing import Any, Callable, Mapping
import numpy as np

from Numerical.experimental.GeneralT3MixedScalarMassBasis import (
    MixedDoubletMassBasis, quadratic_mass_tensor_transport,
    verify_tensor_transport,
)

MODELS = ('B_-1','C_-1')
MATCH_QUARTICS = ('lambdaT3', 'lambdaT3ExtraQ020020I1',
                  'lambdaT3ExtraQ020200I1')


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError('B/C running one-scale: '+message)


def scalar(z: Any, label: str) -> complex:
    if isinstance(z, dict):
        require(set(z) == {'real', 'imag'},label+': expected real/imag')
        result = complex(float(z['real']),float(z['imag']))
    else:
        result = complex(z)
    require(np.isfinite(result.real) and np.isfinite(result.imag),
            label+': nonfinite')
    return result


def real(z: Any, label: str) -> float:
    v=scalar(z,label)
    require(abs(v.imag)<1e-10*max(1.,abs(v.real)),label+': expected real')
    return float(v.real)


def matrix(z: Any, label: str) -> np.ndarray:
    require(isinstance(z,dict) and set(z)=={'real','imag'},
            label+': expected packed real/imag matrix')
    a=np.asarray(z['real'],float);b=np.asarray(z['imag'],float)
    require(a.shape==(3,3) and b.shape==(3,3),label+': expected 3x3')
    require(np.isfinite(a).all() and np.isfinite(b).all(),label+': nonfinite')
    return a+1j*b


def packed(z: np.ndarray) -> dict:
    a=np.asarray(z,complex)
    require(a.shape==(3,3) and np.isfinite(a).all(),'packed: expected 3x3 finite')
    return {'real':a.real.tolist(),'imag':a.imag.tolist()}


def takagi_majorana(M: np.ndarray, *, symmetry_rtol: float=1e-10,
                    gap_rtol: float=1e-7) -> tuple[np.ndarray,np.ndarray,dict]:
    """Takagi U^T M U = diag(m>0), sorted ascending, no heavy-basis guess.

    For distinct singular values, eigenvectors of M^dagger M are Takagi
    vectors modulo a phase. Fix each by the phase of diag(U^T M U).
    Reject degenerate singular values because flavour vectors are ambiguous.
    """
    A=np.asarray(M,complex)
    require(A.shape==(3,3) and np.isfinite(A).all(),'MF must be finite 3x3')
    norm=max(float(np.linalg.norm(A)),1.)
    require(np.linalg.norm(A-A.T)<symmetry_rtol*norm,
            'MF is not complex symmetric in Majorana convention')
    eig, U=np.linalg.eigh(A.conj().T @ A)
    require(eig[0]>0,'nonpositive heavy Majorana mass')
    masses=np.sqrt(np.maximum(eig,0.))
    require(np.all(np.diff(masses)>gap_rtol*np.maximum(masses[1:],1.)),
            'degenerate heavy singular masses; column basis ambiguous')
    middle=U.T@A@U
    phase=np.exp(-.5j*np.angle(np.diag(middle)))
    V=U@np.diag(phase)
    diagonal=V.T@A@V
    residual=float(np.linalg.norm(diagonal-np.diag(masses))/norm)
    recon=float(np.linalg.norm(A-V.conj()@np.diag(masses)@V.conj().T)/norm)
    require(max(residual,recon)<max(5e-10,10*symmetry_rtol),
            f'Majorana Takagi factorization failed residual={residual:g}')
    return masses,V,{'takagi_residual':residual,
                     'reconstruction_residual':recon,
                     'heavy_mass_basis_is_takagi':True}


def uv_boundary(model: str, report: Mapping[str,Any]) -> dict:
    """Extract all matching inputs from the UV endpoint at muF=1e5 GeV.

    Never obtain y or MF from EFT1 (they were integrated out). Never use the
    separately evolved EFT1 scalar masses with UV Yukawas in a claimed
    common-scale amplitude. The grouped F scale is a candidate boundary.
    """
    require(model in MODELS and report.get('model')==model,
            'wrong model')
    require(report.get('status')=='CandidateCoupledUVAndEFT1FiveRealMassesIntegrated'
            and report.get('solver_success') is True,
            'source is not successful coupled B/C five-real-mass run')
    require(report.get('physical_matchete_rgbeta_mapping_certified') is False
            and report.get('physical_neutrino_mass_prediction') is False,
            'unexpected source certification flags')
    start,finish=report.get('uv_range_gev',[None,None])
    require(start==1e7 and finish==1e5,'unsupported grouped fermion endpoint')
    endpoint=report['uv_threshold'];state=endpoint['state']
    require(math.isclose(float(endpoint['mu_gev']),float(finish),rel_tol=1e-12),
            'UV endpoint scale mismatch')
    require(state['representation']['d_f']==(1 if model=='B_-1' else 3)
            and state['representation']['alpha']==-1,
            'unexpected B/C heavy-fermion representation')
    masses,V,check=takagi_majorana(matrix(state['MF'],'MF'))
    # T3Topology.wl builds y1 with CConj[F], y2 with F. Under an
    # ordinary complex Weyl unitary heavy-field redefinition, those two
    # spinor orientations transform with conjugate heavy matrices. Actual
    # Matchete spinor/phase normalization is NOT yet independently certified.
    # Retain the same-V alternative as a diagnostic, not another prediction.
    y1_original=matrix(state['y1'],'y1')
    y2_original=matrix(state['y2'],'y2')
    y1=y1_original@V.conj()
    y2=y2_original@V
    y1_alternative=y1_original@V
    if not (np.isfinite(y1).all() and np.isfinite(y2).all()):
        raise ValueError('nonfinite heavy-eigenstate Yukawa')
    m1=real(state['mS1Sq'],'mS1Sq');m2=real(state['mS2Sq'],'mS2Sq')
    mix=scalar(endpoint['muT3GeneralS1S2_gev2'],'mixed scalar bilinear')
    basis=MixedDoubletMassBasis.build(model,m1,m2,mix)
    require(not basis.tachyonic and not basis.degenerate,
            'tachyonic or degenerate scalar mass eigenbasis')
    transport=quadratic_mass_tensor_transport(basis)
    source=state
    quartics={
        'lambdaT3':scalar(source['lambdaT3'],'lambdaT3'),
        'lambdaT3ExtraQ020020I1':-scalar(source['lambdaHHdagS2S2'],'lambdaHHdagS2S2'),
        'lambdaT3ExtraQ020200I1':-scalar(source['lambdaHHdagS1barS1bar'],'lambdaHHdagS1barS1bar'),
    }
    # The two negative signs are the existing Stage17 candidate mapping.
    # Their full inter-package physical normalization is NOT yet certified.
    words=[('H','H','S1','S2dag'),('H','H','S1','S1'),
           ('Hdag','Hdag','S2','S2')]
    tensor_errors={','.join(w):verify_tensor_transport(w,basis.U,seed=11,trials=4)
                   for w in words}
    require(masses[0]>float(basis.eigenvalues_sq[1]**.5),
            'MF eigenstate not heavier than both scalars; source ordered hierarchy invalid')
    return {'model':model,'mu_gev':float(finish),'mf':masses,'y1':y1,'y2':y2,
            'y1_same_V_alternative':y1_alternative,
            'scalar_basis':basis,'quartics':quartics,'heavy_basis':V,
            'heavy_checks':check,'scalar_tensor_checks':tensor_errors,
            'mass_tensor_checks':transport,
            'quartic_mapping':'Stage17 conditional: Matchete lambdaT3=RGBeta lambdaT3; extra020020=-lambdaHHdagS2S2; extra020200=-lambdaHHdagS1barS1bar',
            'heavy_spinor_phase_mapping':'conditional T3Topology orientation: y1 couples CConj[F] -> y1@V*, y2 couples F -> y2@V; same-V alternative retained',
            'physical_heavy_spinor_phase_mapping_certified':False,
            'physical_coupling_mapping_certified':False}


def assemble_ordered_c5(masses: np.ndarray, y1: np.ndarray, y2: np.ndarray,
                        mx: float, my: float, U: np.ndarray,quartics: Mapping[str,complex],
                        evaluator: Callable[...,tuple[complex,complex,complex]]) -> tuple[np.ndarray, list[dict]]:
    """Reuse source-v22 direct_ordered for each Takagi heavy mass.

    Source v22 convention: A=(a11,a12,a22) for *barred* Yukawa factors;
    C5 = sum_r (conj(y1_r),conj(y2_r)) [[2*a11,a12],[a12,2*a22]]
         (conj(y1_r),conj(y2_r))^T . No conjugation of the outer C5.
    """
    masses=np.asarray(masses,float)
    y1=np.asarray(y1,complex);y2=np.asarray(y2,complex)
    require(masses.shape==(3,) and y1.shape==(3,3) and y2.shape==(3,3),
            'flavour inputs not three heavy 3x3 matrices')
    require(np.isfinite(masses).all() and np.isfinite(y1).all() and np.isfinite(y2).all(),
            'nonfinite matching inputs')
    require(float(masses[0])>my>mx>0,'one-scale source mass ordering invalid')
    require(np.asarray(U).shape==(2,2) and
            np.linalg.norm(U.conj().T@U-np.eye(2))<1e-10,
            'scalar U not unitary')
    require(set(quartics)==set(MATCH_QUARTICS),'incomplete restricted quartic inventory')
    out=np.zeros((3,3),complex)
    pieces=[]
    for i,mf in enumerate(masses):
        a11,a12,a22=map(complex,evaluator(float(mf),float(mx),float(my),U,quartics))
        K=np.array([[2*a11,a12],[a12,2*a22]],complex)
        yy=np.column_stack((y1[:,i].conj(),y2[:,i].conj()))
        term=yy@K@yy.T
        require(np.isfinite(term).all(),'nonfinite direct source coefficient')
        out+=term
        pieces.append({'heavy_column':i,'MF_gev':float(mf),
                       'ordered_A11_A12_A22_gev_inverse':[[v.real,v.imag] for v in (a11,a12,a22)],
                       'contribution_frobenius_norm_gev_inverse':float(np.linalg.norm(term))})
    require(np.linalg.norm(out-out.T)<=1e-11*max(np.linalg.norm(out),1e-40),
            'C5 ordered flavour symmetrization failed')
    return out,pieces


def evaluate_one_scale(boundary:dict,evaluator:Callable,*, vev_gev:float) -> dict:
    b=boundary['scalar_basis'];b.verify()
    require(math.isfinite(vev_gev) and 100<vev_gev<400,
            'EW vev must be explicitly supplied, in GeV')
    mx,my=np.sqrt(b.eigenvalues_sq)
    c5,terms=assemble_ordered_c5(boundary['mf'],boundary['y1'],boundary['y2'],
                        float(mx),float(my),b.U,boundary['quartics'],evaluator)
    alt_c5,_=assemble_ordered_c5(boundary['mf'],boundary['y1_same_V_alternative'],
                        boundary['y2'],float(mx),float(my),b.U,boundary['quartics'],evaluator)
    spinor_alt_residual=float(np.linalg.norm(c5-alt_c5)/max(np.linalg.norm(c5),1e-250))
    mnu=-.5*vev_gev**2*c5*1e9
    masses=np.sort(np.linalg.svd(mnu,compute_uv=False))[::-1]
    require(bool(np.isfinite(masses).all()),'nonfinite reference neutrino masses')
    return {
      'format':'t3-bc-common-uv-scale-restricted-direct-matchete-triangle-v1',
      'status':'RunningCommonScaleRestrictedTriangleEvaluatedNotPhysicalC5',
      'model':boundary['model'],'matching_input_scale_gev':boundary['mu_gev'],
      'input_scale_note':'every coupling/mass here from SAME grouped UV endpoint; no EFT1 quantities are substituted',
      'muX_gev':float(mx),'muY_gev':float(my),
      'MF_takagi_masses_gev':boundary['mf'].tolist(),
      'MF_takagi_heavy_rotation_3x3':packed(boundary['heavy_basis']),
      'running_heavy_eigenbasis_y1_3x3':packed(boundary['y1']),
      'running_heavy_eigenbasis_y2_3x3':packed(boundary['y2']),
      'MF_takagi_checks':boundary['heavy_checks'],
      'heavy_matrix_rotation':'U_F^T MF U_F=diag(M); candidate y1_mass=y1_UV conj(U_F), y2_mass=y2_UV U_F; Matchete spinor phase NOT certified',
      'heavy_spinor_phase_mapping':boundary['heavy_spinor_phase_mapping'],
      'heavy_spinor_phase_mapping_physically_certified':False,
      'alternative_same_V_spinor_rotation_relative_C5_difference':spinor_alt_residual,
      'scalar_mass_basis':b.info(),
      'scalar_tensor_transport_checks':boundary['scalar_tensor_checks'],
      'mass_tensor_transport_checks':boundary['mass_tensor_checks'],
      'conditional_quartic_dictionary':boundary['quartic_mapping'],
      'conditional_matchete_quartics':{k:{'real':float(z.real),'imag':float(z.imag)}
                       for k,z in boundary['quartics'].items()},
      'full_three_heavy_flavor_matrix_ordered_source_convention':True,
      'heavy_terms':terms,
      'candidate_ordered_source_C5_3x3_gev_inverse_NOT_PHYSICAL':packed(c5),
      'mnu_at_matching_scale_using_reference_EW_vev_3x3_ev_NOT_PHYSICAL':packed(mnu),
      'reference_EW_vev_gev':float(vev_gev),
      'reference_mnu_singular_ev_descending_NOT_PHYSICAL':masses.tolist(),
      'source_matching_scheme':'Matchete direct all-heavy restricted O(y^2 lambda) v22 inputform; hbar=1/(16*pi^2)',
      'physical_C5_certified':False,
      'physical_neutrino_masses_certified':False,
      'full_fermion_X_Y_sequential_RG_transport':False,
      'full_finite_field_parameter_EOM_scheme_certified':False,
      'general_extra_mixed_mass_source_matching_certified':False,
      'native_matchete_rgbeta_action_dictionary_certified':False,
    }
