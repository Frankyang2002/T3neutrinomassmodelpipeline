"""Conditional lambdaT3 one-channel response for the six real-split T3 models.

This module uses the NATIVE Matchete Cartesian phase rotation (NOT the
spherical-basis W matrix) documented in RunSixSelfBilinearEigenMatchete.wl:
  S_j = exp(i theta_F) (R_j + i I_j)/sqrt(2).
The legacy fermion-first source orders the split real multiplets MINUS (I),
PLUS (R). We only accept a positive, real F-boundary mixing input; the actual
F source did not certify a complex phase convention. The S-scale mass Hessian
is transported into that fixed F-source basis, including sign reversals of
mu2; all scalar eigenmodes are integrated out at ONE threshold.

The RGBeta->native Matchete quartic and self-mass normalization is a declared
hypothesis, never a physical certification. The full SU(2) projection decides
whether the conditional amplitude can even be called a Weinberg projection.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import isfinite
from pathlib import Path
from typing import Any

import numpy as np

from RGE.running.intermediate.All16PreliminaryNumerics import (
    FORMAT as PREVIOUS_FORMAT, DATA_ROOT, SPLIT_MODELS, _lambda, _mass_c5,
    _real_from_number, tag,
)
from RGE.running.intermediate.B1PreliminaryNumerics import (
    PilotError, cjson, digest, read, require,
)

FORMAT = 't3-six-self-split-lambdaT3-conditional-v1'
NATIVE_ROTATION_SOURCE = 'Lagrangian/RunSixSelfBilinearEigenMatchete.wl::T3SBRotate'
MASS_INPUT_SOURCE = 'output/general_observable/all16_v1'
SPLIT_SHAPES = {
    'A_-2': (1, 3, 2), 'A_+0': (1, 3, 1),
    'D_-2': (3, 1, 2), 'D_+0': (3, 1, 1),
    'E_-2': (3, 3, 2), 'E_+0': (3, 3, 1),
}


@dataclass(frozen=True)
class SplitScalarPoint:
    model: str
    source_to_original_R: np.ndarray
    mass2_boundary_muS: np.ndarray
    muS: float
    lambda_candidate: complex
    inputs: dict[str, Any]
    status: str


def native_cartesian_lift(d1: int, d2: int, split_species: int, phase: complex = 1.) -> np.ndarray:
    """s_original_real = R @ s_Freal, with F splitting I (minus), R (plus).

    Interleaved original Cartesian real components per complex SU(2) component:
    S_j=(x_j+i*y_j)/sqrt(2).  Split F field ordering: [I_0,...,I_d-1,
    R_0,...,R_d-1] for the selected species, while unsplit scalar species
    retains the original real/imaginary interleaved component layout.

    A complex phase rotates both real fields, but the released F source
    has not independently certified non-real phase conventions; production
    strictly accepts positive-real F holomorphic bilinear and phase=+1.
    """
    require(split_species in (1, 2) and d1 in (1, 3) and d2 in (1, 3),
            'split model must be A/D/E self-bilinear representation')
    phase = complex(phase)
    require(np.isfinite([phase.real, phase.imag]).all() and
            abs(abs(phase) - 1.) <= 1e-12, 'native Cartesian phase must be unit modulus')
    dims=(d1, d2)
    n=2*(d1+d2)
    R=np.zeros((n,n),dtype=float)
    for k,d in enumerate(dims):
        start=2*sum(dims[:k])
        if k+1 != split_species:
            R[start:start+2*d,start:start+2*d]=np.eye(2*d)
            continue
        for j in range(d):
            a,b=start+2*j,start+2*j+1
            i_idx,r_idx=start+j,start+d+j
            # z_original = phase*(R_source+i I_source)/sqrt(2),
            # and both sides use identical sqrt(2) real-field normalisation.
            R[a,r_idx], R[b,r_idx] = phase.real,phase.imag
            R[a,i_idx], R[b,i_idx] = -phase.imag,phase.real
    require(np.linalg.norm(R.T@R-np.eye(n),ord=np.inf)<2e-12 and
            np.linalg.norm(R@R.T-np.eye(n),ord=np.inf)<2e-12,
            'native self-bilinear Cartesian lift is not real orthogonal')
    return R


def native_original_mass2(d1: int,d2: int,split_species: int,
                          m1:float,m2:float,mu2:complex)->np.ndarray:
    """Hessian in original interleaved real Cartesian basis, all fields active.

    H_{R,I} = [[m^2+2 Re mu2, 2 Im mu2],
               [2 Im mu2, m^2-2 Re mu2]],
    for V=m^2 Sdag S + (mu2 (Sdag)^T Sdag + h.c.).
    NO extra 1/2 on holomorphic term. This is a conditional Matchete/RGBeta map.
    """
    d1=int(d1);d2=int(d2)
    require((d1,d2,split_species) in {(1,3,2),(1,3,1),(3,1,2),(3,1,1),(3,3,2),(3,3,1)},
            'unsupported real Cartesian scalar representation')
    m1=float(m1);m2=float(m2);mu2=complex(mu2)
    require(all(isfinite(x) for x in (m1,m2,mu2.real,mu2.imag)),
            'nonfinite mass or self-bilinear parameter')
    n=2*(d1+d2)
    mass=np.zeros((n,n),float)
    for k,(d,m) in enumerate(((d1,m1),(d2,m2))):
        start=2*sum((d1,d2)[:k])
        for j in range(d):
            p=start+2*j
            if k+1==split_species:
                mass[p:p+2,p:p+2]=[[m+2*mu2.real,2*mu2.imag],
                                  [2*mu2.imag,m-2*mu2.real]]
            else: mass[p:p+2,p:p+2]=np.eye(2)*m
    eig=np.linalg.eigvalsh(mass)
    require(eig[0]>0., 'tachyonic/nonpositive scalar mass at common muS')
    self_mass=(m1,m2)[split_species-1]
    expected=sorted(([m1]*2*d1 if split_species==2 else
                     [m2]*2*d2) +
                    [self_mass-2*abs(mu2)]*(d1 if split_species==1 else d2)+
                    [self_mass+2*abs(mu2)]*(d1 if split_species==1 else d2))
    require(np.allclose(eig,expected,rtol=1e-11,atol=1e-7),
            'self-bilinear Cartesian Hessian spectrum does not match m2 +/- 2|mu2|')
    return mass


def split_scalar_point(root: Path, model: str, source: Any, ch:Any)->SplitScalarPoint:
    """Read existing illustrative RGBeta trajectory, NOT a physical muS endpoint."""
    require(model in SPLIT_MODELS, 'model is not a split-real self-bilinear case')
    d1,d2,split=SPLIT_SHAPES[model]
    root=Path(root)
    path=root/MASS_INPUT_SOURCE/tag(model)/'running.json'
    j=read(path)
    require(j.get('model')==model and j.get('status')=='NumericalRenormalisableRunningOnly' and
            j.get('illustrative_benchmark') is True and
            j.get('general_matching_performed') is False and
            j.get('physical_mass_parameters_available') is False and
            j.get('eft1_wilson_running_performed') is False,
            f'{model}: illustrative scalar-running status or provenance changed')
    scales=j['scales_gev']
    require(float(scales['F'])==float(source.mu_f_gev) and
            0<float(scales['S'])<float(scales['F']),
            f'{model}: inconsistent F/S matching scales')
    f,s=j['eft1_at_fermion_threshold'],j['eft1_at_scalar_threshold']
    key='muT3GeneralS1Sq' if split==1 else 'muT3GeneralS2Sq'
    def getmu(v):
        value=v[key]
        if isinstance(value,dict):
            require(set(value)>={'real','imag'},'mu2 lacks real/imag fields')
            return complex(float(value['real']),float(value['imag']))
        return complex(value)
    muf,mus=getmu(f),getmu(s)
    require(all(isfinite(x) for x in (muf.real,muf.imag,mus.real,mus.imag)),
            f'{model}: invalid mu2 trajectory')
    # The released F split source explicitly states complex phases unsupported.
    # A nonreal/negative F mu2 would invalidate its minus/plus source labeling.
    require(muf.real>0 and abs(muf.imag)<=1e-12*max(abs(muf.real),1.),
            f'{model}: F boundary minus/plus source not tied to a positive-real mu2; '
            'native phase identification requires independent certification')
    x1=_real_from_number(s['mS1Sq'],f'{model}: mS1Sq')
    x2=_real_from_number(s['mS2Sq'],f'{model}: mS2Sq')
    mass_original=native_original_mass2(d1,d2,split,x1,x2,mus)
    R=native_cartesian_lift(d1,d2,split)
    mass_source=R.T@mass_original@R
    require(tuple(source.labels)==(tuple(('NewScalar1',)*(2*d1) if split==2 else
            ('T3Scalar1Minus',)*d1+('T3Scalar1Plus',)*d1) +
            tuple(('NewScalar2',)*(2*d2) if split==1 else
            ('T3Scalar2Minus',)*d2+('T3Scalar2Plus',)*d2)),
            f'{model}: F-source real field order is not the native-minus/plus assumption')
    require(tuple(ch.source_basis_labels)==('NewScalar1',)*(2*d1)+('NewScalar2',)*(2*d2),
            f'{model}: Matchete quartic original real-component order changed')
    basis_doc=root/DATA_ROOT/'fermion_first_self_fix_v2'/tag(model)/'self_mass_basis.json'
    basis=read(basis_doc)
    require(basis.get('mass_basis_kind')=='RealSelfBilinearSplitLightFields' and
            basis.get('split_species')==split and
            basis.get('complex_self_bilinear_phase_supported') is False and
            basis.get('field_basis_normalization_physically_certified') is False and
            basis.get('mass_relation')=='Mminus^2=MS^2-2 Re(kappa); Mplus^2=MS^2+2 Re(kappa)',
            f'{model}: F-source split field convention/provisional status changed')
    require(int(basis.get('split_species'))==split and
            basis.get('minus_field')==f'T3Scalar{split}Minus' and
            basis.get('plus_field')==f'T3Scalar{split}Plus',
            f'{model}: F-boundary split field names inconsistent')
    return SplitScalarPoint(model,R,mass_source,float(scales['S']),_lambda(s['lambdaT3']),
          {'muS_running_json_path':str(path.relative_to(root)).replace('\\','/'),
           'muS_running_json_sha256':digest(path),
           'F_self_bilinear_parameter_RGBeta':cjson(muf),
           'S_self_bilinear_parameter_RGBeta':cjson(mus),
           'F_to_S_self_bilinear_sign_changed':bool(muf.real*mus.real<0.),
           'muF_source_real_split_metadata':str(basis_doc.relative_to(root)).replace('\\','/'),
           'muF_source_real_split_metadata_sha256':digest(basis_doc),
           'native_cartesian_source_rotation':NATIVE_ROTATION_SOURCE,
           'native_mass_bilinear_hypothesis':'mu2_RGBeta equals native Cartesian potential mu2; 2|mu2| mass splitting',
           'source_split_basis_hypothesis':'native Cartesian I=minus, R=plus at muF with unit phase',
           'muF_matched_to_this_running_benchmark':False,
           'quartic_component_map_certified':False,
           'GammaCC_weyl_normalization_certified':False},
           'IllustrativeSelfBilinearCartesianHypothesis_NotCertified')


def evaluate_split_model(root: Path,model:str)->dict:
    """Source-resolved six-case conditional response; no output on failed Ward test."""
    from Numerical.joint_scalar_matching.fermion_boundary_bridge import load_source
    from Numerical.joint_scalar_matching.source_quartic_v4 import (
        build_source_channels, source_seed_path, evaluate_channels,
    )
    from Numerical.joint_scalar_matching.core import project_weinberg
    root=Path(root)
    require(model in SPLIT_MODELS,f'{model}: not a split-real model')
    sr=root/DATA_ROOT/'lo_fermion_boundary/releases/source_021dcc1e8fe8ba0f'/tag(model)/'boundary.json'
    native=root/DATA_ROOT/'lo_joint_scalar_matching/source_quartic_v4/native_quartic_channels'/f'{tag(model)}.json'
    source=load_source(sr,expected_model=model)
    require(source.mu_f_gev==100000. and source.normalization_certified is False,
            f'{model}: source normalization/status changed')
    seed=source_seed_path(root/'output/full16_symbolic',model)
    ch=build_source_channels(seed,model)
    record=read(native)
    require(record.get('format')=='t3-source-hhss-quartic-channel-export-v4' and
            record.get('model')==model and
            record.get('seed_quartic_sha256')==ch.seed_sha256 and
            record.get('source_boundary_sha256')==source.hash and
            record.get('couplings_required_at_muS')==sorted(ch.by_name) and
            record.get('source_basis_matches_fermion_boundary') is False,
            f'{model}: native original-field quartic source mismatch')
    point=split_scalar_point(root,model,source,ch)
    names=tuple(sorted(ch.by_name))
    require('lambdaT3' in names,f'{model}: source lacks topology lambdaT3')
    # All other quartics excluded as an explicit channel truncation, not zero physics.
    channels,provided=evaluate_channels(ch,{name: cjson(point.lambda_candidate) if name=='lambdaT3'
                                else 0. for name in names})
    require(provided==names and channels.shape==(len(names),4,4,ch.n,ch.n),
            f'{model}: native quartic components malformed')
    qorig=channels[names.index('lambdaT3')]
    R=point.source_to_original_R
    qsource=np.einsum('ABab,ai,bj->ABij',qorig,R,R,optimize=True)
    response=_mass_c5(source.raw,qsource,point.mass2_boundary_muS,source.mu_f_gev,point.muS)
    # Source-level flavour decomposition without erasing suppressed/cancelled terms.
    proj_total=np.zeros((3,3),complex)
    per_source=[]
    for i,idx in enumerate(source.source_indices):
        part=np.einsum('fgijab,ABab->fgijAB',source.by_source[i],qsource,optimize=True)
        proj,residual=project_weinberg(part,strict=False)
        proj_total+=2*proj
        per_source.append({'PL_source_index':int(idx),
            'beta_C5_master_per_N_gev_inverse':cjson(2*proj),
            'source_SU2_projection_residual':float(residual)})
    expected=np.asarray(response['beta_C5_master_16pi2_per_N_gev_inverse']['real'])+1j*\
             np.asarray(response['beta_C5_master_16pi2_per_N_gev_inverse']['imag'])
    require(np.linalg.norm(proj_total-expected)<=2e-9*max(np.linalg.norm(expected),1e-30),
            f'{model}: source-resolved beta does not reproduce total')
    return {'format':FORMAT,'model':model,'status':response['status'],
        'source_point_class':point.status,'fermion_boundary_sha256':source.hash,
        'native_quartic_seed_sha256':ch.seed_sha256,
        'native_quartic_export_sha256':digest(native),'muF_gev':source.mu_f_gev,
        'muS_gev':point.muS,'source_count':len(source.source_indices),
        'n_real_heavy_scalars':len(source.labels),'all_native_quartic_channels':list(names),
        'included_native_quartic_channel':'lambdaT3',
        'other_quartic_channels_excluded_not_zero':[x for x in names if x!='lambdaT3'],
        'lambdaT3_RGBeta_as_Matchete_hypothesis':cjson(point.lambda_candidate),
        'cross_workstream_input_provenance':point.inputs,
        'original_to_F_source_real_R':R.tolist(),
        'scalar_mass_squared_muS_in_F_source_real_basis':point.mass2_boundary_muS.tolist(),
        'per_PL_source_beta':per_source,'response':response,
        'unknown_GammaCC_factor':None,'Matchete_RGBeta_coupling_normalization_certified':False,
        'self_bilinear_component_map_physically_certified':False,
        'physical_C5_computed':False,'full_eft1_seven_term_running_included':False,
        'hard_F_one_loop_matching_included':False,
        'both_scalars_integrated_together_at_muS':True,
        'interpretation':('Conditional self-bilinear Cartesian-rotation hypothesis; RGBeta mu2 and lambdaT3 '
          'identified with native Matchete potential couplings only for a diagnostic. '
          'muF split phase and minus/plus field identification are not independently certified; '
          'one quartic only; no GammaCC absolute factor or complete hard fermion matching. '
          'Not a physical neutrino mass.')}
