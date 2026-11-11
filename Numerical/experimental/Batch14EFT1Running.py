"""Candidate Batch 14 F-threshold handoff and full scalar EFT1 RGBeta running.

Project path: Numerical/experimental/Batch14EFT1Running.py

Uses the saved Batch13 final UV state, existing tree-level F projection and
canonical EFT1 state vector. Both S1 and S2 stay active throughout this step.
This is NOT a physical F-threshold Wilson match or scalar mass mixing solution.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Mapping

import numpy as np
from scipy.integrate import solve_ivp

from Numerical.experimental.Batch13GeneralUVRunning import (
    EXPECTED_CONVENTION, EXPECTED_FORMAT, EXPECTED_STATUS,
    illustrative_initial_state, validate_uv_export, REPRESENTATIONS,
)
from Numerical.core.StateVector import pack_uv_state, unpack_uv_state
from Numerical.core.RGBetaEvaluator import evaluate_inputform_expression
from Numerical.core.BetaVector import beta_values_to_derivative
from Numerical.running.FermionThresholdBoundary import (
    build_intermediate_scalar_boundary, fermion_threshold_masses,
)
from Numerical.running.IntermediateScalarStateVector import (
    pack_intermediate_scalar_state, unpack_intermediate_scalar_state,
)
from Numerical.running.IntermediateScalarRGBetaEvaluator import (
    intermediate_scalar_state_environment,
)

FERMION_KEYS = frozenset({'y1','y2','MF','y1Conj','y2Conj'})
SCALAR_GENERAL = frozenset({'muT3GeneralS1Sq', 'muT3GeneralS2Sq', 'muT3GeneralS1S2'})


def load_saved_uv_endpoint(path: Path | str, uv_payload: Mapping[str,Any], model: str):
    """Deserialize all real/imag parts without reseeding the evolved endpoint."""
    validate_uv_export(uv_payload, model)
    document = json.loads(Path(path).read_text(encoding='utf-8-sig'))
    if document.get('model') != model or document.get('status') != 'CandidateUVIntegrated':
        raise ValueError(f'{model}: invalid Batch13 UV endpoint provenance')
    if document.get('restricted_yukawa') is not True or document.get('matchete_rgbeta_normalisation_certified') is not False:
        raise ValueError(f'{model}: unexpected UV physical-status flags')
    mu = float(document['mu_final_gev'])
    if not math.isfinite(mu) or mu <= 0:
        raise ValueError(f'{model}: invalid final UV scale')
    source = document.get('final')
    if not isinstance(source, dict) or set(source) != set(uv_payload['report_betas']):
        raise ValueError(f'{model}: UV endpoint coupling inventory mismatch')
    # The illustrative state fixes only the canonical layout and shapes;
    # every value is overwritten by the serialised, evolved endpoint.
    template = illustrative_initial_state(uv_payload, model, mu_gev=mu)
    _, layout = pack_uv_state(template)
    vec = np.zeros(layout.size, dtype=float)
    for block in layout.blocks:
        raw = source[block.name]
        if not isinstance(raw, dict) or set(raw) != {'real','imag'}:
            raise ValueError(f'{model}: invalid encoded coupling {block.name}')
        real, imag = np.asarray(raw['real'],dtype=float), np.asarray(raw['imag'],dtype=float)
        if real.shape != block.shape or imag.shape != block.shape:
            raise ValueError(f'{model}: incorrect encoded shape for {block.name}')
        if not np.isfinite(real).all() or not np.isfinite(imag).all():
            raise ValueError(f'{model}: non-finite UV endpoint {block.name}')
        if block.kind == 'real_scalar':
            if float(imag) != 0:
                raise ValueError(f'{model}: real state variable {block.name} has imaginary part')
            vec[block.start] = float(real)
        elif block.kind == 'complex_scalar':
            vec[block.start:block.stop] = [float(real),float(imag)]
        elif block.kind == 'complex_matrix':
            n = real.size
            vec[block.start:block.start+n] = real.ravel()
            vec[block.start+n:block.stop] = imag.ravel()
        else:
            raise ValueError(f'{model}: invalid block kind {block.kind}')
    return unpack_uv_state(vec, layout, mu_gev=mu)


def validate_eft1_export(payload: Mapping[str,Any], uv_payload: Mapping[str,Any], model: str) -> tuple[str,...]:
    """Enforce identical scalar operator inventory across F decoupling."""
    validate_uv_export(uv_payload,model)
    if payload.get('format') != EXPECTED_FORMAT or payload.get('status') != EXPECTED_STATUS:
        raise ValueError(f'{model}: unsupported EFT1 export format/status')
    meta = payload.get('metadata')
    if not isinstance(meta,dict):
        raise ValueError(f'{model}: missing EFT1 metadata')
    family,alpha_string = model.split('_',1)
    d1,d2,df = REPRESENTATIONS[family]
    required = {'GeneralScalarModel':model, 'GeneralScalarStage':'EFT1',
                'dS1':d1,'dS2':d2,'dF':df,'alpha':int(alpha_string),
                'ReportBetaConvention':EXPECTED_CONVENTION,
                'GeneralScalarBasisStatus':'CandidateNativeRGBetaUnmappedToMatchete'}
    for k,v in required.items():
        if meta.get(k) != v:
            raise ValueError(f'{model}: EFT1 metadata {k} mismatch: {meta.get(k)!r} != {v!r}')
    if payload.get('physical_mapping_certified') is not False:
        raise ValueError(f'{model}: EFT1 mapping claimed certified')
    betas = payload.get('report_betas')
    if not isinstance(betas,dict) or not betas or not all(isinstance(x,str) and x.strip() for x in betas.values()):
        raise ValueError(f'{model}: missing EFT1 beta expressions')
    uv_keys = set(uv_payload['report_betas'])
    eft_keys = set(betas)
    if eft_keys & FERMION_KEYS:
        raise ValueError(f'{model}: heavy fermion beta functions present in EFT1')
    if eft_keys != uv_keys - FERMION_KEYS:
        raise ValueError(f'{model}: scalar EFT1/UV inventory mismatch: missing={sorted((uv_keys-FERMION_KEYS)-eft_keys)}, extra={sorted(eft_keys-(uv_keys-FERMION_KEYS))}')
    if set(meta.get('GeneralScalarNewBetaKeys',[])) - eft_keys:
        raise ValueError(f'{model}: missing general EFT1 quartic betas')
    if meta.get('QuadraticKind') != 'none' and meta.get('QuadraticMassSymbol') not in eft_keys:
        raise ValueError(f'{model}: missing additional quadratic beta')
    return tuple(betas)


def run_eft1(uv_state, uv_payload: Mapping[str,Any], eft1_payload: Mapping[str,Any],
             model: str, *, mu_scalar_gev: float=1e3, rtol: float=1e-7,
             atol: float=1e-10, max_nfev: int=50000) -> dict[str,Any]:
    """Tree-level scalar handoff at mu_F, followed by renormalisable EFT1 RGE.

    Does not compute the F-generated Wilson tensor, the scalar mass basis,
    or a neutrino-mass prediction. All threshold corrections are outside scope.
    """
    validate_eft1_export(eft1_payload,uv_payload,model)
    uv_state=uv_state.validated()
    mu_f=float(uv_state.mu_gev)
    mu_s=float(mu_scalar_gev)
    if not (np.isfinite(mu_s) and mu_f>mu_s>0):
        raise ValueError('Require mu_F > mu_S > 0')
    fermions=fermion_threshold_masses(uv_state)
    if fermions.minimum_mass_gev <= mu_f:
        raise ValueError(f'{model}: F-first threshold hierarchy violated')
    start=build_intermediate_scalar_boundary(uv_state,matching_scale_gev=mu_f)
    y0,layout=pack_intermediate_scalar_state(start)
    if set(layout.names)!=set(eft1_payload['report_betas']):
        raise ValueError(f'{model}: threshold projected EFT1 state differs from beta inventory')
    if set(start.general_couplings)!=set(uv_state.general_couplings):
        raise ValueError(f'{model}: F threshold dropped an extra scalar coupling')
    beta_text=eft1_payload['report_betas']
    nfev=0
    def rhs(t,y):
        nonlocal nfev
        nfev+=1
        if nfev>max_nfev:
            raise RuntimeError(f'{model}: EFT1 maximum beta evaluations exceeded')
        state=unpack_intermediate_scalar_state(y,layout,mu_gev=math.exp(t))
        env=intermediate_scalar_state_environment(state)
        values={}
        for k,v in beta_text.items():
            try:
                values[k]=evaluate_inputform_expression(v,env)
            except Exception as exc:
                raise ValueError(f'{model}: cannot evaluate EFT1 beta {k}: {exc}') from exc
        return beta_values_to_derivative(values,layout,divide_by_loop_factor=True)
    rhs(math.log(mu_f),y0)
    result=solve_ivp(rhs,(math.log(mu_f),math.log(mu_s)),y0,
                     method='DOP853',rtol=rtol,atol=atol)
    if not result.success:
        raise RuntimeError(f'{model}: EFT1 integration failed: {result.message}')
    final=unpack_intermediate_scalar_state(result.y[:,-1],layout,mu_gev=mu_s)
    return {'model':model,'status':'CandidateFThresholdAndEFT1Integrated',
            'boundary':start,'state':final,'nfev':int(result.nfev),
            'n_couplings':len(layout.names),'n_parameters_real':int(layout.size),
            'mu_f_gev':mu_f,'mu_s_gev':mu_s,
            'fermion_singular_masses_gev':fermions.masses_gev.tolist(),
            'both_scalars_active_in_eft1':True,
            'restricted_yukawa':True,
            'fermion_threshold_wilson_matching_performed':False,
            'general_scalar_mass_matrix_diagonalised':False,
            'matchete_rgbeta_normalisation_certified':False,
            'physical_neutrino_mass_prediction':False}


def serialise_eft1_state(state) -> dict[str,Any]:
    from Numerical.running.IntermediateScalarRGBetaEvaluator import intermediate_scalar_state_environment
    output={}
    for name, value in intermediate_scalar_state_environment(state).items():
        arr=np.asarray(value,dtype=complex)
        output[name]={'real':arr.real.tolist(), 'imag':arr.imag.tolist()}
    return output
