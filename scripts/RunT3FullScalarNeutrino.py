"""Numerical C5 -> low-energy neutrino continuation from continuous T3 scalar RGEs.

This stage DOES NOT fabricate a general mass-basis matching formula:
  * B(-1), C(-1): existing mass-eigenstate direct one-loop Matchete A5,
    using the actual evolved UV Yukawa/MF state and EFT1 scalar endpoint.
  * Eight models with no extra holomorphic quadratic: existing model-specific
    hierarchical (or designated legacy) matching and low-energy SM+K running.
  * Six models with nonzero self-bilinears: evaluate the old no-mixing C5
    ONLY as an explicitly labelled control. No nonzero-mixing C5 is claimed.

The extra quartics affect all native UV/EFT1 beta functions. For B/C the two
C5-active extra quartics are included using the source-level sign mapping.
No general-scalar fermion-first Wilson matching or arbitrary mass-basis
Matchete/RGBeta normalization is certified. Outputs never modify old fits.

    python -m scripts.RunT3FullScalarNeutrino --all --run-id continuous_v1
    python -m scripts.RunT3FullScalarNeutrino --all --run-id continuous_v1 --resume
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import fields, replace
from hashlib import sha256
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from scripts import RunT3FullScalarContinuous as full

ROOT = full.ROOT
MODELS = full.MODELS
SOURCE_ROOT = full.OUTPUT
FORMAT = 't3-full-scalar-continuous-neutrino-candidate-v1'
LOW_SCALE = 100.0
SCALAR_SCALE = 1000.0
MIXED = frozenset(('B_-1', 'C_-1'))
SELF = frozenset(('A_+0','A_-2','D_+0','D_-2','E_+0','E_-2'))


def require(value: bool, message: str) -> None:
    if not value:
        raise ValueError('Continuous T3 neutrino: ' + message)


def load(path: Path) -> dict:
    require(path.is_file(), 'missing '+str(path))
    data = json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(data, dict), 'not JSON object '+str(path))
    return data


def fingerprint(path: Path) -> str:
    require(path.is_file(), 'missing source '+str(path))
    return sha256(path.read_bytes()).hexdigest()


def atomic_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name+'.tmp')
    tmp.write_text(json.dumps(payload,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    tmp.replace(path)


def complex_pair(raw: Any, label: str) -> complex:
    require(isinstance(raw,dict) and set(raw)=={'real','imag'}, label+': missing complex real/imag')
    z=complex(float(raw['real']),float(raw['imag']))
    require(math.isfinite(z.real) and math.isfinite(z.imag), label+': nonfinite')
    return z


def complex_matrix(raw: Any, label: str) -> np.ndarray:
    require(isinstance(raw,dict) and set(raw)=={'real','imag'},label+': invalid complex matrix')
    matrix=np.asarray(raw['real'],dtype=float)+1j*np.asarray(raw['imag'],dtype=float)
    require(matrix.shape==(3,3) and np.isfinite(matrix.real).all() and
            np.isfinite(matrix.imag).all(),label+': expected finite 3x3')
    return matrix


def matrix_payload(matrix: np.ndarray) -> dict:
    m=np.asarray(matrix,dtype=complex)
    require(m.shape==(3,3) and np.isfinite(m.real).all() and np.isfinite(m.imag).all(),
            'nonfinite 3x3 matrix')
    return {'real':m.real.tolist(),'imag':m.imag.tolist()}


def validated_c5(matrix: Any) -> np.ndarray:
    k=np.asarray(matrix,dtype=complex)
    require(k.shape==(3,3) and np.isfinite(k.real).all() and np.isfinite(k.imag).all(),
            'C5 requires finite 3x3 matrix')
    require(np.allclose(k,k.T,rtol=1e-10,atol=1e-18),
            'C5 must be flavour symmetric')
    return 0.5*(k+k.T)


def decode_uv(raw: dict):
    """Reconstruct the exact saved continuous UV state (not a fresh fitted seed)."""
    from Numerical.core.State import T3Representation, SMNumericalState, T3UVState
    require(isinstance(raw,dict), 'missing UV threshold state')
    rep=T3Representation(**raw['representation']).validated()
    sm=raw['sm']
    smv=SMNumericalState(gY=sm['gY'],g2=sm['g2'],g3=sm['g3'],
        lambdaH=sm['lambdaH'],yu=complex_matrix(sm['yu'],'yu'),
        yd=complex_matrix(sm['yd'],'yd'),ye=complex_matrix(sm['ye'],'ye')).validated()
    values={}
    for f in fields(T3UVState):
        n=f.name
        if n in ('representation','sm'):
            continue
        if n in ('y1','y2','MF'):
            values[n]=complex_matrix(raw[n],n)
        elif n=='lambdaT3' or n.startswith('lambdaHHdag') or n in {
                'lambdaS1bar2S2bar2','lambdaS1barS2S2bar2','lambdaS1S1bar2S2bar'}:
            values[n]=None if raw[n] is None else complex_pair(raw[n],n)
        else:
            values[n]=raw[n]
    state=T3UVState(representation=rep,sm=smv,**values).validated()
    require(math.isclose(state.mu_gev,full.s30.STOP,rel_tol=1e-10),
            'UV endpoint is not fermion threshold')
    return state


def decode_eft(raw: dict):
    """Reuse existing exact EFT decoder (its 1e5 GeV precondition is bypassed only for the saved 1e3 endpoint)."""
    from scripts import IntermediateEFTRunning as s34
    require(isinstance(raw,dict) and math.isclose(float(raw.get('mu_gev',0)),SCALAR_SCALE,rel_tol=1e-9),
            'EFT1 state not saved at scalar endpoint')
    source=dict(raw,mu_gev=float(s34.START))
    state=s34.decode_intermediate_state(source)
    return replace(state,mu_gev=SCALAR_SCALE).validated()


def verify_continuous(model: str, path: Path) -> dict:
    data=load(path)
    require(data.get('format')==full.FORMAT and data.get('status')==full.STATUS and
            data.get('model')==model and
            data.get('physical_matchete_rgbeta_mapping_certified') is False and
            data.get('full_sequential_neutrino_mass_prediction') is False and
            data.get('fermion_threshold_wilson_matching_performed') is False and
            data.get('scalar_threshold_wilson_matching_performed') is False,
            'continuous running source schema or certification changed')
    require(data.get('scalar_endpoint_physically_admissible') is True and
            data.get('scalar_endpoint_spectrum',{}).get('tachyonic') is False,
            'tachyonic scalar endpoint')
    require(data.get('quadratic_kind')==('mixed' if model in MIXED else
                    'none' if model not in SELF else
                    'self_s1' if model in {'A_+0','D_+0','E_+0'} else 'self_s2'),
            'quadratic operator classification changed')
    require(math.isclose(data['settings']['scalar_endpoint_gev'],SCALAR_SCALE) and
            math.isclose(data['settings']['fermion_threshold_gev'],full.s30.STOP),
            'unexpected threshold scales')
    return data


def matched_nonmixed_c5(model: str, uv: Any, eft: Any) -> tuple[np.ndarray,Path,str]:
    """Use existing exact project model-specific C5 kernels; NO mass splitting inserted here."""
    from scripts.RunT3GeneralScalar48Campaign import matching
    kernel, path, kind=matching(model)
    c5=validated_c5(kernel.evaluate(uv,eft))
    return c5,path,kind


def matched_bc_c5(model: str, uv: Any, endpoint_dir: Path):
    """Use saved mass-eigenfield Matchete A5, *evolved* heavy Yukawas/MF, and scalar endpoint."""
    from Numerical.running.FinalC5TrajectoryAdapter import _heavy_mass_basis_data
    from scripts.RunT3QuadraticEndpointC5 import (
        endpoint, numerical_c5, QUARTIC_MAP, DEFAULT_SOURCE,
    )
    import scripts.RunGeneralT3NumericalC5 as reference
    source_run, _raw, basis, base, _extra, rgbeta = endpoint(model,endpoint_dir,'final')
    heavy=_heavy_mass_basis_data(uv,offdiagonal_atol=1e-7,imaginary_atol=1e-9,
                 takagi_residual_rtol=1e-8,takagi_residual_atol=1e-7)
    masses=np.asarray(heavy.masses_gev,dtype=float)
    require(masses.shape==(3,) and np.all(masses>0), 'bad UV-evolved heavy fermion masses')
    y1=np.asarray(heavy.y1,dtype=complex)
    y2=np.asarray(heavy.y2,dtype=complex)
    require(y1.shape==(3,3) and y2.shape==(3,3),'bad evolved heavy Yukawa shapes')
    encode=lambda m: [[{'real':float(z.real),'imag':float(z.imag)} for z in row] for row in m]
    config={'y1':encode(y1),'y2':encode(y2),'fermion_masses_gev':masses.tolist(),
            'vev_gev':246.22}
    stem=full.model_tag(model)+'_final'
    proof_path=DEFAULT_SOURCE/(stem+'_one_loop_probe.json')
    proof=load(proof_path)
    require(proof.get('status')=='OneLoopDirectSymbolicWeinbergProbePassedNumericalC5Pending' and
            proof.get('model')==model and proof.get('one_loop_match_stage_completed') is True and
            proof.get('production_sequential_matching_performed') is False and
            proof.get('physical_matching_performed') is False,
            'saved direct mass-eigenstate matching source unavailable')
    exprpath=DEFAULT_SOURCE/(stem+'_one_loop_probe_weinberg_ordered_coefficient.inputform.txt')
    expression=exprpath.read_text(encoding='utf-8')
    require(set(__import__('re').findall(r'Coupling\[([A-Za-z0-9]+)',expression)) ==
            reference.ALLOWED_COUPLINGS, 'Matchete direct C5 source symbols changed')
    quartics={'lambdaT3':complex_pair(base['lambdaT3'],'lambdaT3')}
    for mat in QUARTIC_MAP: quartics[mat]=0j
    c0, _m0, _sv0=numerical_c5(expression,config,basis,quartics)
    with_extra=dict(quartics)
    for mat,rgb in QUARTIC_MAP.items():
        with_extra[mat]=-rgbeta[rgb] # only two source-level sign-mapped B/C contractions
    c5, _mnu, _singular=numerical_c5(expression,config,basis,with_extra)
    c5=validated_c5(c5)
    c0=validated_c5(c0)
    return c5,exprpath,'mass-eigenstate-one-scale-direct-A5',{
        'source_probe_sha256':fingerprint(proof_path),
        'scalar_endpoint_source_sha256':fingerprint(source_run),
        'heavy_yukawa_source':'UV-evolved y1/y2 rotated to heavy fermion mass basis',
        'heavy_mass_basis_residual':float(max(heavy.takagi_residual,heavy.basis_residual)),
        'mixed_scalar_masses_gev':[float(np.sqrt(x)) for x in basis.eigenvalues_sq],
        'quadratic_mu2_gev2':{'real':float(basis.mu_sq.real),'imag':float(basis.mu_sq.imag)},
        'signed_direct_quartics_rgbeta':{k:{'real':float(v.real),'imag':float(v.imag)}
                for k,v in rgbeta.items()},
        'direct_quartic_source_convention':'lambda_Matchete = - lambda_RGBeta for two B/C(-1) direct structures',
        'relative_extra_quartic_amplitude':(float(np.linalg.norm(c5-c0)/np.linalg.norm(c0))
                    if np.linalg.norm(c0)>0 else None),
        'matching_scope':'simultaneous all-heavy direct one-scale only',
    }


def low_energy(model: str, c5: np.ndarray, eft: Any, *, vev_gev: float=246.22) -> dict:
    from Numerical.running.ScalarThresholdBoundary import (
        build_final_sm_boundary, build_sm_weinberg_initial_conditions)
    from Numerical.running.SMWeinbergEvolution import evolve_sm_weinberg
    from Numerical.fitting.ParameterScan import (
        neutrino_mass_from_c5, charged_lepton_mass_basis_matrix)
    from Numerical.fitting.OscillationFit import prediction_from_mass_matrix
    boundary=build_final_sm_boundary(eft,matching_scale_gev=SCALAR_SCALE)
    init=build_sm_weinberg_initial_conditions(boundary,validated_c5(c5))
    low=evolve_sm_weinberg(init,SCALAR_SCALE,LOW_SCALE,rtol=1e-8,atol=1e-11)
    mass=neutrino_mass_from_c5(low.K,vev_gev=vev_gev)
    flavour=charged_lepton_mass_basis_matrix(mass,low.ye)
    prediction=prediction_from_mass_matrix(flavour,ordering='NO')
    eigen=np.linalg.svd(flavour,compute_uv=False)*1e9
    require(np.isfinite(eigen).all() and np.all(eigen>=0), 'invalid neutrino singular masses')
    return {'c5_at_1tev_gev_inverse':matrix_payload(c5),
            'c5_at_100gev_gev_inverse':matrix_payload(low.K),
            'mnu_at_100gev_ev':matrix_payload(flavour*1e9),
            'takagi_masses_ev_descending':eigen.tolist(),
            'oscillation_observables_normal_ordering':prediction.as_dict(),
            'takagi_residual':float(prediction.takagi_residual),
            'weinberg_rge_solver_nfev':int(low.nfev),
            'sm_weinberg_rge_performed':True}


def compute_one(model: str, *, source_root: Path, destination_root: Path,
                resume: bool) -> dict:
    require(model in MODELS, 'invalid model')
    folder=source_root/full.model_tag(model)
    source=folder/'continuous_trajectory.json'
    resultfile=destination_root/(full.model_tag(model)+'_neutrino.json')
    upstream=verify_continuous(model,source)
    source_hash=fingerprint(source)
    runner_hash=fingerprint(Path(__file__).resolve())
    sidecar=folder/upstream['saved_candidate_eft1_endpoint']
    sidecar_hash=fingerprint(sidecar)
    if resultfile.is_file():
        existing=load(resultfile)
        require(resume and existing.get('format')==FORMAT and
                existing.get('source_continuous_sha256')==source_hash and
                existing.get('source_endpoint_sha256')==sidecar_hash and
                existing.get('analysis_runner_sha256')==runner_hash and
                fingerprint(ROOT/existing['matching_source'])==existing.get('matching_source_sha256'),
                'output exists or sources changed; use a new --output-id')
        return existing
    uv=decode_uv(upstream['uv_fermion_threshold_state'])
    eft=decode_eft(upstream['eft1_scalar_endpoint_state'])
    vev=246.22
    extra={}
    if model in MIXED:
        c5,path,kind,extra=matched_bc_c5(model,uv,folder)
        status='CandidateMassEigenBCEvolvedYukawasSMWeinberg'
        mixing_in_c5=True
    else:
        c5,path,kind=matched_nonmixed_c5(model,uv,eft)
        mixing_in_c5=False
        status=('RestrictedKernelControlSelfBilinearNotMatched'
                if model in SELF else 'CandidateNoAdditionalQuadraticMatchedC5')
    physical=low_energy(model,c5,eft,vev_gev=vev)
    out={'format':FORMAT,'status':status,'model':model,
        'quadratic_kind':upstream['quadratic_kind'],
        'nonzero_bilinear_in_matching':mixing_in_c5,
        'scalar_spectrum_includes_bilinear':True,
        'source_continuous_report':source.relative_to(ROOT).as_posix(),
        'source_continuous_sha256':source_hash,
        'source_endpoint_sha256':sidecar_hash,
        'analysis_runner_sha256':runner_hash,
        'matching_source':path.relative_to(ROOT).as_posix(),
        'matching_source_sha256':fingerprint(path),
        'matching_kernel_kind':kind,
        'extra_mass_basis_matching_details':extra,
        'perturbative_physics_status':(
            'Nonzero holomorphic self mass is omitted from this C5 formula: '
            'this output is ONLY an ordinary reduced-C5 diagnostic on the running background'
            if model in SELF else
            'Conditional one-scale mass-eigenstate C5: full sequential EFT Wilson matching is absent'
            if model in MIXED else
            'Existing original topology C5 on fully coupled quartic running background; '
            'extra quartics have no direct one-vertex C5 in this diagonal basis'),
        'full_general_scalar_physical_c5_certified':False,
        'full_sequential_threshold_matching_performed':False,
        'frozen_reference_yukawas_used':False,
        'uv_evolved_yukawas_used':True,
        'additional_quadratic_c5_omitted':bool(model in SELF),
        **physical}
    atomic_json(resultfile,out)
    return out


def main() -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    select=parser.add_mutually_exclusive_group(required=True)
    select.add_argument('--all',action='store_true')
    select.add_argument('--model',choices=MODELS)
    parser.add_argument('--run-id',required=True,help='Existing continuous running ID')
    parser.add_argument('--source-root',type=Path,default=SOURCE_ROOT)
    parser.add_argument('--output-id',default='neutrino_v1')
    parser.add_argument('--resume',action='store_true')
    args=parser.parse_args()
    require(bool(__import__('re').fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,48}',args.output_id)),
            'invalid output-id')
    source_root=args.source_root/args.run_id
    target=source_root/args.output_id
    require(source_root.is_dir(),'continuous running directory missing: '+str(source_root))
    selected=MODELS if args.all else (args.model,)
    rows=[]
    for model in selected:
        try:
            out=compute_one(model,source_root=source_root,
                            destination_root=target,resume=args.resume)
            masses=out['takagi_masses_ev_descending']
            row={'model':model,'status':out['status'],
                 'largest_neutrino_mass_ev':masses[0],
                 'middle_neutrino_mass_ev':masses[1],
                 'smallest_neutrino_mass_ev':masses[2],
                 'quadratic_kind':out['quadratic_kind'],
                 'matching_includes_nonzero_bilinear':out['nonzero_bilinear_in_matching'],
                 'self_bilinear_c5_omitted':out['additional_quadratic_c5_omitted'],
                 **out['oscillation_observables_normal_ordering']}
            print(f'{model}: {out["status"]}; masses={masses}',flush=True)
        except (ValueError,KeyError,TypeError,OSError,RuntimeError,OverflowError,
                FloatingPointError,np.linalg.LinAlgError) as exc:
            row={'model':model,'status':'Blocked','error':str(exc)[:550]}
            print(f'{model}: BLOCKED {exc}',flush=True)
        rows.append(row)
        summary={'format':FORMAT,'source_run':args.run_id,
                 'neutrino_output_id':args.output_id,
                 'models_requested':list(selected),'models_processed':len(rows),
                 'models_successful':sum(r['status']!='Blocked' for r in rows),
                 'models_blocked':sum(r['status']=='Blocked' for r in rows),
                 'models_with_mixed_bilinear_c5':sum(r.get('matching_includes_nonzero_bilinear',False) for r in rows),
                 'models_with_unmatched_self_bilinear':sum(r.get('self_bilinear_c5_omitted',False) for r in rows),
                 'physical_full_general_scalar_matching_complete':False,
                 'rows':rows}
        atomic_json(target/'summary.json',summary)
        target.mkdir(parents=True,exist_ok=True)
        columns=['model','status','quadratic_kind','matching_includes_nonzero_bilinear',
                 'self_bilinear_c5_omitted','largest_neutrino_mass_ev',
                 'middle_neutrino_mass_ev','smallest_neutrino_mass_ev',
                 'delta_m21_sq_ev2','delta_m3l_sq_ev2','sin2_theta12','sin2_theta13','sin2_theta23',
                 'error']
        with (target/'summary.csv.tmp').open('w',newline='',encoding='utf-8') as f:
            writer=csv.DictWriter(f,fieldnames=columns,extrasaction='ignore')
            writer.writeheader()
            writer.writerows(rows)
        (target/'summary.csv.tmp').replace(target/'summary.csv')
    print('Summary:',target/'summary.json')
    return 2 if any(r['status']=='Blocked' for r in rows) else 0

if __name__=='__main__':
    raise SystemExit(main())
