"""Batch 2 post-export audit of the A(-4) renormalisable EFT1 tensors.

READ-ONLY inputs: the seven-scale NPZ, its basis/couplings JSON, the native
EFT1 RGBeta beta export, and (optionally) the checked-out executable sources.
This does not transform a Matchete Wilson boundary, match scalars, or solve C.
All source-invariant comparisons use *explicit* reference tensor columns.
An algebraic discrepancy is recorded as an obstruction, never 'fixed' by
renaming couplings or changing their normalisation to match the betas.
"""
from __future__ import annotations

from hashlib import sha256
from pathlib import Path
from typing import Any, Mapping
import json
import math

import numpy as np

from Numerical.running import Aminus4EFT1Couplings as eft

F_OUTPUT = 'A_m4_eft1_tensors.npz'
B_OUTPUT = 'A_m4_eft1_tensor_basis.json'
SCALAR_NAMES = (
    'lambdaH','lambdaS1','lambdaS2','lambdaS2Adj',
    'lambdaH1','lambdaH2','lambdaH2Adj','lambda12',
    'lambdaT3.real','lambdaT3.imag',
)


class AuditError(ValueError):
    """A malformed export or a failed numerical source reconstruction."""


def _require(ok: bool, description: str) -> None:
    if not ok:
        raise AuditError(description)


def _finite_array(array: Any, shape: tuple[int,...], label: str) -> np.ndarray:
    a=np.asarray(array)
    _require(a.shape == shape, f'{label}: expected shape {shape}, got {a.shape}')
    _require(np.issubdtype(a.dtype,np.number) and np.isfinite(a).all(),
             f'{label}: nonnumeric or nonfinite data')
    return a


def _residual(a: np.ndarray, b: np.ndarray) -> dict:
    a=np.asarray(a); b=np.asarray(b)
    diff=np.asarray(a-b)
    return {'max_absolute':float(np.max(np.abs(diff))),
            'relative_frobenius':float(np.linalg.norm(diff.ravel())/
                max(1.e-30,float(np.linalg.norm(np.asarray(b).ravel()))))}


def _scalar_coordinate_columns() -> np.ndarray:
    columns=eft.quartic_columns()
    basis=[np.asarray(columns[k],float).ravel() for k in SCALAR_NAMES[:8]]
    basis.extend([2*columns['lambdaT3'].real.ravel(),
                  -2*columns['lambdaT3'].imag.ravel()])
    return np.stack(basis,axis=1)


def scalar_beta_tensor(lam:np.ndarray) -> np.ndarray:
    """Scalar-only 16*pi^2 beta(lambda_abcd), with real-scalar 1/4! potential.

      beta_abcd = lambda_abef lambda_efcd + lambda_acef lambda_efbd
                + lambda_adef lambda_efbc

    This checks the universal quartic-loop contribution. Other gauge/Yukawa
    contributions remain in the full RGBeta ODE and are not discarded there.
    """
    return (np.einsum('abef,efcd->abcd',lam,lam,optimize=True)
            + np.einsum('acef,efbd->abcd',lam,lam,optimize=True)
            + np.einsum('adef,efbc->abcd',lam,lam,optimize=True))


def scalar_beta_comparison(point:Mapping, beta_payload:Mapping,
                           *, evaluator=None) -> dict:
    """Independent tensor-loop versus actual RGBeta scalar beta expressions.

    Set only active gauge and SM Yukawa inputs to zero for this diagnostic;
    this is NOT the trajectory used by the RGE. All quartics, including the
    complex topology coupling, retain their supplied values.
    """
    eft.validate_beta_payload(beta_payload)
    p=eft._point(point)
    q={**p, 'gY':0., 'g2':0., 'g3':0.,
       'yu':np.zeros((3,3),complex), 'yd':np.zeros((3,3),complex),
       'ye':np.zeros((3,3),complex)}
    if evaluator is None:
        from Numerical.core.RGBetaEvaluator import evaluate_inputform_expression as evaluator
    names=SCALAR_NAMES[:8]+('lambdaT3',)
    observed={name:evaluator(beta_payload['report_betas'][name],q)
              for name in names}
    actual=[]
    for name in SCALAR_NAMES:
        if name == 'lambdaT3.real':z=complex(observed['lambdaT3']).real
        elif name == 'lambdaT3.imag':z=complex(observed['lambdaT3']).imag
        else:
            z=complex(observed[name])
            _require(abs(z.imag)<1e-9*max(1.,abs(z.real)),
                     f'{name}: scalar-only native beta is unexpectedly complex')
            z=z.real
        actual.append(float(z))
    A=_scalar_coordinate_columns()
    _require(np.linalg.matrix_rank(A,tol=1e-10)==10,
             'Reference real quartic basis is not rank 10')
    lam=eft.quartic_tensor(q)
    full_beta=scalar_beta_tensor(lam).ravel()
    predicted=np.linalg.lstsq(A,full_beta,rcond=1e-13)[0]
    closure=float(np.linalg.norm(A@predicted-full_beta)/
                  max(1e-30,float(np.linalg.norm(full_beta))))
    rows={}
    for name,a,b in zip(SCALAR_NAMES,predicted,actual):
        abs_err=abs(float(a)-float(b))
        rows[name]={'reference_scalar_loop':float(a),
                    'native_rgbeta_scalar_loop':float(b),
                    'absolute_difference':abs_err,
                    'relative_difference':abs_err/max(1e-14,abs(float(a)),abs(float(b)))}
    mismatches=[name for name,rec in rows.items()
                if rec['absolute_difference']>2e-9*max(1,abs(rec['reference_scalar_loop']),
                                                       abs(rec['native_rgbeta_scalar_loop']))]
    return {'reference_basis_rank':10,
            'scalar_loop_projection_closure_relative':closure,
            'coupling_beta_comparison':rows,
            'mismatched_coupling_betas':mismatches,
            'reference_and_native_scalar_loop_match':closure<1e-9 and not mismatches,
            'native_to_reference_basis_assumed_identical':False,
            'no_coupling_rescaling_fitted':True}


def scalar_unit_channel_audits(beta_payload:Mapping, *, evaluator=None) -> list[dict]:
    """Source-level one-coupling checks, including the complex topology phase.

    Unlike inspecting the beta of one running benchmark, these probes make the
    independent invariant coefficient of each quadratic beta monomial visible.
    No beta coefficient is fitted or rescaled.
    """
    probes=[]
    for name in SCALAR_NAMES:
        p={'gY':0.,'g2':0.,'g3':0.,
           'yu':np.zeros((3,3),complex),'yd':np.zeros((3,3),complex),
           'ye':np.zeros((3,3),complex),'mS1Sq':1e6,'mS2Sq':1e6}
        p.update({n:(0j if n=='lambdaT3' else 0.) for n in eft.SCALAR_COUPLINGS})
        if name.startswith('lambdaT3.'):
            p['lambdaT3']=1j if name.endswith('.imag') else (1+0j)
        else:p[name]=1.
        comparison=scalar_beta_comparison(p,beta_payload,evaluator=evaluator)
        probes.append({'unit_coupling_direction':name,
                       'mismatched_beta_components':comparison['mismatched_coupling_betas'],
                       'closure_relative':comparison['scalar_loop_projection_closure_relative'],
                       'all_agree':comparison['reference_and_native_scalar_loop_match'],
                       'components':{k:v for k,v in comparison['coupling_beta_comparison'].items()
                                     if (abs(v['reference_scalar_loop'])>1e-12 or
                                         abs(v['native_rgbeta_scalar_loop'])>1e-12)}})
    return probes


def _gauge_ward_residual(Y: np.ndarray, lam:np.ndarray, generators:Mapping) -> dict:
    checks={}
    for group in ('U1','SU2','SU3'):
        yf=0.;sv=0.
        for f,s in zip(generators[group]['fermion'],generators[group]['scalar']):
            wy=(np.einsum('ik,kja->ija',f.T,Y,optimize=True)
                +np.einsum('ika,kj->ija',Y,f,optimize=True)
                +np.einsum('ijb,ba->ija',Y,s,optimize=True))
            wl=(np.einsum('ea,ebcd->abcd',s,lam,optimize=True)
                +np.einsum('eb,aecd->abcd',s,lam,optimize=True)
                +np.einsum('ec,abed->abcd',s,lam,optimize=True)
                +np.einsum('ed,abce->abcd',s,lam,optimize=True))
            yf=max(yf,float(np.max(np.abs(wy))))
            sv=max(sv,float(np.max(np.abs(wl))))
        checks[group]={'yukawa_max_absolute':yf,'quartic_max_absolute':sv}
    return checks


def _state_from_export(mu:float,couplings:dict):
    """Reuse the already checked Batch-2 native-state constructor."""
    from scripts.ExportAminus4EFT1Couplings import make_state
    sm_keys=('gY','g2','g3','lambdaH','yu','yd','ye')
    obj={'mu_gev':float(mu),'sm':{k:couplings[k] for k in sm_keys}}
    obj.update({k:v for k,v in couplings.items() if k not in sm_keys})
    return make_state(obj)


def _verify_ode_intervals(scales:np.ndarray,states:list[dict],payload:Mapping,
                          *,rtol:float=2.e-9,atol:float=1.e-11) -> list[dict]:
    from Numerical.running.IntermediateScalarRunner import run_intermediate_scalar_segment
    from Numerical.running.IntermediateScalarStateVector import pack_intermediate_scalar_state
    records=[]
    # One-step re-integration of every saved interval, from its own stored
    # initial state, is independent of interpolation and of source summary flags.
    for k in range(len(scales)-1):
        state=_state_from_export(scales[k],states[k]);target=_state_from_export(scales[k+1],states[k+1])
        result=run_intermediate_scalar_segment(state,payload,float(scales[k+1]),
                        save_scales_gev=[float(scales[k]),float(scales[k+1])],
                        rtol=rtol,atol=atol)
        numeric,layout=pack_intermediate_scalar_state(result.final_state)
        saved,saved_layout=pack_intermediate_scalar_state(target)
        _require(layout.blocks == saved_layout.blocks, 'ODE and exported state layouts differ')
        raw=np.abs(numeric-saved)
        denom=atol+2e-6*np.maximum(np.abs(numeric),np.abs(saved))
        max_scaled=float(np.max(raw/denom))
        records.append({'from_gev':float(scales[k]),
                        'to_gev':float(scales[k+1]),
                        'maximum_absolute_state_difference':float(np.max(raw)),
                        'maximum_scaled_state_error':max_scaled,
                        'passed':bool(max_scaled<1.)})
    return records


def audit_export(folder:Path,beta_file:Path,*,root:Path|None=None,
                 check_sources:bool=True,verify_ode:bool=True) -> dict:
    """Read-only verification. Returns diagnostic results, including blockers.

    Malformed/tampered export => AuditError. A genuine source-basis
    normalization mismatch => detailed diagnostic with certified=False.
    """
    folder=Path(folder);beta_file=Path(beta_file)
    basis_path=folder/B_OUTPUT;npz_path=folder/F_OUTPUT
    _require(npz_path.is_file() and basis_path.is_file(),
             'Missing Batch-2 NPZ or tensor-basis JSON')
    _require(beta_file.is_file(),'Missing native A(-4) EFT1 beta file')
    basis=json.loads(basis_path.read_text(encoding='utf-8-sig'))
    payload=json.loads(beta_file.read_text(encoding='utf-8-sig'))
    eft.validate_beta_payload(payload)
    _require(basis.get('model')==eft.MODEL,'Wrong model in saved basis')
    _require(basis.get('fermion_basis')==list(eft.F_LABELS),
             'Saved Weyl ordering differs from adapter')
    _require(basis.get('real_scalar_basis')==list(eft.S_LABELS),
             'Saved scalar ordering differs from adapter')
    _require(basis.get('matchete_scalar_and_Wilson_basis_map_certified') is False,
             'Unexpected physical-basis certification claim in input')
    _require(basis.get('rge_payload_sha256')==sha256(beta_file.read_bytes()).hexdigest(),
             'Native beta file SHA-256 differs from the one used for the export')
    states=basis.get('running_couplings')
    _require(isinstance(states,list) and len(states)>=3,'Missing multi-scale couplings')
    with np.load(npz_path,allow_pickle=False) as archive:
        arrays={k:np.array(archive[k]) for k in archive.files}
    scales=np.asarray(arrays.get('mu_gev',[]),float)
    _require(scales.ndim==1 and len(scales)==len(states) and
             np.isfinite(scales).all() and np.all(scales>0) and
             np.all(np.diff(scales)<0), 'Invalid EFT1 sampled scale grid')
    _require(np.allclose(scales,basis.get('sampled_scales_gev',[]),rtol=1e-13,atol=0),
             'NPZ and JSON scale grids differ')
    n=len(scales)
    Y=_finite_array(arrays.get('Y',[]),(n,45,45,12),'Y')
    lam=_finite_array(arrays.get('lambda4',[]),(n,12,12,12,12),'lambda4')
    direct_g=eft.gauge_generators()
    for kind,sector in direct_g.items():
        for rep in ('fermion','scalar'):
            key=f'generator_{kind}_{rep}'
            observed=_finite_array(arrays.get(key,[]),
                                   np.stack(sector[rep]).shape,key)
            _require(np.allclose(observed,np.stack(sector[rep]),rtol=0,atol=1e-13),
                     f'{key} differs from generated SU(3)xSU(2)xU(1) representation')
    if check_sources:
        _require(root is not None,'Pass the repository root when check_sources=True')
        provenance=eft.source_contract(Path(root),payload)
        previous=basis.get('provenance',{}).get('source_sha256',{})
        _require(previous==provenance['source_sha256'],
                 'Executable-source hashes changed after tensor export')
    else:provenance={'source_check':'not requested (fixture/offline mode)'}
    sample_records=[]
    for k,(mu,s) in enumerate(zip(scales,states)):
        _require(set(s)==eft.BETA_NAMES,'Saved scale has missing/extra EFT1 couplings')
        pt=eft.tensor_point(float(mu),s)
        dy=_residual(pt.Y,Y[k]);dl=_residual(pt.lambda4,lam[k])
        sy=float(np.max(np.abs(Y[k]-Y[k].transpose(1,0,2))))
        sl=float(max(np.max(np.abs(lam[k]-lam[k].transpose(1,0,2,3))),
                     np.max(np.abs(lam[k]-lam[k].transpose(0,2,1,3))),
                     np.max(np.abs(lam[k]-lam[k].transpose(0,1,3,2)))))
        ward=_gauge_ward_residual(Y[k],lam[k],direct_g)
        passed=dy['max_absolute']<5e-12 and dl['max_absolute']<5e-11 and \
            sy<5e-12 and sl<5e-11 and all(
            d['yukawa_max_absolute']<5e-10 and d['quartic_max_absolute']<5e-9
            for d in ward.values())
        sample_records.append({'mu_gev':float(mu),'tensor_reconstruction':{
             'Y':dy,'lambda4':dl},
             'fermion_exchange_symmetry_max_absolute':sy,
             'quartic_permutation_symmetry_max_absolute':sl,
             'gauge_ward':ward, 'passed':bool(passed)})
    _require(all(s['passed'] for s in sample_records),
             'Saved tensors fail direct-coupling reconstruction or Ward tests')
    beta_records=[]
    for k in sorted(set((0,n//2,n-1))):
        beta=scalar_beta_comparison(states[k],payload)
        beta_records.append({'mu_gev':float(scales[k]),**beta})
    unit_audits=scalar_unit_channel_audits(payload)
    ode=_verify_ode_intervals(scales,states,payload) if verify_ode else []
    if verify_ode:
        _require(all(v['passed'] for v in ode),
                 'Independent interval re-integration does not reproduce exported RG trajectory')
    return {'format':'t3-a-m4-batch2-post-export-audit-v1','model':eft.MODEL,
            'source_revision_expected':eft.REVISION_CHECKED,
            'input_tensors_sha256':sha256(npz_path.read_bytes()).hexdigest(),
            'input_basis_sha256':sha256(basis_path.read_bytes()).hexdigest(),
            'beta_sha256':sha256(beta_file.read_bytes()).hexdigest(),
            'sample_count':n,'scale_range_gev':[float(scales[-1]),float(scales[0])],
            'field_dimensions':{'weyl':45,'real_scalar':12},
            'source_provenance':provenance,
            'tensor_reconstruction_and_gauge_ward_passed':True,
            'samples':sample_records,
            'scalar_beta_diagnostics':beta_records,
            'source_unit_invariant_beta_audits':unit_audits,
            'all_native_reference_scalar_beta_coefficients_agree':all(
                  z['reference_and_native_scalar_loop_match'] for z in beta_records) and
                all(z['all_agree'] for z in unit_audits),
            'independent_ode_interval_reintegration_checked':bool(verify_ode),
            'ode_interval_reintegration':ode,
            'physical_matchete_to_RGBeta_map_certified':False,
            'ready_for_physically_certified_master_Wilson_RGE':False,
            'reason':'A complete native Matchete-to-RGBeta invariant/CG map has not been algebraically certified',
            'note':'No Wilson ODE, scalar matching, GammaCC conversion or post-threshold SMEFT running performed'}
