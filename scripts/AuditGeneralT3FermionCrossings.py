"""Stage 2: read-only full-UV heavy mass-crossing replay for all 16 ordinary T3.

Replays only the SAVED original continuous_v1 UV RGBeta trajectory using its
same fitted UV input, candidate quartics, complex mu2 and beta export. Verifies
original SHA256 input inventory, saved UV scalar samples, full MF/Yukawa endpoint
and grouped-F singular masses BEFORE locating UV mass crossings on a dense ODE
solution. Reuses Stage 1 scalar crossing diagnostics; no EFT RGE is changed.

Example:
  python -m scripts.AuditGeneralT3FermionCrossings --all
  python -m scripts.AuditGeneralT3FermionCrossings --model B_-1

A found crossing of F2 on the unde-coupled UV background is NOT a correct EFT
matching boundary for F1: proper sequential removal and intermediate RGEs are
still Workstreams B/C dependencies. Nothing in the saved pipeline is edited.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from Numerical.experimental.GeneralT3RunningFermionThresholds import (
    locate_fermion_crossings, mass_matrix_from_vector, singular_masses,
)
from Numerical.experimental.GeneralT3ScalarRGState import (
    load_registry, pack_general_scalar_state, unpack_general_scalar_state,
)
from Numerical.experimental.GeneralT3ScalarQuadraticMasses import (
    load_manifest, quadratic_kind,
)
from Numerical.running.FermionThresholdBoundary import fermion_threshold_masses
from scripts import RunT3FullScalarContinuous as continuous
from scripts import AuditGeneralT3MassCrossings as stage1

STATUS='CandidateRunningFermionCrossingAuditNoMatching'
FORMAT='general-t3-running-heavy-mass-crossing-audit-v1'
ROOT=continuous.ROOT
SOURCE_ROOT=stage1.SOURCE_ROOT
DEFAULT_OUT=SOURCE_ROOT/'workstream_d_fermion_audit_v1/audit_summary.json'


def require(ok: bool, message: str) -> None:
    if not ok:
        raise ValueError('Workstream D fermion audit: '+message)


def verify_saved_inputs(record: dict, model: str, paths: dict) -> None:
    """Verify current source file hashes against the original saved run."""
    hashes=record.get('input_sha256')
    require(isinstance(hashes,dict) and set(hashes)==set(paths),
            'saved continuous UV input inventory incomplete or altered')
    for name,path in paths.items():
        require(Path(path).is_file(), 'missing original source '+str(path))
        require(sha256(Path(path).read_bytes()).hexdigest()==hashes[name],
                'saved UV source SHA256 mismatch: '+name+' ('+str(path)+')')


def _complex_matrix(raw: dict, name: str) -> np.ndarray:
    require(isinstance(raw,dict) and set(raw)=={'real','imag'},
            'missing complex saved matrix '+name)
    out=np.asarray(raw['real'],dtype=float)+1j*np.asarray(raw['imag'],dtype=float)
    require(out.shape==(3,3) and np.isfinite(out).all(),
            'invalid saved matrix '+name)
    return out


def replay_and_find(record: dict, model: str, *, scan_points: int = 65) -> dict:
    """Reintegrate the original UV ODE only; never use candidate matching."""
    require(record.get('model')==model and
            record.get('format')==continuous.FORMAT and
            record.get('status')==continuous.STATUS,
            'saved continuous run metadata differs')
    settings=record.get('settings',{})
    require(settings.get('model')==model and
            settings.get('uv_start_gev')==continuous.s30.START and
            settings.get('fermion_threshold_gev')==continuous.s30.STOP,
            'saved UV scale/model setting differs')
    kind=quadratic_kind(load_manifest(continuous.MANIFEST),model)
    require(record.get('quadratic_kind')==kind,'quadratic kind changed')
    registry=load_registry(continuous.REGISTRY)
    full,conf,stage29=continuous.initial_uv(model,registry,continuous.STAGE22)
    beta_file=continuous.BETA_DIR/f'{continuous.model_tag(model)}_uv_quadratic_live_rgbeta.json'
    eft_file=continuous.BETA_DIR/f'{continuous.model_tag(model)}_eft1_quadratic_live_rgbeta.json'
    sources={'registry':continuous.REGISTRY,'manifest':continuous.MANIFEST,
             'stage29_boundary':stage29,'fitted_config':conf,
             'uv_mass_beta':beta_file,'eft1_mass_beta':eft_file,
             'quadratic_rgbeta_source':ROOT/'RGE/running/rgbeta/RunGeneralT3QuadraticRGBeta.wl',
             'continuous_runner':Path(continuous.__file__).resolve()}
    verify_saved_inputs(record,model,sources)
    original_uv=record.get('uv_trajectory')
    saved_points=settings.get('saved_points_per_segment')
    require(isinstance(saved_points,int) and saved_points>=4 and
            isinstance(original_uv,list) and len(original_uv)==saved_points,
            'missing/inconsistent saved UV sampled trajectory')
    y0,layout=pack_general_scalar_state(full,registry)
    mu2=stage1.complex_pair(record['quadratic_mu2_uv_gev2'])
    mu2_fraction=float(settings['mu2_fraction'])
    phase=float(settings['mu2_phase_rad'])
    require(math.isfinite(mu2_fraction) and math.isfinite(phase)
            and 0<=mu2_fraction<=0.45,'bad original quadratic UV settings')
    expected=(mu2_fraction*min(full.baseline.mS1Sq,full.baseline.mS2Sq)
              *complex(math.cos(phase),math.sin(phase)) if kind!='none' else 0j)
    require(abs(expected-mu2)<2e-12*max(1.,abs(expected)),
            'original saved UV mu2 disagrees with benchmark settings')
    if kind!='none':
        y0=np.concatenate((y0,np.array([mu2.real,mu2.imag])))
    betas=continuous.beta_inventory(continuous.read(beta_file),model=model,
        stage='UV',representation=full.baseline.representation,
        layout=layout,registry=registry,kind=kind)
    mu_hi=float(continuous.s30.START); mu_lo=float(continuous.s30.STOP)
    t_hi=math.log(mu_hi);t_lo=math.log(mu_lo)
    # EXACT integration scheme/options as RunT3FullScalarContinuous.integrate_segment.
    rhs=lambda t,y:continuous.coupled_rhs(t,y,layout=layout,registry=registry,
                                           stage='UV',kind=kind,betas=betas)
    rhs(t_hi,y0)
    sol=solve_ivp(rhs,(t_hi,t_lo),y0,method='DOP853',
                  t_eval=np.linspace(t_hi,t_lo,saved_points),dense_output=True,
                  rtol=1e-8,atol=1e-11,max_step=0.35)
    require(sol.success and sol.sol is not None and
            sol.y.shape==(len(y0),saved_points) and np.isfinite(sol.y).all(),
            'UV replay ODE failed: '+str(sol.message))
    max_sample_err=0.0
    for j,row in enumerate(original_uv):
        mu=math.exp(float(sol.t[j]));vec=sol.y[:,j]
        require(math.isclose(mu,float(row['mu_gev']),rel_tol=2e-11,abs_tol=1e-7),
                'UV saved RG sampling scales mismatch')
        st=unpack_general_scalar_state(vec[:layout.size],layout,registry,uv=True,mu_gev=mu)
        for key,actual in (('mS1Sq_gev2',st.baseline.mS1Sq),
                           ('mS2Sq_gev2',st.baseline.mS2Sq)):
            saved=float(row[key]);rel=abs(float(actual)-saved)/max(1.,abs(saved))
            max_sample_err=max(max_sample_err,rel)
            require(rel<=2e-7,'saved scalar RG replay sample mismatch: '+key)
        values=row['extra_quartics']
        require(set(values)==set(st.extra),'saved quartic state inventory changed')
        for name,z in st.extra.items():
            saved=stage1.complex_pair(values[name]);rel=abs(z-saved)/max(1e-10,abs(saved))
            max_sample_err=max(max_sample_err,rel)
            require(rel<=2e-7,'saved UV quartic replay sample mismatch: '+name)
        saved_mu2=stage1.complex_pair(row['mu2_gev2'])
        current_mu2=complex(sol.y[layout.size,j],sol.y[layout.size+1,j]) if kind!='none' else 0j
        rel=abs(current_mu2-saved_mu2)/max(1.,abs(saved_mu2))
        max_sample_err=max(max_sample_err,rel)
        require(rel<=2e-7,'saved UV bilinear replay sample mismatch')
    endpoint=unpack_general_scalar_state(sol.y[:layout.size,-1],layout,registry,
                                         uv=True,mu_gev=mu_lo).baseline
    # Majorana symmetric MF is also validated by unpacking each saved-state sample.
    saved_endpoint=record.get('uv_fermion_threshold_state',{})
    for name in ('MF','y1','y2'):
        current=np.asarray(getattr(endpoint,name),dtype=complex)
        target=_complex_matrix(saved_endpoint.get(name),name)
        rel=float(np.linalg.norm(current-target)/max(np.linalg.norm(target),1.))
        require(rel<2e-7,'UV fermion/Yukawa endpoint parity failed: '+name)
        max_sample_err=max(max_sample_err,rel)
    endpoint_mass=fermion_threshold_masses(endpoint).masses_gev
    saved_mass=np.asarray(record['heavy_fermion_threshold_masses_gev'],dtype=float)
    require(saved_mass.shape==(3,) and np.allclose(endpoint_mass,saved_mass,
                    rtol=2e-7,atol=1e-7),'saved grouped fermion singular masses changed')
    majorana=bool(endpoint.representation.majorana_fermion)
    def matrix_at(t:float)->np.ndarray:
        # The full UV solver appends Re/Im(mu2) beyond StateVectorLayout for
        # self-pair and mixed scalar models. Explicitly permit exactly those
        # two coordinates; MF is still read only from its original named block.
        return mass_matrix_from_vector(
            np.asarray(sol.sol(float(t)), dtype=float), layout,
            trailing_components=(2 if kind != 'none' else 0),
        )
    # Endpoint consistency with direct packed-MF extraction and same singular definition.
    require(np.allclose(singular_masses(matrix_at(t_lo),majorana=majorana),
                        endpoint_mass,rtol=2e-10,atol=1e-8),
            'packed MF endpoint mass extraction disagrees with legacy state')
    result=locate_fermion_crossings(matrix_at,mu_hi,mu_lo,
                                    majorana=majorana,scan_points=scan_points)
    result['saved_uv_sample_max_relative_replay_error']=float(max_sample_err)
    result['saved_grouped_F_mass_relative_parity_error']=float(np.max(
        abs(endpoint_mass-saved_mass)/np.maximum(saved_mass,1.)))
    result['uv_replay_nfev']=int(sol.nfev)
    result['uv_beta_export_sha256']=record['input_sha256']['uv_mass_beta']
    result['original_continuous_source_hashes_verified']=True
    return result


def audit_model(record: dict, model: str, manifest: dict, *, scan_points: int) -> dict:
    original=stage1.audit_record(record,manifest)
    ferm=replay_and_find(record,model,scan_points=scan_points)
    scalar=original['eft1_scalar_crossings']
    # Single event *diagnostic* order; no matching, mass removal or piecewise RGE.
    combined=([{'particle':e['particle'],'kind':'fermion',
                'scale_gev':e['scale_gev'],'status':e['status']}
               for e in ferm['crossings']]
              +[{'particle':e['particle'],'kind':'scalar',
                 'scale_gev':e['scale_gev'],'status':e['status']} for e in scalar])
    combined.sort(key=lambda e:(-e['scale_gev'],e['kind'],e['particle']))
    return {'model':model,'status':STATUS,
            'quadratic_kind':original['quadratic_kind'],
            'grouped_F_spectrum_status':original['grouped_F_spectrum_status'],
            'uv_scalar_crossings':original['uv_scalar_crossings'],
            'eft1_scalar_crossings':scalar,
            'uv_fermion_crossings':ferm['crossings'],
            'unresolved_fermion_crossings':ferm['unresolved_crossings'],
            'all_fermion_crossings_bracketed':ferm['all_three_crossings_bracketed'],
            'candidate_event_order':combined,
            'full_uv_mass_diagnostic':ferm,
            'candidate_basis_not_physically_certified':True,
            'full_uv_background_kept_below_first_crossing':True,
            'sequential_fermion_matching_performed':False,
            'finite_scalar_matching_performed':False,
            'physical_neutrino_mass_prediction':False}


def main()->int:
    ap=argparse.ArgumentParser(description=__doc__)
    which=ap.add_mutually_exclusive_group(required=True)
    which.add_argument('--all',action='store_true')
    which.add_argument('--model')
    ap.add_argument('--source-root',type=Path,default=SOURCE_ROOT)
    ap.add_argument('--output',type=Path,default=None)
    ap.add_argument('--scan-points',type=int,default=65)
    args=ap.parse_args()
    manifest=load_manifest(continuous.MANIFEST)
    selected=tuple(manifest['models']) if args.all else (args.model,)
    if any(m not in manifest['models'] for m in selected):ap.error('Unknown ordinary T3 model')
    out=args.output or (args.source_root/'workstream_d_fermion_audit_v1/audit_summary.json')
    if out.is_file():ap.error('Will not overwrite previous audit: '+str(out))
    rows=[]
    for model in selected:
        filename=args.source_root/model.replace('+','p').replace('-','m')/'continuous_trajectory.json'
        try:
            require(filename.is_file(),'saved continuous UV source missing: '+str(filename))
            record=json.loads(filename.read_text(encoding='utf-8-sig'))
            row=audit_model(record,model,manifest,scan_points=args.scan_points)
        except (ValueError,OSError,RuntimeError,KeyError,TypeError,
                ArithmeticError,np.linalg.LinAlgError) as exc:
            row={'model':model,'status':'Blocked','reason':str(exc)}
        rows.append(row)
        if row['status']==STATUS:
            print(model,'fermion crossings:',
                  [(x['particle'],round(x['scale_gev'],6)) for x in row['uv_fermion_crossings']],
                  'unresolved:',[x['particle'] for x in row['unresolved_fermion_crossings']],flush=True)
        else:print(model,'BLOCKED',row.get('reason'),flush=True)
    payload={'format':FORMAT,'source_root':str(args.source_root),
             'status':'CandidateAuditComplete' if all(x['status']==STATUS for x in rows) else 'Incomplete',
             'numerical_stage':'UndecoupledFullUVThenSavedScalarEFT1',
             'physical_matching_certified':False,'models':rows}
    out.parent.mkdir(parents=True,exist_ok=True)
    tmp=out.with_name(out.name+'.tmp')
    tmp.write_text(json.dumps(payload,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    tmp.replace(out)
    print('Saved Stage 2 threshold audit:',out,flush=True)
    return 2 if payload['status']!='CandidateAuditComplete' else 0


if __name__=='__main__':
    raise SystemExit(main())
