"""B(-1)/C(-1) same-run one-loop Matchete C5 candidate -> SM+Weinberg MZ.

Uses the output of RunGeneralT3BCNormal and the *existing* native Matchete
one-loop ordered coefficient evaluator. This is a direct all-heavy one-scale
kernel evaluated on the current running endpoints; it is NOT a certified
sequential F -> scalar EFT matching, nor complete EFT1 Wilson running.

Run: python -X utf8 -m scripts.RunGeneralT3BCDirectCandidate --all
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np

from Numerical.running.SMWeinbergEvolution import (
    SMWeinbergInitialConditions, evolve_sm_weinberg,
)
from physics.NeutrinoTrajectory import (
    neutrino_mass_from_c5, charged_lepton_mass_basis_matrix,
)
from physics.NeutrinoObservables import calculate_neutrino_observables

ROOT = Path(__file__).resolve().parents[1]
RUNNING_ROOT = ROOT / 'output/general_numerical_BC_normal'
OUTPUT = ROOT / 'output/general_numerical_BC_direct_candidate'
SOURCE = ROOT / 'scripts/DirectVsSequentialOneLoop.experimental.py'
MODELS = ('B_-1', 'C_-1')
QUARTICS = ('lambdaT3', 'lambdaT3ExtraQ020020I1', 'lambdaT3ExtraQ020200I1')
MAP_AUDIT = ROOT / 'output/audits/t3_bc_missing_quartic_map.json'
RGBETA_KEYS = {
    'lambdaT3ExtraQ020020I1': 'lambdaHHdagS2S2',
    'lambdaT3ExtraQ020200I1': 'lambdaHHdagS1barS1bar',
}


def mapped_quartics(model, eft, audit):
    require(audit.get('status') == 'ConditionalNativeTensorMapAudited' and
            audit.get('physical_mapping_certified') is False,
            'expected uncertified native mapping audit')
    require('lambdaT3' in eft, 'missing native lambdaT3 in running output')
    quartics = {'lambdaT3': complex_value(eft['lambdaT3'])}
    rows = audit.get('entries')
    require(isinstance(rows, list), 'invalid audit entry list')
    for native, rgbeta in RGBETA_KEYS.items():
        matches = [r for r in rows if isinstance(r, dict) and
                   r.get('model') == model and r.get('matchete') == native and
                   r.get('rgbeta') == rgbeta]
        require(len(matches) == 1, f'{model}/{native}: unique mapping entry missing')
        entry = matches[0]
        require(entry.get('physical_mapping_certified') is False and
                entry.get('conditional_matchete_from_rgbeta_factor') == '-1',
                f'{model}/{native}: expected conditional -1 mapping')
        require(rgbeta in eft, f'missing RGBeta running coupling {rgbeta}')
        quartics[native] = -complex_value(eft[rgbeta])
    return quartics



def require(test, message):
    if not bool(test):
        raise ValueError('B/C native C5 candidate: ' + message)


def read(path):
    require(path.is_file(), 'missing input '+str(path))
    return json.loads(path.read_text(encoding='utf-8-sig'))


def complex_value(v):
    if isinstance(v, dict) and set(v) == {'real','imag'}:
        return complex(float(v['real']), float(v['imag']))
    if isinstance(v, (int,float)):
        return complex(v)
    if isinstance(v, list) and len(v)==2 and all(isinstance(x,(int,float)) for x in v):
        return complex(*v)
    raise ValueError('unexpected complex value '+repr(v)[:120])


def matrix(raw, n=3):
    if isinstance(raw,dict) and set(raw)=={'real','imag'}:
        a=np.asarray(raw['real'],float)+1j*np.asarray(raw['imag'],float)
    else:
        a=np.asarray([[complex_value(v) for v in row] for row in raw],complex)
    require(a.shape==(n,n) and np.isfinite(a).all(),f'expected finite {n}x{n} matrix')
    return a


def encode(a):
    a=np.asarray(a,complex)
    return {'real':a.real.tolist(),'imag':a.imag.tolist()}


def encode_json(x):
    if isinstance(x,np.ndarray):
        return encode_json(x.tolist())
    if isinstance(x,(tuple,list)):
        return [encode_json(v) for v in x]
    if isinstance(x,(complex,np.complexfloating)):
        return {'real':float(x.real),'imag':float(x.imag)}
    if isinstance(x,np.generic):
        return encode_json(x.item())
    if isinstance(x,dict):
        return {str(k):encode_json(v) for k,v in x.items()}
    return x


def save(path, obj):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def load_evaluator():
    require(SOURCE.is_file(), 'missing source evaluator '+str(SOURCE))
    spec=importlib.util.spec_from_file_location('t3_bc_direct_native_v22',SOURCE)
    require(spec is not None and spec.loader is not None,'could not load native evaluator')
    obj=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    require(tuple(obj.Q_NAMES)==QUARTICS,'native evaluator quartic inventory changed')
    return obj


def model_result(model, src, root, out, mu_low, vev, audit, audit_path):
    tag=model.replace('-','m')
    running_path=root/'running'/f'{tag}.json'
    basis_path=root/f'{tag}_mass_basis_handoff.json'
    running=read(running_path)
    handoff=read(basis_path)
    require(running.get('model')==model and running.get('status')=='NumericalRenormalisableRunningOnly',
            'invalid B/C normal-running JSON')
    require(handoff.get('model')==model and handoff.get('C5_computed') is False,
            'invalid B/C scalar mass-basis handoff')
    require(running.get('uv_solver',{}).get('success') is True and
            running.get('eft1_solver',{}).get('success') is True,'RGE failed')
    uv=running['uv_at_fermion_threshold']
    eft=running['eft1_at_scalar_threshold']
    masses=np.asarray(handoff['basis']['masses_gev'],float)
    require(masses.shape==(2,) and np.isfinite(masses).all() and
            np.all(masses>0) and masses[1]>masses[0], 'invalid mixed-scalar eigenmasses')
    U=matrix(handoff['basis']['unitary_flavour_rotation_U'],2)
    require(np.linalg.norm(U.conj().T@U-np.eye(2))<1e-10,'nonunitary scalar U')
    MF=matrix(uv['MF'])
    require(np.linalg.norm(MF-np.diag(np.diag(MF))) < 1e-9*max(1,np.linalg.norm(MF)),
            'running MF is not diagonal: heavy-fermion Takagi rotation required')
    fm=np.real(np.diag(MF))
    require(np.max(abs(np.imag(np.diag(MF))))<1e-9 and np.all(fm>masses[1]),
            'require real positive diagonal MF_i > upper scalar eigenmass')
    y1,y2=matrix(uv['y1']),matrix(uv['y2'])
    quartics=mapped_quartics(model,eft,audit)
    # Validate this native coefficient and its own Matchete source report.
    # The v19/v21 historical gates refer to an archived zero-mixing scenario;
    # they are not a prerequisite or a certificate for this running endpoint.
    direct_path=ROOT/src.BASE/'one_loop_probe'/f'{tag}_final_one_loop_probe_weinberg_ordered_coefficient.inputform.txt'
    source_report_path=ROOT/src.BASE/'one_loop_probe'/f'{tag}_final_one_loop_probe.json'
    require(direct_path.is_file() and direct_path.stat().st_size>0,
            'native Matchete coefficient missing/empty '+str(direct_path))
    source_report=read(source_report_path)
    require(source_report.get('format')=='t3-matchete-native-heavy-mass-eigenbasis-v1' and
            source_report.get('model')==model and
            source_report.get('one_loop_match_stage_completed') is True and
            source_report.get('weinberg_extraction_status')=='Success' and
            source_report.get('weinberg_term_count')==236 and
            source_report.get('weinberg_holomorphic_term_count')==118 and
            source_report.get('weinberg_conjugate_term_count')==118 and
            source_report.get('ordered_coefficient_is_canonical_flavour_symmetric_c5') is False and
            source_report.get('production_sequential_matching_performed') is False and
            source_report.get('numerical_c5_computed') is False and
            source_report.get('all_17_quartics_present') is True and
            source_report.get('both_yukawa_interactions_present') is True and
            source_report.get('heavy_mass_basis_detailed_check_passed') is True and
            source_report.get('weinberg_ordered_coefficient_inputform_file')==direct_path.name,
            'native one-loop Matchete source report missing, stale or incompatible: '+str(source_report_path))
    compiled,denom,_,_=src.compile_source(direct_path.read_text(encoding='utf-8-sig'),model)
    c5=np.zeros((3,3),complex)
    for i,mf in enumerate(fm):
        a11,a12,a22=src.direct_ordered(compiled,denom,float(mf),*masses,U,quartics)
        kernel=np.array([[2*a11,a12],[a12,2*a22]],complex)
        yi=np.column_stack((y1[:,i].conj(),y2[:,i].conj()))
        c5+=yi@kernel@yi.T
    residual=float(np.linalg.norm(c5-c5.T)/max(np.linalg.norm(c5),1e-250))
    require(residual<1e-11 and np.isfinite(c5).all(),'C5 non-symmetric/nonfinite')
    scales=running['scales_gev']
    mu_s=float(scales['S'])
    initial=SMWeinbergInitialConditions(
        gY=float(eft['gY']),g2=float(eft['g2']),g3=float(eft['g3']),
        lambdaH=float(eft['lambdaH']),
        ye=matrix(eft['ye']),yu=matrix(eft['yu']),yd=matrix(eft['yd']),K=c5)
    evolved=evolve_sm_weinberg(initial,mu_s,mu_low)
    rotated=charged_lepton_mass_basis_matrix(neutrino_mass_from_c5(evolved.K,vev_gev=vev),evolved.ye)
    observables=calculate_neutrino_observables(rotated,ordering='AUTO')
    masses_ev=np.sort(np.linalg.svd(rotated,compute_uv=False))*1e9
    report={'model':model,'status':'SuccessDirectOneScaleCandidate_NotSequentialMatching',
        'C5_source':'native_matchete_direct_one_loop_expression_evaluated_on_current_RGE_endpoints',
        'archived_C5_used':False,'mixed_scalar_diagonalization_used':True,
        'same_benchmark_input':True,'sequential_matching_certified':False,
        'physical_predictions_certified':False,
        'quartic_mapping_scope':'ConditionalNativeTensorMapAudited',
        'quartic_mapping_physical_certified':False,
        'native_one_loop_source_report_checked':True,
        'historical_v19_v21_topology_gates_used':False,
        'boundary_approximation':'MF and Yukawas from mu_F; scalar quartics/masses and SM couplings from mu_S; missing EFT1 Wilson running and finite threshold matching',
        'mu_F_gev':scales['F'],'mu_S_gev':mu_s,'mu_low_gev':mu_low,
        'C5_at_mu_S_GeV_inverse':encode(c5),'C5_at_mu_low_GeV_inverse':encode(evolved.K),
        'masses_ev':masses_ev.tolist(),
        'observables':encode_json(vars(observables)),
        'C5_symmetry_residual':residual,'solver_nfev':evolved.nfev,
        'source_hashes':{str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()
                        for p in (running_path,basis_path,SOURCE,direct_path,source_report_path,audit_path)},
    }
    save(out/f'{tag}_direct_candidate.json',report)
    print(f'{model}: native direct one-scale C5 candidate -> MZ; masses [eV]={report["masses_ev"]}',flush=True)
    return {'model':model,'status':report['status'],'masses_ev':report['masses_ev'],
            'file':str(out/f'{tag}_direct_candidate.json')}


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    selection=p.add_mutually_exclusive_group(required=True)
    selection.add_argument('--all',action='store_true')
    selection.add_argument('--model',choices=MODELS)
    p.add_argument('--running-root',type=Path,default=RUNNING_ROOT)
    p.add_argument('--output',type=Path,default=OUTPUT)
    p.add_argument('--mu-low',type=float,default=91.1876)
    p.add_argument('--vev',type=float,default=246.22)
    p.add_argument('--quartic-map-audit',type=Path,default=MAP_AUDIT)
    args=p.parse_args(argv)
    audit_path=args.quartic_map_audit.resolve()
    audit=read(audit_path)
    rows=[]
    for model in MODELS if args.all else (args.model,):
        try:
            rows.append(model_result(model,load_evaluator(),args.running_root.resolve(),
                                     args.output.resolve(),args.mu_low,args.vev,audit,audit_path))
        except Exception as exc:
            print(f'{model}: BLOCKED {type(exc).__name__}: {exc}',file=sys.stderr,flush=True)
            rows.append({'model':model,'status':'Blocked','error':str(exc)})
    save(args.output/'summary.json',{'scope':'same_run_native_one_scale_candidate',
        'physical_predictions_certified':False,'sequential_matching_certified':False,
        'quartic_mapping_scope':'ConditionalNativeTensorMapAudited',
        'quartic_mapping_physical_certified':False,'models':rows})
    success=sum(r['status'].startswith('Success') for r in rows)
    print(f'B/C current-running native C5 candidate: {success}/{len(rows)}')
    return 0 if success==len(rows) else 1

if __name__=='__main__':
    raise SystemExit(main())
