"""General T3: experimental EFT1 RG flow including all allowed scalar bilinears.

Reuses the Stage33 threshold boundary and the existing general-quartic state.
Regenerates *live* RGBeta EFT1 beta functions with the additional gauge-
invariant holomorphic scalar mass term for eight of the 16 ordinary models.

No production source or matching is changed. The scalar mass convention is a
DECLARED spherical SU(2)/RGBeta potential convention, not a certified physical
Matchete/RGBeta field or coupling transformation.

Examples:
    python -m scripts.RunGeneralT3QuadraticEFT1 --model B_-1 --mu2-fraction 0.01
    python -m scripts.RunGeneralT3QuadraticEFT1 --all --mu2-fraction 0.01
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
import subprocess
from typing import Any

import numpy as np
from scipy.integrate import solve_ivp

from scripts import Stage33FermionThresholdHandoff as s33
from scripts import Stage34IntermediateEFTRunning as s34
from Numerical.experimental.GeneralT3ScalarRGState import (
    GeneralT3ScalarRGState, load_registry, pack_general_scalar_state,
    unpack_general_scalar_state, registry_model,
)
from Numerical.experimental.GeneralT3ScalarQuadraticMasses import (
    load_manifest, quadratic_kind, quadratic_mass_hessian,
)
from Numerical.running.IntermediateScalarRGBetaEvaluator import (
    intermediate_scalar_state_environment,
)
from Numerical.core.RGBetaEvaluator import evaluate_inputform_expression
from Numerical.core.BetaVector import beta_values_to_derivative

ROOT = s34.ROOT
REGISTRY = ROOT / 'RGE/general/generated/general_rgbeta_quartic_registry.json'
MANIFEST = ROOT / 'docs/general_scalar_trial_operator_manifest.json'
WOLFRAM = ROOT / 'RGE/running/rgbeta/RunGeneralT3QuadraticRGBeta.wl'
OUTPUT = ROOT / 'output/general_t3_quadratic_eft1'
FORMAT = 'general-t3-candidate-quadratic-eft1-v1'
MASS = {'self_s1':'muT3GeneralS1Sq', 'self_s2':'muT3GeneralS2Sq',
        'mixed':'muT3GeneralS1S2'}
LOOP = 16.0 * math.pi**2


def fail(msg: str) -> None:
    raise ValueError('Candidate T3 quadratic EFT1: ' + msg)


def fingerprint(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def decode_complex(raw: Any, name: str) -> complex:
    if not isinstance(raw, dict) or set(raw) != {'real', 'imag'}:
        fail(name+': expected a real/imag object')
    val = complex(float(raw['real']),float(raw['imag']))
    if not np.isfinite([val.real,val.imag]).all():
        fail(name+': nonfinite complex scalar')
    return val


def command(model: str, rep: Any, path: Path, *, stage: str,
            exe: str='wolframscript') -> list[str]:
    if stage not in ('UV','EFT1'):
        fail('Invalid RGBeta stage')
    if model not in s34.MODELS:
        fail('Unknown model '+model)
    alpha_token = f'm{-rep.alpha}' if rep.alpha < 0 else str(rep.alpha)
    return [exe,'-file',str(WOLFRAM),stage,str(rep.d_s1),str(rep.d_s2),
            str(rep.d_f),alpha_token,model,str(path)]


def export_live(model: str, rep: Any, path: Path, *, stage: str,
                exe: str='wolframscript') -> None:
    if not WOLFRAM.is_file():
        fail('Missing RGBeta quadratic Wolfram runner: '+str(WOLFRAM))
    path.parent.mkdir(parents=True,exist_ok=True)
    path.unlink(missing_ok=True)
    cmd = command(model,rep,path,stage=stage,exe=exe)
    try:
        proc = subprocess.run(cmd,cwd=ROOT,check=False,capture_output=True,text=True)
    except OSError as e:
        fail(f'{model}/{stage}: wolframscript could not start: {e}')
    if proc.returncode or not path.is_file():
        fail(f'{model}/{stage}: RGBeta export failed ({proc.returncode})\n'
             +proc.stdout[-3500:]+'\n'+proc.stderr[-3500:])


def checked_payload(model: str, rep: Any, kind: str,
                    layout: Any, registry: dict, path: Path) -> dict:
    payload = s34.read_json(path)
    info = registry_model(registry,rep)
    meta = payload.get('metadata')
    if payload.get('format')!='general-t3-quadratic-rgbeta-experimental-v1' or \
            payload.get('status')!='SuccessCandidateMassBeta' or \
            payload.get('physical_mapping_certified') is not False or \
            not isinstance(meta,dict):
        fail(model+': malformed/incomplete live RGBeta export')
    expected = {
        'SharedScalar':False, 'dS1':rep.d_s1, 'dS2':rep.d_s2,
        'dF':rep.d_f, 'alpha':rep.alpha,
        'GeneralScalarModel':model, 'GeneralScalarStage':'EFT1',
        'GeneralScalarBasisStatus':'CandidateNativeRGBetaUnmappedToMatchete',
        'GeneralScalarQuarticBetaBasis':'RGBetaFullQuarticProjection',
        'GeneralScalarQuadraticMassMixing':'IncludedCandidate',
        'GeneralScalarAdditionalQuartics':info['total_additional_quartics'],
        'GeneralScalarExistingStructuralCount':len(info['preexisting_rgbeta_symbols']),
        'QuadraticKind':kind,'QuadraticMassSymbol':MASS.get(kind),
        'QuadraticPrefactorConvention':'V2 = mu2 * invariant + h.c.; no 1/2',
        'PhysicalMatcheteRGBetaMassBasisCertified':False,
        'ReportBetaConvention':'16*pi^2*dX/dln(mu); gauge BetaTerm divided by 2*g',
    }
    for key,val in expected.items():
        if meta.get(key)!=val:
            fail(f'{model}: inconsistent beta metadata {key}: {meta.get(key)!r} != {val!r}')
    wanted_new = {x['name'] for x in info['new_rgbeta_symbols']}
    if set(meta.get('GeneralScalarNewBetaKeys',[])) != wanted_new:
        fail(model+': wrong new-quartic inventory in quadratic beta system')
    wanted=set(layout.names) | ({MASS[kind]} if kind!='none' else set())
    report=payload.get('report_betas')
    if not isinstance(report,dict) or set(report)!=wanted or \
            set(payload.get('beta_symbols',[]))!=wanted:
        fail(f'{model}: incomplete mass/quartic beta inventory, '
             f'missing={sorted(wanted-set(report or {}))}, '
             f'excess={sorted(set(report or {})-wanted)}')
    if any(not isinstance(v,str) or not v or any(b in v for b in
           ('BetaTerm[','Missing[','$Failed','$Aborted')) for v in report.values()):
        fail(model+': unresolved symbolic beta expression')
    return report


def derivative(vector: np.ndarray, *, scale_gev: float, layout: Any,
               kind: str, registry: dict, report: dict) -> np.ndarray:
    size=layout.size
    if vector.shape != (size+(2 if kind!='none' else 0),):
        fail('Incorrect EFT1 extended ODE vector size')
    if not np.isfinite(vector).all():
        fail('Nonfinite EFT1 numerical state')
    base=unpack_general_scalar_state(vector[:size],layout,registry,
                                     uv=False,mu_gev=scale_gev)
    env=intermediate_scalar_state_environment(base.baseline)
    env.update(base.extra)
    if kind!='none':
        env[MASS[kind]]=complex(float(vector[size]),float(vector[size+1]))
    values={name:evaluate_inputform_expression(expr,env)
            for name,expr in report.items()}
    if kind=='none':
        return beta_values_to_derivative(values,layout,divide_by_loop_factor=True)
    mu_beta=values.pop(MASS[kind])
    try:
        mu_beta=complex(mu_beta)
    except (TypeError,ValueError) as e:
        fail('Invalid complex mass beta: '+str(e))
    if not np.isfinite([mu_beta.real,mu_beta.imag]).all():
        fail('Nonfinite complex quadratic mass beta')
    base_beta=beta_values_to_derivative(values,layout,divide_by_loop_factor=True)
    return np.concatenate((base_beta,[mu_beta.real/LOOP,mu_beta.imag/LOOP]))


def run_model(model: str, *, input_dir: Path, outdir: Path,
              mu2_fraction: float, mu2_phase: float,
              exe: str='wolframscript', export_uv: bool=False) -> dict:
    if model not in s34.MODELS:
        fail('Unsupported model: '+model)
    if not (math.isfinite(mu2_fraction) and 0<=mu2_fraction<=0.45):
        fail('--mu2-fraction must lie in [0,0.45]')
    if not math.isfinite(mu2_phase):
        fail('--mu2-phase must be finite')
    reg=load_registry(REGISTRY)
    manifest=load_manifest(MANIFEST)
    single, baseline, sidecar, reused, hashes=s34.verify_boundary(model,input_dir)
    info=registry_model(reg,baseline.representation)
    expected_new={x['name'] for x in info['new_rgbeta_symbols']}
    expected_reused={x['existing_name'] for x in info['preexisting_rgbeta_symbols']}
    if set(sidecar)!=expected_new or set(reused)!=expected_reused:
        fail(model+': Stage33 quartic sidecars are not the Stage29 inventory')
    initial=GeneralT3ScalarRGState(baseline,{
        key:decode_complex(val,key) for key,val in sidecar.items()}).validated(reg)
    kind=quadratic_kind(manifest,model)
    if (kind!='none') != bool(info['quadratics_not_implemented']):
        fail(model+': quadratic manifests disagree')
    mu_in=0j
    if kind!='none':
        scale=min(baseline.mS1Sq,baseline.mS2Sq)
        if scale<=0:
            fail(model+': nonpositive scalar mass-squared at threshold')
        mu_in=mu2_fraction*scale*complex(math.cos(mu2_phase),math.sin(mu2_phase))
    y0,layout=pack_general_scalar_state(initial,reg)
    if kind!='none':
        y0=np.concatenate((y0,[mu_in.real,mu_in.imag]))
    outdir.mkdir(parents=True,exist_ok=True)
    tag=s34.safe_model(model)
    export=outdir/f'{tag}_eft1_quadratic_live_rgbeta.json'
    export_live(model,baseline.representation,export,stage='EFT1',exe=exe)
    report=checked_payload(model,baseline.representation,kind,layout,reg,export)
    first=derivative(y0,scale_gev=baseline.mu_gev,layout=layout,
                     kind=kind,registry=reg,report=report)
    if first.shape!=y0.shape:
        fail(model+': initial derivative has incorrect dimension')
    def rhs(t: float,y:np.ndarray) -> np.ndarray:
        return derivative(y,scale_gev=float(np.exp(t)),layout=layout,
                          kind=kind,registry=reg,report=report)
    result=solve_ivp(rhs,(math.log(baseline.mu_gev),math.log(s34.STOP)),
                     y0,method='DOP853',rtol=1e-8,atol=1e-11)
    if not result.success or not np.isfinite(result.y).all():
        fail(model+': quadratic-plus-quartic integration failed: '+result.message)
    y_final=np.asarray(result.y[:,-1],dtype=float)
    full_final=unpack_general_scalar_state(y_final[:layout.size],layout,reg,
                                         uv=False,mu_gev=s34.STOP)
    mu_final=complex(*y_final[layout.size:]) if kind!='none' else 0j
    initial_spectrum=quadratic_mass_hessian(manifest,model,
         baseline.mS1Sq,baseline.mS2Sq,mu_in)
    final_spectrum=quadratic_mass_hessian(manifest,model,
         full_final.baseline.mS1Sq,full_final.baseline.mS2Sq,mu_final)
    uv_path=None
    if export_uv:
        uv_path=outdir/f'{tag}_uv_quadratic_live_rgbeta.json'
        export_live(model,baseline.representation,uv_path,stage='UV',exe=exe)
        uv=s34.read_json(uv_path)
        uv_meta=uv.get('metadata',{})
        uv_betas=uv.get('report_betas')
        if (uv.get('status')!='SuccessCandidateMassBeta' or
                uv.get('physical_mapping_certified') is not False or
                uv_meta.get('GeneralScalarStage')!='UV' or
                uv_meta.get('GeneralScalarModel')!=model or
                uv_meta.get('QuadraticKind')!=kind or
                uv_meta.get('QuadraticMassSymbol')!=MASS.get(kind) or
                not isinstance(uv_betas,dict) or
                set(uv_meta.get('GeneralScalarNewBetaKeys',[]))!=expected_new or
                (kind!='none' and MASS[kind] not in uv_betas) or
                any(not isinstance(v,str) or not v or
                    any(b in v for b in ('BetaTerm[','Missing[','$Failed','$Aborted'))
                    for v in uv_betas.values())):
            fail(model+': UV symbolic mass-beta export incomplete/malformed')
    data={
        'format':FORMAT,'status':'CandidateQuadraticAndFullQuarticEFT1Integrated',
        'model':model,'quadratic_kind':kind,'mu2_rgbeta_symbol':MASS.get(kind),
        'mu2_fraction_initial':mu2_fraction if kind!='none' else 0,
        'mu2_phase_radians_initial':mu2_phase if kind!='none' else 0,
        'mu2_initial_gev2':s33.json_state(mu_in),
        'mu2_final_gev2':s33.json_state(mu_final),
        'mu2_beta_included':kind!='none',
        'complete_one_loop_quartic_class_projected':True,
        'new_quartics_evolved':len(expected_new),
        'reused_quartics_evolved':len(expected_reused),
        'quadratic_operators_evolved':int(kind!='none'),
        'mu_initial_gev':baseline.mu_gev,'mu_final_gev':s34.STOP,
        'solver_nfev':int(result.nfev),
        'beta_keys':sorted(report),
        'initial_scalar_mass_spectrum':{
            key:value for key,value in initial_spectrum.items()
            if key!='real_scalar_hessian_gev2'},
        'final_scalar_mass_spectrum':{
            key:value for key,value in final_spectrum.items()
            if key!='real_scalar_hessian_gev2'},
        'initial_baseline':s33.json_state(initial.baseline),
        'final_baseline':s33.json_state(full_final.baseline),
        'initial_new_quartics':{key:s33.json_state(val)
                                for key,val in initial.extra.items()},
        'final_new_quartics':{key:s33.json_state(val)
                              for key,val in full_final.extra.items()},
        'input_sha256':{**hashes,str(REGISTRY):fingerprint(REGISTRY),
                        str(MANIFEST):fingerprint(MANIFEST),
                        str(WOLFRAM):fingerprint(WOLFRAM),
                        str(export):fingerprint(export),
                        **({str(uv_path):fingerprint(uv_path)} if uv_path else {})},
        'physical_matchete_rgbeta_mapping_certified':False,
        'triplet_component_basis_certified':False,
        'scalar_mass_threshold_matching_performed':False,
        'wilson_coefficient_matching_performed':False,
        'physical_neutrino_mass_prediction':False,
    }
    outfile=outdir/f'{tag}_eft1_quadratic_running.json'
    outfile.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    return data


def main() -> None:
    ap=argparse.ArgumentParser(description=__doc__)
    selection=ap.add_mutually_exclusive_group(required=True)
    selection.add_argument('--model',choices=s34.MODELS)
    selection.add_argument('--all',action='store_true')
    ap.add_argument('--mu2-fraction',type=float,default=0.0)
    ap.add_argument('--mu2-phase',type=float,default=0.0)
    ap.add_argument('--stage33-dir',type=Path,default=s34.DATA)
    ap.add_argument('--output-dir',type=Path,default=OUTPUT)
    ap.add_argument('--wolframscript',default='wolframscript')
    ap.add_argument('--also-export-uv',action='store_true')
    opts=ap.parse_args()
    rows=[]
    for model in (s34.MODELS if opts.all else (opts.model,)):
        row=run_model(model,input_dir=opts.stage33_dir,outdir=opts.output_dir,
                      mu2_fraction=opts.mu2_fraction,mu2_phase=opts.mu2_phase,
                      exe=opts.wolframscript,export_uv=opts.also_export_uv)
        rows.append(row)
        print(f'{model}: {row["new_quartics_evolved"]} new + '
              f'{row["reused_quartics_evolved"]} reused quartics; '
              f'quadratic={row["quadratic_kind"]}; '
              f'nfev={row["solver_nfev"]}; physical map UNCERTIFIED',flush=True)
    if opts.all:
        totals=[sum(row[key] for row in rows) for key in
            ('new_quartics_evolved','reused_quartics_evolved','quadratic_operators_evolved')]
        if len(rows)!=16 or totals!=[38,12,8]:
            fail(f'16-model coverage mismatch: {totals}')
        summary={
            'status':'Candidate16ModelQuadraticEFT1Integrated',
            'models_checked':16,'new_quartics_evolved':38,
            'reused_quartics_evolved':12,'quadratic_model_slots_evolved':8,
            'mu2_fraction':opts.mu2_fraction,'mu2_phase':opts.mu2_phase,
            'physical_mapping_certified':False,
            'matching_extended':False,
            'models':{row['model']:{'quadratic_kind':row['quadratic_kind'],
                    'nfev':row['solver_nfev'],
                    'tachyonic':row['final_scalar_mass_spectrum']['tachyonic']}
                    for row in rows}}
        p=opts.output_dir/'all_models_quadratic_eft1_summary.json'
        p.write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
        print('16/16 candidate quadratic+quartic EFT1 models; 38 new, 12 reused, 8 bilinears')

if __name__=='__main__':
    main()
