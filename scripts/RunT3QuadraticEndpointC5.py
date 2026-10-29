"""Conditional B/C(-1) EFT1 scalar-endpoint -> one-scale Matchete C5 bridge.

Consumes saved candidate quadratic EFT1 states, the frozen *reference Yukawas*,
and an independently extracted one-scale Matchete expression. This is NOT
sequential EFT matching. Candidate RGBeta quartic-to-Matchete signs are exposed.

python -m scripts.RunT3QuadraticEndpointC5 --all
"""
from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path
import numpy as np
import mpmath as mp

from Numerical.experimental.GeneralT3MixedScalarMassBasis import MixedDoubletMassBasis
from scripts.RunGeneralT3NumericalC5 import (ROOT, DEFAULT_SOURCE, DEFAULT_PARAMETERS,
    compile_coefficient, complex_value, matrix_values, read_json)
from scripts.RunGeneralT3ReferenceC5 import MANIFEST, REFERENCE_ID, make_transfer

MODELS = ('B_-1', 'C_-1')
DEFAULT_RUNNING = ROOT / 'output/general_t3_quadratic_eft1'
DEFAULT_OUTPUT = ROOT / 'output/studies/general_scalar_extended/quadratic_endpoint_c5'
QUARTIC_MAP = {'lambdaT3ExtraQ020020I1':'lambdaHHdagS2S2',
               'lambdaT3ExtraQ020200I1':'lambdaHHdagS1barS1bar'}


def require(condition, message):
    if not condition:
        raise ValueError('Quadratic endpoint C5: ' + message)


def pair(raw):
    require(isinstance(raw, dict) and set(raw) == {'real','imag'}, 'invalid complex pair')
    val = complex(float(raw['real']), float(raw['imag']))
    require(np.isfinite([val.real,val.imag]).all(), 'nonfinite complex pair')
    return val


def cj(z):
    z=complex(z)
    return {'real':float(z.real),'imag':float(z.imag)}


def endpoint(model, folder, which):
    stem=model.replace('-','m').replace('+','p')
    path=folder / f'{stem}_eft1_quadratic_running.json'
    data=read_json(path)
    require(data.get('format')=='general-t3-candidate-quadratic-eft1-v1' and
            data.get('status')=='CandidateQuadraticAndFullQuarticEFT1Integrated' and
            data.get('model')==model and data.get('quadratic_kind')=='mixed' and
            data.get('mu2_beta_included') is True, 'missing compatible candidate quadratic flow')
    require(data.get('physical_matchete_rgbeta_mapping_certified') is False and
            data.get('scalar_mass_threshold_matching_performed') is False and
            data.get('wilson_coefficient_matching_performed') is False,
            'unexpected upstream certification or changed flow contract')
    base=data.get(which+'_baseline')
    require(isinstance(base,dict), 'missing endpoint baseline')
    mu=pair(data.get('mu2_'+which+'_gev2'))
    basis=MixedDoubletMassBasis.build(model,base['mS1Sq'],base['mS2Sq'],mu)
    require(not basis.tachyonic and not basis.degenerate and min(basis.eigenvalues_sq)>0,
            'unphysical or degenerate scalar endpoint')
    extra=data.get(which+'_new_quartics')
    require(isinstance(extra,dict), 'missing endpoint quartics')
    vals={k:pair(base[k]) for k in QUARTIC_MAP.values()}
    require('lambdaT3' in base and isinstance(base['lambdaT3'],dict), 'missing endpoint topology quartic')
    return path, data, basis, base, extra, vals


def numerical_c5(expression, config, basis, quartics):
    mp.mp.dps=55
    code=compile_coefficient(expression)
    u=basis.U
    env={f'T3MEU{i+1}{j+1}':mp.mpc(str(u[i,j].real),str(u[i,j].imag))
         for i in range(2) for j in range(2)}
    env.update({'MSX':mp.mpf(str(np.sqrt(basis.eigenvalues_sq[0]))),
                'MSY':mp.mpf(str(np.sqrt(basis.eigenvalues_sq[1]))),
                'hbar':1/(16*mp.pi**2), 'conj':mp.conj,'sqrt':mp.sqrt,'log':mp.log})
    for name,z in quartics.items():
        env[name]=mp.mpc(str(z.real),str(z.imag))
    y1=matrix_values(config['y1'],'y1',rows=3)
    y2=matrix_values(config['y2'],'y2',rows=3)
    require(len(y1[0])==len(y2[0])>=1,'Yukawa columns inconsistent')
    masses=config.get('fermion_masses_gev')
    require(isinstance(masses,list) and len(masses)==len(y1[0]), 'missing reference fermion masses')
    mass=[mp.mpf(str(v)) for v in masses]
    require(all(mp.isfinite(v) and v>0 for v in mass),'invalid fermion masses')
    ordered=[[mp.mpc(0) for _ in range(3)] for _ in range(3)]
    for i in range(3):
        for j in range(3):
            for k,mf in enumerate(mass):
                local=dict(env,MF=mf,y1a=y1[i][k],y1b=y1[j][k],
                           y2a=y2[i][k],y2b=y2[j][k])
                val=eval(code,{'__builtins__':{}},local)
                require(mp.isfinite(mp.re(val)) and mp.isfinite(mp.im(val)), 'nonfinite coefficient')
                ordered[i][j]+=val
    c5=np.array([[complex(ordered[i][j]+ordered[j][i]) for j in range(3)]
                 for i in range(3)],dtype=complex)
    require(np.isfinite(c5).all() and np.allclose(c5,c5.T,rtol=1e-12,atol=1e-24),
            'non-symmetric or nonfinite C5')
    vev=float(config['vev_gev'])
    require(np.isfinite(vev) and vev>0,'invalid vev')
    mnu=-vev**2/2*c5*1e9
    masses_ev=np.linalg.svd(mnu,compute_uv=False)
    return c5,mnu,masses_ev


def calculate_model(model, *, which, running_dir, source_dir, parameters_dir,
                    frozen_dir, manifest_path, output_dir, include_candidate_extras):
    require(model in MODELS,'unsupported model')
    runpath,data,basis,base,extra,rgbeta=endpoint(model,running_dir,which)
    model_id=REFERENCE_ID[model][0]
    frozen=frozen_dir/(model_id+'.json')
    source_params=parameters_dir/(model.replace('-','m')+'_final_parameters.json')
    transfer=make_transfer(model,read_json(manifest_path),frozen,source_params)
    stem=model.replace('-','m')+'_final'
    probe=read_json(source_dir/(stem+'_one_loop_probe.json'))
    require(probe.get('status')=='OneLoopDirectSymbolicWeinbergProbePassedNumericalC5Pending' and
            probe.get('model')==model and probe.get('one_loop_match_stage_completed') is True and
            probe.get('physical_matching_performed') is False and
            probe.get('production_sequential_matching_performed') is False,
            'one-scale Matchete source certificate missing')
    source=source_dir/(stem+'_one_loop_probe_weinberg_ordered_coefficient.inputform.txt')
    expression=source.read_text(encoding='utf-8')
    topology=pair(base['lambdaT3'])
    common={'lambdaT3':topology}
    for mat in QUARTIC_MAP:
        common[mat]=0j
    c0,m0,s0=numerical_c5(expression,transfer,basis,common)
    mapped=dict(common)
    if include_candidate_extras:
        for mat,rgb in QUARTIC_MAP.items():
            mapped[mat]=-rgbeta[rgb]  # Source-derived only for these two B/C(-1) field words.
    c1,m1,s1=numerical_c5(expression,transfer,basis,mapped)
    delta=m1-m0
    denom=float(np.linalg.norm(m0))
    result={
        'status':'ConditionalEFT1ScalarEndpointOneScaleC5NotPhysicalSequential',
        'model':model,'endpoint':which,'mu_gev':data['mu_'+which+'_gev'],
        'running_report':str(runpath.relative_to(ROOT)),
        'running_sha256':sha256(runpath.read_bytes()).hexdigest(),
        'one_scale_expression_sha256':sha256(source.read_bytes()).hexdigest(),
        'frozen_yukawa_source_sha256':sha256(frozen.read_bytes()).hexdigest(),
        'quadratic_mu2_gev2':cj(basis.mu_sq),
        'scalar_masses_gev':[float(np.sqrt(x)) for x in basis.eigenvalues_sq],
        'mixing_sin2theta':basis.mixing_sin_2theta,
        'candidate_rgbeta_endpoint_quartics':{k:cj(v) for k,v in rgbeta.items()},
        'matchete_quartics_used':{k:cj(v) for k,v in mapped.items()},
        'extra_quartic_map':'lambda_Matchete=-lambda_RGBeta (two source-level B/C structures only)' if include_candidate_extras else 'extra quartics set to zero for control',
        'baseline_one_scale_mnu_ev':[[cj(v) for v in row] for row in m0],
        'candidate_one_scale_mnu_ev':[[cj(v) for v in row] for row in m1],
        'delta_mnu_frobenius_ev':float(np.linalg.norm(delta)),
        'relative_delta_mnu_frobenius':float(np.linalg.norm(delta)/denom) if denom>0 else None,
        'baseline_takagi_masses_ev_descending':s0.tolist(),
        'candidate_takagi_masses_ev_descending':s1.tolist(),
        'physical_matchete_rgbeta_mapping_certified':False,
        'scalar_threshold_matching_performed':False,
        'sequential_EFT1_wilson_running_performed':False,
        'yukawas_evolved_with_candidate_scalar_background':False,
        'physical_neutrino_mass_prediction':False,
        'note':'Frozen reference Yukawas + candidate EFT1 scalar endpoint + original direct one-scale Matchete expression. No fermion/scalar sequential matching, no full basis certification.',
    }
    output_dir.mkdir(parents=True,exist_ok=True)
    path=output_dir/(model.replace('-','m')+'_'+which+'_endpoint_one_scale.json')
    path.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    return result,path


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    sel=ap.add_mutually_exclusive_group(required=True)
    sel.add_argument('--all',action='store_true')
    sel.add_argument('--model',choices=MODELS)
    ap.add_argument('--which',choices=('initial','final'),default='final')
    ap.add_argument('--running-dir',type=Path,default=DEFAULT_RUNNING)
    ap.add_argument('--source-dir',type=Path,default=DEFAULT_SOURCE)
    ap.add_argument('--parameters-dir',type=Path,default=DEFAULT_PARAMETERS)
    ap.add_argument('--frozen-dir',type=Path,default=MANIFEST.parent)
    ap.add_argument('--manifest',type=Path,default=MANIFEST)
    ap.add_argument('--output-dir',type=Path,default=DEFAULT_OUTPUT)
    ap.add_argument('--zero-extra-quartics',action='store_true')
    args=ap.parse_args()
    rows=[]
    for model in MODELS if args.all else (args.model,):
        result,path=calculate_model(model,which=args.which,running_dir=args.running_dir,
            source_dir=args.source_dir,parameters_dir=args.parameters_dir,
            frozen_dir=args.frozen_dir,manifest_path=args.manifest,
            output_dir=args.output_dir,include_candidate_extras=not args.zero_extra_quartics)
        rows.append({'model':model,'file':str(path.relative_to(ROOT)),
            'masses_eV':result['candidate_takagi_masses_ev_descending'],
            'relative_delta':result['relative_delta_mnu_frobenius']})
        print(model,'candidate masses (eV):',result['candidate_takagi_masses_ev_descending'],
              'relative extra-quartic delta:',result['relative_delta_mnu_frobenius'],
              'PHYSICAL SEQUENTIAL C5: NOT CERTIFIED',flush=True)
    summary=args.output_dir/'summary.json'
    summary.write_text(json.dumps({'status':'ConditionalOneScaleEndpointOnly',
        'models':rows,'physical_neutrino_mass_prediction':False},indent=2)+'\n',encoding='utf-8')
    print('Summary:',summary)

if __name__=='__main__':
    main()
