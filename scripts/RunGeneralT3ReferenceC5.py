"""Transfer frozen fitted-reference Yukawas to the verified X/Y mass-basis benchmark.

This is an explicitly HYBRID one-scale comparison: the frozen reference scalar
mass parameters differ from the independent experimental eigenstate inputs.
It neither reproduces the original fitted oscillation spectrum nor performs
fermion-first EFT1 matching. Nothing is written to the frozen reference files.

Run from the repository root:
    python -m scripts.RunGeneralT3ReferenceC5 --all
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import numpy as np

from scripts.RunGeneralT3NumericalC5 import (
    ROOT, DEFAULT_SOURCE, DEFAULT_PARAMETERS, MODELS, REFERENCE_FORMAT, REFERENCE_STATUS,
    fail, read_json, numerical_model,
)

MANIFEST = ROOT / 'Numerical/benchmarks/fitted_reference/benchmark_manifest.json'
DEFAULT_OUTPUT = ROOT / 'output/general_t3_numerical_c5_reference_transfer'
REFERENCE_ID = {
    'B_-1': ('T3_dS1_2_dS2_2_dF_1_alpha_m1', (2, 2, 1, -1)),
    'C_-1': ('T3_dS1_2_dS2_2_dF_3_alpha_m1', (2, 2, 3, -1)),
}
STATUS = REFERENCE_STATUS
REQUIRED_ZERO_QUARTICS = {
    'lambdaHHdagS2S2': 'lambdaT3ExtraQ020020I1',
    'lambdaHHdagS1barS1bar': 'lambdaT3ExtraQ020200I1',
}


def _matrix(data: dict, label: str, n: int = 3) -> list[list[dict]]:
    if not isinstance(data, dict) or set(data) != {'real', 'imag'}:
        fail(label+': expected complex real/imag 3x3 matrix')
    a=np.asarray(data['real'],dtype=float)
    b=np.asarray(data['imag'],dtype=float)
    if a.shape != (n,n) or b.shape != (n,n) or not (np.isfinite(a).all() and np.isfinite(b).all()):
        fail(label+': expected finite complex 3x3 matrix')
    return [[{'real':float(a[i,j]),'imag':float(b[i,j])} for j in range(n)] for i in range(n)]


def _diagonal_mass(data: dict) -> list[float]:
    if not isinstance(data,dict) or set(data)!={'real','imag'}:
        fail('reference MF must be a complex 3x3 diagonal matrix')
    a=np.asarray(data['real'],dtype=float)
    b=np.asarray(data['imag'],dtype=float)
    if a.shape!=(3,3) or b.shape!=(3,3) or not np.isfinite(a).all() or not np.isfinite(b).all():
        fail('invalid 3x3 frozen reference MF')
    if np.max(np.abs(b))>1e-12 or np.max(np.abs(a-np.diag(np.diag(a))))>1e-12:
        fail('reference fermion MF must be diagonal and real (no silent rotation)')
    masses=np.diag(a)
    if np.min(masses)<=0:fail('reference fermion MF must be positive')
    return masses.tolist()


def make_transfer(model: str, manifest: dict, frozen_path: Path, params_path: Path) -> dict:
    if model not in REFERENCE_ID:fail('unsupported frozen reference transfer model')
    model_id, dimensions=REFERENCE_ID[model]
    rows=manifest.get('models')
    if not isinstance(rows,list) or len(rows)!=16:fail('expected frozen 16-model reference manifest')
    matches=[row for row in rows if row.get('model')==model_id]
    if len(matches)!=1:fail(f'{model}: missing or duplicated frozen manifest model')
    row=matches[0]
    expected_path=Path('Numerical/benchmarks/fitted_reference')/(model_id+'.json')
    if row.get('config_path')!=expected_path.as_posix():fail('unexpected frozen reference config path')
    if not isinstance(row.get('sha256'),str) or len(row['sha256'])!=64:
        fail('missing reference SHA256')
    if sha256(frozen_path.read_bytes()).hexdigest()!=row['sha256']:
        fail(f'{model}: frozen reference hash changed; refuse mismatched inputs')
    frozen=read_json(frozen_path)
    r=frozen.get('representation')
    if not isinstance(r,dict) or tuple(r.get(x) for x in ('d_s1','d_s2','d_f','alpha'))!=dimensions or r.get('shared_scalar') is not False:
        fail('incorrect model representation in frozen reference')
    if frozen.get('benchmark_search',{}).get('enabled') is not False or frozen.get('sensitivity',{}).get('enabled') is not False:
        fail('reference unexpectedly enables a fit or scan')
    orig=frozen.get('base_state',{}).get('ordinary')
    if not isinstance(orig,dict):fail('reference missing original couplings')
    yuk1=_matrix(orig.get('y1'),'frozen y1')
    yuk2=_matrix(orig.get('y2'),'frozen y2')
    fermion_masses=_diagonal_mass(orig.get('MF'))
    topology=orig.get('lambdaT3')
    if not isinstance(topology,dict) or set(topology)!={'real','imag'}:
        fail('missing reference lambdaT3 complex input')
    try:
        topology_real=float(topology['real']); topology_imag=float(topology['imag'])
    except (TypeError,ValueError):fail('invalid reference lambdaT3')
    if not np.isfinite([topology_real,topology_imag]).all():fail('nonfinite reference topology coupling')
    quartics={'lambdaT3':{'real':topology_real,'imag':topology_imag}}
    for freeze,matchete in REQUIRED_ZERO_QUARTICS.items():
        raw=orig.get(freeze)
        if not isinstance(raw,dict) or set(raw)!={'real','imag'}:
            fail(f'reference missing complex quartic {freeze}')
        if complex(float(raw['real']),float(raw['imag']))!=0:
            fail(f'{freeze} is nonzero: Matchete/RGBeta normalization not certified')
        quartics[matchete]={'real':0.,'imag':0.}
    params=read_json(params_path)
    if params.get('model')!=model or not params.get('verified_eight_real_mode_hessian'):
        fail('experimental mass-basis parameters missing/mismatched')
    a=params.get('original_mass_squared_A_gev2',{})
    try:
        m2=np.asarray(a['real'],dtype=float)
        m2imag=np.asarray(a['imag'],dtype=float)
        ref2=np.array([float(orig['mS1Sq']),float(orig['mS2Sq'])])
    except (KeyError,ValueError,TypeError):fail('missing reference/mass-basis scalar masses')
    if m2.shape!=(2,2) or m2imag.shape!=(2,2) or ref2.shape!=(2,) or not np.isfinite(ref2).all() or not np.isfinite(m2).all():
        fail('invalid scalar squared masses')
    diff=(m2.diagonal()-ref2)
    v=float(frozen.get('numerical',{}).get('vev_gev',0))
    if not np.isfinite(v) or v<=0:fail('frozen reference missing valid vev')
    return {
        'format':REFERENCE_FORMAT,
        'benchmark_label':'hybrid_fitted_reference_yukawas_on_experimental_mass_basis',
        'reference_transfer_only':True,
        'mass_basis_compatible_with_frozen_reference':False,
        'physical_input_certified':False,
        'sequential_EFT1':False,
        'fermion_masses_gev':fermion_masses,
        'vev_gev':v,
        'y1':yuk1,'y2':yuk2,'quartics':quartics,
        'reference_provenance':{
            'frozen_config_path':expected_path.as_posix(),
            'frozen_sha256':row['sha256'],
            'original_reference_mS1Sq_gev2':float(ref2[0]),
            'original_reference_mS2Sq_gev2':float(ref2[1]),
            'experimental_mass_basis_original_A_diagonal_gev2':m2.diagonal().tolist(),
            'mass_squared_difference_gev2':diff.tolist(),
            'mass_squared_relative_difference':[float(diff[i]/ref2[i]) for i in range(2)],
            'frozen_yukawa_scale':'saved base_state values; no RG running applied',
            'matchete_extra_couplings_set_zero_only_because_frozen_inputs_are_zero':True,
            'matchete_rgbeta_nonzero_coupling_map_certified':False,
            'reproduces_original_fit':False,
        },
    }


def run_model(model: str, manifest_path: Path, benchmark_dir: Path,
              source_dir: Path, parameters_dir: Path, output_dir: Path):
    manifest=read_json(manifest_path)
    slug=REFERENCE_ID[model][0]
    reference_path=benchmark_dir/(slug+'.json')
    params_path=parameters_dir/(MODELS[model]+'_final_parameters.json')
    config=make_transfer(model,manifest,reference_path,params_path)
    result,dest=numerical_model(model,config,source_dir,parameters_dir,output_dir)
    if (result.get('n_fermion_flavours')!=3 or result.get('physical_c5_certified') is not False or
        result.get('status')!=STATUS):
        fail('unexpected transfer result flags')
    return result,dest


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--model',choices=tuple(REFERENCE_ID))
    group.add_argument('--all',action='store_true')
    parser.add_argument('--manifest',type=Path,default=MANIFEST)
    parser.add_argument('--benchmark-dir',type=Path,default=MANIFEST.parent)
    parser.add_argument('--source-dir',type=Path,default=DEFAULT_SOURCE)
    parser.add_argument('--parameters-dir',type=Path,default=DEFAULT_PARAMETERS)
    parser.add_argument('--output-dir',type=Path,default=DEFAULT_OUTPUT)
    args=parser.parse_args()
    items=[]
    for model in (tuple(REFERENCE_ID) if args.all else (args.model,)):
        result,dest=run_model(model,args.manifest,args.benchmark_dir,
                              args.source_dir,args.parameters_dir,args.output_dir)
        deviations=result['benchmark_config']['reference_provenance']['mass_squared_relative_difference']
        print(f'{model}: scalar source mismatch relative to frozen input: '
              f'[{deviations[0]:+.3%}, {deviations[1]:+.3%}] (squared masses)',flush=True)
        items.append({'model':model,'file':str(dest),'masses_eV':result['takagi_singular_values_eV_descending'],
                      'max_abs_entry_eV':result['max_abs_mnu_eV']})
        print(f'{model}: three-flavour HYBRID transfer; largest Takagi mass '
              f'{result["largest_takagi_mass_eV"]:.7g} eV; '
              f'masses={result["takagi_singular_values_eV_descending"]} eV; '
              'NOT the fitted reference result or EFT1 matching',flush=True)
    args.output_dir.mkdir(parents=True,exist_ok=True)
    dest=args.output_dir/'summary_final.json'
    dest.write_text(json.dumps({'status':STATUS,'models':items,
        'same_scalar_inputs_as_frozen_fit':False,'physical_c5_certified':False,
        'sequential_EFT1_performed':False},indent=2)+'\n',encoding='utf-8')
    print('Summary:',dest)


if __name__=='__main__':main()
