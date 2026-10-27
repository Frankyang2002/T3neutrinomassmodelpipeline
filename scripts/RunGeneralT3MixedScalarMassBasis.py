"""Experimental B_-1/C_-1 scalar mass diagonalization and full vertex transport.

Uses the *already computed* Stage quadratic+quartic EFT1 JSON; no Wolfram
recalculation, Matchete matching or production changes. All RGBeta quartic
source operators, including separate same-field contraction tensors, are
transported without equating Matchete and RGBeta normalizations.

Examples:
  python -m scripts.RunGeneralT3MixedScalarMassBasis --model B_-1
  python -m scripts.RunGeneralT3MixedScalarMassBasis --all
"""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
from typing import Any

import numpy as np

from Numerical.experimental.GeneralT3MixedScalarMassBasis import (
    MIXED_MODELS, SPECIES, EPS, MixedDoubletMassBasis, cjson,
    transport_operator, verify_su2_covariance, verify_tensor_transport,
    quadratic_mass_tensor_transport,
)
from Numerical.experimental.GeneralT3ScalarQuadraticMasses import (
    quadratic_mass_hessian,
)

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / 'output/general_t3_quadratic_eft1'
OUTPUT = ROOT / 'output/general_t3_mass_basis'
MANIFEST = ROOT / 'docs/general_scalar_trial_operator_manifest.json'
REGISTRY = ROOT / 'RGE/general/generated/general_rgbeta_quartic_registry.json'
STATUS = 'ExperimentalMassBasisTensorTransportNoPhysicalMatching'

# Exact field order from T3RGBetaModel.wl (both B/C alpha=-1):
BASE_QUARTICS = {
    'lambdaH':('Hdag','H','Hdag','H'),
    'lambdaS1':('S1dag','S1','S1dag','S1'),
    'lambdaS2':('S2dag','S2','S2dag','S2'),
    'lambdaH1':('Hdag','H','S1dag','S1'),
    'lambdaH2':('Hdag','H','S2dag','S2'),
    'lambda12':('S1dag','S1','S2dag','S2'),
    'lambdaT3':('H','H','S1','S2dag'),
    'lambdaH1Adj':('Hdag','H','S1dag','S1'),
    'lambdaH2Adj':('Hdag','H','S2dag','S2'),
    'lambda12Adj':('S1dag','S1','S2dag','S2'),
}
# The B/C(-1) reused operators' RGBeta order is verified from T3RGBetaModel.wl.
REUSED_QUARTICS = {
    'lambdaHHdagS2S2':('Hdag','Hdag','S2','S2'),
    'lambdaHHdagS1barS1bar':('Hdag','Hdag','S1dag','S1dag'),
    'lambdaS1bar2S2bar2':('S1dag','S1dag','S2dag','S2dag'),
    'lambdaS1barS2S2bar2':('S1dag','S2','S2dag','S2dag'),
    'lambdaS1S1bar2S2bar':('S1','S1dag','S1dag','S2dag'),
    'lambdaHHdagS1barS2barCross':('H','Hdag','S1dag','S2dag'),
}
# Yukawa scalar legs from T3RGBetaAddT3Yukawas: other CG slots (l,F)
# are not altered numerically; their SU(2) tensor must be supplied separately.
YUKAWA_SCALAR_LEGS = {'y1':'S1dag', 'y2':'S2'}
REAL_QUARTICS = frozenset(BASE_QUARTICS)-{'lambdaT3'}


def fail(message: str) -> None:
    raise ValueError('T3 mixed scalar mass-basis check: '+message)


def _read_json(path: Path) -> dict:
    if not path.is_file():
        fail(f'Missing {path}. Run RunGeneralT3QuadraticEFT1 first.')
    obj=json.loads(path.read_text(encoding='utf-8'))
    if not isinstance(obj,dict):
        fail('Expected JSON object in '+str(path))
    return obj


def _pair(obj: Any, label: str) -> complex:
    if not isinstance(obj,dict) or set(obj)!={'real','imag'}:
        fail(label+': expected {real,imag}')
    try:
        z=complex(float(obj['real']),float(obj['imag']))
    except (TypeError,ValueError) as ex:
        fail(label+': invalid complex pair: '+str(ex))
    if not np.isfinite([z.real,z.imag]).all():
        fail(label+': nonfinite complex value')
    return z


def _coupling(raw: Any, name: str, *, real: bool) -> complex:
    z=(_pair(raw,name) if isinstance(raw,dict) else complex(raw))
    if not np.isfinite([z.real,z.imag]).all():
        fail(name+': nonfinite coupling')
    if real and abs(z.imag)>1e-12*max(1,abs(z)):
        fail(name+': source real scalar quartic has nonzero imaginary component')
    return z


def _fields_from_counts(counts: list[int]) -> list[str]:
    if not isinstance(counts,list) or len(counts)!=6 or any(type(x)!=int or x<0 for x in counts):
        fail('Malformed registry/manifest field counts')
    return [f for f,k in zip(SPECIES,counts) for _ in range(k)]


def _field_counts(fields) -> tuple[int,...]:
    counter=Counter(fields)
    return tuple(counter.get(k,0) for k in SPECIES)


def _rgbeta_fields(fields: list[str]) -> list[str]:
    normal={'H':'H','Bar @ H':'Hdag','S1':'S1','Bar @ S1':'S1dag',
            'S2':'S2','Bar @ S2':'S2dag'}
    try:
        return [normal[x] for x in fields]
    except KeyError as exc:
        fail('Unrecognized RGBeta field: '+repr(exc))


def quartic_operators(manifest:dict, registry:dict, model:str) -> dict[str,tuple[str,...]]:
    if model not in MIXED_MODELS or model not in manifest.get('models',{}) or \
            model not in registry.get('models',{}):
        fail('Missing B/C(-1) source operator manifest')
    m=manifest['models'][model]
    r=registry['models'][model]
    if registry.get('totals',{}).get('matchete_extra_quartics')!=50 or \
            registry.get('totals',{}).get('rgbeta_new_vertices')!=38 or \
            registry.get('totals',{}).get('rgbeta_preexisting_structural')!=12:
        fail('General quartic registry global inventory changed')
    if m.get('extra_quartic_complex_couplings')!=7 or \
            m.get('quadratic_complex_couplings_not_enabled')!=1 or \
            r.get('total_additional_quartics')!=7 or \
            r.get('quadratics_not_implemented')!=1 or \
            r.get('dimensions')!=[2,2,1 if model.startswith('B') else 3] or \
            r.get('alpha')!=-1:
        fail(model+': unexpected manifest/registry model structure')
    if len(r.get('new_rgbeta_symbols',[]))!=1 or \
            set(x['existing_name'] for x in r.get('preexisting_rgbeta_symbols',[])) != set(REUSED_QUARTICS):
        fail(model+': missing/changed reused or new quartic registry')
    if Counter({tuple(x['counts']):x['invariants'] for x in m['extra_quartic_structures']}) != \
            Counter({tuple(x['counts']):1 for x in r['preexisting_rgbeta_symbols']} |
                    {tuple(x['counts']):2 for x in r['new_rgbeta_symbols']}):
        fail(model+': independent Matchete operator multiplicities differ from RGBeta slots')
    operators=dict(BASE_QUARTICS)
    operators.update(REUSED_QUARTICS)
    for item in r['preexisting_rgbeta_symbols']:
        name=item['existing_name']
        if _field_counts(operators[name])!=tuple(item['counts']):
            fail(name+': registered field-content mismatch')
    for item in r['new_rgbeta_symbols']:
        name=item['name']
        fields=tuple(_rgbeta_fields(item['fields']))
        if fields!=('H','Hdag','S1dag','S2dag') or \
                _field_counts(fields)!=tuple(item['counts']) or \
                name!='lambdaT3GeneralQ110101I1' or name in operators:
            fail(model+': unexpected new independent quartic structure')
        operators[name]=fields
    if len(operators)!=17 or len([k for k in operators if k.startswith('lambdaT3General')])!=1:
        fail('Expected all 10 ordinary + 6 reused + 1 new quartic operators')
    return operators


def _state_value(state: dict, name: str, sidecar: dict) -> complex:
    if name in sidecar:
        return _coupling(sidecar[name],name,real=False)
    if name=='lambdaH':
        raw=state.get('sm')
        if not isinstance(raw,dict): fail('Missing nested SM state')
        if 'lambdaH' not in raw: fail('Missing SM Higgs quartic')
        return _coupling(raw['lambdaH'],name,real=True)
    if name not in state or state[name] is None:
        fail('Missing EFT1 source quartic: '+name)
    return _coupling(state[name],name,real=(name in REAL_QUARTICS))


def run_model(model: str, *, folder: Path=INPUT, output: Path=OUTPUT,
              which: str='final', trials: int=4) -> dict:
    if model not in MIXED_MODELS:
        fail('Only B_-1 and C_-1 contain off-diagonal doublet bilinears')
    if which not in ('initial','final') or not (1<=trials<=100):
        fail('Invalid stage or trial count')
    input_file=folder/f'{model.replace("-","m").replace("+","p")}_eft1_quadratic_running.json'
    data=_read_json(input_file)
    if data.get('format')!='general-t3-candidate-quadratic-eft1-v1' or \
            data.get('status')!='CandidateQuadraticAndFullQuarticEFT1Integrated' or \
            data.get('model')!=model or data.get('quadratic_kind')!='mixed' or \
            data.get('mu2_beta_included') is not True or \
            data.get('physical_matchete_rgbeta_mapping_certified') is not False or \
            data.get('scalar_mass_threshold_matching_performed') is not False:
        fail('Incompatible candidate EFT1 report or physical certification flags')
    manifest=_read_json(MANIFEST)
    registry=_read_json(REGISTRY)
    operators=quartic_operators(manifest,registry,model)
    state=data.get(which+'_baseline')
    extra=data.get(which+'_new_quartics')
    if not isinstance(state,dict) or not isinstance(extra,dict):
        fail('Missing complete quadratic EFT1 endpoint state')
    mu=_pair(data.get('mu2_'+which+'_gev2'), 'mu2_'+which)
    basis=MixedDoubletMassBasis.build(model,state['mS1Sq'],state['mS2Sq'],mu)
    prototype=quadratic_mass_hessian(manifest,model,basis.m_s1_sq,basis.m_s2_sq,mu)
    expected=np.repeat(basis.eigenvalues_sq,4)
    masses=np.asarray(prototype['mass_squared_eigenvalues_gev2'],dtype=float)
    if masses.shape!=(8,) or not np.allclose(masses,expected,rtol=2e-11,atol=1e-5):
        fail(model+': 8-real-mode Hessian disagrees with complex-doublet eigenvalues')
    saved=data.get(which+'_scalar_mass_spectrum',{})
    old=np.asarray(saved.get('mass_squared_eigenvalues_gev2',[]),dtype=float)
    if old.shape!=(8,) or not np.allclose(old,expected,rtol=1e-9,atol=1e-5):
        fail(model+': previously saved numerical Hessian spectrum disagrees')
    U=basis.U
    gauge_error=verify_su2_covariance(basis)
    result={'format':'t3-mixed-doublet-mass-basis-transport-v1',
            'status':STATUS,'model':model,'endpoint':which,
            'scale_gev':data['mu_'+('final' if which=='final' else 'initial')+'_gev'],
            'basis':basis.info(),'quadratic_hessian_agreement':True,
            'exact_quadratic_tensor_transport':quadratic_mass_tensor_transport(basis),
            'max_su2_covariance_absolute_error':gauge_error,
            'source_quartic_count':17,'source_new_quartics':1,
            'source_reused_quartics':6,'source_baseline_quartics':10,
            'quartics':{},'yukawa_scalar_leg_transforms':{},
            'physical_mass_matching_performed':False,
            'MatcheteRGBetaPhysicalCouplingMapVerified':False,
            'MatcheteHeavyMassBasisCheckOnRedefinedLagrangian':False,
            'source_eft1_dirac_yukawas_active':False,
            'source_sha256':{
                p.relative_to(ROOT).as_posix():sha256(p.read_bytes()).hexdigest()
                for p in (input_file,MANIFEST,REGISTRY) if p.is_relative_to(ROOT)},
           }
    max_error=0.
    for name,fields in operators.items():
        z=_state_value(state,name,extra)
        residual=verify_tensor_transport(fields,U,trials=trials)
        max_error=max(max_error,residual)
        result['quartics'][name]={
            'source_rgbeta_fields':list(fields),'source_rgbeta_coupling':cjson(z),
            'independent_tensor_preserved':True,
            'source_invariant_component_tensor_supplied':False,
            'coefficient_normalization_certified':False,
            'rotation_terms':transport_operator(fields,U),
            'transport_identity_max_relative_error':residual,
        }
    for name,field in YUKAWA_SCALAR_LEGS.items():
        max_error=max(max_error,verify_tensor_transport((field,),U,trials=trials))
        result['yukawa_scalar_leg_transforms'][name]={
            'rgbeta_source_scalar_field':field,
            'non_scalar_indices_and_original_cg_unchanged':True,
            'yukawa_numerical_value_not_defined_in_eft1':True,
            'rotation_terms':transport_operator((field,),U),
        }
    if max_error>3e-12: fail('Tensor transport failed numerical identity')
    result['max_source_tensor_transport_relative_error']=max_error
    path=output/f'{model.replace("-","m").replace("+","p")}_{which}_mass_basis.json'
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    return result


def main() -> None:
    p=argparse.ArgumentParser(description=__doc__)
    sel=p.add_mutually_exclusive_group(required=True)
    sel.add_argument('--model',choices=MIXED_MODELS)
    sel.add_argument('--all',action='store_true')
    p.add_argument('--which',choices=['initial','final'],default='final')
    p.add_argument('--input-dir',type=Path,default=INPUT)
    p.add_argument('--output-dir',type=Path,default=OUTPUT)
    p.add_argument('--trials',type=int,default=4)
    a=p.parse_args()
    for name in MIXED_MODELS if a.all else (a.model,):
        result=run_model(name,folder=a.input_dir,output=a.output_dir,
                         which=a.which,trials=a.trials)
        masses=result['basis']['masses_gev']
        print(f'{name}: rotated all 17 quartic tensors + 2 Yukawa legs; '
              f'mass eigenvalues [GeV]={masses}; '
              f'identity error={result["max_source_tensor_transport_relative_error"]:.2e}; '
              'physical matching NOT performed',flush=True)
    print('Mass-basis tensor transport: experimental only; Matchete/RGBeta map UNCERTIFIED')

if __name__=='__main__':
    main()
