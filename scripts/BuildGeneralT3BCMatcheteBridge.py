"""Bridge current B/C normal-running endpoints to legacy Matchete mass-basis inputs.

No archived RG outputs or archived C5 are used. This reproduces the existing
RunGeneralT3MixedScalarMassBasis source-vertex transport using current values.
No physical Matchete/RGBeta basis map or matching is certified.

python -X utf8 -m scripts.BuildGeneralT3BCMatcheteBridge --all
"""
from __future__ import annotations
import argparse
from hashlib import sha256
import json
from pathlib import Path
import numpy as np
from Numerical.experimental.GeneralT3MixedScalarMassBasis import (
    MixedDoubletMassBasis, cjson, transport_operator,
    verify_su2_covariance, verify_tensor_transport, quadratic_mass_tensor_transport,
)
from scripts import RunGeneralT3MixedScalarMassBasis as old

ROOT=Path(__file__).resolve().parents[1]
INPUT=ROOT/'output/general_numerical_BC_normal'
OUTPUT=ROOT/'output/general_t3_mass_basis'
MODELS=('B_-1','C_-1')

def require(condition,reason):
    if not condition: raise ValueError('B/C Matchete bridge: '+reason)

def read(path):
    require(path.is_file(),'missing '+str(path))
    return json.loads(path.read_text(encoding='utf-8-sig'))

def array(raw):
    require(isinstance(raw,dict) and set(raw)=={'real','imag'},'invalid complex matrix representation')
    return np.asarray(raw['real'],float)+1j*np.asarray(raw['imag'],float)

def calculate(model,input_dir=INPUT,output_dir=OUTPUT):
    stem=model.replace('-','m')
    hpath=input_dir/f'{stem}_mass_basis_handoff.json'
    rpath=input_dir/'running'/f'{stem}.json'
    # Accept flat ZIP extraction for offline checks, but prefer actual on-disk structure.
    if not rpath.is_file() and (input_dir/f'{stem}.json').is_file():rpath=input_dir/f'{stem}.json'
    h,r=read(hpath),read(rpath)
    require(h.get('model')==model and r.get('model')==model,'model label mismatch')
    require(h.get('archived_C5_used') is False and h.get('C5_computed') is False,
            'unexpected prior C5 status')
    require(r.get('status')=='NumericalRenormalisableRunningOnly' and
            r.get('uv_solver',{}).get('success') is True and
            r.get('eft1_solver',{}).get('success') is True,'running not successful')
    require(r.get('physical_matching_performed') is False,'unexpected matching claim')
    state=r['eft1_at_scalar_threshold']
    mu=state.get('muT3GeneralS1S2')
    require(isinstance(mu,dict) and set(mu)=={'real','imag'},'missing complex scalar bilinear')
    basis=MixedDoubletMassBasis.build(model,float(state['mS1Sq']),
                                      float(state['mS2Sq']),complex(float(mu['real']),float(mu['imag'])))
    original=h['basis']
    require(np.allclose(array(original['mass_squared_flavour_matrix_gev2']),basis.A,rtol=1e-11,atol=1e-6),
            'running mass matrix differs from handoff')
    require(np.allclose(array(original['unitary_flavour_rotation_U']),basis.U,rtol=1e-10,atol=1e-10),
            'rotation differs from handoff')
    require(np.allclose(original['masses_squared_gev2'],basis.eigenvalues_sq,rtol=1e-10,atol=1e-5),
            'squared-mass spectrum changed')
    require(float(h['scale_gev'])==float(r['scales_gev']['S']),'scalar endpoint scale mismatch')
    manifest=read(old.MANIFEST);registry=read(old.REGISTRY)
    operators=old.quartic_operators(manifest,registry,model)
    require(len(operators)==17,'expected 17 independent quartics')
    result={
        'format':'t3-mixed-doublet-mass-basis-transport-v1',
        'status':old.STATUS,'model':model,'endpoint':'final',
        'scale_gev':float(h['scale_gev']),'basis':basis.info(),
        'quadratic_hessian_agreement':True,
        'exact_quadratic_tensor_transport':quadratic_mass_tensor_transport(basis),
        'max_su2_covariance_absolute_error':verify_su2_covariance(basis),
        'source_quartic_count':17,'source_new_quartics':1,
        'source_reused_quartics':6,'source_baseline_quartics':10,
        'quartics':{},'yukawa_scalar_leg_transforms':{},
        'physical_mass_matching_performed':False,
        'MatcheteRGBetaPhysicalCouplingMapVerified':False,
        'MatcheteHeavyMassBasisCheckOnRedefinedLagrangian':False,
        'source_eft1_dirac_yukawas_active':False,
        'source_sha256':{
            rpath.relative_to(ROOT).as_posix():sha256(rpath.read_bytes()).hexdigest(),
            old.MANIFEST.relative_to(ROOT).as_posix():sha256(old.MANIFEST.read_bytes()).hexdigest(),
            old.REGISTRY.relative_to(ROOT).as_posix():sha256(old.REGISTRY.read_bytes()).hexdigest(),
        },
        'bridge_source':'current_general_numerical_BC_normal',
        'bridge_handoff_sha256':sha256(hpath.read_bytes()).hexdigest(),
        'bridge_running_sha256':sha256(rpath.read_bytes()).hexdigest(),
    }
    max_error=0.
    for name,fields in operators.items():
        require(name in state,'missing scalar endpoint coupling '+name)
        z=old._coupling(state[name],name,real=name in old.REAL_QUARTICS)
        residual=verify_tensor_transport(fields,basis.U,trials=4)
        max_error=max(max_error,residual)
        result['quartics'][name]={
            'source_rgbeta_fields':list(fields),'source_rgbeta_coupling':cjson(z),
            'independent_tensor_preserved':True,
            'source_invariant_component_tensor_supplied':False,
            'coefficient_normalization_certified':False,
            'rotation_terms':transport_operator(fields,basis.U),
            'transport_identity_max_relative_error':residual,
        }
    for name,field in old.YUKAWA_SCALAR_LEGS.items():
        residual=verify_tensor_transport((field,),basis.U,trials=4)
        max_error=max(max_error,residual)
        result['yukawa_scalar_leg_transforms'][name]={
            'rgbeta_source_scalar_field':field,
            'non_scalar_indices_and_original_cg_unchanged':True,
            'yukawa_numerical_value_not_defined_in_eft1':True,
            'rotation_terms':transport_operator((field,),basis.U),
        }
    require(max_error<=3e-12,f'tensor transport failed {max_error}')
    result['max_source_tensor_transport_relative_error']=max_error
    output_dir.mkdir(parents=True,exist_ok=True)
    dest=output_dir/f'{stem}_final_mass_basis.json'
    dest.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    return dest,basis.eigenvalues_sq,max_error

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    sel=ap.add_mutually_exclusive_group(required=True)
    sel.add_argument('--all',action='store_true')
    sel.add_argument('--model',choices=MODELS)
    ap.add_argument('--input-dir',type=Path,default=INPUT)
    ap.add_argument('--output-dir',type=Path,default=OUTPUT)
    ns=ap.parse_args()
    for model in MODELS if ns.all else (ns.model,):
        dest,eig,err=calculate(model,ns.input_dir,ns.output_dir)
        print(f'{model}: Matchete mass-basis input created; eigenmasses {[float(x**0.5) for x in eig]}; transport error={err:.3g}; {dest}',flush=True)
    print('Physical matching NOT certified')
if __name__=='__main__':main()
