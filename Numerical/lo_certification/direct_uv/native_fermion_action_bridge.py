"""Read-only Workstream 4 native fermion-action/conditional component bridge.

This deliberately does NOT identify native Matchete GammaCC with an absolute
four-by-four C phase or RGBeta's symbolic SU(2) epsilon with concrete axes.
The frozen benchmark tensors are checked in the all-left-handed candidate
basis, not promoted to a physically source-normalised UV result.
"""
from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
from typing import Any

import numpy as np
import sympy as sp

CLASSES = {
    'A': (1, 3, 2, (-4, -2, 0)),
    'B': (2, 2, 1, (-3, -1, 1)),
    'C': (2, 2, 3, (-3, -1, 1)),
    'D': (3, 1, 2, (-2, 0, 2)),
    'E': (3, 3, 2, (-4, -2, 0, 2)),
}
REPRESENTATIVES = {'A_p0': 'A_+0', 'B_m1': 'B_-1', 'C_m1': 'C_-1'}
SRC = {
    'rgbeta': 'RGE/running/rgbeta/T3RGBetaModel.wl',
    'topology': 'Lagrangian/interactions/T3Topology.wl',
    'model': 'Lagrangian/model/T3Fields.wl',
    'builder': 'scripts/BuildT3DirectOneLoopInputs.py',
}
COMPONENT_AUDITS = Path('output/general_observable/lo_direct_one_loop/source_audits')
FROZEN = Path('output/general_observable/lo_direct_one_loop/provisional_inputs')
BENCHMARK = Path('Numerical/benchmarks/fitted_reference')


class FermionBridgeError(ValueError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise FermionBridgeError(message)


def load(path: Path) -> dict:
    require(path.is_file(), 'Missing evidence: '+str(path))
    data = json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(data, dict), 'Expected JSON object: '+str(path))
    return data


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def compact(source: str) -> str:
    # Source-level contract only; this isn't a Mathematica parser.
    return re.sub(r'\s+', '', re.sub(r'\(\*.*?\*\)', '', source, flags=re.S))


def guard_sources(root: Path) -> dict:
    data = {}
    for name, relative in SRC.items():
        p = root / relative
        require(p.is_file(), 'Missing original source '+str(p))
        data[name] = compact(p.read_text(encoding='utf-8-sig'))
    guards = {
        'rgbeta': (
            'AddYukawa[y1,{Bar@S1,l,F},',
            'AddYukawa[y2,{S2,l,secondHeavy},',
            'secondHeavy=If[yF===0,F,Fc];',
            'AddFermionMass[MF,{F,F},GroupInvariant->fermionInvariant,',
            'AddFermionMass[MF,{F,Fc},GroupInvariant->fermionInvariant,',
            'dF===2,(eps[T3RGBetaRep[dF],#1,#2]&)',
            'dF===3,(del[T3RGBetaRep[dF],#1,#2]&)',
        ),
        'topology': (
            'useCConj=which===1;',
            'scalarBar=which===2;',
            'y[p,r]NCM[Bar[l[i,p]],f]scalarCG[cg,labels]',
            'lambdaT3[]H[i]H[j]s1s2CG[T3MixCG,labels]',
        ),
        'model': (
            'SelfConjugate->True',
            'Mass->T3MassSpec["F",Lookup[f,"MassSymbol",MF]]',
        ),
        'builder': (
            "metric=np.asarray(us[df].T*us[df],complex)",
            "y1dag[:,index,s]+=np.conjugate(y1[:,r])*factor",
            "y2plain[:,r*df+f,ds1+s]+=np.conjugate(y2[:,r])*factor",
        ),
    }
    for name, needles in guards.items():
        for needle in needles:
            require(needle in data[name], 'Source contract changed '+name+': '+needle)
    return {name: {'path': rel, 'sha256': sha(root / rel)}
            for name, rel in SRC.items()}


def validate_native(rgb: dict, mat: dict, label: str) -> dict:
    require(rgb.get('format') == 't3-ws4-simultaneous-native-rgbeta-yukawa-mass-v1'
            and rgb.get('status') == 'NativeRGBetaYukawaAndMassRegistryVerified'
            and rgb.get('model') == label, label+': RGBeta native evidence failed')
    require(mat.get('format') == 't3-ws4-restricted-native-matchete-fermion-action-v1'
            and mat.get('status') == 'NativeRestrictedFermionActionCaptured'
            and mat.get('model') == label
            and mat.get('native_lagrangian_validation') is True,
            label+': Matchete native evidence failed')
    require(mat.get('checks') and all(v is True for v in mat['checks'].values()),
            label+': failed Matchete action checks')
    require(rgb.get('physical_Matchete_RGBeta_dictionary_certified') is False
            and mat.get('physical_Matchete_RGBeta_dictionary_certified') is False
            and rgb.get('four_spinor_to_weyl_phase_measured') is False
            and rgb.get('native_su2_numeric_component_axes_measured') is False
            and mat.get('physical_charge_conjugation_phase_measured') is False,
            label+': unsupported native physical certificate claim')
    expected_dims = {'A_p0':[1,3,2,0], 'B_m1':[2,2,1,-1], 'C_m1':[2,2,3,-1]}[label]
    require(rgb.get('dimensions_alpha') == expected_dims, label+': model dimensions changed')
    maj = label != 'A_p0'
    require(rgb.get('majorana') is maj, label+': wrong heavy fermion type')
    ys = rgb.get('yukawas')
    require(isinstance(ys,list) and len(ys)==2
            and [y.get('name') for y in ys]==['y1','y2'], label+': missing Yukawas')
    for y in ys:
        require(y.get('chirality')=='Right' and
                y.get('coupling_internal_factor') in ('2^Rational[1, 2]', 'Sqrt[2]', '2^(1/2)')
                and y.get('couplingbar_internal_factor') in ('2^Rational[1, 2]', 'Sqrt[2]', '2^(1/2)')
                and y.get('physical_action_coefficient_certified') is False,
                label+': unexpected Yukawa prefactors')
    require(ys[0]['fields_inputform'].replace(' ','') == '{Bar[S1],l,F}'
            and ys[1]['fields_inputform'].replace(' ','') == ('{S2,l,F}' if maj else '{S2,l,Fc}'),
            label+': Yukawa field orientation inconsistent')
    mass = rgb.get('mass')
    require(isinstance(mass,dict) and mass.get('fields_inputform').replace(' ','') ==
            ('{F,F}' if maj else '{F,Fc}')
            and mass.get('chirality') == 'Right'
            and mass.get('invariant_kind') == ('epsilon' if label=='A_p0' else 'singlet' if label=='B_m1' else 'delta')
            and mass.get('factor_is_physical_Lagrangian_coefficient') is False,
            label+': mass registry incompatible')
    factor = sp.sympify(mass['mass_internal_tensor_factor'].replace('^','**'),
                        locals={'Rational':sp.Rational,'Sqrt':sp.sqrt})
    require(bool(factor.is_number and factor.is_finite and factor!=0),
            label+': invalid internal MF tensor factor')
    lt = mat['yukawa_and_mixing_interactions_inputform']
    lf = mat['mass_free_lagrangian_inputform']
    require(all(x in lt for x in ('Coupling[y1,','Coupling[y2,','NewFermion','GammaCC')),
            label+': native assembled Yukawa expression incomplete')
    require('MF' in lf and 'NewFermion' in lf,
            label+': native assembled free mass expression incomplete')
    return {'model':REPRESENTATIVES[label],
            'native_rgbeta_yukawa_mass_simultaneously_registered':True,
            'native_matchete_action_extracted':True,
            'mf_internal_factor':str(factor),
            'mass_invariant':mass['invariant_kind'],
            'yukawa_internal_factors':'sqrt(2) for both registered Yukawas',
            'rgbeta_y2_heavy_field':'F' if maj else 'Fc',
            'matchete_gamma_cc_numerical_phase_fixed':False,
            'physical_cross_package_normalisation_certified':False}


def arr(obj: Any, label: str) -> np.ndarray:
    if isinstance(obj, dict):
        require(set(obj)=={'real','imag'}, label+': malformed complex array')
        a=np.asarray(obj['real'],dtype=float)+1j*np.asarray(obj['imag'],dtype=float)
    else:
        a=np.asarray(obj,dtype=complex)
    require(np.all(np.isfinite(a)), label+': nonfinite matrix')
    return a


def exact(value: str) -> complex:
    val=sp.sympify(value,locals={'I':sp.I,'Sqrt':sp.sqrt,'Rational':sp.Rational})
    require(bool(val.is_number and val.is_finite), 'Non-numeric exact CG value '+str(value))
    return complex(val.evalf(18))


def rotated_tensor(vertex: dict, flags: tuple[bool,...], us: dict[int,np.ndarray]) -> np.ndarray:
    dims = vertex['dimensions']
    require(len(flags)==len(dims),'Bad native tensor axes')
    maps=[np.conjugate(us[d]) if conj else us[d] for d,conj in zip(dims,flags)]
    result=np.zeros(tuple(dims),dtype=complex)
    for component in vertex['components']:
        idx=tuple(int(z)-1 for z in component['indices'])
        require(len(idx)==len(dims) and all(0<=v<d for v,d in zip(idx,dims)),
                'Native CG support index out of range')
        r=exact(component['cg'])
        outer=np.asarray(r,complex)
        for axis,old in enumerate(idx):
            outer=np.multiply.outer(outer,maps[axis][old])
        result +=outer
    return result


def rotations(analysis:dict) -> dict[int,np.ndarray]:
    us={1:np.ones((1,1),complex)}
    for rec in analysis['representations']:
        d=int(rec['dimension'])
        raw=rec['U_columns_are_weight_eigenvectors']
        u=np.array([[exact(z) for z in row] for row in raw],complex)
        require(u.shape==(d,d) and np.linalg.norm(u.conj().T@u-np.eye(d))<1e-10,
                'Generator weight transform not unitary')
        require(rec.get('status')=='generator_algebra_verified','Generator audit status changed')
        us[d]=u
    require(set(us)=={1,2,3},'Incomplete generator audit representations')
    return us


def benchmark_name(short: str) -> str:
    label=short[0];alpha=int(short[2:]);d1,d2,df,_=CLASSES[label]
    return f'T3_dS1_{d1}_dS2_{d2}_dF_{df}_alpha_{"m" if alpha<0 else "p"}{abs(alpha)}'


def charge(d:int,y:float,i:int)->float:
    return (d-1)/2-i+y


def conditional_mass_and_vertices(short:str, case:dict, bench:dict, inp:dict,
                                  us:dict[int,np.ndarray]) -> dict:
    label=short[0];alpha=int(short[2:]);d1,d2,df,_=CLASSES[label]
    cc=case['model'];require((cc['dS1'],cc['dS2'],cc['dF'],cc['alpha'])==
        (d1,d2,df,alpha),short+': CG case mismatch')
    top={k:rotated_tensor(case['vertices'][k],flags,us) for k,flags in {
        'T3Y1CG':(True,True,False), 'T3Y2CG':(True,False,True)}.items()}
    base=bench['base_state']['ordinary'];y1=arr(base['y1'],'y1');y2=arr(base['y2'],'y2')
    mf=arr(base['MF'],'MF'); require(all(x.shape==(3,3) for x in (y1,y2,mf)),
                              short+': unexpected flavour dimensions')
    maj=bool(df%2 and alpha==-1); nc=3*df; nf=nc if maj else 2*nc;ns=d1+d2
    metric=us[df].T@us[df]
    m=np.zeros((nf,nf),complex)
    for r in range(3):
        for t in range(3):
            if maj:
                m[r*df:(r+1)*df,t*df:(t+1)*df]=mf[r,t]*metric
            else:
                for k in range(df):
                    m[r*df+k,nc+t*df+k]=mf[r,t]
                    m[nc+t*df+k,r*df+k]=mf[r,t]
    require(np.linalg.norm(m-m.T)<1e-9*max(np.linalg.norm(m),1),
            short+': heavy Weyl mass symmetry broken')
    y=np.zeros((3,nf,ns),complex);z=np.zeros((3,nf,ns),complex)
    for r in range(3):
        for k in range(df):
            q=charge(df,(alpha+1)/2,k)
            for s in range(d1):
                fac=np.conjugate(top['T3Y1CG'][0,k,s]);
                if abs(fac)<1e-13:continue
                require(abs(q-charge(d1,alpha/2,s))<1e-10,
                        short+': y1 charge mismatch')
                ix=r*df+k if maj else nc+r*df+k
                z[:,ix,s]+=np.conjugate(y1[:,r])*fac
            for s in range(d2):
                fac=np.conjugate(top['T3Y2CG'][0,k,s]);
                if abs(fac)<1e-13:continue
                require(abs(q-charge(d2,(alpha+2)/2,s))<1e-10,
                        short+': y2 charge mismatch')
                if maj:
                    for j in range(df):
                        y[:,r*df+j,d1+s]+=np.conjugate(y2[:,r])*fac*metric[k,j]
                else:
                    y[:,r*df+k,d1+s]+=np.conjugate(y2[:,r])*fac
    im=arr(inp['fermion_mass_gev'],'frozen mass')
    yy=next(v for v in inp['vertices'] if v['name']=='y2_HC')
    zz=next(v for v in inp['vertices'] if v['name']=='y1_HC')
    yf=arr(yy['field'],'frozen y2');zf=arr(zz['conjugate'],'frozen y1')
    require(im.shape==m.shape and yf.shape==y.shape and zf.shape==z.shape,
            short+': frozen tensor shapes changed')
    def rel(a,b):return float(np.linalg.norm(a-b)/max(np.linalg.norm(b),1.e-25))
    residuals={'mass':rel(m,im),'y1_conjugate_scalar':rel(z,zf),'y2_field_scalar':rel(y,yf)}
    require(max(residuals.values())<2e-10,
            short+': frozen Weyl candidate differs from original CG source: '+str(residuals))
    return {'model':short,'fermion_type':'Majorana' if maj else 'Dirac',
            'n_heavy_weyl_components':nf, 'recomputed_candidate_vertex_and_mass_residuals':residuals,
            'native_absolute_spinor_phase_certified':False,
            'original_conditional_frozen_components_reproduced':True}


def evaluate(root:Path,native:Path, *, source:Path|None=None) -> dict:
    source=source or (root / COMPONENT_AUDITS)
    guard=guard_sources(root)
    for name in ('component_support_audit.json','matchete_generator_analysis.json',
                 'cg_covariance_audit.json'):
        require((source/name).is_file(),'Missing canonical source audit: '+str(source/name))
    cases=load(source/'component_support_audit.json')
    generators=load(source/'matchete_generator_analysis.json')
    covariance=load(source/'cg_covariance_audit.json')
    require(not cases.get('errors') and len(cases['cases'])==16,
            'Incomplete component-source coverage')
    require(covariance.get('status')=='cg_generator_covariance_verified',
            'CG covariance audit not verified')
    us=rotations(generators)
    lookup={(int(c['model']['dS1']),int(c['model']['dS2']),int(c['model']['dF']),int(c['model']['alpha'])):c
            for c in cases['cases']}
    cov={(int(c['model']['dS1']),int(c['model']['dS2']),int(c['model']['dF']),int(c['model']['alpha'])):c['all_covariant']
            for c in covariance['models']}
    native_rows=[]
    for label in REPRESENTATIVES:
        nr=load(native/(label+'_rgbeta.json'))
        nm=load(native/(label+'_matchete.json'))
        native_rows.append(validate_native(nr,nm,label))
    records=[]
    for label,(d1,d2,df,alphas) in CLASSES.items():
        for alpha in alphas:
            short=f'{label}_{alpha:+d}';key=(d1,d2,df,alpha)
            require(key in lookup and cov.get(key) is True,short+': no covariant native CG case')
            frozen_file=root/FROZEN/(short+'.json')
            if not frozen_file.is_file():frozen_file=root/'input/lo_direct_one_loop'/(short+'.json')
            inp=load(frozen_file)
            prov=inp['provenance']
            require(prov.get('cg_support_sha256')==sha(source/'component_support_audit.json')
                    and prov.get('generator_analysis_sha256')==sha(source/'matchete_generator_analysis.json'),
                    short+': frozen provenance hash differs from current component source')
            benchmark_file=root/BENCHMARK/(benchmark_name(short)+'.json')
            require(prov.get('benchmark_sha256')==sha(benchmark_file),
                    short+': frozen benchmark provenance hash mismatch')
            require(inp.get('model')==short,short+': wrong frozen model id')
            records.append(conditional_mass_and_vertices(short,lookup[key],load(benchmark_file),inp,us))
    require(len(records)==16,'Wrong all-16 model coverage')
    return {'format':'t3-ws4-native-fermion-spinor-bridge-v1',
            'status':'native_registry_and_candidate_weyl_components_verified_physical_phase_pending',
            'source_guards':guard,
            'native_simultaneous_yukawa_and_mass_representatives':native_rows,
            'native_representative_count':len(native_rows),
            'original_frozen_mass_and_yukawa_tensors_recomputed':len(records),
            'model_count':len(records),'models':records,
            'physical_source_to_weyl_normalisation_certified_models':0,
            'not_measured':[
                'Matchete installed numerical GammaCC phase in a chosen 4x4 representation',
                'RGBeta native SU2 epsilon/triplet axes aligned to Matchete',
                'cross-package signed and complex-conjugation-aware y1,y2,MF mapping',
                'physical relative lambdaT3 coupling across independent packages'],
            'EFT1_heavy_fermion_matching_performed':False,
            'SMEFT_scalar_matching_performed':False,
            'order':'UV -> integrate out all fermions -> integrate out both scalars -> SMEFT'}
