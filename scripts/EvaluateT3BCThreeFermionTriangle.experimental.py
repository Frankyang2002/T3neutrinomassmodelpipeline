#!/usr/bin/env python3
"""Exploratory B/C three-heavy-flavour neutrino mass topology, using v17.

NOT physical matching. Uses three frozen UV-reference Yukawa columns and MF
masses with current scalar mass-basis parameters. The signed triangle master
and ordered Warsaw projection are *source audited* in v16/v17. This script
only performs the remaining linear sum over the three NFlavor columns.

Reference Yukawa entries belong to a 1e7 GeV base state, while the heavy
fermions are 100/110/120 TeV. No running or parameter matching is performed.
The two auxiliary T3 quartics are deliberately ZERO for the illustrative
frozen-reference scenario, not inferred from a Matchete beta-function map.

Never label the output as renormalized physical C5, neutrino masses, or fit.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

MODELS = {
    'B_m1': ('B_-1', 1, 'T3_dS1_2_dS2_2_dF_1_alpha_m1.json'),
    'C_m1': ('C_-1', 3, 'T3_dS1_2_dS2_2_dF_3_alpha_m1.json'),
}
SQUARE_ROOT = Path('output/general_t3_mass_eigen_matchete')
AUDITS = Path('output/audits')
Q_NAMES=('lambdaT3','lambdaT3ExtraQ020020I1','lambdaT3ExtraQ020200I1')

class ReferenceError(ValueError):
    pass

def require(cond, message):
    if not cond:
        raise ReferenceError(message)

def load_json(path,fmt=None):
    require(path.is_file(), 'Required input missing: '+str(path))
    data=json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(data,dict),'Invalid source JSON object: '+str(path))
    if fmt is not None:
        require(data.get('format')==fmt,'Evidence schema changed: '+str(path))
    return data

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def as_matrix(obj, size, name):
    import numpy as np
    require(isinstance(obj,dict) and set(('real','imag'))<=set(obj),
            name+' needs explicit real and imag arrays')
    x=np.asarray(obj['real'],dtype=float)+1j*np.asarray(obj['imag'],dtype=float)
    require(x.shape==size and np.isfinite(x).all(),name+' invalid matrix shape or values')
    return x

def as_complex(obj,name):
    require(isinstance(obj,dict) and 'real' in obj and 'imag' in obj,
            name+' must be explicit complex object')
    z=complex(float(obj['real']),float(obj['imag']))
    require(math.isfinite(z.real) and math.isfinite(z.imag),name+' nonfinite')
    return z

def load_frozen_reference(root,stem):
    import numpy as np
    model,df,filename=MODELS[stem]
    path=root/'Numerical/benchmarks/fitted_reference'/filename
    data=load_json(path)
    rep=data.get('representation',{})
    require(tuple(rep.get(x) for x in ('d_s1','d_s2','d_f','alpha'))==(2,2,df,-1),
            model+': fitted reference representation changed')
    require(data.get('kind')=='t3_pipeline_numerical_results_v1',
            model+': unknown reference configuration')
    state=data.get('base_state',{})
    require(float(state.get('mu_gev',0))==1e7,
            model+': changed reference Yukawa definition scale')
    src=state.get('ordinary',{})
    y1=as_matrix(src['y1'],(3,3),'y1')
    y2=as_matrix(src['y2'],(3,3),'y2')
    M=as_matrix(src['MF'],(3,3),'MF')
    require(np.linalg.norm(M-np.diag(np.diag(M)))<1e-10,
            model+': fermion mass matrix no longer diagonal (Takagi needed)')
    mf=np.diag(M)
    require(np.max(np.abs(mf.imag))<1e-12 and np.all(mf.real>0),
            model+': MF not positive real masses')
    require(np.all(np.diff(mf.real)>0),model+': unexpected MF ordering')
    lam=as_complex(src['lambdaT3'],'lambdaT3')
    # This scenario assumes these two independent native quartics are zero.
    # Do not equate them to arbitrary RGBeta couplings without a verified map.
    quartics={Q_NAMES[0]:lam,Q_NAMES[1]:0j,Q_NAMES[2]:0j}
    vev=float(data['numerical']['vev_gev'])
    require(math.isfinite(vev) and vev>0,'Invalid vev')
    return {'model':model,'file':str(path.relative_to(root)),'file_sha256':sha(path),
            'y1':y1,'y2':y2,'mf':mf.real,'quartics':quartics,
            'vev_gev':vev,'frozen_source_mu_gev':float(state['mu_gev'])}

def polynomial_values(vec,substitutions):
    import numpy as np
    import sympy as sp
    result=[]
    for expression in vec:
        require(not (expression.free_symbols-set(substitutions)),
                'Unassigned native spurions: '+str(expression.free_symbols-set(substitutions)))
        result.append(complex(sp.N(expression.subs(substitutions),18)))
    z=np.asarray(result,dtype=complex)
    require(z.shape==(3,) and np.isfinite(z).all(),
            'Invalid three-Yukawa-structure coefficients')
    return z

def sum_fermions(weights,reference,masses,mode,v17):
    """Return a 3x3 source-normalized *reference topology*, NOT physical C5."""
    import numpy as np
    mf=reference['mf'];y1=reference['y1'];y2=reference['y2']
    MX,MY=masses
    out=np.zeros((3,3),complex)
    for col in range(3):
        M=float(mf[col]); A=M*M;B=MY*MY;C=MX*MX
        if mode=='strict':
            phi_xx=1+math.log(C/A)
            phi_yy=1+math.log(B/A)
            phi_xy=math.log(B/A)
        elif mode=='finite_XY':
            phi_xx=1+math.log(C/A)
            phi_yy=1+math.log(B/A)
            phi_xy=math.log(B/A)-v17.soft_factor(C/B)
        elif mode=='full_mass':
            phi_xx=float(v17.signed_triangle_master(A,C,C))
            phi_yy=float(v17.signed_triangle_master(A,B,B))
            phi_xy=float(v17.signed_triangle_master(A,B,C))
        else:
            raise ReferenceError('Unrecognized triangle comparison mode')
        Wxx,Wyy,Wxy=weights
        k=Wxx*phi_xx+Wyy*phi_yy+Wxy*phi_xy
        ker=np.array([[2*k[0],k[1]],[k[1],2*k[2]]],complex)
        # PL Wilson orientation determined in v11.2, not a Hermitian product.
        L=np.column_stack((y1[:,col].conj(),y2[:,col].conj()))
        out+=(L@ker@L.T)/(16*math.pi**2*M)
    require(np.all(np.isfinite(out)),'Nonfinite topology matrix')
    require(np.linalg.norm(out-out.T)<=2e-12*max(1e-30,np.linalg.norm(out)),
            'C5 reference must be complex symmetric (not Hermitian)')
    return out

def export_matrix(mat):
    return [[[float(z.real),float(z.imag)] for z in row] for row in mat]

def measure_topology(C5,vev):
    import numpy as np
    # 1 GeV = 1e9 eV; 1/2-normalized Weinberg convention.
    mnu=-0.5*vev**2*C5*1e9
    masses=np.sort(np.linalg.svd(mnu,compute_uv=False))
    return {
        'reference_C5_3x3_GeV_inverse_NOT_PHYSICAL':export_matrix(C5),
        'reference_mnu_3x3_eV_NOT_PHYSICAL':export_matrix(mnu),
        'reference_Takagi_singular_masses_eV_ascending_NOT_PHYSICAL':[float(v) for v in masses],
        'illustrative_delta_m21_squared_eV2':float(masses[1]**2-masses[0]**2),
        'illustrative_delta_m31_squared_eV2':float(masses[2]**2-masses[0]**2),
        'physical_C5':None,'physical_neutrino_masses':None,
    }

def self_test():
    import numpy as np
    import sympy as sp
    class Mini:
        @staticmethod
        def signed_triangle_master(a,b,c):
            # Synthetic test values are all deliberately equal/positive.
            return -0.5 if a==b==c else -1.0
        @staticmethod
        def soft_factor(r):
            return r*math.log(r)/(1-r)
    source={
        'mf':np.array([100000.,110000.,120000.]),
        'y1':np.array([[1,0,0],[0,1,0],[0,0,1]],complex)*0.001,
        'y2':np.zeros((3,3),complex), 'vev_gev':246.22,
    }
    weights=(np.array([0.,0.,0.]),np.array([1.,0.,0.]),np.array([0.,0.,0.]))
    C=sum_fermions(weights,source,(1000.,1150.),'strict',Mini)
    assert C.shape==(3,3) and np.max(np.abs(C-C.T))<1e-14
    # Each heavy column supplies a separate mass-denominator contribution.
    expected=np.array([2e-6*(1+math.log(1150.**2/M**2))/(16*math.pi**2*M)
                       for M in source['mf']])
    assert np.max(np.abs(np.diag(C)-expected))<1e-25
    assert np.max(np.abs(C-np.diag(np.diag(C))))==0
    stats=measure_topology(C,246.22)
    assert len(stats['reference_Takagi_singular_masses_eV_ascending_NOT_PHYSICAL'])==3
    u1,u2=sp.symbols('U1 U2')
    w=polynomial_values((u1+u2,u1-u2,u1*u2),{u1:1+1j,u2:2-1j})
    assert np.max(np.abs(w-np.array([3.,-1+2j,3+1j])))<1e-12
    try:
        polynomial_values((u1,u2,sp.Symbol('unknown')),{u1:1,u2:2})
    except ReferenceError:
        pass
    else:
        raise ReferenceError('Unassigned symbolic quartic accepted')
    print('PASS: three independent MF columns, symmetric Majorana projection, eV unit conversion')
    print('PASS: source polynomial substitution and unsupported-symbol rejection')
    print('PASS: physical C5 and neutrino masses remain explicitly uncertified')

def run(root,out):
    import numpy as np
    import sympy as sp
    from importlib.util import spec_from_file_location,module_from_spec
    def load_script(path,name):
        require(path.is_file(),'Missing upstream v16/v17 '+str(path))
        spec=spec_from_file_location(name,path)
        require(spec is not None and spec.loader is not None,'Cannot import '+name)
        module=module_from_spec(spec);spec.loader.exec_module(module)
        return module
    v16=load_script(root/'scripts/AssembleT3BCFormalFiniteScaleEnvelope.experimental.py','t3bc_v18_v16')
    v17=load_script(root/'scripts/CompareT3BCFiniteEnvelopeToTriangle.experimental.py','t3bc_v18_v17')
    v14=v16.read_json(root/AUDITS/'t3_bc_sequential_batch_v14/summary.json',
                      't3-bc-sequential-matching-fixed-order-batch-v14')
    v17audit=v16.read_json(root/AUDITS/'t3_bc_triangle_master_v17/summary.json',
                      't3-bc-independent-three-propagator-topology-comparison-v17')
    extract=v16.read_json(root/AUDITS/'t3_bc_scalar_finite_term_extraction_v1/summary.json',
                      't3-bc-native-local-matchete-pole-log-split-v1')
    require(v17audit.get('source_renormalized_matching_certified') is False and
            v17audit.get('physical_C5') is None,
            'v17 certification state unexpectedly changed')
    deps=v16.import_deps(root)
    _,kernel,v12,v13,v14script=deps
    result={'format':'t3-bc-frozen-reference-three-fermion-triangle-v18',
         'status':'IllustrativeFrozenReferenceThreeFermionTopology_NotPhysicalMatching',
         'scope':'B/C three heavy-flavour columns, source-audited scalar-loop topology only',
         'UV_yukawa_scale_frozen_GeV':1e7,
         'quartic_scenario':'lambdaT3 from fitted_reference, both lambdaT3ExtraQ020020I1 and lambdaT3ExtraQ020200I1 set explicitly to 0',
         'running_performed':False,'full_matching_scheme_certified':False,
         'full_uv_one_loop_vertex_normalization_certified':False,
         'physical_C5':None,'physical_neutrino_masses':None,
         'models':{}}
    for stem,(model,_,filename) in MODELS.items():
        prior=v17audit['models'][stem]
        for name in ('source_verified_F_finite_equals_F_pole',
                     'source_verified_K_X_constant_zero',
                     'source_verified_P_F_plus_P_Y_plus_P_X_zero',
                     'strict_pairwise_triangle_master_equals_v16_envelope'):
            require(prior.get(name) is True,model+': missing v17 gate '+name)
        sectors,poly=v16.compare_finite_and_poles(root,stem,v14,extract,deps)
        source=root/SQUARE_ROOT/'eft1_fermion_first_oneloop'/f'{stem}_final_eft1_fermion_first_oneloop_raw_hard_llhh_sector.inputform.txt'
        fsec=v13.split_sectors(source.read_text(encoding='utf-8-sig'))
        kf=v16.component_polynomials(fsec['constant'],model,kernel,v12,68)
        pf,py,px=poly['F'],poly['Y'],poly['X']
        require(all(sp.expand(kf[i]-pf[i])==0 for i in range(3)),model+': F finite/pole changed')
        ky=tuple(sp.expand(poly['constant'][i]-kf[i]) for i in range(3))
        wx=tuple(-p for p in px)
        wy=tuple(sp.expand(-py[i]+ky[i]) for i in range(3))
        wxy=tuple(-p for p in ky)
        params=load_json(root/SQUARE_ROOT/f'{stem}_final_parameters.json')
        require(params.get('model')==model,model+': mismatched mass basis')
        mx,my=map(float,params['heavy_masses_gev'])
        require(0<mx<my and math.isfinite(mx) and math.isfinite(my),model+': scalar spectrum invalid')
        U=as_matrix(params['unitary_U'],(2,2),'unitary_U')
        require(np.linalg.norm(U.conj().T@U-np.eye(2))<1e-9,model+': nonunitary mixing')
        reference=load_frozen_reference(root,stem)
        sub={}
        for i in range(2):
            for j in range(2):
                name=f'T3MEU{i+1}{j+1}'
                sub[sp.Symbol(name)]=U[i,j]
                sub[sp.Symbol('bar_'+name)]=U[i,j].conjugate()
        for name,val in reference['quartics'].items():
            sub[sp.Symbol(name)]=val
            sub[sp.Symbol('bar_'+name)]=val.conjugate()
        weights=tuple(polynomial_values(w,sub) for w in (wx,wy,wxy))
        cases={}
        matrices={}
        for mode in ('strict','finite_XY','full_mass'):
            c5=sum_fermions(weights,reference,(mx,my),mode,v17)
            matrices[mode]=c5
            cases[mode]=measure_topology(c5,reference['vev_gev'])
        full=matrices['full_mass']
        norm=float(np.linalg.norm(full))
        rels={mode:float(np.linalg.norm(matrices[mode]-full)/max(norm,1e-300))
              for mode in ('strict','finite_XY')}
        result['models'][stem]={
            'model':model,
            'fitted_reference_path':reference['file'],
            'fitted_reference_sha256':reference['file_sha256'],
            'matchete_mass_basis_path':str((SQUARE_ROOT/f'{stem}_final_parameters.json')),
            'mass_basis_sha256':sha(root/SQUARE_ROOT/f'{stem}_final_parameters.json'),
            'MF_GeV':[float(x) for x in reference['mf']],
            'MX_GeV':mx,'MY_GeV':my,
            'r_MX2_over_MY2':(mx/my)**2,
            'quartics_explicit_scenario':{n:[float(c.real),float(c.imag)] for n,c in reference['quartics'].items()},
            'reference_vev_GeV':reference['vev_gev'],
            'source_polynomial_gate':'v17 audited WXX/WYY/WXY with v16 live sources rechecked',
            'strict_and_full_triangle':cases,
            'relative_frobenius_deviation_from_full_mass':rels,
            'physical_C5':None,'physical_neutrino_masses':None,
        }
        vals=cases['full_mass']['reference_Takagi_singular_masses_eV_ascending_NOT_PHYSICAL']
        print(f'{model}: frozen 3-MF benchmark; M_X={mx:.4f} M_Y={my:.4f} GeV; '
              f'full-mass reference singular masses [eV]: '
              f'{vals[0]:.4g}, {vals[1]:.4g}, {vals[2]:.4g}; '
              f'strict/full relative difference={rels["strict"]:.4g}')
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print('PASS: source-audited three-heavy-flavour frozen-reference topology; NOT physical C5')
    print('Saved: '+str(out))

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path.cwd())
    p.add_argument('--output',type=Path,default=Path('output/audits/t3_bc_three_fermion_triangle_v18/summary.json'))
    p.add_argument('--self-test',action='store_true')
    p.add_argument('--require-physical',action='store_true')
    args=p.parse_args()
    try:
        if args.self_test:
            self_test();return 0
        if args.require_physical:
            raise ReferenceError('Physical C5 not established: missing matched running, full EFT finite counterterms and UV vertex sign')
        root=args.root.resolve()
        out=args.output if args.output.is_absolute() else root/args.output
        run(root,out)
        return 0
    except (ReferenceError,KeyError,ValueError,OSError,ImportError,TypeError,IndexError) as exc:
        print('FAIL (no physical C5): '+str(exc),file=sys.stderr)
        return 1

if __name__=='__main__':sys.exit(main())
