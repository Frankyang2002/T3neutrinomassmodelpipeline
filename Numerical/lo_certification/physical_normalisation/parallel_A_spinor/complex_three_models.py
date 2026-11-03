#!/usr/bin/env python3
"""Track A: three-model chosen-Weyl complex source experiment.

No frozen UV neutrino matrix, provisional engine, matching or RGE imports.
`--mode native` verifies *existing* installed-Matchete witnesses and source hashes;
`--mode audited-constants` uses independently transcribed, published native CGs,
NOT a new installed-native test. Both modes leave source certification FALSE.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np
import mpmath as mp

REL = Path('Numerical/lo_certification/physical_normalisation/parallel_A_spinor')
NATIVE = {
    'A_+0': ('A_p0', 'output/general_observable/physical_normalisation/a0_native_quadratic_20261103_113729/A_+0_native_action.json', 't3-b1-a0-live-native-scalar-quadratic-v1'),
    'B_-1': ('B_m1', 'output/general_observable/physical_normalisation/bm1_native_v2_20261103_121426/B_-1_native_action_v2.json', 't3-b1-bm1-live-native-scalar-quadratic-v2'),
    'C_-1': ('C_m1', 'output/general_observable/physical_normalisation/cm1_native_20261103_122426/C_-1_native_action.json', 't3-b1-cm1-live-native-scalar-quadratic-v1'),
}
CG_DECLARED = {
    'A_+0': (1+0j, -1j*math.sqrt(2)*3**(-.25), 1j*math.sqrt(2)*3**(-.25)),
    'B_-1': (1+0j, 1+0j, -2/math.sqrt(3)+0j),
    'C_-1': (3**(-.25)+0j, -3**(-.25)+0j, -2/math.sqrt(3)+0j),
}

def require(t, message):
    if not t:
        raise ValueError(message)

def rel(x, y):
    return float(np.linalg.norm(x-y)/max(np.linalg.norm(x), np.linalg.norm(y), 1e-90))

def enc(z):
    a=np.asarray(z)
    return {'real':a.real.tolist(),'imag':a.imag.tolist()}

def scalar(z):
    z=complex(z)
    return {'real':z.real,'imag':z.imag}

def fingerprint(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def read_witness(path, root, expected_format, expected_model):
    require(path.is_file(), f'Missing already-produced native witness: {path}')
    row=json.loads(path.read_text(encoding='utf-8-sig'))
    require(row.get('format')==expected_format and row.get('model')==expected_model, 'Native model/format mismatch')
    require(row.get('physical_source_certified') is False,'Witness incorrectly declares physical certification')
    checks=row.get('native_component_checks',row.get('native_checks',{}))
    require(checks and all(x is True for x in checks.values()), 'Native CG checks absent/failed')
    for file, sha in row['source_sha256'].items():
        p=root/file
        require(p.is_file() and fingerprint(p).lower()==sha.lower(), f'Native source SHA256 drift: {file}')
    return row

def native_tensor(row, key, dims):
    data=row['native_tensors'][key]
    require(tuple(data['shape'])==dims, f'{key}: wrong native tensor shape')
    z=np.array(data['real_flat'],float)+1j*np.array(data['imag_flat'],float)
    require(z.size==math.prod(dims) and np.isfinite(z).all(), f'{key}: invalid native CG')
    return z.reshape(dims)

def witness_cg(label,row):
    c=3**(-.25)
    if label=='A_+0':
        a=native_tensor(row,'T3Y1CG',(2,2))
        b=native_tensor(row,'T3Y2CG',(2,2,3))
        h=native_tensor(row,'T3MixCG',(2,2,3))
        u=np.array([1j,0,1],complex)/math.sqrt(2)
        vals=(a[0,1], np.dot(b[0,1,:],u.conj()), np.dot(h[1,1,:],u.conj()))
    elif label=='B_-1':
        a=native_tensor(row,'T3Y1CG',(2,2))
        b=native_tensor(row,'T3Y2CG',(2,2))
        h=native_tensor(row,'T3MixCG',(2,2,2,2))
        vals=(a[0,0],b[0,1],h[1,1,0,1])
    else:
        a=native_tensor(row,'T3Y1CG',(2,3,2))
        b=native_tensor(row,'T3Y2CG',(2,3,2))
        h=native_tensor(row,'T3MixCG',(2,2,2,2))
        require(np.count_nonzero(abs(a[0,:,0])>1e-12)==1 and np.count_nonzero(abs(b[0,:,1])>1e-12)==1,
                'Triplet neutral charge channel not unique')
        vals=(a[0,1,0],b[0,1,1],h[1,1,0,1])
    require(max(abs(x-y) for x,y in zip(vals,CG_DECLARED[label]))<2e-13,
            f'{label}: native charge-projected CG mismatch')
    return tuple(complex(x) for x in vals)

def native_mode(root,rule_witness):
    cg={}; hashes={}
    for label,(short,file,fmt) in NATIVE.items():
        path=root/file
        row=read_witness(path,root,fmt,label)
        cg[label]=witness_cg(label,row)
        hashes[file]=fingerprint(path)
    # This existing witness is source-native, but only abstract GammaCC.
    file='output/general_observable/physical_normalisation/b1_complex_spinor_20261103_124038/native_spinor_action.json'
    path=root/file
    require(path.is_file(),f'Missing original B1 spinor witness: {file}')
    row=json.loads(path.read_text(encoding='utf-8-sig'))
    require(row.get('format')=='t3-b1-native-spinor-convention-v1' and
            row.get('physical_source_certified') is False, 'Unacceptable spinor native witness')
    require(row.get('installed_gammaCC',{}).get('explicit_numeric_4x4_exposed') is False,
            'Changed installed spinor observation requires manual re-audit')
    for file2,sha in row['source_sha256'].items():
        require(fingerprint(root/file2).lower()==sha.lower(),f'Spinor witness source drift: {file2}')
    for model in ('A_p0','B_m1','C_m1'):
        m=row['models'][model]
        require(m['native_action_passed'] is True and all(m['checks'].values()),f'{model}: native spinor contract failed')
        require(m['mass']['coefficients']==(['-1'] if model=='A_p0' else ['-1/2']),f'{model}: mass factor changed')
        for key in ('y1_original','y1_hc','y2_original','y2_hc'):
            require(m[key]['coefficients']==['1'],f'{model}: Yukawa/HC coefficient changed')
    hashes[file]=fingerprint(path)
    require(rule_witness is not None and rule_witness.is_file(),
            'Missing new Matchete rule-level witness; run A2 first')
    require(rule_witness.resolve().is_relative_to(root), 'Rule witness outside local project')
    rules=json.loads(rule_witness.read_text(encoding='utf-8-sig'))
    require(rules.get('format')=='t3-trackA-native-feynman-rule-probe-v1' and
            rules.get('status')=='NativeMatcheteRulesCaptured_RelativeWeylPhaseNotEstablished'
            and rules.get('physical_source_certified') is False and
            rules.get('native_feynman_rules_invoked') is True and
            rules.get('isolated_wrong_HC_sign_negative_control_detected') is True,
            'New native Matchete rule output missing or unverified')
    require(set(rules.get('models',{}))=={'A_p0','B_m1','C_m1'},'Incomplete A2 native rule suite')
    for model,record in rules['models'].items():
        require(record['negative_control'] and all(record['negative_control'].values()) and
                record.get('relative_numeric_GammaCC_to_Weyl_phase_measured') is False,
                f'{model}: rule-level negative control not validated or phase status changed')
    for path,sha in rules['source_sha256'].items():
        require(fingerprint(root/path).lower()==sha.lower(), f'A2 native source changed: {path}')
    hashes[str(rule_witness.relative_to(root))]=fingerprint(rule_witness)
    return cg,hashes

def group_metrics():
    eps=np.array([[0,1],[-1,0]],complex)
    # Matchete source SU(2) convention: epsilon_12=+1, not delta for two doublets.
    R=np.array([[.8,.6],[-.6,.8]],complex)
    epsres=rel(R.T@eps@R,eps)
    T3=np.array([[0,0,-1j],[0,0,0],[1j,0,0]],complex)
    up=np.array([-1j,0,1],complex)/math.sqrt(2)
    zero=np.array([0,1,0],complex)
    um=np.array([1j,0,1],complex)/math.sqrt(2)
    U=np.column_stack([up,zero,um]);G=U.T@U
    want=np.array([[0,0,1],[0,1,0],[1,0,0]],complex)
    weights=np.diag([1,0,-1])
    metricres=rel(G,want)
    charge_res=float(np.linalg.norm(weights.T@G+G@weights))
    invalid_metric_control=float(np.linalg.norm(weights.T@np.eye(3)+np.eye(3)@weights))
    require(epsres<1e-14 and metricres<1e-14 and charge_res<1e-14 and invalid_metric_control>2,
            'Doublet epsilon or triplet Majorana metric inconsistent')
    return {'doublet_epsilon_invariance':epsres,'triplet_G_equals_UtU':metricres,
            'triplet_UG_charge_residual':charge_res,'wrong_triplet_identity_charge_violation':invalid_metric_control,
            'triplet_metric':enc(G)}

def draw_source(seed):
    rng=np.random.default_rng(seed)
    # Deliberately chosen nontrivial complex couplings, independent of benchmark fits.
    y1=(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))*8.5e-5
    y2=(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))*9.3e-5
    masses=np.array([100000.,111000.,123000.])
    lam=0.13*np.exp(0.29j)
    return y1,y2,masses.astype(complex),lam

def mpz(z):
    return mp.mpc(repr(float(complex(z).real)),repr(float(complex(z).imag)))

def two_methods(M2,s_lo,s_hi):
    delta=s_hi-s_lo
    require(s_lo>0 and delta>0,'Invalid scalar eigenvalues')
    # Independent integral and separate closed form; no provisional engine imported.
    quad=-mp.quad(lambda x:mp.log1p((1-x)*delta/(x*M2+(1-x)*s_lo)),[0,1])/delta
    def b0(s): return 1+s*mp.log(s/M2)/(M2-s)
    analytic=(b0(s_hi)-b0(s_lo))/delta
    return quad,analytic

def evaluate(y1,y2,masses,lam,cg,wrong_y1=False,wrong_lam=False,wrong_mass=False,
             bad_symfactor=False):
    mp.mp.dps=55
    c1,c2,cm=cg
    a=mp.mpf('1306122.448979592');b=mp.mpf('1000000');v=mp.mpf('246.22')
    A21=-mpz(lam.conjugate() if wrong_lam else lam)*v*v*mpz(cm)/2
    spread=mp.sqrt(((a-b)/2)**2+abs(A21)**2)
    lo=(a+b)/2-spread;hi=(a+b)/2+spread
    # Source y1,y2 are coefficients of original 4-spinor Lagrangian. The
    # selected left-Weyl external holomorphic vertices use the H.c. sources.
    h1=(-np.conj((y1.conj() if wrong_y1 else y1)*c1))
    h2=(-np.conj(y2*c2))
    out_quad=mp.matrix(3);out_b0=mp.matrix(3)
    for k,M in enumerate(masses):
        physical_sq=mp.mpf(repr(float(abs(M)**2)))
        q,closed=two_methods(physical_sq,lo,hi)
        for x, dd in ((out_quad,q),(out_b0,closed)):
            # For L_mass=-M psi psi / 2 + H.c., <psi psi> has M*
            # in its numerator. Wrong control deliberately uses M instead.
            numerator=mpz(M if wrong_mass else M.conjugate())
            factor=numerator*A21*dd/(16*mp.pi**2)
            if bad_symfactor:factor/=2
            for i in range(3):
                for j in range(3):
                    x[i,j]+=factor*(mpz(h1[i,k])*mpz(h2[j,k])+mpz(h2[i,k])*mpz(h1[j,k]))
    x=np.array([[complex(out_quad[i,j])*1e9 for j in range(3)] for i in range(3)])
    z=np.array([[complex(out_b0[i,j])*1e9 for j in range(3)] for i in range(3)])
    return x,rel(x,z)

def cp_odd(m):
    H=m@m.conj().T
    return float(np.imag(H[0,1]*H[1,2]*H[2,0])/(np.trace(H).real/3)**3)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mode',choices=('audited-constants','native'),default='audited-constants')
    p.add_argument('--repo-root',type=Path,default=Path('.'))
    p.add_argument('--rule-witness',type=Path,default=None)
    p.add_argument('--output-dir',type=Path,required=True)
    args=p.parse_args()
    root=args.repo_root.resolve()
    out=args.output_dir.resolve()
    require(not out.exists(),f'Non-overwrite policy: output already exists: {out}')
    if args.mode=='native':
        require((root/'Lagrangian/model/T3ModelCatalog.wl').is_file(),'Not a mounted project root')
        cgs,hashes=native_mode(root,args.rule_witness)
    else:
        cgs=CG_DECLARED;hashes={}
    group=group_metrics()
    results={}
    for label,seed in (('A_+0',704),('B_-1',805),('C_-1',906)):
        y1,y2,masses,lam=draw_source(seed)
        original,diff=evaluate(y1,y2,masses,lam,cgs[label]);require(diff<5e-12,'Feynman / B0 disagreement')
        # Positive controls: source global phases AND independent heavy and
        # scalar field rephasings of every source monomial, not global mnu sign.
        a,b,c,d=(.23,-.31,.17,-.41)
        yf1=y1*np.exp(1j*(-a+c));yf2=y2*np.exp(1j*(-b-d))
        Mf=masses*np.exp(1j*(a+b));lf=lam*np.exp(1j*(c-d))
        transported,dr=evaluate(yf1,yf2,Mf,lf,cgs[label]); require(dr<5e-12,'transport loop method mismatch')
        phaseres=rel(original,transported);require(phaseres<1e-12,'Failed Dirac/fermion/scalar rephase invariance')
        # In Majorana models a == b; impose the same Weyl-field rephase.
        if label!='A_+0':
            alpha=.337
            mm,mr=evaluate(y1*np.exp(-1j*alpha),y2*np.exp(-1j*alpha),
                           masses*np.exp(2j*alpha),lam,cgs[label]);require(mr<5e-12,'Majorana loop method mismatch')
            majoranares=rel(original,mm);require(majoranares<1e-12,'Majorana rephase covariance failure')
        else:majoranares=None
        # Global spurion source phases e^(i(phi_lambda-phi_y1-phi_y2)).
        x,y,z=.21,-.17,.43
        glob,_=evaluate(y1*np.exp(1j*x),y2*np.exp(1j*y),masses,lam*np.exp(1j*z),cgs[label])
        globalres=rel(glob,original*np.exp(1j*(z-x-y)));require(globalres<1e-12,'Wrong global phase transport')
        bad_y1,_=evaluate(y1,y2,masses,lam,cgs[label],wrong_y1=True)
        bad_lam,_=evaluate(y1,y2,masses,lam,cgs[label],wrong_lam=True)
        bad_mass,_=evaluate(y1,y2,masses,lam,cgs[label],wrong_mass=True)
        # Complex heavy-mass control uses a genuinely nonreal mass point.
        mm=masses*np.exp(1j*np.array([.24,-.39,.52]))
        exact,_=evaluate(y1,y2,mm,lam,cgs[label]);wrong,_=evaluate(y1,y2,mm,lam,cgs[label],wrong_mass=True)
        bad_factor,_=evaluate(y1,y2,masses,lam,cgs[label],bad_symfactor=True)
        zero_mix,_=evaluate(y1,y2,masses,0.+0j,cgs[label])
        require(np.linalg.norm(zero_mix)==0, 'lambdaT3=0 did not decouple loop')
        ctr={'wrong_y1_conjugation':rel(original,bad_y1),'wrong_lambda_conjugation':rel(original,bad_lam),
             'wrong_complex_mass_conjugation':rel(exact,wrong),'extra_half_factor':rel(original,bad_factor)}
        require(min(ctr.values())>1e-3,'Negative control did not fail')
        require(rel(original,original.T)<1e-14,'Neutrino Majorana symmetry failed')
        # A complete -mnu is not a physical sign/phase discriminator.
        require(np.allclose(np.linalg.svd(original,compute_uv=False),np.linalg.svd(-original,compute_uv=False)),
                'Global unphysical sign changed Takagi masses')
        results[label]={
            'cg_assumption':{'c1':scalar(cgs[label][0]),'c2':scalar(cgs[label][1]),'cm':scalar(cgs[label][2])},
            'source_complex_point':{'y1':enc(y1),'y2':enc(y2),'mass_gev':enc(masses),'lambda':scalar(lam)},
            'mnu_ev':enc(original),'takagi_singular_masses_ev_descending':np.linalg.svd(original,compute_uv=False).tolist(),
            'cp_odd_rephasing_invariant':cp_odd(original),'matrix_symmetry_residual':rel(original,original.T),
            'feynman_quadrature_vs_closed_B0':diff,'source_rephasing_covariance_residual':phaseres,
            'global_spurion_phase_residual':globalres,'majorana_field_rephase_residual':majoranares,
            'negative_controls':ctr,'wrong_overall_sign_not_physical_discriminator':True,
            'lambda_zero_neutrino_mass_norm_ev':float(np.linalg.norm(zero_mix)),
            'relative_source_spinor_phase_certified':False,'physical_source_certified':False}
    data={'format':'t3-trackA-complex-three-models-v1','mode':args.mode,'github_head_audited':'d2163e6bb95109c0fac419f7e41f335f0d347b62',
          'source_native_witness_sha256':hashes,'group_tests':group,'models':results,
          'one_loop_computed_without_frozen_uv_or_provisional_engine':True,
          'native_GammaCC_relative_phase_fixed':False,'physical_source_certified':False,'certified_count_of_16':0}
    out.mkdir(parents=True,exist_ok=False)
    target=out/'three_model_complex_results.json'
    target.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
    for label,r in results.items():
        print(f'{label}: B0={r["feynman_quadrature_vs_closed_B0"]:.3e}; '
              f'sym={r["matrix_symmetry_residual"]:.3e}; '
              f'field={r["source_rephasing_covariance_residual"]:.3e}; '
              f'spurion={r["global_spurion_phase_residual"]:.3e}; '
              f'wrong y1={r["negative_controls"]["wrong_y1_conjugation"]:.3e}')
    print(f'Native mode={args.mode=="native"}; physical certifications 0/16; output={target}')
    return 0
if __name__=='__main__':
    try:sys.exit(main())
    except (ValueError,KeyError,TypeError,IndexError,ArithmeticError) as e:
        print('BLOCKED:',e,file=sys.stderr);sys.exit(2)
