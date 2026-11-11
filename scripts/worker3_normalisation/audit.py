"""Worker 3, independent read-only convention and Wilson scaling audit.

The tests here certify mathematics and consistency, not the installed Matchete
GammaCC -> canonical Weyl physical coefficient. No production file is modified.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
import numpy as np

OUTPUT = Path('output/audits/worker3_normalisation')
SOURCES = {
 'master': 'RGE/running/intermediate/GeneralEFT1C5Mixing.py',
 'general': 'RGE/general/WilsonTensorRGE.py',
 'polynomial': 'Numerical/running/GeneralT3WilsonNormalizationAudit.py',
 'b1': 'Numerical/lo_certification/physical_normalisation/B1_MATCHETE_WEYL_DERIVATION.md',
 'b1py': 'Numerical/lo_certification/physical_normalisation/b1_matchete_weyl.py',
 'tree': 'Lagrangian/RunGeneralT3FermionFirstEFT1Tree.wl',
 'hierarchy': 'Lagrangian/RunT3EFT1WilsonHierarchyNative.wl',
 'coupled': 'Numerical/orchestration/T3BIntegratedCoupledPipeline.py',
}
# Source blob hashes observed using the read-only GitHub connector, 2026-11-12.
BLOBS = {
 'master':'2ceafa9afd09375a96b2a54a3eec30468ebe93a5',
 'general':'bc02723d4ffd6176f4c8ac1b1132ae4376f55a0c',
 'polynomial':'3099ad2c321d049c6a757247ae7bb03e53cc4493',
 'b1':'187c310f26f81f465cb62da040159a54f2d1a69c',
 'b1py':'93c3f8cd833dcb8141b71e8fac7f9bbe0ef2850a',
 'tree':'76064f0ee1639402f252ee28583b9a7ce4ef3e44',
 'hierarchy':'9f736cc288cd7336efdfb19a4ce150592f17488f',
 'coupled':'ee7525cf482e8bdea0cc4ee2eb6a28ef8db3313a',
}
CONVENTIONS = [
 ('tree LLSS', 'Matchete tree EFT1', 'LLSS in fermion-first intermediate EFT, not C5', 'verified source semantics; physical conversion open', 'tree:4-6,1330-1392'),
 ('direct LLHH', 'Matchete one-loop/direct matching', 'separate operator and matching channel; cannot reuse LLSS factor', 'conversion unresolved', 'master:1-21; hierarchy:4-6'),
 ('Matchete PL', 'native GammaCC Proj[-1]', 'holomorphic PL; Proj[1] is barred conjugate channel', 'source verified, physical gamma phase open', 'polynomial:identify_chirality,verify_native_convention'),
 ('Matchete four-spinor to Weyl', 'GammaCC, CConj, PlusHc', 'chosen C=i gamma2 gamma0 proves conditional identities; installed Matchete C not measured', 'unresolved physical factor/sign', 'b1:sections 1-3,7'),
 ('Majorana free mass', '-(1/2) MF bar F F', 'canonical Weyl -(1/2) MF chi chi + h.c.; no extra factor two', 'chosen convention analytic; installed basis conditional', 'b1:section 2'),
 ('fermion-pair exchange', 'C_ijab symmetric in i,j', 'independent 2! combinatoric; not automatic Matchete multiplier', 'master internally verified', 'master:_read_wilson'),
 ('scalar-pair exchange', 'C_ijab symmetric in a,b', 'independent 2! combinatoric, real components', 'master internally verified', 'master:_read_wilson'),
 ('scalar potential', 'V=lambda_abcd phi_a phi_b phi_c phi_d / 4!', 'fully symmetric quartic; source portal CG normalization separate', 'master internally verified; physical CG open', 'master:1-12,_read_quartic'),
 ('charged and neutral Higgs', 'H+=(h+R+i h+I)/sqrt2; H0=(h0R+i h0I)/sqrt2', 'complex-to-real 1/sqrt2 per scalar insertion', 'master internally verified', 'master:_weinberg_template'),
 ('Weinberg', 'L5=(1/2)C5_pq(Lp.H)(Lq.H)+h.c.', 'C[nu_p,nu_q,h0R,h0R]=C5_pq/4; mnu=-v^2 C5/2', 'certified conditional on declared basis', 'master:417-454,558'),
 ('triplet weight basis', 'G=U.T@U', 'opposite T3 weight pairs and same CG phase in vertex', 'chosen basis analytic; native U unresolved', 'b1:section 4'),
 ('all-operator source factor', 'unknown map per operator a', 'must not equate LLSS, LLHH, and conjugate channels without derivation', 'unresolved', 'master:CERTIFICATES'),
 ('complex Wilson scaling', 'beta=A(C)+B(conj(C))', 'real scalar r commutes; complex z generally fails unless B=0 or z real', 'analytic', 'master:conjugate_coefficient_yukawa'),
]

def complex_matrix(x):
    if isinstance(x,dict) and set(x)>={'real','imag'}:
        z=np.asarray(x['real'],float)+1j*np.asarray(x['imag'],float)
    else:z=np.asarray(x,complex)
    if z.shape!=(3,3) or not np.isfinite(z).all():raise ValueError('Expected finite 3x3 complex matrix')
    return z

def encode(x):
    a=np.asarray(x,complex)
    return {'real':a.real.tolist(),'imag':a.imag.tolist()}

def obs(m):
    """Charged-lepton diagonal basis, nondegenerate masses; eigenvector ordering by mass.

    Eigenvectors of m^dagger m are U of Takagi factorisation up to conjugation;
    magnitudes give identical PMNS angles. CP phases are deliberately not inferred.
    """
    h=m.conj().T@m
    eig,u=np.linalg.eigh(h)
    masses=np.sqrt(np.clip(eig,0,None))
    if len(set(np.round(masses / max(np.max(masses),1e-100),9)))!=3:
        raise ValueError('Degenerate Takagi spectrum; mixing not unique')
    a=np.abs(u)
    s13=a[0,2]**2
    s12=a[0,1]**2/(1-s13)
    s23=a[1,2]**2/(1-s13)
    return {'masses':masses.tolist(),'dm21':float(masses[1]**2-masses[0]**2),
      'dm31':float(masses[2]**2-masses[0]**2),
      'ratios':[float(masses[0]/masses[2]),float(masses[1]/masses[2])],
      'sin2_angles':[float(s12),float(s13),float(s23)],
      'pmns_abs':a.tolist()}

def invariant_checks(m, r, z):
    base=obs(m); real=obs(r*m); complex_o=obs(z*m)
    def maximum(a,b):return float(np.max(np.abs(np.asarray(a)-np.asarray(b))))
    return {'real_factor':r,'complex_factor':encode(np.array([z])),
      'mass_real_error':maximum(real['masses'],np.asarray(base['masses'])*abs(r)),
      'dm_real_error':maximum([real['dm21'],real['dm31']],np.array([base['dm21'],base['dm31']])*r*r),
      'mass_ratio_real_error':maximum(real['ratios'],base['ratios']),
      'angle_real_error':maximum(real['sin2_angles'],base['sin2_angles']),
      'pmns_abs_real_error':maximum(real['pmns_abs'],base['pmns_abs']),
      'angle_complex_global_error':maximum(complex_o['sin2_angles'],base['sin2_angles']),
      'mass_complex_global_error':maximum(complex_o['masses'],np.asarray(base['masses'])*abs(z)),
      'base':base}

def synthetic_matrix(seed):
    rng=np.random.default_rng(seed)
    t=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
    return (t+t.T)/2

def phase_obstruction(c, a, b, z):
    """Real-linear beta of C: A*C+B*conj(C); A,B act elementwise as illustration."""
    beta=lambda x:a*x+b*x.conj()
    mismatch=beta(z*c)-z*beta(c)
    expected=b*(z.conjugate()-z)*c.conj()
    np.testing.assert_allclose(mismatch,expected,atol=1e-12)
    return float(np.linalg.norm(mismatch)),float(np.linalg.norm(beta(c)))

def complex_channel_checks(a,b,r):
    full=a+b
    alternative=r*a+b
    o=obs(full); p=obs(alternative)
    return {'relative_channel_factor':r,
     'frobenius_shape_shift':float(np.linalg.norm(alternative/np.linalg.norm(alternative)-full/np.linalg.norm(full))),
     'angles_shift':float(np.max(np.abs(np.array(o['sin2_angles'])-p['sin2_angles']))),
     'mass_ratio_shift':float(np.max(np.abs(np.array(o['ratios'])-p['ratios'])))}

def stage_ratio_check(start,end,r):
    """Ratios of any homogeneous mass norm at two stages, not entrywise if zero."""
    n=lambda a:np.linalg.norm(a)
    if n(start)==0 or n(end)==0:raise ValueError('Zero stage mass norm')
    initial=n(end)/n(start)
    scaled=n(r*end)/n(r*start)
    return {'baseline':float(initial),'scaled':float(scaled),'error':float(abs(initial-scaled))}

def source_status(root):
    data={}
    for name,rel in SOURCES.items():
        p=root/rel
        data[name]={'path':rel,'github_blob_sha':BLOBS[name], 'local_exists':p.is_file()}
        if p.is_file():
            b=p.read_bytes();data[name]['local_sha256']=hashlib.sha256(b).hexdigest()
    # Check the exact repo contracts only when files are locally available.
    p=root/SOURCES['master']
    if p.is_file():
        s=p.read_text(encoding='utf-8')
        for needle in ['conjugate_coefficient_yukawa','def project_c5','matrix = 4 *', '0.25 * c5[p, q]','synthetic']:
            if needle not in s: raise RuntimeError('Source convention drift: '+needle)
    return data

def run(root,out,seed=3403, fixtures=None):
    out.mkdir(parents=True,exist_ok=True)
    status=source_status(root)
    with (out/'convention_table.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f);w.writerow(['channel','source_convention','derived_factor_or_warning','status','reference']);w.writerows(CONVENTIONS)
    rng=np.random.default_rng(seed)
    records=[]
    cases=[('synthetic_'+str(i),synthetic_matrix(seed+i),synthetic_matrix(seed+30+i),synthetic_matrix(seed+60+i)) for i in range(32)]
    # Optional user-supplied or production-exported complex 3x3 matrices, not silently interpreted.
    if fixtures:
        j=json.loads(Path(fixtures).read_text(encoding='utf-8'))
        if j.get('format')!='worker3-c5-trajectories-v1':raise ValueError('Wrong fixture format')
        for item in j['cases']:
            cases.append((item['name'],complex_matrix(item['C5_start']),complex_matrix(item['channel_A']),complex_matrix(item['channel_B'])))
    for name,m,a,b in cases:
        result=invariant_checks(m,2.7,1.3+0.7j)
        channel=complex_channel_checks(a,b,1.9)
        stage=stage_ratio_check(m,m+0.02*a,2.7)
        obstruction,norm=phase_obstruction(m,0.17-0.3j,0.13+0.2j,1.3+0.7j)
        real_obstruction,_=phase_obstruction(m,0.17-0.3j,0.13+0.2j,2.7)
        assert result['mass_real_error'] < 1e-10
        assert result['dm_real_error'] < 1e-9
        assert result['angle_real_error'] < 1e-10
        assert result['angle_complex_global_error'] < 1e-10
        assert result['mass_ratio_real_error'] < 1e-10
        assert stage['error'] < 1e-12
        assert real_obstruction < 1e-12
        assert obstruction > 1e-8
        records.append({'name':name,'source':'fixture' if name.startswith('fixture_') else 'synthetic_or_fixture',
                        'invariants':{k:v for k,v in result.items() if k!='base'},
                        'channel_sensitivity':channel,'stage':stage,
                        'complex_beta_obstruction':obstruction,'complex_beta_norm':norm})
    summary={'format':'worker3-normalisation-audit-v1','source_status':status,
       'cases_passed':len(records),'fixture_path':str(fixtures) if fixtures else None,
       'physical_matchete_to_weyl_conversion_certified':False,
       'source_to_master_wilson_conversion_common_factor_certified':False,
       'tests_do_not_claim_full_coupled_ivp_replay':True,
       'complex_rescaling_in_beta_generally_invalid':True,
       'notes':'Synthetic matrix tests verify scaling identities, not a source-to-Weyl coefficient.'}
    (out/'regression_cases.json').write_text(json.dumps(records,indent=2),encoding='utf-8')
    (out/'source_manifest.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    return summary

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--root',default='.')
    ap.add_argument('--output',default=str(OUTPUT))
    ap.add_argument('--fixtures')
    a=ap.parse_args()
    root=Path(a.root).resolve();out=Path(a.output)
    report=run(root,out,fixtures=a.fixtures)
    print(json.dumps({'cases_passed':report['cases_passed'],'report_dir':str(out),'physical_certificate':False},indent=2))
if __name__=='__main__':main()
