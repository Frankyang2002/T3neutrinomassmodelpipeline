#!/usr/bin/env python3
"""Independent B_-1 Majorana UV self-energy from a LIVE native Matchete witness.

Does not import production core.py, frozen provisional tensors or their B0.
No EFT/matching/running. Requires fresh output directory; refuses overwrites.

Run: python scripts/ReviewPhysicalBm1NativeLoop.py --repo-root . \
         --native-json output/.../B_-1_native_action.json \
         --output-dir output/general_observable/physical_normalisation/bm1_.../commonMF
"""
from __future__ import annotations
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys
import mpmath as mp
import numpy as np

FROZEN = 'Numerical/benchmarks/fitted_reference/T3_dS1_2_dS2_2_dF_1_alpha_m1.json'
PROVISIONAL = 'output/general_observable/lo_direct_one_loop/provisional/B_-1_direct_one_loop.json'
OUTPUT_BASE = 'output/general_observable/physical_normalisation'
EXPECTED_FORMAT = 't3-b1-bm1-live-native-scalar-quadratic-v1'


def require(ok: bool, msg: str) -> None:
    if not ok:
        raise ValueError(msg)


def read_json(path: Path) -> dict:
    require(path.is_file(), f'Missing source: {path}')
    value = json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(value, dict), f'Not JSON object: {path}')
    return value


def hash_file(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def decode(obj: object) -> np.ndarray:
    if isinstance(obj, dict) and set(obj) == {'real', 'imag'}:
        a = np.asarray(obj['real'], dtype=float) + 1j*np.asarray(obj['imag'], dtype=float)
    else:
        a = np.asarray(obj, dtype=complex)
    require(np.isfinite(a.real).all() and np.isfinite(a.imag).all(), 'Nonfinite tensor')
    return a


def native_tensor(native: dict, key: str, shape: tuple[int, ...]) -> np.ndarray:
    row = native['native_tensors'][key]
    require(tuple(row['shape']) == shape, f'{key}: shape mismatch')
    arr = np.asarray(row['real_flat'], dtype=float) + 1j*np.asarray(row['imag_flat'], dtype=float)
    require(arr.size == int(np.prod(shape)) and np.isfinite(arr).all(), f'{key}: bad data')
    return arr.reshape(shape)


def native_project(native: dict, key: str, shape: tuple[int, ...]) -> np.ndarray:
    row = native['native_H0_projection'][key]
    require(tuple(row['shape']) == shape, f'{key}: projected shape mismatch')
    return (np.asarray(row['real_flat'],float)+1j*np.asarray(row['imag_flat'],float)).reshape(shape)


def verify_source(native: dict, root: Path) -> tuple[complex, complex, complex, dict]:
    require(native.get('format') == EXPECTED_FORMAT and native.get('model') == 'B_-1',
            'Wrong native witness format/model')
    require(native.get('status') == 'NativeBm1ScalarAndMajoranaActionExtracted_PhysicalDictionaryPending'
            and native.get('native_lagrangian_validated') is True,
            'Live Matchete validation missing')
    require(native.get('physical_source_certified') is False, 'Witness overclaims certification')
    checks=native.get('native_checks',{})
    require(len(checks) >= 8 and all(v is True for v in checks.values()), 'Failed native CG/action checks')
    require(native['native_mass_term']['coefficients'] == ['-1/2'], 'Majorana mass 1/2 missing')
    for name in ('S1Inv1','S1Inv2','S2Inv1','S2Inv2','mixing','mixing_hc','y1','y2','y1_hc','y2_hc'):
        require(native['native_action_terms'][name]['coefficients'] == ['1'],
                f'{name}: unexpected source interaction coefficient')
    require(all(native['restricted_environment']['effective'][key]=='0' for key in
            ('T3_EXPERIMENTAL_GENERAL_FERMIONS','T3_EXPERIMENTAL_GENERAL_QUADRATICS',
             'T3_EXPERIMENTAL_GENERAL_QUARTICS','T3_EXPERIMENTAL_FLAVOR_MASSES')),
            'Nonrestricted source environment')
    for rel,expected in native['source_sha256'].items():
        require(hash_file(root/rel).lower() == expected.lower(), f'Native source changed after extraction: {rel}')
    y1=native_tensor(native,'T3Y1CG',(2,2))
    y2=native_tensor(native,'T3Y2CG',(2,2))
    mix=native_tensor(native,'T3MixCG',(2,2,2,2))
    p11=native_tensor(native,'T3HiggsPortal1CGInv1',(2,2,2,2))
    p21=native_tensor(native,'T3HiggsPortal2CGInv1',(2,2,2,2))
    pr1=native_project(native,'S1Inv1',(2,2))
    pr2=native_project(native,'S2Inv1',(2,2))
    require(np.allclose(pr1,np.eye(2),rtol=0,atol=2e-13) and
            np.allclose(pr2,np.eye(2),rtol=0,atol=2e-13), 'Identity portals are not native Inv1')
    require(np.allclose(p11[1,1],pr1,rtol=0,atol=2e-13) and
            np.allclose(p21[1,1],pr2,rtol=0,atol=2e-13), 'Inconsistent native portal projection')
    mu=mix[1,1]
    require(np.count_nonzero(abs(mu)>1e-12)==1, 'More than one Higgs-neutral mixing channel')
    c1=complex(y1[0,0]);c2=complex(y2[0,1]);cm=complex(mu[0,1])
    require(abs(c1-1)<2e-13 and abs(c2-1)<2e-13 and
            abs(cm+2/np.sqrt(3))<2e-13, 'Unexpected native neutral component values')
    require(abs(native_project(native,'mixing_S1_S2bar',(2,2))[0,1]-cm)<2e-13,
            'Inconsistent source mixing projection')
    return c1,c2,cm,{'y1':c1,'y2':c2,'mix':cm}


def mp_complex(z: complex) -> mp.mpc:
    z=complex(z)
    return mp.mpc(str(z.real),str(z.imag))


def finite_divided_difference(mf2: mp.mpf, small: mp.mpf, large: mp.mpf) -> mp.mpf:
    """Finite scalar propagator divided difference, no production B0 used."""
    delta=large-small
    if delta==0:
        return -mp.quad(lambda x: (1-x)/(x*mf2+(1-x)*small),[0,1])
    return -mp.quad(lambda x: mp.log1p((1-x)*delta/(x*mf2+(1-x)*small))/delta,[0,1])


def evaluate(source: dict, coeff: tuple[complex,complex,complex], *,
             portal: str, masses: str) -> dict:
    mp.mp.dps=60
    ordinary=source['base_state']['ordinary']
    v=mp.mpf(str(source['numerical']['vev_gev']))
    lam=complex(ordinary['lambdaT3']['real'],ordinary['lambdaT3']['imag'])
    m1=mp.mpf(str(ordinary['mS1Sq']))
    m2=mp.mpf(str(ordinary['mS2Sq']))
    if portal=='benchmark-identity':
        m1-=mp.mpf(str(ordinary['lambdaH1']))*v*v/2
        m2-=mp.mpf(str(ordinary['lambdaH2']))*v*v/2
    else:
        require(portal=='zero', 'Unknown portal assignment')
    # Ensure no additional source terms are silently treated as zero.
    for key in ('lambdaH1Adj','lambdaH2Adj','lambda12Adj',
                'lambdaHHdagS2S2','lambdaHHdagS1barS1bar',
                'lambdaHHdagS1barS2barCross'):
        val=ordinary.get(key,0)
        if isinstance(val,dict):val=complex(val['real'],val['imag'])
        require(abs(complex(val)) < 1e-16, f'Nonzero extra scalar coupling {key}: unsupported')
    mf=decode(ordinary['MF'])
    require(mf.shape==(3,3) and np.linalg.norm(mf-np.diag(np.diag(mf)))<1e-12,
            'Native source requires diagonal heavy masses')
    raw=np.diag(mf)
    require(np.max(abs(raw.imag))<1e-11 and (raw.real>0).all(), 'Nonpositive/complex masses')
    if masses=='native-common':
        heavy=np.repeat(float(raw.real[0]),3)
    else:
        require(masses=='benchmark-diagonal','Unknown fermion mass assumption')
        heavy=raw.real
    c1,c2,cm=coeff
    A21=-(mp_complex(lam)*v*v/2)*mp_complex(cm)
    diff=mp.sqrt((m1-m2)**2+4*abs(A21)**2)
    smin=(m1+m2-diff)/2;smax=(m1+m2+diff)/2
    require(smin>0 and smax>0, 'Nonpositive scalar eigenmass squared')
    y1=decode(ordinary['y1']);y2=decode(ordinary['y2'])
    require(y1.shape==(3,3) and y2.shape==(3,3),'Wrong Yukawa flavour shape')
    g1=-np.conjugate(y1*c1);g2=-np.conjugate(y2*c2)
    mass_eV=np.zeros((3,3),complex)
    loop_per_flavour=[]
    for r,Mi in enumerate(heavy):
        mm=mp.mpf(str(Mi))
        D=finite_divided_difference(mm*mm,smin,smax)
        pref=complex(mm*A21*D/(16*mp.pi**2)*mp.mpf('1e9'))
        mass_eV+=pref*(np.outer(g1[:,r],g2[:,r])+np.outer(g2[:,r],g1[:,r]))
        loop_per_flavour.append({'heavy_M_gev':float(mm),'divided_difference_gev_minus2':float(D),
                                 'prefactor_eV_before_Yukawas':{'real':pref.real,'imag':pref.imag}})
    require(np.linalg.norm(mass_eV-mass_eV.T)<1e-12*max(np.linalg.norm(mass_eV),1e-100),
            'Unexpected nonsymmetric Majorana mass matrix')
    masses_eV=np.linalg.svd(mass_eV,compute_uv=False)
    diag=[m1,m2]
    A=np.diag([float(m1),float(m1),float(m2),float(m2)]).astype(complex)
    A[3,0]=complex(A21);A[0,3]=complex(A21).conjugate()
    return {
      'source_point':{'portal_assignment':portal,'fermion_mass_assignment':masses,
           'source_F_Matchete_one_symbol':'MF',
           'three_distinct_frozen_flavour_masses_directly_native':False},
      'physical_neutral_channel':{'y1_CG':{'real':c1.real,'imag':c1.imag},
           'y2_CG':{'real':c2.real,'imag':c2.imag},
           'mix_CG':{'real':cm.real,'imag':cm.imag},
           'neutral_scalar_A21_gev2':{'real':complex(A21).real,'imag':complex(A21).imag}},
      'scalar_A_gev2':{'real':A.real.tolist(),'imag':A.imag.tolist()},
      'neutral_scalar_eigenmass_squared_gev2':[float(smin),float(smax)],
      'loop_terms':loop_per_flavour,
      'neutrino_mass_matrix_eV':{'real':mass_eV.real.tolist(),'imag':mass_eV.imag.tolist()},
      'takagi_singular_values_eV_descending':masses_eV.tolist(),
      'neutrino_matrix_symmetric':True,
      'integral':'-int_0^1 dx log(1+(1-x)*(splus-sminus)/(x*M^2+(1-x)*sminus))/(splus-sminus)',
      'uv_selfenergy':'sum_r (M_r*A21*D_r/(16*pi^2)) * (g1_ar*g2_br+g2_ar*g1_br)',
      'additional_majorana_loop_multiplier':1,
      'source_real_C_matrix_phase_measurement':False,
      'physical_source_certified':False,
    }


def main() -> int:
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--repo-root',type=Path,default=Path(__file__).resolve().parents[1])
    ap.add_argument('--native-json',type=Path,required=True)
    ap.add_argument('--output-dir',type=Path,required=True)
    ap.add_argument('--portal-assignment',choices=['zero','benchmark-identity'],default='zero')
    ap.add_argument('--fermion-masses',choices=['native-common','benchmark-diagonal'],default='native-common')
    ap.add_argument('--compare-frozen',action='store_true')
    args=ap.parse_args()
    root=args.repo_root.resolve()
    path=lambda p: p.resolve() if p.is_absolute() else (root/p).resolve()
    out=path(args.output_dir)
    try:
        require(out.is_relative_to((root/OUTPUT_BASE).resolve()), 'Output must be in owned evidence directory')
        require(not out.exists(),'Refusing to overwrite existing evidence directory')
        if args.compare_frozen:
            require(args.portal_assignment=='benchmark-identity' and
                    args.fermion_masses=='benchmark-diagonal',
                    'Frozen comparison only at the EXACT frozen portal and three-mass point')
        native_path=path(args.native_json)
        native=read_json(native_path)
        c1,c2,cm,_=verify_source(native,root)
        benchmark_path=root/FROZEN
        benchmark=read_json(benchmark_path)
        rep=benchmark.get('representation',{})
        require(tuple(rep.get(x) for x in ('d_s1','d_s2','d_f','alpha')) == (2,2,1,-1),
                'Frozen benchmark is not B_-1')
        result=evaluate(benchmark,(c1,c2,cm),portal=args.portal_assignment,masses=args.fermion_masses)
        print('B_-1 independent LIVE-NATIVE-source scalar amplitude (eV):')
        print(np.array(result['neutrino_mass_matrix_eV']['real'])+
              1j*np.array(result['neutrino_mass_matrix_eV']['imag']))
        print('Takagi masses descending (eV):',result['takagi_singular_values_eV_descending'])
        compare=None
        if args.compare_frozen:  # deliberately read only AFTER independent evaluation
            provisional_path=root/PROVISIONAL
            reference=read_json(provisional_path)
            ref=decode(reference['total']['complex_symmetric_matrix_gev'])*1e9
            current=decode(result['neutrino_mass_matrix_eV'])
            residual=float(np.linalg.norm(ref-current)/max(np.linalg.norm(ref),1e-100))
            compare={'comparison_is_conditional':True,'frozen_path':PROVISIONAL,
                     'frozen_sha256':hash_file(provisional_path),
                     'relative_frobenius_matrix_residual':residual}
            print('CONDITIONAL frozen matrix residual:',residual)
        report={'format':'t3-b1-bm1-independent-native-loop-v1',
            'status':'IndependentMajoranaLoopEvaluated_NoFullPhysicalSourceCertificate',
            'model':'B_-1','method':'direct Feynman parameter integral, not production B0',
            'native_witness':{'path':str(native_path.relative_to(root)) if native_path.is_relative_to(root) else str(native_path),
                              'sha256':hash_file(native_path),
                              'native_source_git_head':native.get('source_git_head','')},
            'benchmark_input':{'path':FROZEN,'sha256':hash_file(benchmark_path)},
            'independent':result,'conditional_comparison':compare,
            'source_assumptions':['single MF native declaration; benchmark 3 diagonal masses requires an extension',
                                  'only explicit identity Higgs portals are activated in the comparison',
                                  'installed numeric GammaCC phase not measured; native action Majorana 1/2 measured'],
            'physical_certified_model_count_increment':0,'eft_or_rge_performed':False}
        out.mkdir(parents=True,exist_ok=False)
        output=out/'B_-1_independent_native_loop.json'
        output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
        print('Evidence:',output)
        print('Physical source-certified models: 0/16')
        return 0
    except (ValueError,KeyError,TypeError,OSError,OverflowError) as exc:
        print('BLOCKED: '+str(exc),file=sys.stderr)
        return 2

if __name__=='__main__':
    raise SystemExit(main())
