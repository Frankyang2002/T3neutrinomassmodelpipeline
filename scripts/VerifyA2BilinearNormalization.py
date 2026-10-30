"""A(-2) Matchete/RGBeta self-bilinear mass-convention and loop-kernel audit.

Reads the successful 16-model symbolic export. Obtains a NEW native Matchete
quadratic tensor in its native triplet basis. Checks the bilinear coefficient
against the independent RGBeta spherical-basis real-mass decomposition.
Separately tests the exact one-scale Matchete C5 mass dependence against a
three-propagator Euclidean one-loop integral, with no fitted loop factor.

This verifies a restricted scalar mass-convention mapping and kernel SHAPE;
it is NOT a complete physical one-loop matching or a certified full RG map.

PowerShell:
  python -m scripts.VerifyA2BilinearNormalization --symbolic-run-id all16_v1 --run-id a2_norm_v1
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
import os
from pathlib import Path
import re
import subprocess

import mpmath as mp
import numpy as np
import sympy as sp

from Numerical.experimental.T3SelfBilinearMassBasis import (
    SelfBilinearMassBasis, invariant_metric, parse_matchete_constant,
)
from RGE.matching.MatcheteC5Parsing import parse_matchete_c5
from RGE.matching.WeinbergFlavorMatching import extract_one_generation_loop_kernel

ROOT = Path(__file__).resolve().parents[1]
SYM_ROOT = ROOT / 'output/studies/general_scalar_extended/complete_scalar_symbolic'
OUTPUT = ROOT / 'output/studies/general_scalar_extended/a2_bilinear_normalization'
NATIVE = ROOT / 'Lagrangian/RunA2NativeBilinearTensor.wl'
MASS = 'muT3GeneralS2Sq'
MODEL = 'A_-2'
FORMAT = 't3-a-minus2-bilinear-convention-and-loop-shape-v1'


def require(test: bool, text: str) -> None:
    if not test:
        raise ValueError('A(-2) bilinear verification: ' + text)


def read(path: Path) -> dict:
    require(path.is_file(), 'missing file ' + str(path))
    data = json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(data, dict), 'expected JSON object '+str(path))
    return data


def sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def verify_symbolic_export(folder: Path) -> dict:
    status = read(folder / 'status.json')
    lag = read(folder / 'lagrangian.json')
    require(status.get('status') == 'SymbolicRGExportedMatchingUnverified' and
            not status.get('errors'), 'symbolic stage not successful')
    require(lag.get('status') == 'Success' and
            lag.get('quadratic_included') == 1 and
            lag.get('quadratic_expected') == 1 and
            lag.get('matching_performed') is False,
            'full A(-2) Lagrangian lacks the extra bilinear')
    result = {}
    for stage in ('uv', 'eft1'):
        p = folder / (stage + '_rgbeta.json')
        data = read(p)
        meta = data.get('metadata', {})
        betas = data.get('report_betas', {})
        require(data.get('status') == 'SuccessCandidateMassBeta' and
                data.get('physical_mapping_certified') is False and
                meta.get('GeneralScalarModel') == MODEL and
                meta.get('GeneralScalarStage') == stage.upper() and
                meta.get('QuadraticKind') == 'self_s2' and
                meta.get('QuadraticMassSymbol') == MASS and
                isinstance(betas, dict) and
                {'mS1Sq','mS2Sq',MASS} <= set(betas) and
                not any('BetaTerm[' in x for x in betas.values()),
                stage + ' does not have complete candidate mass betas')
        require(MASS in betas[MASS],
                stage + ': holomorphic mass beta unexpectedly independent of its coupling')
        result[stage] = {'beta_count': len(betas), 'mass_beta': betas[MASS],
                         'sha256': sha(p)}
    return result


def parse_native_metric(native: dict) -> complex:
    require(native.get('format')=='t3-a-minus2-native-bilinear-tensor-v1' and
            native.get('model')==MODEL and
            native.get('status')=='NativeMatcheteBilinearMetricExtracted' and
            native.get('native_metric_is_factor_times_identity') is True and
            native.get('mass_operator_in_full_lagrangian') is True and
            native.get('species')==['S2dag','S2dag'],
            'unexpected native Matchete tensor certificate')
    raw=native.get('native_metric_inputform')
    require(isinstance(raw,list) and len(raw)==3 and
            all(isinstance(row,list) and len(row)==3 for row in raw),
            'native metric shape must be 3x3')
    m=np.array([[parse_matchete_constant(v) for v in row] for row in raw],dtype=complex)
    k=parse_matchete_constant(native['native_metric_factor_inputform'])
    require(abs(k)>1e-14 and np.allclose(m,k*np.eye(3),rtol=2e-13,atol=2e-14),
            'native CG normalization disagrees with actual component tensor')
    return k


def validate_source_conventions() -> dict:
    mat = ROOT/'Lagrangian/interactions/GeneralScalarOperators.wl'
    rgb = ROOT/'RGE/running/rgbeta/RunGeneralT3QuadraticRGBeta.wl'
    a,b=mat.read_text(encoding='utf-8'),rgb.read_text(encoding='utf-8')
    # The UV interaction and RGBeta scalar potential have opposite overall sign.
    require('"Expression" -> PlusHc[contribution]' in a and
            'contribution = coupling[] product' in a and
            re.search(r'AddScalarMass\[muT3GeneralS2Sq,\{Bar @ S2,Bar @ S2\}',b),
            'source sign/field-word contracts changed')
    require('V2=mu2*(Bar[S])^2+h.c.' in b,
            'RGBeta declared scalar potential convention changed')
    return {'matchete_extra_bilinear_is_plus_in_lagrangian':True,
            'rgbeta_bilinear_is_plus_in_potential':True,
            'relative_sign_if_same_canonical_field_basis':-1,
            'matchete_source_sha256':sha(mat),'rgbeta_source_sha256':sha(rgb)}


def verify_mass_map(k: complex, m2: float, mu_rgb: complex) -> dict:
    """Match V_M=-mu_M*k*(Sdag)^2+h.c. to V_R=mu_R*C*(Sdag)^2+h.c.

    The SU(2) metric basis W obeys W.T C W = I, independently of any
    Matchete-to-RGBeta complete vertex intertwiners.
    """
    require(m2>0 and np.isfinite(m2) and np.isfinite([mu_rgb.real,mu_rgb.imag]).all(),
            'invalid reference scalar mass')
    mu_mat = -mu_rgb/k
    basis=SelfBilinearMassBasis.build('S2',3,m2,mu_rgb)
    require(not basis.tachyonic,'reference mass eigenstate tachyonic')
    C=invariant_metric(3)
    rng=np.random.default_rng(2602)
    max_relative=0.0
    for _ in range(20):
        r,i=rng.standard_normal((2,3))
        scalar=basis.fields(r,i)
        # In Matchete Cartesian coordinates, the field is W^dag S up to the
        # fixed intertwiners. This is a metric test, not a gauge-vertex map.
        native = basis.intertwiner.conj().T @ scalar
        native_v=(m2*np.vdot(native,native).real -
                  2*np.real(mu_mat*k*np.vdot(native, native.conjugate())))
        rgb_v=(m2*np.vdot(scalar,scalar).real+
               2*np.real(mu_rgb*(scalar.conj().T@C@scalar.conj())))
        target=basis.potential_mass_basis(r,i)
        err=max(abs(native_v-rgb_v),abs(rgb_v-target))/max(1,abs(target))
        max_relative=max(max_relative,float(err))
        require(err<2e-12,'Lagrangian/potential bilinear tensor map mismatch')
    return {'matchete_bilinear_cg_factor': {'real':k.real,'imag':k.imag},
            'mass_map_in_matchete_native_metric_convention':
                'mu_RGBeta = -k_native * mu_Matchete',
            'matchete_mass_parameter_gev2':{'real':mu_mat.real,'imag':mu_mat.imag},
            'rgbeta_mass_parameter_gev2':{'real':mu_rgb.real,'imag':mu_rgb.imag},
            'masses_squared_gev2':list(basis.eigen_masses_squared),
            'mass_splitting_squared_gev2':4*abs(mu_rgb),
            'potential_max_relative_error':max_relative,
            'spherical_triplet_metric_matches_declared_real_basis':True,
            'full_matchete_rgbeta_gauge_vertex_map_certified':False}


def integral_triangle(mf: float, ms1: float, ms2: float) -> mp.mpf:
    """Four-dimensional Euclidean scalar triangle ∫_0∞ t dt/∏(t+m_i²).

    No universal numerical prefactor: this is compared by ratios only.
    """
    with mp.workdps(48):
        s=mp.mpf(str(mf))**2
        b=(mp.mpf(str(ms1))/mp.mpf(str(mf)))**2
        c=(mp.mpf(str(ms2))/mp.mpf(str(mf)))**2
        f=lambda t: t/((t+1)*(t+b)*(t+c))
        return +(mp.quad(f,[0,mp.mpf('0.001'),mp.mpf('0.1'),1,10,mp.inf])/s)


def check_one_scale_loop_shape(source: Path, mf: float, ms1: float,
                               masssq_plus: float, masssq_minus: float) -> dict:
    """An independent mass-ratio check using the saved Matchete source itself."""
    require(source.is_file(),'missing one-scale Matchete C5 coefficient '+str(source))
    raw=source.read_text(encoding='utf-8-sig')
    require(all(token in raw for token in ('Coupling[MF,','Coupling[MS1,',
            'Coupling[MS2,','Coupling[lambdaT3,','Coupling[y1,','Coupling[y2,')),
            'incorrect one-scale source')
    expr=parse_matchete_c5(raw)
    k=extract_one_generation_loop_kernel(expr)
    MF,MS1,MS2,lam,hbar=(sp.Symbol(n) for n in
                           ('MF','MS1','MS2','lambdaT3','hbar'))
    require(k.free_symbols <= {MF,MS1,MS2,lam,hbar},
            'unsupported Matchete loop kernel symbols '+str(k.free_symbols))
    f=sp.lambdify((MF,MS1,MS2,lam,hbar),k,'mpmath')
    with mp.workdps(45):
        def eval_source(ms2):
            return mp.mpc(f(mp.mpf(str(mf)),mp.mpf(str(ms1)),
                            mp.mpf(str(ms2)),1,1))
        vals=[eval_source(math.sqrt(masssq_plus)),
              eval_source(math.sqrt(masssq_minus)),
              eval_source(math.sqrt((masssq_plus+masssq_minus)/2))]
        ints=[integral_triangle(mf,ms1,math.sqrt(masssq_plus)),
              integral_triangle(mf,ms1,math.sqrt(masssq_minus)),
              integral_triangle(mf,ms1,math.sqrt((masssq_plus+masssq_minus)/2))]
        require(all(abs(z)>mp.mpf('1e-80') for z in vals) and
                all(abs(z)>mp.mpf('1e-80') for z in ints),
                'zero source or independent loop integral')
        ratios=[vals[i]/vals[2] for i in (0,1)]
        exact=[ints[i]/ints[2] for i in (0,1)]
        errors=[float(abs(ratios[i]-exact[i])) for i in (0,1)]
        require(max(errors)<2e-9,
                'Matchete one-loop source scalar-mass response does not match independent loop integral')
        return {'status':'OneScaleLoopMassDependenceRatioChecked',
                'mass_response_ratios_matchete':[
                    {'real':float(z.real),'imag':float(z.imag)} for z in ratios],
                'mass_response_ratios_integral':[float(z) for z in exact],
                'maximum_relative_ratio_error':max(errors),
                'unit_normal_kernel_over_center':{'real':float(((ratios[0]+ratios[1])/2).real),
                                                   'imag':float(((ratios[0]+ratios[1])/2).imag)},
                'unit_anomalous_spectral_difference_over_center':{
                    'real':float(((ratios[0]-ratios[1])/2).real),
                    'imag':float(((ratios[0]-ratios[1])/2).imag)},
                'source_sha256':sha(source),
                'full_one_loop_matchete_mass_eigen_matching_performed':False}


def main() -> int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--symbolic-run-id',default='all16_v1')
    p.add_argument('--run-id',required=True)
    p.add_argument('--output-root',type=Path,default=OUTPUT)
    p.add_argument('--wolframscript',default='wolframscript')
    p.add_argument('--mu2-fraction',type=float,default=.01)
    p.add_argument('--mu2-phase',type=float,default=.35)
    args=p.parse_args()
    require(bool(re.fullmatch(r'[a-zA-Z0-9][a-zA-Z0-9_-]{0,63}',args.run_id)),
            'unsafe run id')
    require(math.isfinite(args.mu2_fraction) and 0<args.mu2_fraction<.45 and
            math.isfinite(args.mu2_phase), 'invalid mass-splitting parameters')
    out=args.output_root/args.run_id
    require(not out.exists(),'output already exists: use a fresh --run-id')
    folder=SYM_ROOT/args.symbolic_run_id/'A_m2'
    beta=verify_symbolic_export(folder)
    conventions=validate_source_conventions()
    config=read(ROOT/'Numerical/benchmarks/fitted_reference/T3_dS1_1_dS2_3_dF_2_alpha_m2.json')
    ordinary=config['base_state']['ordinary']
    m2=float(ordinary['mS2Sq'])
    mu_rgb=complex(args.mu2_fraction*m2*math.cos(args.mu2_phase),
                   args.mu2_fraction*m2*math.sin(args.mu2_phase))
    source=ROOT/'output/hypercharge/T3_A_alpha_m2/data/c5_coefficient.txt'
    # Before invoking Wolfram, check loop source is nontrivial and mass sensitive.
    loop=check_one_scale_loop_shape(source,1e5,math.sqrt(float(ordinary['mS1Sq'])),
                                    m2+2*abs(mu_rgb),m2-2*abs(mu_rgb))
    require(NATIVE.is_file(),'missing native tensor Wolfram helper')
    out.mkdir(parents=True)
    native_path=out/'native_matchete_bilinear.json'
    env=os.environ.copy()
    env['T3_EXPERIMENTAL_GENERAL_QUADRATICS']='1'
    env['T3_EXPERIMENTAL_GENERAL_QUARTICS']='1'
    proc=subprocess.run([args.wolframscript,'-file',str(NATIVE),str(native_path)],
                        cwd=ROOT,env=env,capture_output=True,text=True,check=False)
    (out/'native_matchete.log').write_text(proc.stdout+'\n'+proc.stderr,encoding='utf-8')
    require(proc.returncode==0 and native_path.is_file(),
            'Matchete native CG extraction failed: '+proc.stdout[-1200:]+proc.stderr[-1200:])
    native=read(native_path)
    k=parse_native_metric(native)
    mass=verify_mass_map(k,m2,mu_rgb)
    report={'format':FORMAT,
            'status':'A2NativeBilinearMassMapAndOneScaleKernelShapeChecked',
            'symbolic_run_id':args.symbolic_run_id,
            'model':MODEL,
            'input_reference_fitted_mass_squared_gev2':m2,
            'beta_exports':beta,
            'source_conventions':conventions,
            'mass_basis':mass,
            'mass_dependent_one_scale_kernel':loop,
            'native_matchete_certificate_sha256':sha(native_path),
            'native_tensor_source':native_path.as_posix(),
            'scope':'Native Matchete self-pair metric, declared RGBeta spherical mass convention, source-level opposite Lagrangian/potential sign, and one-scale scalar-mass response. No complete physical CG intertwiners or matched bilinear C5 derivation.',
            'physical_matchete_rgbeta_entire_vertex_map_certified':False,
            'independent_one_loop_mass_eigenstate_matching_performed':False,
            'complete_sequential_EFT_matching_performed':False,
            'existing_benchmarks_changed':False}
    (out/'certificate.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('A(-2) native Matchete metric factor:',k)
    print('Mass-squared eigenvalues (GeV^2):',mass['masses_squared_gev2'])
    print('Independent one-scale C5 loop ratio error:',loop['maximum_relative_ratio_error'])
    print('Report:',out/'certificate.json')
    print('NOTE: full mass-eigenstate C5 matching remains uncertified.')
    return 0


if __name__=='__main__':
    raise SystemExit(main())
