"""A(-2) native-tensor and absolute one-scale anomalous-C5 audit.

No modification to production matching, RGBeta exports or benchmark inputs.
Requires the previously passed VerifyA2BilinearNormalization certificate.
The Matchete<->RGBeta triplet-field basis mapping and sequential-EFT matching
remain explicitly uncertified.

PowerShell (repository root):
 python -u -X utf8 -m scripts.VerifyA2ExtraQuarticC5 \
   --normalization-run-id a2_norm_v1 --run-id a2_extra_c5_v1
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
import re
from typing import Any

import mpmath as mp
import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
INPUT_ROOT = ROOT / 'output/studies/general_scalar_extended/a2_bilinear_normalization'
OUTPUT_ROOT = ROOT / 'output/studies/general_scalar_extended/a2_extra_quartic_c5'
MODEL = 'A_-2'
MAT_QUARTIC = 'lambdaT3ExtraQ020101I1'
RGB_QUARTIC = 'lambdaT3GeneralQ020101I1'
FORMAT = 't3-a-minus2-native-anomalous-c5-one-scale-audit-v1'


def require(ok: bool, reason: str) -> None:
    if not ok:
        raise ValueError('A(-2) anomalous C5 audit: ' + reason)


def read(path: Path) -> dict:
    require(path.is_file(), 'missing '+str(path))
    obj = json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(obj, dict), 'not an object: '+str(path))
    return obj


def digest(path: Path) -> str:
    return sha256(Path(path).read_bytes()).hexdigest()


def complex_object(z: complex) -> dict:
    z = complex(z)
    require(np.isfinite([z.real, z.imag]).all(), 'nonfinite complex output')
    return {'real': float(z.real), 'imag': float(z.imag)}


def matrix_object(mat: np.ndarray) -> dict:
    mat = np.asarray(mat, dtype=complex)
    require(mat.shape == (3, 3) and np.isfinite(mat).all(), 'nonfinite or non-3x3 C5')
    return {'real': mat.real.tolist(), 'imag': mat.imag.tolist(),
            'frobenius_norm': float(np.linalg.norm(mat))}


def decode_complex(raw: dict, shape: tuple[int, ...]) -> np.ndarray:
    require(isinstance(raw, dict) and set(raw) == {'real', 'imag'}, 'complex matrix serialization')
    a = np.asarray(raw['real'], dtype=float) + 1j*np.asarray(raw['imag'], dtype=float)
    require(a.shape == shape and np.isfinite(a).all(), 'complex matrix shape/finite')
    return a


def normalized_error(x: np.ndarray, y: np.ndarray) -> float:
    a, b = np.asarray(x, dtype=complex), np.asarray(y, dtype=complex)
    require(a.shape == b.shape, 'comparison shape mismatch')
    scale = max(float(np.linalg.norm(a)),float(np.linalg.norm(b)))
    require(scale > 0 and np.isfinite(scale), 'comparison has zero/nonfinite scale')
    return float(np.linalg.norm(a-b)/scale)


def matching_channels(plus: np.ndarray, minus: np.ndarray, *,
                      lambda_t3: complex, lambda_rgbeta: complex,
                      mu_rgbeta: complex, tensor_ratio: complex) -> tuple[np.ndarray,np.ndarray]:
    """One-scale normal and anomalous C5 from exact mass eigenpropagators.

    Q_mat = Hdag Hdag S1dag S2dag, so its H H conjugate has lambda_mat*.
    With L_mat=+lambda_mat Q+h.c. and L_RGBeta=-lambda_rgb Q+h.c.,
    lambda_mat=-lambda_rgb ONLY after identifying the scalar components.
    The A(-2) anomalous <S2 S2> propagator carries mu/|mu|.
    """
    p, m = np.asarray(plus,dtype=complex),np.asarray(minus,dtype=complex)
    require(p.shape == m.shape == (3,3) and np.isfinite(p).all() and np.isfinite(m).all(),
            'invalid split kernels')
    for name,z in [('lambdaT3',lambda_t3),('extra quartic',lambda_rgbeta),('mu2',mu_rgbeta),('ratio',tensor_ratio)]:
        require(np.isfinite([complex(z).real,complex(z).imag]).all(), name+' nonfinite')
    phase = complex(mu_rgbeta)/abs(mu_rgbeta) if mu_rgbeta else 0j
    lambda_mat = -complex(lambda_rgbeta)
    normal = complex(lambda_t3)*0.5*(p+m)
    anomalous = lambda_mat.conjugate()*phase*complex(tensor_ratio)*0.5*(p-m)
    for name,x in [('normal',normal),('anomalous',anomalous)]:
        require(np.allclose(x,x.T,rtol=3e-12,atol=1e-28),name+' C5 not symmetric')
    return normal,anomalous


def direct_spectral_propagators(p2: float, m_plus_sq: float, m_minus_sq: float,
                                mu2: complex) -> tuple[complex,complex]:
    """Trace-normal and metric-stripped anomalous propagators (A(-2))."""
    require(min(p2,m_plus_sq,m_minus_sq)>=0 and p2+m_plus_sq>0 and p2+m_minus_sq>0,
            'invalid Euclidean propagator masses')
    q1,q2=1/(p2+m_plus_sq),1/(p2+m_minus_sq)
    phase=complex(mu2)/abs(mu2) if mu2 else 0j
    return 0.5*(q1+q2),0.5*phase*(q1-q2)


def scalar_triangle(mf: float, ms1_sq: float, ms2_sq: float) -> mp.mpf:
    """Independent positive Euclidean loop ∫dt t/[(t+MF²)(t+m1²)(t+m2²)]."""
    require(all(math.isfinite(float(x)) and float(x)>0 for x in (mf,ms1_sq,ms2_sq)),
            'nonpositive integral masses')
    with mp.workdps(65):
        m=mp.mpf(str(mf)); a=mp.mpf(str(ms1_sq))/m**2; b=mp.mpf(str(ms2_sq))/m**2
        f=lambda t: t/((t+1)*(t+a)*(t+b))
        return +(mp.quad(f,[0,mp.mpf('0.00001'),mp.mpf('0.001'),mp.mpf('0.1'),1,10,mp.inf])/m**2)


def test_source_kernel_absolute(source: Path, mf: np.ndarray, m1sq: float,
                                eigen_sq: tuple[float,float]) -> tuple[np.ndarray,np.ndarray,dict]:
    """Compare the *absolute* Matchete one-loop kernel to a finite integral.

    For A(-2), stripping both Yukawas and setting lambdaT3=1 gives
      K(M,m1,m2) = -(2/sqrt(3))*hbar*M*J(M,m1,m2).
    The -2/sqrt(3) is from the source CG algebra, not a fitted factor.
    """
    from RGE.matching.MatcheteC5Parsing import parse_matchete_c5
    from RGE.matching.WeinbergFlavorMatching import extract_one_generation_loop_kernel
    raw=Path(source).read_text(encoding='utf-8-sig')
    require(all(x in raw for x in ('Coupling[MF,','Coupling[MS1,','Coupling[MS2,',
               'Coupling[lambdaT3,','Coupling[y1,','Coupling[y2,')), 'invalid mass-dependent one-scale source')
    expr=extract_one_generation_loop_kernel(parse_matchete_c5(raw))
    names=('MF','MS1','MS2','lambdaT3','hbar')
    syms=tuple(sp.Symbol(s) for s in names)
    require(expr.free_symbols <= set(syms) and all(s in expr.free_symbols for s in
            (sp.Symbol('MF'),sp.Symbol('MS1'),sp.Symbol('MS2'),sp.Symbol('lambdaT3'))),
            'wrong symbolic Matchete kernel')
    fn=sp.lambdify(syms,expr,modules='mpmath')
    plus,minus=[],[]
    checks=[]
    with mp.workdps(65):
        hbar=1/(16*mp.pi**2)
        ms1=mp.sqrt(mp.mpf(str(m1sq)))
        for r,mass in enumerate(mf):
            for state,m2sq,arr in [('plus',eigen_sq[0],plus),('minus',eigen_sq[1],minus)]:
                m=mp.mpf(str(float(mass)))
                ms2=mp.sqrt(mp.mpf(str(m2sq)))
                got=mp.mpc(fn(m,ms1,ms2,1,hbar))
                expected=-2/mp.sqrt(3)*hbar*m*scalar_triangle(float(mass),m1sq,m2sq)
                err=abs(got-expected)/max(abs(got),abs(expected))
                require(err<mp.mpf('2e-9'),f'absolute C5 source/integral disagrees at heavy {r}, {state}: {mp.nstr(err,6)}')
                require(abs(got.imag)<mp.mpf('1e-28')*max(abs(got),mp.mpf('1e-90')),
                        'unexpected imaginary one-scale source with real masses')
                arr.append(float(got.real))
                checks.append({'heavy_generation':r+1,'mass_gev':float(mass),'split_state':state,
                               'source_coefficient_gev_inv':float(got.real),
                               'integral_coefficient_gev_inv':float(expected),
                               'relative_error':float(err)})
    return np.asarray(plus),np.asarray(minus),{
        'identity':'K=-2/sqrt(3) * hbar * MF * integral_0^infty[t dt/((t+MF^2)(t+MS1^2)(t+MS2^2))]',
        'full_absolute_coefficient_not_fitted':True,
        'maximum_relative_error':max(x['relative_error'] for x in checks),
        'checks':checks,
        'source_sha256':digest(source)}


def lift_flavour(unit_weights: np.ndarray, y1: np.ndarray, y2: np.ndarray) -> np.ndarray:
    """The source one-loop coefficient lifted with three nondegenerate MF_r."""
    w=np.asarray(unit_weights,dtype=complex)
    require(w.shape==(3,) and y1.shape==y2.shape==(3,3),'heavy flavour dimensions')
    a,b=y1.conjugate(),y2.conjugate()
    result=(a*w)@b.T+(b*w)@a.T
    require(np.all(np.isfinite(result)) and np.allclose(result,result.T,rtol=1e-12,atol=1e-26),
            'flavour-lifted C5 not finite/symmetric')
    return result


def verify_native_tensors(root: Path) -> dict:
    from Numerical.experimental.T3SelfBilinearAnomalousKernel import tensor_ratio
    stage=root/'output/general_t3_stage22/matchete_tensors.json'
    group=root/'output/group_factors/T3_A_uv_cg_registry.json'
    obj=tensor_ratio(read(stage),read(group),MODEL,1,3)
    ratio=complex(obj['ratio'])
    require(abs(ratio-1)<2e-12 and obj['relative_tensor_residual']<2e-12 and
            obj['native_self_metric_residual']<2e-12 and
            obj['matchete_extra_symbol']==MAT_QUARTIC and
            obj['rgbeta_extra_symbol']==RGB_QUARTIC,
            'native Matchete anomalous tensor is not exactly the T3 CG')
    return {'native_ratio':complex_object(ratio),
            'native_relative_tensor_residual':obj['relative_tensor_residual'],
            'native_metric_residual':obj['native_self_metric_residual'],
            'matchete_extra_quartic':MAT_QUARTIC, 'rgbeta_extra_quartic':RGB_QUARTIC,
            'source_sha256':{'stage22':digest(stage),'group_CG':digest(group)},
            'matchete_to_rgbeta_triplet_component_intertwiner_certified':False}


def verify_quartic_source_contract(root: Path) -> dict:
    registry=root/'RGE/general/generated/general_rgbeta_quartic_registry.json'
    extra=(read(registry)['models'][MODEL]['new_rgbeta_symbols'])
    rows=[x for x in extra if x.get('name')==RGB_QUARTIC]
    require(len(rows)==1 and rows[0].get('fields')==[
            'Bar @ H','Bar @ H','Bar @ S1','Bar @ S2'] and
            rows[0].get('matchete_basis_name')=='lambdaT3ExtraQ020101I#' and
            rows[0].get('mapping_status')=='NotNormalised',
            'RGBeta quartic field word/registry changed')
    mat=root/'Lagrangian/interactions/GeneralScalarOperators.wl'
    a=mat.read_text(encoding='utf-8')
    require('"Expression" -> PlusHc[contribution]' in a and
            'contribution = coupling[] product' in a,
            'native Matchete quartic sign changed')
    rgb_logic=root/'Numerical/experimental/GeneralT3LagrangianConvention.py'
    b=rgb_logic.read_text(encoding='utf-8')
    require('RGBETA_MBL_EQUATION' in b and 'L contains - sum_couplings' in b,
            'RGBeta potential/Lagrangian convention provenance unavailable')
    return {
        'matchete_field_monomial':'Hdag Hdag S1dag S2dag',
        'rgbeta_field_monomial':'Hdag Hdag S1dag S2dag',
        'native_matchete_lagrangian':'+ lambdaMat * Q + h.c.',
        'rgbeta_lagrangian_convention':'- lambdaRGB * Q + h.c.',
        'conditional_coefficient_mapping':'lambdaMat = -lambdaRGB (when component fields/intertwiners and tensor normalization align)',
        'hermitian_conjugation':'HH channel uses conjugate(lambdaMat)',
        'relative_sign_in_same_component_basis':-1,
        'triplet_component_basis_map_certified':False,
        'sha256':{'registry':digest(registry),'matchete_builder':digest(mat),
                  'rgbeta_convention_source':digest(rgb_logic)}}


def verify_prerequisite(path: Path) -> tuple[dict, complex,tuple[float,float]]:
    o=read(path)
    require(o.get('status')=='A2NativeBilinearMassMapAndOneScaleKernelShapeChecked' and
            o.get('model')==MODEL and
            o.get('physical_matchete_rgbeta_entire_vertex_map_certified') is False and
            o.get('independent_one_loop_mass_eigenstate_matching_performed') is False,
            'A(-2) normalization prerequisite not valid')
    mass=o.get('mass_basis',{})
    mu=mass.get('rgbeta_mass_parameter_gev2',{})
    require(set(mu)=={'real','imag'},'missing RGBeta mu2 in certificate')
    z=complex(float(mu['real']),float(mu['imag']))
    v=mass.get('masses_squared_gev2',[])
    require(len(v)==2,'missing split eigenvalues')
    eigen=tuple(float(x) for x in v)
    m2=float(o['input_reference_fitted_mass_squared_gev2'])
    require(z!=0 and min(eigen)>0 and
            np.allclose(eigen,[m2+2*abs(z),m2-2*abs(z)],rtol=2e-12,atol=1e-8) and
            o['source_conventions']['relative_sign_if_same_canonical_field_basis']==-1 and
            o['mass_dependent_one_scale_kernel']['maximum_relative_ratio_error']<2e-9,
            'prerequisite mass split/convention/loop check changed')
    return o,z,eigen


def main() -> int:
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--normalization-run-id',default='a2_norm_v1')
    ap.add_argument('--run-id',required=True)
    ap.add_argument('--output-root',type=Path,default=OUTPUT_ROOT)
    ap.add_argument('--extra-rgbeta-real',type=float,default=0.02,
                    help='Illustrative dimensionless extra quartic; NOT a fitted value')
    ap.add_argument('--extra-rgbeta-imag',type=float,default=0.01)
    args=ap.parse_args()
    require(bool(re.fullmatch(r'[a-zA-Z0-9][a-zA-Z0-9_-]{0,63}',args.run_id)), 'unsafe run ID')
    require(all(math.isfinite(x) and abs(x)<=4*math.pi for x in
                (args.extra_rgbeta_real,args.extra_rgbeta_imag)), 'invalid extra quartic')
    output=args.output_root/args.run_id
    require(not output.exists(), 'output exists; use a new --run-id: '+str(output))
    prerequisite=INPUT_ROOT/args.normalization_run_id/'certificate.json'
    old,mu,eigen=verify_prerequisite(prerequisite)
    tensor=verify_native_tensors(ROOT)
    contract=verify_quartic_source_contract(ROOT)
    cfg=ROOT/'Numerical/benchmarks/fitted_reference/T3_dS1_1_dS2_3_dF_2_alpha_m2.json'
    ordinary=read(cfg)['base_state']['ordinary']
    MF=decode_complex(ordinary['MF'],(3,3))
    require(np.allclose(MF,np.diag(np.diag(MF)),rtol=0,atol=1e-10) and
            np.allclose(np.diag(MF).imag,0,rtol=0,atol=1e-10) and
            np.all(np.diag(MF).real>0),'benchmark MF not positive and diagonal')
    masses=np.diag(MF).real
    y1=decode_complex(ordinary['y1'],(3,3))
    y2=decode_complex(ordinary['y2'],(3,3))
    require(np.isclose(float(ordinary['mS2Sq']),sum(eigen)/2,rtol=1e-12),
            'benchmark mass differs from bilinear prerequisite')
    lam=ordinary['lambdaT3'];lam=complex(float(lam['real']),float(lam['imag']))
    lam_extra=complex(args.extra_rgbeta_real,args.extra_rgbeta_imag)
    source=ROOT/'output/hypercharge/T3_A_alpha_m2/data/c5_coefficient.txt'
    wp,wm,integrals=test_source_kernel_absolute(source,masses,float(ordinary['mS1Sq']),eigen)
    kp=lift_flavour(wp,y1,y2)
    km=lift_flavour(wm,y1,y2)
    normal,anomalous=matching_channels(kp,km,lambda_t3=lam,
              lambda_rgbeta=lam_extra,mu_rgbeta=mu,tensor_ratio=1+0j)
    # Construct a completely independent direct Euclidean three-propagator
    # evaluation with the exact A(-2) -2/sqrt(3) group prefactor.
    prefactor=-2/math.sqrt(3)/(16*math.pi**2)
    jplus=np.array([prefactor*float(m)*float(scalar_triangle(float(m),float(ordinary['mS1Sq']),eigen[0]))
                    for m in masses])
    jminus=np.array([prefactor*float(m)*float(scalar_triangle(float(m),float(ordinary['mS1Sq']),eigen[1]))
                     for m in masses])
    n2,a2=matching_channels(lift_flavour(jplus,y1,y2),lift_flavour(jminus,y1,y2),
                        lambda_t3=lam,lambda_rgbeta=lam_extra,mu_rgbeta=mu,tensor_ratio=1)
    errors={'normal_C5':normalized_error(normal,n2),'anomalous_C5':normalized_error(anomalous,a2)} if abs(lam_extra)>0 else {
            'normal_C5':normalized_error(normal,n2),'anomalous_C5':None}
    require(errors['normal_C5']<2e-9 and
            (errors['anomalous_C5'] is None or errors['anomalous_C5']<2e-9),
            'absolute normal/anomalous C5 does not agree with independent integral')
    # Cross-check expected zero limits and actual scalar propagator phase.
    n0,a0=matching_channels(kp,km,lambda_t3=lam,
              lambda_rgbeta=lam_extra,mu_rgbeta=0,tensor_ratio=1)
    require(np.linalg.norm(a0)==0 and np.array_equal(n0,normal), 'zero-mu limit failed')
    _,p_anomalous=direct_spectral_propagators(3e6,eigen[0],eigen[1],mu)
    expected_phase=mu/abs(mu)
    require(abs(p_anomalous-0.5*expected_phase*(1/(3e6+eigen[0])-1/(3e6+eigen[1])))<1e-19,
            'anomalous propagator phase wrong')
    ratio=float(np.linalg.norm(anomalous)/np.linalg.norm(normal)) if np.linalg.norm(normal)>0 else None
    output.mkdir(parents=True)
    result={
        'format':FORMAT,
        'status':'A2NativeTensorAndConditionalOneScaleAnomalousC5Checked',
        'model':MODEL,
        'scope':'Native Matchete tensor contraction + source convention + full absolute one-scale mass-dependent loop coefficient and three-flavour C5. Matchete/RGBeta triplet component map and sequential Wilson matching remain unverified.',
        'input_sha256':{'normalization_certificate':digest(prerequisite),'reference_benchmark':digest(cfg),'one_scale_C5':digest(source),'verification_script':digest(Path(__file__))},
        'native_tensor':tensor,'quartic_convention':contract,
        'mu2_rgbeta_gev2':complex_object(mu),
        'extra_quartic_rgbeta_illustrative':complex_object(lam_extra),
        'extra_quartic_matchete_conditional':complex_object(-lam_extra),
        'lambdaT3_benchmark':complex_object(lam),
        'heavy_masses_benchmark_gev':masses.tolist(),
        'mass_squared_plus_minus_gev2':list(eigen),
        'absolute_one_scale_loop':integrals,
        'unit_c5_plus':matrix_object(kp),
        'unit_c5_minus':matrix_object(km),
        'normal_one_scale_c5':matrix_object(normal),
        'anomalous_one_scale_c5':matrix_object(anomalous),
        'combined_one_scale_c5':matrix_object(normal+anomalous),
        'anomalous_over_normal_norm':ratio,
        'independent_integral_relative_errors':errors,
        'matching_phase_formula':'deltaC5 = -conj(lambdaRGBeta_extra) * (mu2/abs(mu2)) * (K_plus-K_minus)/2, conditional on a unit RGBeta<->Matchete tensor map',
        'full_native_matchete_mass_eigenstate_matching_performed':False,
        'full_rgbeta_matchete_triplet_component_map_certified':False,
        'full_sequential_eft_wilson_matching_performed':False,
        'physical_full_general_scalar_prediction':False,
        'benchmark_was_modified':False,
    }
    target=output/'certificate.json'
    target.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('A(-2) native Matchete tensor ratio:',tensor['native_ratio'])
    print('Source-vs-independent absolute one-scale C5 max relative error:',integrals['maximum_relative_error'])
    print('Anomalous/normal C5 norm (illustrative extra quartic):',ratio)
    print('Report:',target)
    print('NOT CERTIFIED: RGBeta triplet basis mapping and full sequential Wilson matching.')
    return 0


if __name__=='__main__':
    raise SystemExit(main())
