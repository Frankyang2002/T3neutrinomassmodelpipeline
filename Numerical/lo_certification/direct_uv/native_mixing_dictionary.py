"""Read-only live Matchete lambdaT3 CG -> restricted scalar-A witness.

Actual native CGs are input data, NEVER reconstructed from SU(2) guesses.
A successful audit validates the conditional scalar action in the numerical
builder, NOT the physical cross-package identification of lambdaT3 or y1/y2.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import sympy as sp

FORMAT = 't3-ws4-native-mixing-to-frozen-scalar-A-v1'
NATIVE_FORMAT = 't3-ws4-native-mixing-tensor-five-classes-v1'
SOURCE_FILES = (
    'Lagrangian/model/T3ModelCatalog.wl',
    'Lagrangian/model/LagrangianBuilder.wl',
    'Lagrangian/model/SU2Invariants.wl',
    'Lagrangian/model/T3Fields.wl',
    'Lagrangian/interactions/T3Topology.wl',
    'scripts/BuildT3DirectOneLoopInputs.py',
)
CLASSES = {
    'A': (1, 3, 2, (-4, -2, 0)),
    'B': (2, 2, 1, (-3, -1, 1)),
    'C': (2, 2, 3, (-3, -1, 1)),
    'D': (3, 1, 2, (-2, 0, 2)),
    'E': (3, 3, 2, (-4, -2, 0, 2)),
}
REPRESENTATIVES = {'A': 0, 'B': -1, 'C': -1, 'D': 0, 'E': 0}
FLAGS = (False, False, False, True)  # source T3MixCG field order H,H,S1,Bar[S2]
# The actual committed audit tree used by the frozen input builder.
AUDIT_DIR = Path('output/general_observable/lo_direct_one_loop/source_audits')


class MixingAuditError(ValueError):
    pass


def require(ok: bool, msg: str) -> None:
    if not ok:
        raise MixingAuditError(msg)


def read(path: Path) -> dict:
    require(path.is_file(), f'Missing required file: {path}')
    data = json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(data, dict), f'Expected JSON object: {path}')
    return data


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def complex_number(value: Any, *, real: bool = False) -> complex:
    if isinstance(value, dict):
        value = complex(value['real'], value['imag'])
    z = complex(value)
    require(np.isfinite(z.real) and np.isfinite(z.imag), 'Non-finite numerical coupling')
    if real:
        require(abs(z.imag) < 1e-12, 'Expected real scalar parameter')
    return z


def complex_matrix(value: Any, shape: tuple[int, ...] | None = None) -> np.ndarray:
    if isinstance(value, dict):
        value = np.asarray(value['real'], dtype=float) + 1j*np.asarray(value['imag'], dtype=float)
    arr = np.asarray(value, dtype=complex)
    if shape is not None:
        require(arr.shape == shape, f'Tensor shape {arr.shape} rather than {shape}')
    require(bool(np.isfinite(arr).all()), 'Non-finite input tensor')
    return arr


def source_guard(root: Path, report: dict) -> dict:
    native_hashes = report.get('native_source_sha256')
    require(isinstance(native_hashes, dict), 'Native source digest map missing')
    sha = {}
    for rel in SOURCE_FILES:
        path = root / rel
        require(path.is_file(), f'Native guarded source missing: {rel}')
        actual = digest(path)
        require(actual == native_hashes.get(rel), f'Source changed since native export: {rel}')
        sha[rel] = actual
    topology = (root / 'Lagrangian/interactions/T3Topology.wl').read_text(encoding='utf-8')
    builder = (root / 'scripts/BuildT3DirectOneLoopInputs.py').read_text(encoding='utf-8')
    group = (root / 'Lagrangian/model/SU2Invariants.wl').read_text(encoding='utf-8')
    for token in (
        'lambdaT3[] H[i] H[j] s1 s2 CG[T3MixCG, labels]',
        '"Expression" -> PlusHc[',
    ):
        require(token in topology, f'Mixing source operator changed: {token}')
    for token in ('{2, 2, d1, d2}', '{False, False, False, True}', '{1, 2}'):
        require(token in group, f'Mixing group CG flags changed: {token}')
    for token in (
        "cg['T3MixCG'][1,1,s,t]",
        "Vcoeff=-lam*v*v*cg['T3MixCG'][1,1,s,t]/2",
        'A[ds1+t,s]+=Vcoeff',
        'A[s,ds1+t]+=Vcoeff.conjugate()',
    ):
        require(token in builder, f'Frozen source mass-construction convention changed: {token}')
    return {'source_sha256': sha, 'source_forms_verified': True,
            'native_export_sha_matches_current_sources': True}


def models():
    for cls, (d1, d2, df, alphas) in CLASSES.items():
        for alpha in alphas:
            name = f'{cls}_{alpha:+d}'
            key = f'T3_dS1_{d1}_dS2_{d2}_dF_{df}_alpha_{"m" if alpha < 0 else "p"}{abs(alpha)}'
            yield name, key, (d1, d2, df, alpha)


def native_tensor(row: dict, cls: str) -> np.ndarray:
    d1, d2, df, _alphas = CLASSES[cls]
    alpha = REPRESENTATIVES[cls]
    require(row.get('class') == cls and row.get('model') == f'{cls}_{alpha:+d}' and
            row.get('alpha') == alpha and
            (row.get('dS1'), row.get('dS2'), row.get('dF')) == (d1, d2, df),
            f'{cls}: incorrect native representative metadata')
    require(row.get('native_build_status') == 'Success' and
            row.get('native_full_lagrangian_validation') is True and
            row.get('mixing_operator_name') == 'T3ScalarMix_GenericCG' and
            row.get('mixing_coupling') == 'lambdaT3' and
            row.get('native_nonhc_operator_present') is True and
            row.get('native_hc_operator_present') is True and
            row.get('Higgs_fields_in_nonhc_operator') == 'H H' and
            row.get('scalar_fields_in_nonhc_operator') == 'S1 Bar[S2]' and
            row.get('native_nonhc_LInt_coefficient') == '+lambdaT3 * T3MixCG (no 1/2)' and
            row.get('physical_source_to_Weyl_certified') is False,
            f'{cls}: native LInt operator orientation/metadata not verified')
    rec = row.get('native_cg') or {}
    native_shape = tuple(d for d in (2, 2, d1, d2) if d > 1)
    require(rec.get('tensor_accessor') == 'Matchete`GetCGTensor' and
            rec.get('registered_in_assembled_interaction') is True and
            rec.get('symbol') == 'T3MixCG' and
            rec.get('dimensions') == list(native_shape) and
            rec.get('full_dimensions_including_singlets') == [2, 2, d1, d2],
            f'{cls}: tensor is not an expected registered native mixing CG')
    exact = rec.get('components_inputform')
    re = rec.get('components_real')
    im = rec.get('components_imag')
    count = int(np.prod(native_shape))
    require(all(isinstance(v, list) and len(v) == count for v in (exact, re, im)),
            f'{cls}: native mixing CG component count/fields mismatch')
    require(all(isinstance(x, str) and x for x in exact), 'Missing exact native CG source expressions')
    flat = np.asarray(re, dtype=float) + 1j*np.asarray(im, dtype=float)
    tensor = complex_matrix(flat.reshape((2, 2, d1, d2)))
    scale = np.linalg.norm(tensor)
    require(scale > 1e-12, f'{cls}: identically zero registered mixing CG')
    require(np.linalg.norm(tensor - tensor.swapaxes(0, 1)) < 1e-12 * scale,
            f'{cls}: Higgs Bose exchange symmetry broken')
    return tensor


def native_models(report: dict) -> dict[str, np.ndarray]:
    require(report.get('format') == NATIVE_FORMAT and
            report.get('status') == 'five_classes_native_registered_mixing_tensors_extracted' and
            report.get('matching_or_running_performed') is False and
            report.get('native_rgbeta_mixing_coupling_dictionary_certified') is False and
            report.get('native_spinor_and_mass_dictionary_certified') is False,
            'Native mixing export format/status/certification inconsistent')
    env = report.get('restricted_environment') or {}
    require(env.get('isolated_to_wolframscript_process') is True and
            all(env.get('effective', {}).get(key) == '0' for key in (
                'T3_EXPERIMENTAL_GENERAL_FERMIONS',
                'T3_EXPERIMENTAL_GENERAL_QUADRATICS',
                'T3_EXPERIMENTAL_GENERAL_QUARTICS',
                'T3_EXPERIMENTAL_FLAVOR_MASSES',
            )), 'Native mixing export used nonrestricted Lagrangian')
    rows = report.get('models')
    require(isinstance(rows, list) and len(rows) == 5, 'Expected 5 native mixing rows')
    output = {}
    for row in rows:
        cls = row.get('class')
        require(cls in CLASSES and cls not in output, f'Unknown/repeated native class: {cls}')
        output[cls] = native_tensor(row, cls)
    require(set(output) == set(CLASSES), 'Not all five classes exported')
    return output


def source_tensor(case: dict, dims: tuple[int, int, int, int]) -> np.ndarray:
    d1, d2, _df, _a = dims
    vertex = case['vertices']['T3MixCG']
    full = (2, 2, d1, d2)
    require(tuple(vertex['dimensions']) == full, 'Legacy CG audit dimensions differ')
    out = np.zeros(full, dtype=complex)
    for component in vertex['components']:
        indices = tuple(int(i)-1 for i in component['indices'])
        require(len(indices) == 4 and all(0 <= v < d for v, d in zip(indices, full)),
                'Legacy CG component indices invalid')
        expr = sp.sympify(component['cg'], locals={'I': sp.I})
        require(bool(expr.is_number), 'Legacy source CG contains non-numeric symbolic component')
        out[indices] += complex(sp.N(expr, 18))
    require(np.linalg.norm(out) > 1e-12, 'Empty legacy native mixing CG')
    return out


def generator_rotations(data: dict) -> dict[int, np.ndarray]:
    rotations: dict[int, np.ndarray] = {1: np.ones((1, 1), dtype=complex)}
    require(data.get('status') in (
        'generator_algebra_verified_CG_basis_unverified', 'generator_algebra_verified',
        'all_representations_verified', None), 'Unknown generator-audit status')
    for row in data['representations']:
        require(row.get('status') == 'generator_algebra_verified', 'Generator representation unverified')
        d = int(row['dimension'])
        matrix = row['U_columns_are_weight_eigenvectors']
        u = np.asarray([[complex(sp.N(sp.sympify(z, locals={'I': sp.I}), 18))
                         for z in line] for line in matrix], dtype=complex)
        require(u.shape == (d, d) and np.linalg.norm(u.conj().T@u-np.eye(d)) < 1e-10,
                f'Nonunitary SU2 generator rotation d={d}')
        rotations[d] = u
    require(set(rotations) == {1, 2, 3}, 'Expected weight-axis rotations for d=1,2,3')
    return rotations


def rotated_tensor(source: np.ndarray, U: dict[int, np.ndarray], d1: int, d2: int) -> np.ndarray:
    # Independent of the production builder's _tensor() implementation.
    return np.einsum('abcd,ai,bj,ck,dl->ijkl',source,
                     U[2], U[2], U[d1], U[d2].conj(), optimize=True)


def charge(d: int, y: float, index: int) -> float:
    return (d-1)/2-index+y


def finite_hessian(A: np.ndarray) -> tuple[np.ndarray, float]:
    """Independent Cartesian finite differences of the scalar action V2=phi†A phi.

    Tests the -LInt -> V conversion and the complex off-diagonal Hessian;
    numerical differentiation is independent of the direct_one_loop engine.
    """
    n = len(A)
    def v(x):
        phi = (x[:n] + 1j*x[n:])/np.sqrt(2.)
        return np.vdot(phi,A@phi).real
    H = np.empty((2*n,2*n))
    eye = np.eye(2*n)
    for i in range(2*n):
        for j in range(2*n):
            H[i,j] = (v(eye[i]+eye[j])-v(eye[i]-eye[j]))/2.
    target = np.block([[A.real,-A.imag],[A.imag,A.real]])
    residual = np.linalg.norm(H-target)/max(1.,np.linalg.norm(target))
    return H, float(residual)


def mixing_block(tensor: np.ndarray, lam: complex, vev: float) -> np.ndarray:
    require(vev > 0 and np.isfinite(vev), 'Invalid Higgs vev')
    return -lam * vev**2/2. * tensor[1, 1, :, :].T


def all_sixteen(root: Path, native: dict[str, np.ndarray]) -> dict:
    audit_path = root / AUDIT_DIR / 'component_support_audit.json'
    gen_path = root / AUDIT_DIR / 'matchete_generator_analysis.json'
    cov_path = root / AUDIT_DIR / 'cg_covariance_audit.json'
    cases = read(audit_path)
    gen = read(gen_path)
    cov = read(cov_path)
    require(not cases.get('errors') and len(cases.get('cases', [])) == 16,
            'Legacy native component support audit is not complete for 16 models')
    require(cov.get('status') == 'cg_generator_covariance_verified',
            'CG covariance audit is not verified')
    covariance = {(int(r['model']['dS1']),int(r['model']['dS2']),
                   int(r['model']['dF']),int(r['model']['alpha'])):r['all_covariant']
                  for r in cov['models']}
    lookup = {(int(row['model']['dS1']),int(row['model']['dS2']),
               int(row['model']['dF']),int(row['model']['alpha'])):row
              for row in cases['cases']}
    # These JSONs are not arbitrary replacements for the original CG audit:
    # the frozen inputs record the exact SHA256 bytes that generated them.
    audit_sha = digest(audit_path)
    generator_sha = digest(gen_path)
    rows=[]
    U = generator_rotations(gen)
    for model,key,dims in models():
        d1,d2,df,alpha = dims
        cls = model[0]
        require(covariance.get(dims) is True and dims in lookup,
                f'{model}: missing verified source CG / covariance audit')
        old = source_tensor(lookup[dims],dims)
        live = native[cls]
        cg_err = np.linalg.norm(old-live)/max(np.linalg.norm(live),1e-30)
        require(cg_err < 2e-11,
                f'{model}: newly registered Matchete mixing CG differs from old source component audit ({cg_err:.4g})')
        weight = rotated_tensor(live,U,d1,d2)
        weight_old = rotated_tensor(old,U,d1,d2)
        require(np.linalg.norm(weight-weight_old)/max(np.linalg.norm(weight),1e-30) < 2e-11,
                f'{model}: weight basis transport disagrees')
        # For neutral H index 1, the scalar electric charges must match.
        for s in range(d1):
            for t in range(d2):
                if abs(weight[1,1,s,t]) < 1e-12:continue
                require(abs(charge(d1,alpha/2,s)-charge(d2,(alpha+2)/2,t)) < 1e-10,
                        f'{model}: native mixing violates U(1)em charge')
        benchmark_path = root / 'Numerical/benchmarks/fitted_reference'/(key+'.json')
        inp_path = root/'output/general_observable/lo_direct_one_loop/provisional_inputs'/(model+'.json')
        bench, frozen = read(benchmark_path), read(inp_path)
        require(frozen.get('model') == model, f'{model}: frozen input label mismatch')
        provenance = frozen.get('provenance') or {}
        def canonical_source_path(raw: Any) -> str:
            # Input provenance was generated on Windows; normalise separators.
            require(isinstance(raw, str), f'{model}: missing source-audit provenance path')
            return raw.replace('\\', '/')
        require(canonical_source_path(provenance.get('cg_support_path')) ==
                (AUDIT_DIR / 'component_support_audit.json').as_posix(),
                f'{model}: frozen CG provenance points to a different source-audit directory')
        require(canonical_source_path(provenance.get('generator_analysis_path')) ==
                (AUDIT_DIR / 'matchete_generator_analysis.json').as_posix(),
                f'{model}: frozen generator provenance points to a different source-audit directory')
        require(provenance.get('cg_support_sha256') == audit_sha and
                provenance.get('generator_analysis_sha256') == generator_sha,
                f'{model}: source CG/generator audit SHA256 differs from frozen input provenance')
        require(provenance.get('benchmark_sha256') == digest(benchmark_path),
                f'{model}: frozen benchmark SHA256 does not match provenance')
        rep = bench['representation']
        require((rep['d_s1'],rep['d_s2'],rep['d_f'],rep['alpha']) == dims,
                f'{model}: benchmark metadata mismatch')
        ordinary = bench['base_state']['ordinary']
        vev = complex_number(bench['numerical']['vev_gev'],real=True).real
        lam = complex_number(ordinary['lambdaT3'])
        A = complex_matrix(frozen['scalar_A_gev2'], (d1+d2,d1+d2))
        B = complex_matrix(frozen['scalar_B_gev2'], (d1+d2,d1+d2))
        require(np.linalg.norm(A-A.conj().T) < 1e-11 * max(1.,np.linalg.norm(A)),
                f'{model}: frozen scalar-A not Hermitian')
        require(np.linalg.norm(B) < 1e-13 * max(1.,np.linalg.norm(A)),
                f'{model}: nonzero scalar B not supported in restricted mixing audit')
        actual = A[d1:, :d1]
        expected = mixing_block(weight,lam,vev)
        scale = max(float(np.linalg.norm(expected)),1e-50)
        err = float(np.linalg.norm(actual-expected)/scale)
        require(err < 3e-10, f'{model}: native mixing-to-scalar-A mismatch: relative {err:.5g}')
        expected_A = np.zeros_like(A)
        expected_A[d1:, :d1] = expected
        expected_A[:d1,d1:] = expected.conj().T
        _, cartesian_err = finite_hessian(expected_A)
        require(cartesian_err < 1e-11, f'{model}: independent real scalar Hessian mismatch')
        rows.append({
            'model':model,'native_class':cls,
            'benchmark_sha256':digest(benchmark_path),
            'frozen_input_sha256':digest(inp_path),
            'original_source_cg_relative_residual':float(cg_err),
            'relative_complex_scalar_A_mixing_block_residual':err,
            'independent_Cartesian_Hessian_relative_residual':cartesian_err,
            'scalar_mixing_block_gev2':{'real':actual.real.tolist(),'imag':actual.imag.tolist()},
            'mass_mixing_normalization_under_identified_lambdaT3':'-lambdaT3*v^2/2 * CG[H0,H0,S1,Bar[S2]]',
            'native_RGBeta_lambdaT3_normalization_certified':False,
            'physical_source_to_Weyl_certified':False,
        })
    return {'models_checked':len(rows),
        'component_support_audit_path': AUDIT_DIR.as_posix() + '/component_support_audit.json',
        'component_support_audit_sha256': audit_sha,
        'generator_analysis_audit_sha256': generator_sha,
        'frozen_input_provenance_sha256_verified': True,
        'all_native_CG_components_match_legacy_source_audit':True,
        'all_frozen_mixing_blocks_match_conditional_native_action':True,
        'max_native_CG_relative_residual':max(r['original_source_cg_relative_residual'] for r in rows),
        'max_scalar_A_mixing_relative_residual':max(r['relative_complex_scalar_A_mixing_block_residual'] for r in rows),
        'max_independent_Cartesian_Hessian_residual':max(r['independent_Cartesian_Hessian_relative_residual'] for r in rows),
        'models':rows}


def audit(root: Path, native_json: Path, *, check_frozen: bool = True) -> dict:
    native_report = read(native_json)
    tensors = native_models(native_report)
    source = source_guard(root,native_report)
    frozen = all_sixteen(root,tensors) if check_frozen else None
    return {
        'format':FORMAT,
        'status':'native_lambdaT3_mixing_CG_and_restricted_scalar_A_checked_conditional_physical_dictionary',
        'native_report_sha256':digest(native_json),
        'source_contract':source,
        'native_registered_mixing_classes_checked':5,
        'native_CG_tensors_symmetric_under_H_exchange':True,
        'source_LInt_has_lambdaT3_HH_S1_BarS2_and_HC':True,
        'Higgs_neutral_weight_index':1,
        'Higgs_vev_convention':'H0=v/sqrt(2)',
        'conditional_scalar_A_lower_left_block':'A[S2_t,S1_s]=-lambdaT3*v^2/2*C[H0,H0,S1_s,BarS2_t]',
        'no_identical_H_factor_two':'H0*H0=v^2/2; no extra factor two in the potential',
        'actual_frozen_inputs':frozen,
        'physical_Matchete_RGBeta_lambdaT3_coupling_map_certified':False,
        'physical_y1_y2_MF_GammaCC_and_SU2_mass_dictionary_certified':False,
        'full_physical_source_to_Weyl_certified_models':0,
        'remaining_native_tasks':[
            'Fix RGBeta delS2 topology mixing invariant components and coupling normalization in common Matchete weight axes',
            'Derive signed MF and y1,y2 Feynman-rule dictionary using native GammaCC and epsilon in one left-Weyl basis',
            'Measure scalar free-mass sign and normalization against full native quadratic action',
            'Treat six optional symmetry-allowed conjugate-scalar Yukawas as zero only in the restricted benchmark',
            'One-loop hard matching: integrate ALL heavy BSM fermions at muF=1e5 GeV while keeping both BSM scalars and SM fermions',
            'Integrate BOTH BSM scalars at muS=1e3 GeV into SMEFT',
        ],
        'EFT_order':'UV -> integrate all fermions -> EFT1 (both BSM scalars + SM fermions) -> integrate both scalars -> SMEFT',
        'one_loop_or_EFT_matching_performed':False,
    }
