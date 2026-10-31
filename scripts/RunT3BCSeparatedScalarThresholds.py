"""Locate B/C(-1) running X,Y thresholds from the SAME EFT1 benchmark as
RunT3QuadraticEndpointC5; prepare native mass-eigenstate Yukawa boundaries.

Numerically solves mu = m_Y(mu), mu = m_X(mu) using the already-exported
candidate one-loop full-quartic + holomorphic-mass beta system. Above Y this
is an EFT1 calculation; X obtained by continuing the *un-decoupled* two-scalar
beta system below Y is PROVISIONAL until an X-only EFT2 is implemented.

No Wilson coefficient matching, finite scalar threshold, or physical
neutrino-mass prediction is performed. Original files are read-only.

  python -u -X utf8 -m scripts.RunT3BCSeparatedScalarThresholds --all \
      --run-id bc_separated_thresholds_v1
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
from typing import Any, Callable

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from Numerical.experimental.GeneralT3MixedScalarMassBasis import (
    MixedDoubletMassBasis, matrix_json, cjson,
)
from Numerical.experimental.GeneralT3ScalarRGState import (
    GeneralT3ScalarRGState, load_registry, pack_general_scalar_state,
    unpack_general_scalar_state,
)
from scripts import IntermediateEFTRunning as s34
from scripts import RunGeneralT3QuadraticEFT1 as qeft
from scripts import RunGeneralT3MassEigenMatchete as mass_eigen
from scripts import RunT3QuadraticEndpointC5 as one_scale
from scripts import RunT3FullScalarNeutrino as neutrino
from scripts import RunT3FullScalarContinuous as continuous

ROOT = continuous.ROOT
MODELS = ('B_-1', 'C_-1')
DEFAULT_EFT1 = ROOT / 'output/general_t3_quadratic_eft1'
DEFAULT_EIGEN = ROOT / 'output/general_t3_mass_eigen_matchete'
DEFAULT_ONESCALE = ROOT / 'output/studies/general_scalar_extended/quadratic_endpoint_c5'
DEFAULT_OUTPUT = ROOT / 'output/audits/t3_bc_separated_scalar_thresholds'
FORMAT = 't3-bc-running-separated-scalar-threshold-preparation-v1'
STATUS = 'RunningYThresholdLocatedXEstimateUnDecoupledMatchingNotPerformed'
PHASE_Y2 = {'B_-1': -1, 'C_-1': 1}


def require(test: bool, msg: str) -> None:
    if not bool(test):
        raise ValueError('B/C separated scalar thresholds: ' + msg)


def read(path: Path) -> dict:
    require(path.is_file(), f'missing input {path}')
    obj = json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(obj, dict), 'invalid JSON object ' + str(path))
    return obj


def sha(path: Path) -> str:
    require(path.is_file(), 'missing hash source ' + str(path))
    return sha256(path.read_bytes()).hexdigest()


def checked_complex(value: Any, name: str) -> complex:
    return neutrino.complex_pair(value, name)


def exact_scale_crossing(basis_at: Callable[[float], MixedDoubletMassBasis],
                         t_low: float, t_high: float, which: int) -> float:
    """Find m_which(mu)=mu on a descending RG interval, in log(mu)."""
    require(which in (0, 1) and math.isfinite(t_low) and math.isfinite(t_high)
            and t_low < t_high, 'invalid eigenstate or RG interval')
    def residual(t: float) -> float:
        basis = basis_at(float(t))
        require(not basis.tachyonic and basis.eigenvalues_sq[which] > 0,
                'tachyonic/zero scalar squared mass along RG trajectory')
        return float(np.sqrt(basis.eigenvalues_sq[which]) - math.exp(t))
    low, high = residual(t_low), residual(t_high)
    require(low * high < 0, f'no sign-changing eigenstate {which} mass crossing '
            f'in [{math.exp(t_low):g}, {math.exp(t_high):g}] GeV')
    root = brentq(residual, t_low, t_high, xtol=1e-13, rtol=9e-15)
    require(abs(residual(root)) < 1e-8 * math.exp(root),
            'mass crossing equation not satisfied')
    return float(root)


def yukawa_weights(model: str, U: np.ndarray, alpha1: int, alpha2: int) -> dict:
    """Weights multiplying distinct native Matchete y1/y2 spinor/CG vertices.

    y1: CConj[F], scalar S1. y2: F, scalar Bar[S2]. Do not add vertices!
    """
    require(model in MODELS and (alpha1, alpha2) == (1, PHASE_Y2[model]),
            'wrong native Matchete Yukawa CG phase')
    U = np.asarray(U, dtype=complex)
    require(U.shape == (2, 2) and np.isfinite(U.real).all() and
            np.isfinite(U.imag).all() and
            np.allclose(U.conj().T @ U, np.eye(2), rtol=1e-12, atol=1e-12),
            'mixing matrix invalid')
    w1, w2 = alpha1 * U[0, :], -alpha2 * U[1, :]
    require(np.allclose(np.linalg.norm(w1), 1) and
            np.allclose(np.linalg.norm(w2), 1) and
            abs(np.vdot(w1, w2)) < 1e-11,
            'Yukawa mixing weights inconsistent with unitarity')
    return {
        'y1_CConj_fermion_S1_source': [cjson(z) for z in w1],
        'y2_unconjugated_fermion_BarS2_source': [cjson(z) for z in w2],
        'native_y1_cg_projection': alpha1,
        'native_y2_epsilon_cg_projection': alpha2,
        'y2_epsilon_field_replacement_sign': -1,
    }


def boundary(model: str, t: float, basis: MixedDoubletMassBasis,
             alpha1: int, alpha2: int, eigenstate: str) -> dict:
    require(eigenstate in ('X', 'Y'), 'unsupported eigenstate')
    a = 0 if eigenstate == 'X' else 1
    mass = float(np.sqrt(basis.eigenvalues_sq[a]))
    require(abs(mass / math.exp(t) - 1) < 1e-8, 'off-shell threshold scale')
    p = np.outer(basis.U[:, a], basis.U[:, a].conjugate())
    require(np.allclose(p @ p, p, atol=2e-12) and
            np.allclose(p, p.conj().T, atol=2e-12) and
            abs(np.trace(p) - 1) < 2e-12, 'bad scalar flavour projector')
    return {
        'eigenstate': eigenstate,
        'mu_gev': math.exp(t),
        'running_mass_gev': mass,
        'other_eigenmass_gev': float(np.sqrt(basis.eigenvalues_sq[1-a])),
        'mass_squared_matrix_gev2': matrix_json(basis.A),
        'mu2_gev2': cjson(basis.mu_sq),
        'unitary_U': matrix_json(basis.U),
        'projector_in_phi_multiplet': matrix_json(p),
        'projector_rank': 1,
        'yukawa_weights': yukawa_weights(model, basis.U, alpha1, alpha2),
        'projection_convention': 'Phi=(S1, epsilon*conjugate(S2)); columns U=(X,Y)',
    }


def check_native_mixing(model: str, eigen: dict, source: dict,
                        endpoint: MixedDoubletMassBasis) -> tuple[int, int]:
    require(eigen.get('status') == mass_eigen.STATUS_PASSED and
            eigen.get('model') == model and
            eigen.get('mass_basis_ordinary_check_passed') is True and
            eigen.get('heavy_mass_basis_detailed_check_passed') is True and
            eigen.get('both_yukawa_interactions_present') is True and
            eigen.get('all_17_quartics_present') is True and
            eigen.get('physical_matching_performed') is False,
            'saved live Matchete mass-eigenbasis certificate invalid')
    y1 = eigen.get('native_yukawa1_tensor', {})
    y2 = eigen.get('native_yukawa2_epsilon_tensor', {})
    alpha1 = y1.get('OldToNewCGFactor')
    alpha2 = y2.get('OldTimesEpsilonCGFactor')
    require(y1.get('ExactTensorProjectionVerified') is True and
            y2.get('ExactTensorProjectionVerified') is True and
            alpha1 == '1' and alpha2 == str(PHASE_Y2[model]),
            'native Yukawa CG projection factors changed')
    spurions = eigen.get('rotation_spurion_values_from_numerical_U', {})
    require(set(spurions) == {f'T3MEU{i}{j}' for i in (1, 2) for j in (1, 2)},
            'Matchete numerical U inventory changed')
    saved_U = np.array([[checked_complex(spurions[f'T3MEU{i}{j}'], 'U')
                         for j in (1, 2)] for i in (1, 2)], dtype=complex)
    require(np.allclose(saved_U, endpoint.U, atol=5e-9, rtol=5e-9),
            'native Matchete mass-eigenstate U differs from EFT1 endpoint')
    prepared = mass_eigen._complex_matrix(source.get('unitary_U'), 'prepared U', 2)
    require(np.allclose(prepared, endpoint.U, rtol=5e-9, atol=5e-9),
            'saved prepared mass eigenbasis differs from quadratic EFT1 endpoint')
    return int(alpha1), int(alpha2)


def validate_inputs(model: str, q: dict, reference: dict) -> None:
    require(model in MODELS and q.get('model') == model and
            q.get('status') == 'CandidateQuadraticAndFullQuarticEFT1Integrated' and
            q.get('format') == qeft.FORMAT and
            q.get('quadratic_kind') == 'mixed' and
            q.get('mu2_beta_included') is True and
            q.get('new_quartics_evolved') == 1 and
            q.get('reused_quartics_evolved') == 6 and
            q.get('physical_matchete_rgbeta_mapping_certified') is False and
            q.get('scalar_mass_threshold_matching_performed') is False and
            q.get('wilson_coefficient_matching_performed') is False,
            'quadratic EFT1 run incompatible with this candidate audit')
    require(reference.get('status') == 'ConditionalEFT1ScalarEndpointOneScaleC5NotPhysicalSequential' and
            reference.get('model') == model and
            reference.get('endpoint') == 'final' and
            reference.get('running_sha256') is not None and
            reference.get('scalar_threshold_matching_performed') is False,
            'one-scale reference invalid')


def solve(model: str, *, input_dir: Path, eigen_dir: Path,
          reference_dir: Path, out_dir: Path) -> dict:
    tag = model.replace('-', 'm').replace('+', 'p')
    source_path = input_dir / f'{tag}_eft1_quadratic_running.json'
    beta_path = input_dir / f'{tag}_eft1_quadratic_live_rgbeta.json'
    eigen_path = eigen_dir / f'{tag}_final_matchete.json'
    prepared_path = eigen_dir / f'{tag}_final_parameters.json'
    reference_path = reference_dir / f'{tag}_final_endpoint_one_scale.json'
    q, ref = read(source_path), read(reference_path)
    validate_inputs(model, q, ref)
    require(ref['running_sha256'] == sha(source_path),
            'one-scale replay used a different quadratic EFT1 source')
    registry = load_registry(continuous.REGISTRY)
    baseline = s34.decode_intermediate_state(q['initial_baseline'])
    initial = GeneralT3ScalarRGState(baseline, {
        name: checked_complex(raw, name)
        for name, raw in q['initial_new_quartics'].items()
    }).validated(registry)
    y0, layout = pack_general_scalar_state(initial, registry)
    mu2_init = checked_complex(q['mu2_initial_gev2'], 'mu2_initial_gev2')
    y0 = np.concatenate((y0, [mu2_init.real, mu2_init.imag]))
    require(q['mu_initial_gev'] == s34.START and
            q['mu_final_gev'] == s34.STOP and
            set(q.get('beta_keys', [])) == set(layout.names) | {'muT3GeneralS1S2'},
            'quadratic EFT1 scale/coupling inventory changed')
    betas = qeft.checked_payload(model, baseline.representation, 'mixed', layout,
                                 registry, beta_path)
    # The original qeft runner recorded the hashed beta file path using the
    # native OS spelling. Match normalized paths, but do NOT waive the hash.
    hashes = q.get('input_sha256', {})
    normalized = {str(k).replace('\\', '/'): v for k, v in hashes.items()}
    key = str(beta_path).replace('\\', '/')
    require(key in normalized and normalized[key] == sha(beta_path),
            'saved qeft beta source differs from current live RGBeta beta')
    for file in (continuous.REGISTRY, continuous.MANIFEST):
        source_key = str(file).replace('\\', '/')
        require(source_key in normalized and normalized[source_key] == sha(file),
                'quadratic EFT1 registry/manifest source changed: ' + str(file))

    def rhs(t: float, z: np.ndarray) -> np.ndarray:
        return qeft.derivative(z, scale_gev=math.exp(t), layout=layout,
                               kind='mixed', registry=registry, report=betas)
    t_hi = math.log(s34.START)
    t_lo = math.log(s34.STOP)
    # Same equation, solver and tolerances as RunGeneralT3QuadraticEFT1.
    result = solve_ivp(rhs, (t_hi, t_lo), y0, method='DOP853',
                       rtol=1e-8, atol=1e-11, dense_output=True)
    require(result.success and result.sol is not None and
            np.isfinite(result.y).all(), 'candidate quadratic EFT1 replay failed')
    final = unpack_general_scalar_state(result.y[:layout.size, -1],
                                        layout, registry, uv=False, mu_gev=s34.STOP)
    final_mass = np.array([final.baseline.mS1Sq, final.baseline.mS2Sq], dtype=float)
    stored_mass = np.array([q['final_baseline']['mS1Sq'], q['final_baseline']['mS2Sq']], dtype=float)
    require(np.allclose(final_mass, stored_mass, rtol=4e-7, atol=2e-3),
            'EFT1 running scalar mass endpoint mismatch')
    final_mu2 = complex(float(result.y[layout.size, -1]),
                        float(result.y[layout.size+1, -1]))
    require(abs(final_mu2 - checked_complex(q['mu2_final_gev2'], 'mu2_final')) <
            4e-7 * max(1., abs(final_mu2)), 'EFT1 running holomorphic mass endpoint mismatch')
    require(set(final.extra) == set(q['final_new_quartics']),
            'new quartic endpoint keys mismatch')
    for name, value in final.extra.items():
        z = checked_complex(q['final_new_quartics'][name], name)
        require(abs(value-z) < 4e-7 * max(1., abs(z)),
                'new quartic endpoint mismatch: ' + name)

    def basis_at(t: float) -> MixedDoubletMassBasis:
        z = np.asarray(result.sol(t), dtype=float)
        st = unpack_general_scalar_state(z[:layout.size], layout, registry,
                                         uv=False, mu_gev=math.exp(t))
        mu2 = complex(z[layout.size], z[layout.size+1])
        b = st.baseline
        return MixedDoubletMassBasis.build(model, float(b.mS1Sq),
                                          float(b.mS2Sq), mu2)

    endpoint = basis_at(t_lo)
    masses = np.sqrt(endpoint.eigenvalues_sq)
    require(np.allclose(masses, ref.get('scalar_masses_gev', []),
                        rtol=1e-8, atol=1e-5) and
            abs(endpoint.mu_sq-checked_complex(ref['quadratic_mu2_gev2'], 'reference mu2')) < 1e-7,
            'saved one-scale benchmark differs from reconstructed mass-eigenstate endpoint')
    eigen = read(eigen_path)
    prepared = read(prepared_path)
    require(prepared.get('model') == model and prepared.get('endpoint') == 'final' and
            prepared.get('matchete_mass_basis_verified') is False,
            'prepared Matchete mass source invalid')
    alpha1, alpha2 = check_native_mixing(model, eigen, prepared, endpoint)
    t_y = exact_scale_crossing(basis_at, t_lo, t_hi, 1)
    t_x = exact_scale_crossing(basis_at, t_lo, t_hi, 0)
    require(t_hi > t_y > t_x > t_lo,
            'F>Y>X>scalar-endpoint hierarchy not realized')
    by, bx = basis_at(t_y), basis_at(t_x)
    y_boundary = boundary(model, t_y, by, alpha1, alpha2, 'Y')
    x_estimate = boundary(model, t_x, bx, alpha1, alpha2, 'X')
    py = np.outer(by.U[:, 1], by.U[:, 1].conj())
    px = np.outer(by.U[:, 0], by.U[:, 0].conj())
    require(np.allclose(py+px, np.eye(2), atol=2e-12) and
            np.allclose(px@py, np.zeros((2, 2)), atol=2e-12),
            'X/Y projectors do not resolve the scalar flavour space at muY')

    paths = {'quadratic_eft1': source_path, 'rgbeta_eft1_betas': beta_path,
             'live_matchete_eigenbasis': eigen_path,
             'prepared_mass_eigenstate': prepared_path,
             'conditional_one_scale_reference': reference_path}
    payload = {
        'format': FORMAT, 'status': STATUS, 'model': model,
        'input_sha256': {k: sha(v) for k, v in paths.items()},
        'source_stage': 'un-decoupled candidate full-quartic EFT1 RGEs',
        'numerical_rge_solver': 'DOP853, same 1e-8/1e-11 tolerances as quadratic EFT1 source',
        'eft1_nfev': int(result.nfev),
        'rge_interval_gev': [float(s34.START), float(s34.STOP)],
        'fermion_threshold_reference_gev': float(s34.START),
        'running_Y_threshold': y_boundary,
        'running_X_threshold_un_decoupled_estimate': x_estimate,
        'log_muY_over_muX_full_eft1_estimate': t_y - t_x,
        'X_Y_projectors_at_Y_threshold_complete': True,
        'quadratic_eft1_endpoint_replayed': True,
        'one_scale_reference_mass_basis_parity': True,
        'one_scale_reference_neutrino_masses_ev_descending':
            ref['candidate_takagi_masses_ev_descending'],
        'one_scale_neutrino_masses_recomputed': False,
        'full_theory_running_above_Y_completed': True,
        'hard_fermion_Wilson_boundary_evaluated_here': False,
        'Y_threshold_Wilson_matching_performed': False,
        'X_only_EFT2_beta_functions_available_here': False,
        'X_threshold_Wilson_matching_performed': False,
        'finite_scalar_threshold_matching_performed': False,
        'physical_sequential_C5_certified': False,
        'physical_neutrino_mass_prediction': False,
        'notes': (
            'The running Y crossing is obtained within EFT1, above the first scalar '
            'decoupling. The X crossing continues the same two-scalar beta functions '
            'below Y and is a PROVISIONAL scale estimate, not a matched EFT2 result. '
            'No one-loop scalar matching or Wilson evolution across Y/X is inferred.'
        ),
    }
    out_dir.mkdir(parents=True, exist_ok=True)
    target = out_dir / f'{tag}_separated_thresholds.json'
    require(not target.exists(), 'output exists; choose a new --run-id')
    tmp = target.with_name(target.name + '.tmp')
    tmp.write_text(json.dumps(payload, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    tmp.replace(target)
    return payload


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--all', action='store_true')
    group.add_argument('--model', choices=MODELS)
    parser.add_argument('--run-id', required=True, type=continuous.run_id)
    parser.add_argument('--eft1-dir', type=Path, default=DEFAULT_EFT1)
    parser.add_argument('--eigen-dir', type=Path, default=DEFAULT_EIGEN)
    parser.add_argument('--one-scale-dir', type=Path, default=DEFAULT_ONESCALE)
    parser.add_argument('--output-root', type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    out_dir = args.output_root / args.run_id
    errors = []
    for model in MODELS if args.all else (args.model,):
        try:
            payload = solve(model, input_dir=args.eft1_dir,
                            eigen_dir=args.eigen_dir, reference_dir=args.one_scale_dir,
                            out_dir=out_dir)
            my = payload['running_Y_threshold']['mu_gev']
            mx = payload['running_X_threshold_un_decoupled_estimate']['mu_gev']
            print(f'{model}: running muY={my:.9g} GeV; '
                  f'un-decoupled muX estimate={mx:.9g} GeV; '
                  f'log(Y/X)={payload["log_muY_over_muX_full_eft1_estimate"]:.6g}',
                  flush=True)
        except (ValueError, KeyError, TypeError, RuntimeError, OSError,
                ArithmeticError, np.linalg.LinAlgError) as exc:
            errors.append(model)
            print(f'{model}: BLOCKED {type(exc).__name__}: {exc}', flush=True)
    if errors:
        return 2
    print('Y-threshold running and X estimate complete; '
          'finite Y/X Wilson matching and physical C5 NOT performed.', flush=True)
    print('Results:', out_dir)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
