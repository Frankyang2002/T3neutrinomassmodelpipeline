"""Continuous, *candidate-basis* 16-model T3 scalar UV -> EFT1 evolution.

Integrates the existing complete one-loop RGBeta UV/EFT1 quartic and quadratic
beta exports in the SAME state, with a tree-level fermion-threshold handoff.
For B_-1/C_-1, optionally evaluates the EXISTING direct all-heavy one-scale
Matchete C5 at the resulting scalar endpoint (with frozen-reference Yukawas).
This is NOT a physical sequential Wilson matching calculation.

Usage from repository root:
    python -m scripts.RunT3FullScalarContinuous --all --run-id continuous_v1 --mu2-fraction 0.01
    python -m scripts.RunT3FullScalarContinuous --all --run-id continuous_v1 --mu2-fraction 0.01 --resume
    python -m scripts.RunT3FullScalarContinuous --model B_-1 --run-id single_v1 --mu2-fraction 0.01

No historical inputs, UV reports, or fitted benchmark files are overwritten.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
import re
from typing import Any, Mapping

import numpy as np
from scipy.integrate import solve_ivp

from scripts import CoupledA0UVRunning as s30
from scripts import CoupledUVRunning as s31
from scripts import FermionThresholdHandoff as s33
from scripts import IntermediateEFTRunning as s34
from scripts import RunGeneralT3QuadraticEFT1 as qeft
from Numerical.experimental.GeneralT3ScalarRGState import (
    GeneralT3ScalarRGState, load_registry, pack_general_scalar_state,
    unpack_general_scalar_state,
)
from Numerical.experimental.GeneralT3ScalarQuadraticMasses import (
    load_manifest, quadratic_kind, quadratic_mass_hessian,
)
from Numerical.core.StateVector import unpack_uv_state
from Numerical.core.RGBetaEvaluator import evaluate_inputform_expression, state_environment
from Numerical.core.BetaVector import beta_values_to_derivative
from Numerical.running.IntermediateScalarRGBetaEvaluator import intermediate_scalar_state_environment
from Numerical.running.FermionThresholdBoundary import (
    fermion_threshold_masses, build_intermediate_scalar_boundary,
)

ROOT = s30.ROOT
STAGE22 = ROOT / 'output/general_t3_stage22'
BETA_DIR = ROOT / 'output/general_t3_quadratic_eft1'
REGISTRY = ROOT / 'RGE/general/generated/general_rgbeta_quartic_registry.json'
MANIFEST = ROOT / 'docs/general_scalar_trial_operator_manifest.json'
OUTPUT = ROOT / 'output/studies/general_scalar_extended/continuous_full_scalar'
MODELS = tuple(s33.MODELS)
FORMAT = 't3-continuous-general-scalar-candidate-v1'
STATUS = 'CandidateContinuousUVFermionThresholdEFT1Integrated'
LOOP = 16.0 * math.pi**2
MODEL_MASS = {'self_s1': 'muT3GeneralS1Sq', 'self_s2': 'muT3GeneralS2Sq',
              'mixed': 'muT3GeneralS1S2', 'none': None}
EXPECTED_BETA_CONVENTION = '16*pi^2*dX/dln(mu); gauge BetaTerm divided by 2*g'


def require(ok: bool, message: str) -> None:
    if not ok:
        raise ValueError('Continuous T3 scalar pipeline: ' + message)


def digest(path: Path) -> str:
    require(path.is_file(), 'missing source: ' + str(path))
    return sha256(path.read_bytes()).hexdigest()


def read(path: Path) -> dict:
    require(path.is_file(), 'missing JSON input: ' + str(path))
    raw = json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(raw, dict), 'expected JSON object: ' + str(path))
    return raw


def write_json_atomic(path: Path, payload: Mapping[str, Any]) -> None:
    """Publish completed JSON only; interrupted writes leave the checkpoint intact."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + '.tmp')
    temporary.write_text(json.dumps(payload, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    temporary.replace(path)


def cjson(value: complex) -> dict[str, float]:
    z = complex(value)
    require(math.isfinite(z.real) and math.isfinite(z.imag), 'non-finite complex value')
    return {'real': float(z.real), 'imag': float(z.imag)}


def model_tag(model: str) -> str:
    require(model in MODELS, 'unknown model ' + str(model))
    return model.replace('+', 'p').replace('-', 'm')


def run_id(value: str) -> str:
    require(bool(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,48}', value)),
            'run-id must be 1..49 letters/digits/underscore/hyphen')
    return value


def beta_inventory(payload: Mapping[str, Any], *, model: str, stage: str,
                   representation: Any, layout: Any, registry: Mapping[str, Any],
                   kind: str) -> dict[str, str]:
    """Fail closed on coupling inventories and source one-loop conventions."""
    require(payload.get('format') == 'general-t3-quadratic-rgbeta-experimental-v1',
            model + '/' + stage + ': wrong beta export format')
    require(payload.get('status') == 'SuccessCandidateMassBeta' and
            payload.get('physical_mapping_certified') is False and
            payload.get('mass_eigenstate_threshold_matching_performed') is False,
            model + '/' + stage + ': source status/certification changed')
    meta = payload.get('metadata')
    require(isinstance(meta, Mapping), model + '/' + stage + ': metadata absent')
    entry = registry['models'][model]
    expected = {
        'GeneralScalarModel': model, 'GeneralScalarStage': stage,
        'GeneralScalarBasisStatus': 'CandidateNativeRGBetaUnmappedToMatchete',
        'GeneralScalarQuarticBetaBasis': 'RGBetaFullQuarticProjection',
        'GeneralScalarQuadraticMassMixing': 'IncludedCandidate',
        'GeneralScalarAdditionalQuartics': entry['total_additional_quartics'],
        'GeneralScalarExistingStructuralCount': len(entry['preexisting_rgbeta_symbols']),
        'QuadraticKind': kind, 'QuadraticMassSymbol': MODEL_MASS[kind],
        'QuadraticPrefactorConvention': 'V2 = mu2 * invariant + h.c.; no 1/2',
        'PhysicalMatcheteRGBetaMassBasisCertified': False,
        'ReportBetaConvention': EXPECTED_BETA_CONVENTION,
        'SharedScalar': False, 'dS1': representation.d_s1,
        'dS2': representation.d_s2, 'dF': representation.d_f,
        'alpha': representation.alpha,
    }
    for k, v in expected.items():
        require(meta.get(k) == v, f'{model}/{stage}: inconsistent {k}: {meta.get(k)!r} != {v!r}')
    new = {item['name'] for item in entry['new_rgbeta_symbols']}
    require(set(meta.get('GeneralScalarNewBetaKeys', [])) == new,
            model + '/' + stage + ': additional quartic beta basis differs')
    mass = MODEL_MASS[kind]
    expected_keys = set(layout.names) | ({mass} if mass else set())
    betas = payload.get('report_betas')
    require(isinstance(betas, dict) and set(betas) == expected_keys and
            set(payload.get('beta_symbols', [])) == expected_keys,
            model + '/' + stage + ': incomplete/mismatched full beta inventory; '
            f'missing={sorted(expected_keys - set(betas or {}))}')
    require(all(isinstance(expr, str) and expr.strip() and
                not any(bad in expr for bad in ('BetaTerm[', 'Missing[', '$Failed', '$Aborted'))
                for expr in betas.values()),
            model + '/' + stage + ': unresolved beta expression')
    return betas


def coupled_rhs(t: float, y: np.ndarray, *, layout: Any, registry: dict,
                stage: str, kind: str, betas: Mapping[str, str]) -> np.ndarray:
    """All base/new-quartic/holomorphic-mass beta functions backreact."""
    require(stage in ('UV', 'EFT1'), 'invalid running stage')
    size = layout.size
    require(y.shape == (size + (2 if kind != 'none' else 0),) and np.isfinite(y).all(),
            stage + ': invalid ODE state')
    mu = float(math.exp(t))
    full = unpack_general_scalar_state(y[:size], layout, registry,
                                       uv=(stage == 'UV'), mu_gev=mu)
    base = full.baseline
    env = (state_environment(base) if stage == 'UV'
           else intermediate_scalar_state_environment(base))
    env.update(full.extra)
    symbol = MODEL_MASS[kind]
    if symbol:
        env[symbol] = complex(float(y[size]), float(y[size+1]))
    values = {name: evaluate_inputform_expression(expr, env)
              for name, expr in betas.items()}
    if symbol:
        mass_beta = values.pop(symbol)
        require(np.ndim(mass_beta) == 0, 'mass beta is not a complex scalar')
        zm = complex(mass_beta)
        require(math.isfinite(zm.real) and math.isfinite(zm.imag), 'mass beta is nonfinite')
    out = beta_values_to_derivative(values, layout, divide_by_loop_factor=True)
    if symbol:
        out = np.concatenate((out, np.array([zm.real/LOOP, zm.imag/LOOP])))
    require(out.shape == y.shape and np.isfinite(out).all(), 'bad coupled derivative')
    return out


def integrate_segment(y0: np.ndarray, *, mu_start: float, mu_end: float,
                      layout: Any, registry: dict, stage: str, kind: str,
                      betas: Mapping[str, str], points: int,
                      rtol: float = 1e-8, atol: float = 1e-11):
    require(4 <= points <= 200, 'points must be between 4 and 200')
    require(math.isfinite(mu_start) and math.isfinite(mu_end) and
            mu_start > mu_end > 0, 'RG scales must decrease and be positive')
    def rhs(t, y):
        return coupled_rhs(t, y, layout=layout, registry=registry,
                           stage=stage, kind=kind, betas=betas)
    rhs(math.log(mu_start), np.asarray(y0, dtype=float))
    result = solve_ivp(rhs, (math.log(mu_start), math.log(mu_end)),
                       np.asarray(y0, dtype=float), method='DOP853',
                       t_eval=np.linspace(math.log(mu_start), math.log(mu_end), points),
                       rtol=rtol, atol=atol, max_step=0.35)
    require(result.success and result.y.shape == (len(y0), points) and
            np.isfinite(result.y).all(), stage + ': ODE failed: ' + result.message)
    return result


def initial_uv(model: str, registry: dict, stage22: Path):
    """Reuse exact Stage29 illustrative quartic boundaries, not guessed values."""
    stage29_path = stage22 / 'stage29_convention_uv_handoff.json'
    stage29 = read(stage29_path)
    s30.verify_stage29_registry_provenance(stage29)
    if model == s31.BASELINE:
        names = s30.quartic_names(registry)
        vals = s30.uv_quartics(stage29, registry)
        conf = s30.CONFIG
        st, _, _ = s30.benchmark_state(read(conf))
        reused = ()
    else:
        entry, names, reused = s31.inventories(registry, model)
        vals = s31.stage29_values(stage29, entry, names, reused)
        conf = s31.config_path(entry)
        base_vector, base_layout = s31.prepare_state(read(conf), entry, reused, vals)
        st = unpack_uv_state(base_vector, base_layout, mu_gev=s30.START).validated()
    require(tuple(sorted(names)) == tuple(sorted(x['name'] for x in
                    registry['models'][model]['new_rgbeta_symbols'])),
            'UV additional quartic inventory changed')
    full = GeneralT3ScalarRGState(st, {name: vals[name] for name in names}).validated(registry)
    return full, conf, stage29_path


def trajectory_rows(result: Any, *, layout: Any, registry: dict,
                    stage: str, kind: str) -> list[dict]:
    rows = []
    for i, t in enumerate(result.t):
        mu = math.exp(float(t))
        y = result.y[:, i]
        st = unpack_general_scalar_state(y[:layout.size], layout, registry,
                                         uv=(stage == 'UV'), mu_gev=mu)
        base = st.baseline
        item = {'mu_gev': mu, 'mS1Sq_gev2': float(base.mS1Sq),
                'mS2Sq_gev2': float(base.mS2Sq),
                'extra_quartics': {k: cjson(v) for k, v in st.extra.items()},
                'mu2_gev2': cjson(complex(y[layout.size], y[layout.size+1])
                                if kind != 'none' else 0j)}
        rows.append(item)
    return rows


def provenance(paths: Mapping[str, Path]) -> dict[str, str]:
    return {name: digest(path) for name, path in paths.items()}


def result_is_resumable(path: Path, *, expected_inputs: dict, expected_settings: dict) -> bool:
    if not path.is_file():
        return False
    old = read(path)
    require(old.get('format') == FORMAT and old.get('status') == STATUS,
            'existing output is not a successful continuous integration: ' + str(path))
    require(old.get('input_sha256') == expected_inputs and
            old.get('settings') == expected_settings,
            'checkpoint inputs/settings changed; choose a NEW --run-id')
    return True


def compute_model(model: str, *, stage22: Path, beta_dir: Path,
                  dest: Path, mu2_fraction: float, phase: float,
                  points: int, resume: bool, c5: bool) -> dict[str, Any]:
    require(model in MODELS, 'unsupported model')
    registry = load_registry(REGISTRY)
    manifest = load_manifest(MANIFEST)
    kind = quadratic_kind(manifest, model)
    full_uv, conf_path, stage29_path = initial_uv(model, registry, stage22)
    uv_beta_path = beta_dir / f'{model_tag(model)}_uv_quadratic_live_rgbeta.json'
    eft_beta_path = beta_dir / f'{model_tag(model)}_eft1_quadratic_live_rgbeta.json'
    common_paths = {'registry': REGISTRY, 'manifest': MANIFEST,
                    'stage29_boundary': stage29_path, 'fitted_config': conf_path,
                    'uv_mass_beta': uv_beta_path, 'eft1_mass_beta': eft_beta_path,
                    'quadratic_rgbeta_source': ROOT / 'RGE/running/rgbeta/RunGeneralT3QuadraticRGBeta.wl',
                    'continuous_runner': Path(__file__).resolve()}
    hashes = provenance(common_paths)
    settings = {'model': model, 'mu2_fraction': mu2_fraction,
                'mu2_phase_rad': phase, 'saved_points_per_segment': points,
                'uv_start_gev': s30.START, 'fermion_threshold_gev': s30.STOP,
                'scalar_endpoint_gev': s34.STOP}
    target = dest / 'continuous_trajectory.json'
    if resume and result_is_resumable(target, expected_inputs=hashes,
                                      expected_settings=settings):
        payload = read(target)
        if c5 and model in ('B_-1', 'C_-1'):
            run_c5(model, dest)
        return payload
    require(not target.exists(), 'output exists; use --resume or a new --run-id: ' + str(target))
    uv_y, uv_layout = pack_general_scalar_state(full_uv, registry)
    require(s30.START == full_uv.baseline.mu_gev and s30.STOP == s34.START,
            'the UV/EFT1 threshold scales do not meet')
    require(0 <= mu2_fraction <= 0.45 and math.isfinite(mu2_fraction),
            'mu2 fraction must be within [0,0.45]')
    require(math.isfinite(phase), 'nonfinite quadratic phase')
    mu2_initial = (mu2_fraction * min(full_uv.baseline.mS1Sq, full_uv.baseline.mS2Sq) *
                   complex(math.cos(phase), math.sin(phase)) if kind != 'none' else 0j)
    require(min(full_uv.baseline.mS1Sq, full_uv.baseline.mS2Sq) > 0,
            'nonpositive UV scalar mass input')
    if kind != 'none':
        uv_y = np.concatenate((uv_y, [mu2_initial.real, mu2_initial.imag]))
    uv_betas = beta_inventory(read(uv_beta_path), model=model, stage='UV',
               representation=full_uv.baseline.representation, layout=uv_layout,
               registry=registry, kind=kind)
    uv = integrate_segment(uv_y, mu_start=s30.START, mu_end=s30.STOP,
                           layout=uv_layout, registry=registry, stage='UV',
                           kind=kind, betas=uv_betas, points=points)
    uv_end = unpack_general_scalar_state(uv.y[:uv_layout.size, -1], uv_layout,
                                        registry, uv=True, mu_gev=s30.STOP)
    mu2_at_f = (complex(float(uv.y[uv_layout.size, -1]),
                        float(uv.y[uv_layout.size+1, -1])) if kind != 'none' else 0j)
    fermion_masses = fermion_threshold_masses(uv_end.baseline)
    intermediate = build_intermediate_scalar_boundary(
        uv_end.baseline, matching_scale_gev=s30.STOP)
    eft_initial = GeneralT3ScalarRGState(intermediate, uv_end.extra).validated(registry)
    eft_y, eft_layout = pack_general_scalar_state(eft_initial, registry)
    if kind != 'none':
        eft_y = np.concatenate((eft_y, [mu2_at_f.real, mu2_at_f.imag]))
    eft_betas = beta_inventory(read(eft_beta_path), model=model, stage='EFT1',
               representation=intermediate.representation, layout=eft_layout,
               registry=registry, kind=kind)
    eft = integrate_segment(eft_y, mu_start=s34.START, mu_end=s34.STOP,
                            layout=eft_layout, registry=registry, stage='EFT1',
                            kind=kind, betas=eft_betas, points=points)
    eft_end = unpack_general_scalar_state(eft.y[:eft_layout.size, -1], eft_layout,
                                         registry, uv=False, mu_gev=s34.STOP)
    mu2_final = (complex(float(eft.y[eft_layout.size, -1]),
                         float(eft.y[eft_layout.size+1, -1])) if kind != 'none' else 0j)
    spectrum = quadratic_mass_hessian(manifest, model, eft_end.baseline.mS1Sq,
                                      eft_end.baseline.mS2Sq, mu2_final)
    require(spectrum['analytic_spectrum_verified'] is True,
            'scalar eigenvalue check not verified')
    # Retain the numerical RG result even when a benchmark develops a tachyon;
    # the point is then not an admissible physical scalar spectrum.
    dest.mkdir(parents=True, exist_ok=True)
    # Preserve the established Stage37 schema for the existing B/C Matchete C5
    # numerical bridge. This output is a NEW continuous source, not Stage37 data.
    sidecar = dest / f'{model_tag(model)}_eft1_quadratic_running.json'
    endpoint_payload = {
        'format': qeft.FORMAT,
        'status': 'CandidateQuadraticAndFullQuarticEFT1Integrated',
        'model': model, 'quadratic_kind': kind,
        'mu_initial_gev': float(s34.START),
        'mu_final_gev': float(s34.STOP),
        'mu2_beta_included': kind != 'none',
        'mu2_initial_gev2': cjson(mu2_at_f),
        'mu2_final_gev2': cjson(mu2_final),
        'initial_baseline': s33.json_state(eft_initial.baseline),
        'final_baseline': s33.json_state(eft_end.baseline),
        'initial_new_quartics': {k: cjson(v) for k, v in eft_initial.extra.items()},
        'final_new_quartics': {k: cjson(v) for k, v in eft_end.extra.items()},
        'physical_matchete_rgbeta_mapping_certified': False,
        'scalar_mass_threshold_matching_performed': False,
        'wilson_coefficient_matching_performed': False,
        'physical_neutrino_mass_prediction': False,
        'continuous_uv_quadratic_running_performed': True,
        'continuous_uv_fermion_threshold_eft1': True,
        'input_sha256': hashes, 'settings': settings,
    }
    # Checkpoint after BOTH ODE segments are complete, before C5. Never overwrite
    # a pre-existing historical Stage37 report.
    if sidecar.exists():
        previous = read(sidecar)
        require(resume and previous.get('input_sha256') == hashes and
                previous.get('settings') == settings and
                previous.get('format') == qeft.FORMAT,
                'incomplete endpoint checkpoint exists with incompatible inputs; use a new --run-id')
    write_json_atomic(sidecar, endpoint_payload)
    payload = {
        'format': FORMAT, 'status': STATUS,
        'model': model, 'settings': settings, 'input_sha256': hashes,
        'quadratic_kind': kind, 'quadratic_mu2_uv_gev2': cjson(mu2_initial),
        'quadratic_mu2_at_fermion_threshold_gev2': cjson(mu2_at_f),
        'quadratic_mu2_at_scalar_endpoint_gev2': cjson(mu2_final),
        'new_extra_quartics_evolved': len(full_uv.extra),
        'reused_extra_quartics_evolved': len(registry['models'][model]['preexisting_rgbeta_symbols']),
        'extra_quadratic_terms_evolved': int(kind != 'none'),
        'uv_ode_evaluations': int(uv.nfev),
        'eft1_ode_evaluations': int(eft.nfev),
        'heavy_fermion_threshold_masses_gev': fermion_masses.masses_gev.tolist(),
        'uv_trajectory': trajectory_rows(uv, layout=uv_layout, registry=registry,
                                         stage='UV', kind=kind),
        'eft1_trajectory': trajectory_rows(eft, layout=eft_layout, registry=registry,
                                           stage='EFT1', kind=kind),
        'uv_fermion_threshold_state': s33.json_state(uv_end.baseline),
        'eft1_scalar_endpoint_state': s33.json_state(eft_end.baseline),
        'scalar_endpoint_physically_admissible': not spectrum['tachyonic'],
        'scalar_endpoint_spectrum': {k: v for k, v in spectrum.items()
                                     if k != 'real_scalar_hessian_gev2'},
        'saved_candidate_eft1_endpoint': sidecar.name,
        'mass_bilinear_tree_threshold_prescription': 'continuous coefficient at fermion threshold; finite threshold corrections omitted',
        'physical_matchete_rgbeta_mapping_certified': False,
        'fermion_threshold_wilson_matching_performed': False,
        'scalar_threshold_wilson_matching_performed': False,
        'full_sequential_neutrino_mass_prediction': False,
        'limitations': [
            'Candidate RGBeta scalar-field/basis normalization not fully certified',
            'Additional UV quartic boundary is the existing illustrative Stage29 assignment',
            'Heavy-fermion threshold is renormalisable tree-level projection only',
            'No full fermion-first C12/C5 Wilson matching, running, or scalar-threshold matching',
            'B/C optional C5 is direct simultaneous all-heavy one-scale with frozen Yukawas',
            'No A/D/E self-bilinear mass-eigenfield Matchete C5 evaluator',
        ],
    }
    write_json_atomic(target, payload)
    if c5 and model in ('B_-1', 'C_-1'):
        run_c5(model, dest)
    return payload


def run_c5(model: str, dest: Path) -> dict:
    """Use pre-existing direct one-scale B/C Matchete expression, not substitute a guessed loop function."""
    from scripts.RunT3QuadraticEndpointC5 import (
        calculate_model, DEFAULT_SOURCE, DEFAULT_PARAMETERS, MANIFEST as FROZEN_MANIFEST,
    )
    result, path = calculate_model(
        model, which='final', running_dir=dest,
        source_dir=DEFAULT_SOURCE, parameters_dir=DEFAULT_PARAMETERS,
        frozen_dir=FROZEN_MANIFEST.parent, manifest_path=FROZEN_MANIFEST,
        output_dir=dest / 'one_scale_c5', include_candidate_extras=True)
    require(result.get('physical_neutrino_mass_prediction') is False and
            result.get('sequential_EFT1_wilson_running_performed') is False,
            'upstream one-scale routine incorrectly claims physical matching')
    return {'status': result['status'],
            'candidate_takagi_masses_ev_descending': result['candidate_takagi_masses_ev_descending'],
            'relative_extra_quartic_delta': result['relative_delta_mnu_frobenius'],
            'report': str(path.relative_to(ROOT))}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    select = parser.add_mutually_exclusive_group(required=True)
    select.add_argument('--all', action='store_true')
    select.add_argument('--model', choices=MODELS)
    parser.add_argument('--run-id', required=True, type=run_id)
    parser.add_argument('--stage22', type=Path, default=STAGE22)
    parser.add_argument('--quadratic-betas', type=Path, default=BETA_DIR)
    parser.add_argument('--output-root', type=Path, default=OUTPUT)
    parser.add_argument('--mu2-fraction', type=float, default=0.01)
    parser.add_argument('--mu2-phase', type=float, default=0.)
    parser.add_argument('--saved-points', type=int, default=16)
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--skip-c5', action='store_true')
    parser.add_argument('--strict-c5', action='store_true',
                        help='Stop if an optional B/C direct one-scale C5 calculation is unavailable.')
    args = parser.parse_args()
    base = args.output_root / args.run_id
    if base.exists() and not args.resume:
        parser.error('Run ID directory exists; use --resume or choose a new --run-id: ' + str(base))
    selected = MODELS if args.all else (args.model,)
    rows = []
    for model in selected:
        location = base / model_tag(model)
        print(f'{model}: starting coupled candidate UV+EFT1 integration', flush=True)
        try:
            result = compute_model(model, stage22=args.stage22,
                beta_dir=args.quadratic_betas, dest=location,
                mu2_fraction=args.mu2_fraction, phase=args.mu2_phase,
                points=args.saved_points, resume=args.resume,
                c5=False)
            masses_sq = result['scalar_endpoint_spectrum']['mass_squared_eigenvalues_gev2']
            row = {'model': model, 'status': result['status'],
                   'quadratic_kind': result['quadratic_kind'],
                   'extra_quartic_slots_evolved': result['new_extra_quartics_evolved'] + result['reused_extra_quartics_evolved'],
                   'extra_quadratic_slots_evolved': result['extra_quadratic_terms_evolved'],
                   'scalar_mass_squared_min_gev2': min(masses_sq),
                   'scalar_endpoint_tachyonic': min(masses_sq) < 0,
                   'c5_status': 'NotApplicableOneScaleExtractorNotImplemented',
                   'full_sequential_neutrino_mass_prediction': False}
            if model in ('B_-1', 'C_-1') and not args.skip_c5 and min(masses_sq) > 0:
                try:
                    c5 = run_c5(model, location)
                    row.update({'c5_status': c5['status'],
                                'one_scale_candidate_takagi_ev': c5['candidate_takagi_masses_ev_descending'],
                                'one_scale_relative_extra_delta': c5['relative_extra_quartic_delta']})
                except (ValueError, KeyError, RuntimeError, OSError, TypeError) as error:
                    if args.strict_c5:
                        raise
                    row.update({'c5_status': 'BlockedExistingOneScaleInputOrConvention',
                                'c5_error': str(error)})
            print(f'{model}: coupled UV+EFT1 OK; {result["quadratic_kind"]}, '
                  f'min(mS^2)={min(masses_sq):.6g} GeV^2; C5={row["c5_status"]}', flush=True)
        except (ValueError, KeyError, RuntimeError, OSError, TypeError, OverflowError) as error:
            row = {'model': model, 'status': 'Blocked', 'error': str(error),
                   'full_sequential_neutrino_mass_prediction': False}
            print(f'{model}: BLOCKED: {error}', flush=True)
        rows.append(row)
        base.mkdir(parents=True, exist_ok=True)
        summary = {'format': FORMAT, 'status': 'Incomplete' if any(r['status']=='Blocked' for r in rows) else 'CandidateRunningCompleted',
                   'models_requested': list(selected), 'models_processed': len(rows),
                   'models_succeeded': sum(r['status']==STATUS for r in rows),
                   'models_blocked': sum(r['status']=='Blocked' for r in rows),
                   'quartic_slots_completed': sum(r.get('extra_quartic_slots_evolved',0) for r in rows),
                   'quadratic_slots_completed': sum(r.get('extra_quadratic_slots_evolved',0) for r in rows),
                   'physical_full_scalar_sequential_matching_completed': False,
                   'rows': rows}
        write_json_atomic(base / 'summary.json', summary)
    print('Summary:', base / 'summary.json')
    return 2 if summary['models_blocked'] else 0


if __name__ == '__main__':
    raise SystemExit(main())
