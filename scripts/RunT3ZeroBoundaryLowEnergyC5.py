"""Batch 8: B/C zero-extra UV benchmark -> existing Stage36 low-energy C5.

Replays the verified Batch7 full native RGBeta UV/EFT1 one-loop ODE from an
exactly zero additional-coupling UV boundary, then evaluates each model's
*existing* hierarchical fixed-order C5 matching artifact on the resulting
baseline heavy/scalar numerical states and evolves SM+Weinberg to 100 GeV.

This is a controlled CONDITIONAL comparison to saved Stage36 (whose UV
additional-quartic boundary differs). It is not new general-scalar matching,
and cannot isolate the effect of a radiatively generated extra C5 quartic.
Neither Matchete/RGBeta coupling map nor physical neutrino fit is certified.

    python -m scripts.RunT3ZeroBoundaryLowEnergyC5 --all
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np

from scripts import RunT3ZeroExtraBoundaryBenchmark as b7

ROOT = b7.ROOT
B7_DIR = b7.OUT
STAGE36_DIR = ROOT / 'output/general_t3_stage22'
OUT = ROOT / 'output/scalar_relevance/zero_boundary_low_energy'
MODELS = b7.MODELS
FORMAT = 't3-bc-zero-boundary-existing-hierarchical-c5-low-energy-v1'
STATUS = 'ConditionalZeroBoundaryExistingHierarchicalC5Run'


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError('Batch 8 low-energy C5: ' + message)


def safe_tag(model: str) -> str:
    require(model in MODELS, 'unsupported B/C model')
    return model.replace('-', 'm').replace('+', 'p')


def _positive(value: Any, label: str) -> float:
    try:
        v = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f'Batch 8 low-energy C5: {label} is not numeric') from exc
    require(math.isfinite(v) and v > 0, label + ' must be positive finite')
    return v


def validate_saved_batch7(report: dict, model: str, *, root: Path = ROOT) -> None:
    """Require a real local Batch7 run, with intact SHA256 source provenance."""
    require(model in MODELS, 'unsupported model')
    require(report.get('format') == 't3-zero-extra-boundary-bc-renormalisable-running-v1'
            and report.get('status') == 'CandidateBCZeroBoundaryTwoSegmentRunning'
            and report.get('model') == model,
            model + ': missing/wrong Batch7 result')
    require(report.get('boundary_scale_gev') == 1e7
            and report.get('fermion_threshold_gev') == 1e5
            and report.get('scalar_threshold_gev') == 1e3
            and report.get('boundary_all_additional_couplings_zero') is True
            and report.get('uv_extra_quartic_slots') == 7
            and report.get('uv_quadratic_slots') == 1
            and report.get('physical_matchete_rgbeta_map_certified') is False
            and report.get('full_sequential_C5_computed') is False
            and report.get('physical_neutrino_mass_error_bound') is False,
            model + ': Batch7 zero-boundary/scope flags inconsistent')
    hashes = report.get('source_sha256')
    require(isinstance(hashes, dict) and len(hashes) >= 5,
            model + ': Batch7 source hashes missing')
    root = Path(root).resolve()
    for relative, recorded in hashes.items():
        require(isinstance(relative, str) and isinstance(recorded, str)
                and len(recorded) == 64,
                model + ': malformed Batch7 source hash')
        relative = relative.replace('\\', '/')
        p = (root / relative).resolve()
        require(p.is_relative_to(root) and p.is_file(),
                model + ': missing or unsafe Batch7 source ' + relative)
        require(sha256(p.read_bytes()).hexdigest() == recorded,
                model + ': Batch7 source has changed: ' + relative)
    for key in ('uv_initial', 'uv_final', 'eft1_initial', 'eft1_final'):
        require(isinstance(report.get(key), dict), model + ': missing ' + key)
    start = report['uv_initial']
    reg = b7.read_json(root / 'RGE/general/generated/general_rgbeta_quartic_registry.json')
    entry = reg.get('models', {}).get(model)
    require(isinstance(entry, dict), model + ': registry missing')
    new = tuple(q['name'] for q in entry['new_rgbeta_symbols'])
    reused = tuple(q['existing_name'] for q in entry['preexisting_rgbeta_symbols'])
    require(len(new) == 1 and len(reused) == 6,
            model + ': unexpected quartic inventory')
    b7.check_zero_boundary(start, new=new, reused=reused)


def _near_number(a: Any, b: Any, label: str) -> None:
    av = float(a); bv = float(b)
    require(math.isfinite(av) and math.isfinite(bv)
            and math.isclose(av, bv, rel_tol=3e-8, abs_tol=2e-10),
            'Batch7 replay differs: ' + label)


def validate_replay(saved: dict, replay: dict) -> None:
    """Compare FULL saved Batch7 snapshots with the reconstructed endpoints."""
    for label in ('uv_initial', 'uv_final', 'eft1_initial', 'eft1_final'):
        expected = saved[label]
        actual = replay[label]
        for key in ('mu_gev', 'mS1Sq_gev2', 'mS2Sq_gev2'):
            _near_number(actual[key], expected[key], f'{label}/{key}')
        for key in ('lambdaT3', 'quadratic_mu2_gev2'):
            for part in ('real', 'imag'):
                _near_number(actual[key][part], expected[key][part], f'{label}/{key}.{part}')
        for key in ('new_quartics', 'reused_quartics'):
            require(set(actual[key]) == set(expected[key]),
                    label + '/' + key + ' registry changed')
            for coupling in expected[key]:
                for part in ('real', 'imag'):
                    _near_number(actual[key][coupling][part], expected[key][coupling][part],
                                 f'{label}/{coupling}.{part}')


def replay_full_states(model: str, saved: dict, *, n_points: int = 20):
    """Use Batch7 checked integrator; return UV(MF) and EFT1(MS) full states."""
    from scripts import CoupledUVRunning as s31
    from scripts import CoupledA0UVRunning as s30
    from Numerical.core.StateVector import unpack_uv_state
    from Numerical.experimental.GeneralT3ScalarRGState import (
        load_registry, GeneralT3ScalarRGState, pack_general_scalar_state,
        unpack_general_scalar_state, project_general_after_fermion,
    )

    require(model in MODELS, 'invalid model')
    registry = load_registry(b7.REGISTRY)
    entry, new, reused = s31.inventories(registry, model)
    require(len(new) == 1 and len(reused) == 6, 'wrong B/C quartic basis')
    config = b7.read_json(s31.config_path(entry))
    vec_base, layout_base = s31.prepare_state(config, entry, reused,
                                               {n: 0j for n in reused})
    base = unpack_uv_state(vec_base, layout_base, mu_gev=s30.START)
    initial = GeneralT3ScalarRGState(base, {n: 0j for n in new}).validated(registry)
    vec, layout = pack_general_scalar_state(initial, registry)
    y0 = np.concatenate((vec, np.zeros(2)))

    uv_data = b7.read_json(b7.beta_path(model, 'UV'))
    uv_betas = b7.check_payload(uv_data, model=model, stage='UV',
                               rep=base.representation,
                               base_names=tuple(layout.names[:-len(new)]),
                               new_names=new, reused_names=reused)
    original_file = s31.DATA / f'stage31_{safe_tag(model)}_full_uv_betas.json'
    original = b7.read_json(original_file)
    old_betas = s31.validate_export(original, entry, new, reused,
                                    layout_base.names)
    old_rhs = s30.combined_rhs(math.log(1e7), y0[:-2], layout=layout_base,
                               names=new, betas=old_betas)
    fresh_rhs = b7.rhs_at(math.log(1e7), y0, layout=layout,
                          registry=registry, betas=uv_betas, stage='UV')[:-2]
    require(np.max(np.abs(old_rhs - fresh_rhs) / np.maximum(1.0, np.abs(old_rhs))) < 1e-9,
            'quadratic UV beta system not consistent with independent Stage31')
    uv_rhs = lambda t, y: b7.rhs_at(t, y, layout=layout, registry=registry,
                                   betas=uv_betas, stage='UV')
    uv = b7.integrate_segment(y0, uv_rhs, 1e7, 1e5, n_points=n_points)
    uv_state = unpack_general_scalar_state(uv.y[:layout.size, -1], layout,
                                            registry, uv=True, mu_gev=1e5)
    boundary = project_general_after_fermion(uv_state, registry)
    y1, eft_layout = pack_general_scalar_state(boundary, registry)
    y1 = np.concatenate((y1, uv.y[layout.size:, -1]))
    eft_data = b7.read_json(b7.beta_path(model, 'EFT1'))
    eft_betas = b7.check_payload(eft_data, model=model, stage='EFT1',
                                rep=boundary.baseline.representation,
                                base_names=tuple(eft_layout.names[:-len(new)]),
                                new_names=new, reused_names=reused)
    eft_rhs = lambda t, y: b7.rhs_at(t, y, layout=eft_layout, registry=registry,
                                    betas=eft_betas, stage='EFT1')
    eft = b7.integrate_segment(y1, eft_rhs, 1e5, 1e3, n_points=n_points)
    final = unpack_general_scalar_state(eft.y[:eft_layout.size, -1],
                                         eft_layout, registry, uv=False, mu_gev=1e3)
    snapshots = {
        'uv_initial': b7.snapshot(uv.y[:, 0], layout, registry, stage='UV',
                                  mu_gev=1e7, new=new, reused=reused),
        'uv_final': b7.snapshot(uv.y[:, -1], layout, registry, stage='UV',
                                mu_gev=1e5, new=new, reused=reused),
        'eft1_initial': b7.snapshot(eft.y[:, 0], eft_layout, registry, stage='EFT1',
                                    mu_gev=1e5, new=new, reused=reused),
        'eft1_final': b7.snapshot(eft.y[:, -1], eft_layout, registry, stage='EFT1',
                                  mu_gev=1e3, new=new, reused=reused),
    }
    validate_replay(saved, snapshots)
    return uv_state.baseline, final.baseline, snapshots, int(uv.nfev), int(eft.nfev)


def verify_prior_stage36(prior: dict, model: str, artifact: Path,
                         *, root: Path = ROOT) -> None:
    """Stage36 is a benchmark comparison, never an arbitrary C5 matrix input."""
    from scripts import FinalWeinbergRunning as s36
    require(prior.get('format') == s36.FORMAT and prior.get('status') == s36.STATUS
            and prior.get('model') == model
            and prior.get('source_existing_c5_used') is True
            and prior.get('matching_source_is_hierarchical') is True
            and prior.get('new_general_quartics_included_in_matching') is False
            and prior.get('sm_weinberg_running_performed') is True
            and prior.get('physical_coupling_mapping_certified') is False
            and prior.get('mu_scalar_threshold_gev') == 1e3
            and prior.get('mu_low_gev') == 100.0,
            'saved Stage36 comparison has unsupported provenance or scope')
    require(prior.get('source_matching_artifact_sha256') == b7.digest(artifact),
            'Stage36 source C5 artifact has changed')
    location = (Path(root) / prior.get('source_matching_artifact', '')).resolve()
    require(location == artifact.resolve(),
            'Stage36 model-specific source path is inconsistent')


def _matrix_from_json(value: dict, label: str) -> np.ndarray:
    require(isinstance(value, dict) and set(value) == {'real', 'imag'},
            label + ': missing complex matrix')
    r = np.asarray(value['real'], dtype=float)
    i = np.asarray(value['imag'], dtype=float)
    require(r.shape == (3, 3) and i.shape == (3, 3)
            and np.isfinite(r).all() and np.isfinite(i).all(),
            label + ': invalid matrix')
    mat = r + 1j * i
    require(np.allclose(mat, mat.T, rtol=1e-10, atol=1e-20),
            label + ': nonsymmetric')
    return mat


def fractional_matrix_change(reference: np.ndarray, candidate: np.ndarray) -> float:
    reference = np.asarray(reference, dtype=complex)
    candidate = np.asarray(candidate, dtype=complex)
    require(reference.shape == candidate.shape == (3, 3), 'C5 matrix shapes must be 3x3')
    den = float(np.linalg.norm(reference, ord='fro'))
    require(den > 0 and math.isfinite(den), 'zero reference C5 Frobenius norm')
    value = float(np.linalg.norm(candidate - reference, ord='fro') / den)
    require(math.isfinite(value), 'nonfinite C5 comparison')
    return value


def run_model(model: str, *, batch7_dir: Path = B7_DIR,
              stage36_dir: Path = STAGE36_DIR, outdir: Path = OUT) -> dict:
    from scripts import FinalWeinbergRunning as s36
    from Numerical.running.FinalC5TrajectoryAdapter import extract_final_c5_inputs
    from Numerical.running.SMWeinbergStage import evaluate_final_weinberg_json
    from Numerical.running.ScalarThresholdBoundary import (
        build_final_sm_boundary, build_sm_weinberg_initial_conditions)
    from Numerical.running.SMWeinbergEvolution import evolve_sm_weinberg
    from physics.NeutrinoMass import numerical_neutrino_mass_matrix
    from physics.NeutrinoTrajectory import charged_lepton_mass_basis_matrix
    from physics.NeutrinoObservables import calculate_neutrino_observables

    tag = safe_tag(model)
    saved_path = Path(batch7_dir) / tag / 'result.json'
    saved = b7.read_json(saved_path)
    validate_saved_batch7(saved, model)
    artifact, kind = s36.source_path(model, root=ROOT, allow_legacy=False)
    require(kind == 'hierarchical-fixed-order-hard-plus-direct',
            'expected model-specific hierarchical C5 artifact')
    s36.verify_hierarchical_artifact(artifact)
    prior_path = Path(stage36_dir) / f'stage36_{tag}_low_energy_neutrinos.json'
    prior = b7.read_json(prior_path)
    verify_prior_stage36(prior, model, artifact)
    uv_state, eft_state, snapshots, uv_nfev, eft_nfev = replay_full_states(model, saved)
    for state in (uv_state, eft_state):
        require(state.representation.alpha == -1, 'incorrect model hypercharge')
    ratio = saved['generation']['quadratic_fraction_of_min_diagonal_m2']
    require(math.isfinite(float(ratio)) and float(ratio) <= 1e-11,
            'nonzero scalar quadratic mixing needs a new mass-basis C5 treatment')
    require(min(eft_state.mS1Sq, eft_state.mS2Sq) > 0,
            'non-positive scalar mass squared')
    view = SimpleNamespace(uv_threshold_state=uv_state,
                           eft1_threshold_state=eft_state,
                           mu_scalar_threshold_gev=s36.MU_S)
    matching = extract_final_c5_inputs(view)
    c5_muS = s36._valid_c5(evaluate_final_weinberg_json(artifact, matching.as_config()),
                           'zero-extra existing hierarchical C5')
    sm = build_final_sm_boundary(eft_state, matching_scale_gev=s36.MU_S)
    initials = build_sm_weinberg_initial_conditions(sm, c5_muS)
    low = evolve_sm_weinberg(initials, s36.MU_S, s36.MU_LOW,
                             rtol=1e-8, atol=1e-11)
    require(low.solver_success and low.mu_initial == s36.MU_S
            and low.mu_final == s36.MU_LOW,
            'SM+Weinberg evolution did not reach the low scale')
    c5_low = s36._valid_c5(low.K, 'zero-extra low-scale C5')
    mass = numerical_neutrino_mass_matrix(c5_low, vev_gev=s36.VEV)
    charged = charged_lepton_mass_basis_matrix(mass, low.ye)
    obs = calculate_neutrino_observables(charged, ordering='NO')
    require(np.isfinite(obs.masses_ev).all() and obs.takagi_residual < 1e-7,
            'invalid Takagi masses')
    ref_s = _matrix_from_json(prior['c5_at_scalar_threshold'], 'Stage36 C5 at 1 TeV')
    ref_low = _matrix_from_json(prior['c5_at_low_scale'], 'Stage36 C5 at 100 GeV')
    delta_s = fractional_matrix_change(ref_s, c5_muS)
    delta_low = fractional_matrix_change(ref_low, c5_low)
    ref_masses = np.asarray(prior['neutrino_observables']['masses_ev'], dtype=float)
    require(ref_masses.shape == (3,) and np.isfinite(ref_masses).all()
            and ref_masses[-1] > 0, 'Stage36 reference masses malformed')
    require(np.all(np.diff(obs.masses_ev) >= 0), 'new masses not in ascending order')
    model_out = Path(outdir) / tag
    result = {
        'format': FORMAT, 'status': STATUS, 'model': model,
        'uv_boundary_scale_gev': 1e7, 'mu_f_gev': 1e5,
        'mu_s_gev': 1e3, 'mu_low_gev': 100.0,
        'boundary_all_additional_scalar_couplings_zero': True,
        'candidate_uv_and_eft1_native_rgbeta_running_performed': True,
        'quadratic_mixing_negligible_checked': True,
        'matching_source': 'existing-model-specific-hierarchical-fixed-order-hard-plus-direct',
        'existing_matched_C5_reused_on_new_numerical_background': True,
        'new_general_scalar_C5_matching_performed': False,
        'RG_beta_matchete_physical_mapping_certified': False,
        'full_sequential_wilson_EFT1_rge_performed': False,
        'new_numerical_fitting_performed': False,
        'UV_solver_nfev': uv_nfev, 'EFT1_solver_nfev': eft_nfev,
        'SM_weinberg_solver_nfev': int(low.nfev),
        'matching_input_heavy_masses_gev': matching.heavy_masses_gev.tolist(),
        'matching_input_scalar_masses_gev': [matching.scalar_mass1_gev,
                                              matching.scalar_mass2_gev],
        'matching_input_lambdaT3': b7.complex_json(matching.lambda_t3),
        'c5_at_scalar_threshold': s36.matrix_json(c5_muS),
        'c5_at_low_scale': s36.matrix_json(c5_low),
        'neutrino_mass_matrix_gev': s36.matrix_json(mass),
        'neutrino_observables': {
            'ordering': 'NO', 'masses_ev': obs.masses_ev.tolist(),
            'delta_m21_sq_ev2': float(obs.delta_m21_sq_ev2),
            'delta_m31_sq_ev2': float(obs.delta_m31_sq_ev2),
            'pmns_abs': np.asarray(obs.pmns_abs).tolist(),
            'takagi_residual': float(obs.takagi_residual),
        },
        'comparison_to_saved_stage36_nonzero_stage29_boundary': {
            'fractional_C5_fro_change_at_muS': delta_s,
            'fractional_C5_fro_change_at_100gev': delta_low,
            'previous_masses_ev': ref_masses.tolist(),
            'new_masses_ev': obs.masses_ev.tolist(),
            'largest_mass_ratio_new_over_old': float(obs.masses_ev[-1]/ref_masses[-1]),
            'comparison_not_direct_radiative_extra_C5_error': True,
            'the_UV_additional_quartic_boundaries_differ': True,
        },
        'source_sha256': {
            b7.relpath(saved_path): b7.digest(saved_path),
            b7.relpath(prior_path): b7.digest(prior_path),
            b7.relpath(artifact): b7.digest(artifact),
        },
        'limitations': [
            'Two-stage one-loop RG of candidate RGBeta renormalisable interactions only',
            'Uses existing hierarchical fixed-order C5: no new general-scalar Wilson matching',
            'Scalar quadratic mixing near zero, not solved in physical Matchete mass basis',
            'Stage36 reference used nonzero Stage29 additional-quartic UV boundary',
            'Difference to Stage36 includes changed renormalisable RGE backgrounds; not solely direct extra C5',
            'Native RGBeta to physical Matchete quartic normalisation uncertified',
            'No C12 Wilson evolution, finite threshold matching or refitting',
        ],
    }
    b7.save_json(model_out / 'result.json', result)
    return result


def make_summary_md(rows: list[dict]) -> str:
    lines = [
        '# T3 zero-extra-boundary low-energy neutrinos — Batch 8', '',
        'B(-1)/C(-1), candidate native RGBeta UV+EFT1 10^7 -> 10^5 -> 10^3 GeV,',
        'then existing hierarchical fixed-order C5 and SM+Weinberg 10^3 -> 100 GeV.',
        'All eight model-specific extra couplings (7 quartics + 1 quadratic) are zero at 10^7 GeV.',
        'Matching is **reused**, not recomputed; physical coupling mapping remains uncertified.', '',
        '| Model | m1 (eV) | m2 (eV) | m3 (eV) | Delta m21^2 (eV^2) | Delta m31^2 (eV^2) | C5 change vs Stage36 @ 100 GeV |',
        '|---|---:|---:|---:|---:|---:|---:|',
    ]
    for row in rows:
        nu = row['neutrino_observables']
        a, b, c = nu['masses_ev']
        cmp = row['comparison_to_saved_stage36_nonzero_stage29_boundary']
        lines.append(f"| {row['model']} | {a:.8g} | {b:.8g} | {c:.8g} | "
                     f"{nu['delta_m21_sq_ev2']:.8g} | {nu['delta_m31_sq_ev2']:.8g} | "
                     f"{cmp['fractional_C5_fro_change_at_100gev']:.6g} |")
    lines += [
        '', 'The Stage36 comparison is **between distinct UV boundary scenarios**, not a physical',
        'estimate of the radiatively generated extra-C5 insertion. The latter remains uncertified.',
        'These neutrino spectra reuse the fitted-reference Yukawa starting point; they are not independent fits.',
        'No production files or GitHub files have been modified.', ''
    ]
    return '\n'.join(lines)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    select = p.add_mutually_exclusive_group(required=True)
    select.add_argument('--all', action='store_true')
    select.add_argument('--model', choices=MODELS)
    p.add_argument('--batch7-dir', type=Path, default=B7_DIR)
    p.add_argument('--stage36-dir', type=Path, default=STAGE36_DIR)
    p.add_argument('--outdir', type=Path, default=OUT)
    options = p.parse_args()
    rows = []
    for model in (MODELS if options.all else (options.model,)):
        item = run_model(model, batch7_dir=options.batch7_dir,
                         stage36_dir=options.stage36_dir, outdir=options.outdir)
        rows.append(item)
        masses = item['neutrino_observables']['masses_ev']
        delta = item['comparison_to_saved_stage36_nonzero_stage29_boundary']['fractional_C5_fro_change_at_100gev']
        print(f'{model}: neutrino masses (eV) = {masses}; '
              f'C5 boundary-scenario difference vs Stage36 = {delta:.6g}; '
              'not a direct extra-quartic C5 correction', flush=True)
    dest = Path(options.outdir)
    dest.mkdir(parents=True, exist_ok=True)
    summary = {
        'format': 't3-bc-zero-boundary-low-energy-summary-v1',
        'status': STATUS, 'models': [r['model'] for r in rows],
        'source_zero_extra_boundary_1e7_gev': True,
        'new_wilson_matching_performed': False,
        'neutrino_results_are_conditional': True,
        'model_outputs': {r['model']: {
            'path': f"{safe_tag(r['model'])}/result.json",
            'sha256': b7.digest(dest / safe_tag(r['model']) / 'result.json')
        } for r in rows}
    }
    b7.save_json(dest / 'summary.json', summary)
    (dest / 'summary.md').write_text(make_summary_md(rows), encoding='utf-8')
    print('Batch 8 output:', dest / 'summary.md')


if __name__ == '__main__':
    main()
