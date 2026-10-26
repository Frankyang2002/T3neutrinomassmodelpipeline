"""Experimental grouped F-first B/C(-1) threshold trajectory, with honest Wilson boundary.

Implements the renormalisable UV -> scalar-only -> SM boundary using the
already-implemented RGBeta solvers, then the AUTHORITATIVE baseline hierarchical
C5 and SM+Weinberg running.  Extra quartics are evolved on both sides of F,
but the extra EFT1 Wilson matching/mixing has NOT been derived or implemented.
No extra-quartic C5 is inserted into the physical final-C5 boundary.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from Numerical.experimental.GeneralT3CoupledUV import (
    ALL_COMPLEX, SPECIAL, construct_uv_state,
)
RGBETA_BY_MAT = {
    "lambdaT3ExtraQ020020I1": "lambdaHHdagS2S2",
    "lambdaT3ExtraQ020200I1": "lambdaHHdagS1barS1bar",
}


@dataclass(frozen=True)
class ThresholdScales:
    mu_f: float
    mu_s: float
    mu_low: float


def _positive(label: str, value: Any) -> float:
    v = float(value)
    if not np.isfinite(v) or v <= 0:
        raise ValueError(f'{label}: finite positive energy scale required')
    return v


def _complex_json(value: complex) -> dict[str, float]:
    z = complex(value)
    if not np.isfinite([z.real, z.imag]).all():
        raise ValueError('Nonfinite complex quantity')
    return {'real': float(z.real), 'imag': float(z.imag)}


def _matrix_json(matrix: np.ndarray) -> dict[str, Any]:
    m = np.asarray(matrix, dtype=complex)
    if m.shape != (3, 3) or not np.isfinite(m).all():
        raise ValueError('Nonfinite or non-3x3 C5/mass matrix')
    return {
        'real': m.real.tolist(), 'imag': m.imag.tolist(),
        'frobenius_norm': float(np.linalg.norm(m, ord='fro')),
    }


def _threshold_scales(config: Mapping[str, Any], initial_state) -> ThresholdScales:
    params = config.get('thresholds')
    if not isinstance(params, Mapping):
        raise ValueError('thresholds object is required (mu_F, mu_S, mu_low)')
    mu_f = _positive('mu_F', params['mu_F'])
    mu_s = _positive('mu_S', params['mu_S'])
    mu_low = _positive('mu_low', params['mu_low'])
    mu_uv = float(initial_state.mu_gev)
    if not mu_uv > mu_f > mu_s > mu_low:
        raise ValueError('Expected mu_UV > mu_F > mu_S > mu_low')

    # This *grouped* F-first approximation is only valid for a clean hierarchy.
    # Requiring a gap between the boundary and mass eigenvalues makes the
    # omitted multi-threshold refinements explicit.  A mass crossed in running
    # is checked again at the evolved threshold below.
    heavy_masses = np.linalg.svd(initial_state.MF, compute_uv=False)
    scalar_masses = np.sqrt([initial_state.mS1Sq, initial_state.mS2Sq])
    if not (min(heavy_masses) > mu_f > max(scalar_masses)
            and min(scalar_masses) > mu_s):
        raise ValueError(
            'Invalid grouped fermion-first mass hierarchy: require '
            'min(MF) > mu_F > max(MS1,MS2) > min(MS1,MS2) > mu_S; '
            'the Stage 11-15 MF=[1500,2100,2800], MS2=1780 toy is invalid here'
        )
    return ThresholdScales(mu_f, mu_s, mu_low)


def _guard_running_hierarchy(trajectory, scales: ThresholdScales) -> dict:
    heavy = np.asarray(trajectory.fermion_threshold.masses_gev, dtype=float)
    scalars = np.asarray(trajectory.scalar_threshold.masses_gev, dtype=float)
    if heavy.size != 3 or scalars.size != 2 or not np.isfinite(heavy).all() or not np.isfinite(scalars).all():
        raise ValueError('Expected three running fermion masses and two running scalar masses')
    if not (np.min(heavy) > scales.mu_f > np.max(scalars)
            and np.min(scalars) > scales.mu_s):
        raise ValueError(
            'Running masses violate grouped fermion-first ordering; '
            'the model requires individually ordered thresholds or a new benchmark'
        )
    return {'heavy_at_mu_f_gev': sorted(float(x) for x in heavy),
            'scalar_at_mu_s_gev': sorted(float(x) for x in scalars)}


def _validate_payload(payload: Mapping[str, Any], label: str, model: str) -> None:
    if payload.get('status') != 'Success' or not isinstance(payload.get('report_betas'), Mapping):
        raise ValueError(f'{label}: successful RGBeta report_betas is required')
    for name in ALL_COMPLEX:
        if name not in payload['report_betas']:
            raise ValueError(f'{label}: missing B/C(-1) quartic beta {name}')
    meta = payload.get('metadata', {})
    if not isinstance(meta, Mapping):
        raise ValueError(f'{label}: invalid metadata')
    if meta.get('ReportBetaConvention') != '16*pi^2*dX/dln(mu); gauge BetaTerm divided by 2*g':
        raise ValueError(f'{label}: unrecognised beta convention')
    d_f = 1 if model == 'B_m1' else 3
    expected = {'dS1': 2, 'dS2': 2, 'dF': d_f, 'alpha': -1}
    if any(meta.get(key) != value for key, value in expected.items()):
        raise ValueError(f'{label}: representation metadata disagrees with {model}')
    if label == 'EFT1':
        if meta.get('IntegratedField') != 'F' or meta.get('ActiveBSMFields') != ['S1', 'S2']:
            raise ValueError('EFT1 must contain only S1/S2 after F decoupling')
    elif label == 'UV':
        if meta.get('IntegratedField') is not None or meta.get('MajoranaFermion') is not True:
            raise ValueError('UV payload must contain active Majorana F')
    else:
        raise ValueError('Unknown RGBeta stage label')
    # Detailed state-block beta coverage is also checked by production evaluators.


def _transport_summary(trajectory, baseline_c5: np.ndarray | None, low_mnu: np.ndarray | None,
                       model: str, config: Mapping[str, Any], scales: ThresholdScales) -> dict:
    uv = trajectory.uv_initial_state
    at_f = trajectory.uv_threshold_state
    after_f = trajectory.intermediate_initial_state
    at_s = trajectory.eft1_threshold_state
    by_name = {}
    for name in ALL_COMPLEX:
        f_rgb = complex(getattr(at_f, name))
        f_eft = complex(getattr(after_f, name))
        if not np.isclose(f_rgb, f_eft, rtol=1e-10, atol=1e-25):
            raise ValueError(f'{name}: unaccounted discontinuity at F matching')
        by_name[name] = {
            'uv_initial': _complex_json(getattr(uv, name)),
            'uv_at_mu_f': _complex_json(f_rgb),
            'eft1_at_mu_f': _complex_json(f_eft),
            'eft1_at_mu_s': _complex_json(getattr(at_s, name)),
            'uv_generated': _complex_json(f_rgb - complex(getattr(uv, name))),
            'intermediate_running': _complex_json(complex(getattr(at_s, name)) - f_eft),
        }
    masses = _guard_running_hierarchy(trajectory, scales)
    output = {
        'status': 'RenormalisableThresholdPathAndLegacyC5' if baseline_c5 is not None
                  else 'RenormalisableThresholdPathOnly',
        'model': model,
        'scales_gev': {'mu_uv': float(uv.mu_gev), 'mu_F': scales.mu_f,
                       'mu_S': scales.mu_s, 'mu_low': scales.mu_low},
        'running_masses_gev': masses,
        'uv_nfev': int(trajectory.uv.nfev),
        'eft1_nfev': int(trajectory.eft1.nfev),
        'quartics_rgbeta': by_name,
        'c5_active_coupling_mapping': {mat: {'rgbeta': rg, 'matchete_over_rgbeta': -1}
                                       for mat, rg in RGBETA_BY_MAT.items()},
        'source_mapping_same_field_basis_assumed': True,
        'source_mapping_checked_this_run': False,
        'fermion_threshold_renormalisable_boundary': 'continuous tree-level parameters',
        'scalar_threshold_renormalisable_boundary': 'project to SM',
        'baseline_hierarchical_c5_included': baseline_c5 is not None,
        'legacy_baseline_c5_at_mu_s_gev_inverse': _matrix_json(baseline_c5) if baseline_c5 is not None else None,
        'legacy_baseline_neutrino_mass_at_mu_low_ev': _matrix_json(low_mnu * 1e9) if low_mnu is not None else None,
        'extra_c5_status': 'BlockedMissingExtraEFT1WilsonMatchingAndAnomalousDimensions',
        'extra_c5_at_mu_s_gev_inverse': None,
        'extra_neutrino_mass_low_ev': None,
        'complete_two_loop_prediction': False,
        'strict_per_mass_thresholds': False,
        'grouped_threshold_approximation': True,
        'note': (
            'Actual UV and intermediate scalar RGBeta flow with continuous renormalisable F boundary. '
            'Optional BASELINE uses pre-existing validated hierarchical C5 input, NOT a newly matched '
            'extended theory. The extra quartics need consistent LLSS Wilson matching at F, '
            'operator mixing and matching at S. They are NOT inserted into the baseline C5. '
            'Grouped F/S thresholds require a suitable mass hierarchy and are not individual thresholds.'
        ),
    }
    return output


def run_hierarchical_stage16(model: str, config: Mapping[str, Any],
                              uv_payload: Mapping[str, Any], eft_payload: Mapping[str, Any],
                              *, final_weinberg: Path | None = None) -> dict:
    """Evolve both real renormalisable theories and, optionally, existing baseline C5."""
    if model not in ('B_m1', 'C_m1'):
        raise ValueError('Stage16 is restricted to B/C alpha=-1 ordinary models')
    from Numerical.running.T3Trajectory import (
        run_renormalisable_t3_trajectory, continue_with_weinberg_running,
    )
    from Numerical.running.FinalC5TrajectoryAdapter import evaluate_final_c5_on_trajectory

    uv_state = construct_uv_state(model, config)
    scales = _threshold_scales(config, uv_state)
    _validate_payload(uv_payload, 'UV', model)
    _validate_payload(eft_payload, 'EFT1', model)
    path = run_renormalisable_t3_trajectory(
        uv_state, uv_payload, eft_payload,
        mu_fermion_threshold_gev=scales.mu_f,
        mu_scalar_threshold_gev=scales.mu_s,
        uv_save_scales_gev=np.geomspace(float(uv_state.mu_gev), scales.mu_f, 9),
        eft1_save_scales_gev=np.geomspace(scales.mu_f, scales.mu_s, 9),
        rtol=3e-8, atol=3e-12,
    )
    # Guard against crossing an unintended running mass threshold at the
    # endpoints. The UV stage itself need not have a mass-threshold event.
    _guard_running_hierarchy(path, scales)
    c5 = None
    low_mnu = None
    if final_weinberg is not None:
        input_path = Path(final_weinberg)
        if not input_path.is_file():
            raise FileNotFoundError(f'Authoritative final C5 not found: {input_path}')
        c5 = evaluate_final_c5_on_trajectory(input_path, path)
        if not np.all(np.isfinite(c5)) or not np.allclose(c5, c5.T, rtol=1e-10, atol=1e-20):
            raise ValueError('Authoritative baseline C5 is not finite and symmetric')
        final_path = continue_with_weinberg_running(path, c5, scales.mu_low)
        low_mnu = -(246.22 ** 2 / 2) * np.asarray(final_path.weinberg.K, dtype=complex)
    return _transport_summary(path, c5, low_mnu, model, config, scales)
