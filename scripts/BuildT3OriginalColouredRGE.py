"""Full-general shared-A2 results in the actual original coloured T3 RGE dashboard.

Reads the installed Numerical/plotting/templates/T3OptimisedReferenceRGE.html
(the original three-colour UV/EFT1/SMEFT presentation) unchanged; supplies
validated existing UV/EFT1 solver samples and new SMEFT ODE samples.

Optional VISUAL-ONLY EFT1 C5 proxy uses the actual saved model C5(MS) as its
endpoint and a logarithmic interpolation, NOT an EFT1 beta-function solution.
EFT1 leading-log source mixing is available only with an EXTERNALLY CERTIFIED
master-Weyl boundary. The repository's existing seven-term EFT1 source RGE is
then used unchanged. Synthetic/legacy ordered tensors are never plotted as
physical C5. Full LLSS self-running and finite joint-scalar matching remain
outside this diagnostic; the SMEFT threshold stays conditionally one-scale.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
STUDY_ROOT = ROOT / 'output' / 'studies' / 'general_scalar_extended' / 'shared_a2_general16_reports'
SMEFT_FORMAT = 't3-general-shared-a2-sm-weinberg-sampled-v1'
DEFAULT_RUN = 'shared_A2_general16_gamma08_12_v1'
CERTIFIED_ROOT = ROOT / 'output/general_observable/lo_eft1_mixing/certified_boundaries'
EFT1_SOURCE_FORMAT = 't3-eft1-one-loop-c5-source-v1'
EFT1_PHYSICAL_STATUS = 'CertifiedLeadingOneLoopLog_RequiresFiniteJointScalarMatching'
EFT1_VISUAL_FORMAT = 't3-eft1-c5-endpoint-anchored-visual-proxy-v1'
# Source-defined independent gates, reported even before a boundary file exists.
EFT1_CERTIFICATE_GATES = (
    'matchete_gamma_cc_to_weyl_factor_and_sign',
    'canonical_C_ijab_operator_normalization',
    'complete_tree_dim5_operator_inventory',
    'other_dim5_sources_audited_at_one_loop',
    'complete_EFT1_yukawa_and_gauge_dictionary',
    'complete_scalar_quartic_counterterms',
    'all_heavy_scalar_components_and_mixed_masses',
    'weinberg_C5_symmetric_convention',
    'MSbar_subtraction_and_matching_log_partition',
)
# Only the additive files from the two superseded dashboard packages are removed.
# Do not touch source trajectories, fitted inputs, C5 matching, or the original renderer.
OBSOLETE_LOCAL_FILES = (
    'scripts/BuildT3SharedA2ThreeStageRGE.py',
    'scripts/BuildT3SharedA2OriginalRGE.py',
    'tests/python/test_t3_three_stage_rge.py',
    'tests/python/test_t3_original_rge_smeft.py',
    'README_T3_THREE_STAGE_RGE.md',
    'README_T3_ORIGINAL_RGE.md',
    'README_T3_ORIGINAL_COLOURED_RGE.md',
    'README_T3_EFT1_C5_RUNNING.md',
)

SMEFT_KEYS = {
    'c5_fnorm': '||C5||F [GeV^-1]',
    'c5_ee_abs': '|C5_ee| [GeV^-1]',
    'c5_emu_abs': '|C5_eμ| [GeV^-1]',
    'c5_etau_abs': '|C5_eτ| [GeV^-1]',
    'c5_mumu_abs': '|C5_μμ| [GeV^-1]',
    'c5_mutau_abs': '|C5_μτ| [GeV^-1]',
    'c5_tautau_abs': '|C5_ττ| [GeV^-1]',
    'm1_ev': 'm1 [eV]', 'm2_ev': 'm2 [eV]', 'm3_ev': 'm3 [eV]',
    'delta_m21_sq_ev2': 'Δm21² [eV²]',
    'delta_m3l_sq_ev2': 'Δm31² [eV²], NO',
    'delta_m32_sq_ev2': 'Δm32² [eV²], NO',
    'sin2_theta12': 'sin²θ12', 'sin2_theta13': 'sin²θ13',
    'sin2_theta23': 'sin²θ23',
    'gY': 'gY', 'g2': 'g2', 'g3': 'g3',
    'lambdaH': 'λH',
}


# Six independent entries of the complex symmetric SMEFT Weinberg matrix.
_C5_INDEPENDENT = {
    'ee': (0, 0), 'emu': (0, 1), 'etau': (0, 2),
    'mumu': (1, 1), 'mutau': (1, 2), 'tautau': (2, 2),
}
for _name in _C5_INDEPENDENT:
    SMEFT_KEYS['c5_' + _name + '_real'] = 'Re C5_' + _name + ' [GeV^-1]'
    SMEFT_KEYS['c5_' + _name + '_imag'] = 'Im C5_' + _name + ' [GeV^-1]'


def load(path: Path) -> dict:
    data = json.loads(path.read_text(encoding='utf-8-sig'))
    if not isinstance(data, dict):
        raise ValueError('Expected JSON object: ' + str(path))
    return data


def save(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + '\n',
                    encoding='utf-8')


def matrix(raw: Any) -> np.ndarray:
    """Decode compact {real,imag} or saved per-entry complex 3x3 matrices."""
    if isinstance(raw, dict) and 'real' in raw and 'imag' in raw:
        value = np.asarray(raw['real'], dtype=float) + 1j * np.asarray(raw['imag'], dtype=float)
    elif isinstance(raw, list):
        value = np.asarray([[complex(float(x['real']), float(x['imag']))
                             if isinstance(x, dict) else complex(x) for x in row]
                            for row in raw], dtype=complex)
    else:
        raise ValueError('Unexpected matrix representation')
    if value.shape != (3, 3) or not np.isfinite(value).all():
        raise ValueError('Expected finite complex 3x3 matrix')
    return value


def real(raw: Any) -> float:
    if isinstance(raw, dict):
        z = complex(float(raw['real']), float(raw['imag']))
    else:
        z = complex(raw)
    if not math.isfinite(z.real) or not math.isfinite(z.imag) or abs(z.imag) > 1e-11 * max(abs(z.real), 1.):
        raise ValueError('Expected a finite real-valued renormalisable parameter')
    return float(z.real)


def sample_smeft(centre: dict, n: int, *, rtol: float = 1e-8, atol: float = 1e-11) -> dict:
    """Run the already-installed full-flavour one-loop SMEFT ODE at actual t_eval points."""
    from scipy.integrate import solve_ivp
    from Numerical.running.SMWeinbergEvolution import (
        SMWeinbergInitialConditions, _pack, _unpack, _beta,
    )
    from physics.NeutrinoTrajectory import (
        neutrino_mass_from_c5, charged_lepton_mass_basis_matrix,
    )
    from physics.NeutrinoObservables import calculate_neutrino_observables

    if n < 3 or n > 501:
        raise ValueError('--points must be 3..501')
    if centre.get('status') != 'Success':
        raise ValueError('Centre not marked successful')
    detail = centre.get('numerical_details')
    if not isinstance(detail, dict):
        raise ValueError('No saved gamma=1 detailed boundary')
    if detail.get('matching_scope') != 'full-general Matchete quadratic Yukawa channels, one-scale':
        raise ValueError('Only nonmixed general one-scale C5 centres supported')

    eft = detail['eft1_at_scalar_threshold']
    initial = SMWeinbergInitialConditions(
        gY=real(eft['gY']), g2=real(eft['g2']), g3=real(eft['g3']),
        lambdaH=real(eft['lambdaH']),
        ye=matrix(eft['ye']), yu=matrix(eft['yu']), yd=matrix(eft['yd']),
        K=matrix(detail['c5_at_scalar_threshold']),
    ).validated()
    mu_s, mu_low = 1e3, 1e2  # fixed in verified shared-A2 general fitter
    t = np.linspace(math.log(mu_s), math.log(mu_low), n)
    solution = solve_ivp(_beta, (float(t[0]), float(t[-1])), _pack(initial),
                         method='DOP853', t_eval=t, rtol=rtol, atol=atol)
    if not solution.success or solution.y.shape != (76, n):
        raise RuntimeError('SMEFT ODE failed: ' + str(solution.message))

    all_rows = []
    for k in range(n):
        gY, g2, g3, lam, ye, yu, yd, K = _unpack(solution.y[:, k])
        K = 0.5*(K + K.T)
        if not np.isfinite(K).all():
            raise ValueError('Nonfinite SMEFT C5 sample')
        mass = neutrino_mass_from_c5(K, vev_gev=246.22)
        flavour = charged_lepton_mass_basis_matrix(mass, ye)
        o = calculate_neutrino_observables(flavour, ordering='NO')
        if not np.isfinite(o.masses_ev).all() or not math.isfinite(o.takagi_residual):
            raise ValueError('Invalid Takagi output')
        if o.takagi_residual > 1e-7:
            raise ValueError('Takagi residual too large at SMEFT sample')
        pmns = o.pmns_abs
        den = 1. - float(pmns[0, 2]**2)
        if den <= 0.:
            raise ValueError('Invalid PMNS denominator')
        row = {
            'mu_gev': float(math.exp(t[k])),
            'c5_fnorm': float(np.linalg.norm(K)),
            'c5_ee_abs': float(abs(K[0, 0])), 'c5_emu_abs': float(abs(K[0, 1])),
            'c5_etau_abs': float(abs(K[0, 2])),
            'c5_mumu_abs': float(abs(K[1, 1])),
            'c5_mutau_abs': float(abs(K[1, 2])),
            'c5_tautau_abs': float(abs(K[2, 2])),
            'm1_ev': float(o.masses_ev[0]), 'm2_ev': float(o.masses_ev[1]),
            'm3_ev': float(o.masses_ev[2]),
            'delta_m21_sq_ev2': float(o.delta_m21_sq_ev2),
            'delta_m3l_sq_ev2': float(o.delta_m31_sq_ev2),
            'delta_m32_sq_ev2': float(o.delta_m32_sq_ev2),
            'sin2_theta12': float(pmns[0, 1]**2/den),
            'sin2_theta13': float(pmns[0, 2]**2),
            'sin2_theta23': float(pmns[1, 2]**2/den),
            'gY': float(gY), 'g2': float(g2), 'g3': float(g3), 'lambdaH': float(lam),
        }
        for name, (i, j) in _C5_INDEPENDENT.items():
            row['c5_' + name + '_real'] = float(K[i, j].real)
            row['c5_' + name + '_imag'] = float(K[i, j].imag)
        if any(not math.isfinite(value) for value in row.values()):
            raise ValueError('Nonfinite computed SMEFT observable')
        all_rows.append(row)

    saved_low = matrix(detail['c5_at_low_scale'])
    K_out = _unpack(solution.y[:, -1])[-1]
    mismatch = float(np.linalg.norm(K_out - saved_low) / max(np.linalg.norm(saved_low), 1e-99))
    if mismatch > 3e-6:
        raise ValueError(f'SMEFT C5 endpoint fails saved centre parity: {mismatch:.3e}')
    pred = centre.get('prediction') or {}
    for key in ('delta_m21_sq_ev2', 'delta_m3l_sq_ev2', 'sin2_theta12', 'sin2_theta13', 'sin2_theta23'):
        if key in pred and not math.isclose(all_rows[-1][key], float(pred[key]),
                                             rel_tol=3e-5, abs_tol=2e-10):
            raise ValueError(f'SMEFT observable parity mismatch for {key}')

    return {
        'format': SMEFT_FORMAT, 'model': centre['model'], 'status': 'Success',
        'scheme': 'full general one-scale C5 boundary + existing full-flavour SMEFT RGEs',
        'mu_gev': [r.pop('mu_gev') for r in all_rows],
        'series': {key: [r[key] for r in all_rows] for key in SMEFT_KEYS},
        'nfev': int(solution.nfev), 'relative_endpoint_c5_residual': mismatch,
        'normal_ordering': True, 'vev_gev_fixed': 246.22,
        'notes': [
            'All plotted samples are numerical ODE samples, not endpoint interpolation.',
            'mν(μ) = -v² C5(μ)/2 with v held fixed, as in the existing fitted pipeline.',
            'Saved C5 boundary originates from conditional full-general one-scale Matchete matching.',
            'This stage does not certify missing physical EFT1 Wilson mixing or sequential matching.',
            'The Higgs quadratic mass-squared mH² is not part of the current SMWeinbergEvolution state.',
        ]
    }




def visual_eft1_c5_proxy(centre: dict, mu: list[float], smeft: dict) -> dict:
    """Illustrative log-interpolation, NOT EFT1 physical/Wilson RGE evolution.

    K_proxy(mu) = K_saved(MS) * ln(MF/mu)/ln(MF/MS), MF=1e5, MS=1e3 GeV.
    Anchors the REAL saved one-scale C5(MS) at the lower end and chooses
    zero at MF solely to make a visible numerical display. The zero is NOT a
    UV/EFT1 matching result. No synthetic/external C5 values are imported.
    """
    if centre.get('status') != 'Success' or not isinstance(centre.get('numerical_details'), dict):
        raise ValueError('Cannot construct display proxy without successful saved model centre')
    detail = centre['numerical_details']
    if detail.get('matching_scope') != 'full-general Matchete quadratic Yukawa channels, one-scale':
        raise ValueError('Display proxy requires the nonmixed full-general matched C5 boundary')
    k_ms = matrix(detail['c5_at_scalar_threshold'])
    asym = np.linalg.norm(k_ms-k_ms.T) / max(np.linalg.norm(k_ms),1e-90)
    if asym > 1e-7:
        raise ValueError('Saved threshold C5 is not flavour symmetric')
    k_ms = 0.5 * (k_ms + k_ms.T)
    x = np.asarray(mu, dtype=float)
    if x.ndim != 1 or x.size < 3 or not np.isfinite(x).all() or np.any(x <= 0):
        raise ValueError('Invalid EFT1 scale samples')
    if not np.isclose(x[0], 1e5, rtol=1e-9) or not np.isclose(x[-1], 1e3, rtol=1e-9):
        raise ValueError('Visual EFT1 proxy requires the saved 1e5 -> 1e3 GeV interval')
    if np.any(x[:-1] <= x[1:]):
        raise ValueError('Visual EFT1 scales must descend')
    w = np.log(1e5/x) / math.log(1e5/1e3)
    w[0], w[-1] = 0.0, 1.0
    if np.any(w < -1e-8) or np.any(w > 1+1e-8):
        raise ValueError('EFT1 proxy interpolation escaped [0,1]')
    mats = w[:,None,None] * k_ms[None,:,:]
    series = {'c5_fnorm': np.linalg.norm(mats, axis=(1,2)).astype(float).tolist()}
    for name, (i,j) in _C5_INDEPENDENT.items():
        vals=mats[:,i,j]
        series['c5_'+name+'_abs'] = np.abs(vals).astype(float).tolist()
        series['c5_'+name+'_real'] = vals.real.astype(float).tolist()
        series['c5_'+name+'_imag'] = vals.imag.astype(float).tolist()
    # Every proxy channel must join the ACTUAL SMEFT solver boundary; don't
    # silently draw a discontinuity or invent a threshold matching factor.
    for key in SMEFT_C5:
        green = smeft.get('series', {}).get(key)
        if not isinstance(green, list) or not green:
            raise ValueError('Missing SMEFT C5 sample for proxy continuity: '+key)
        if not math.isclose(series[key][-1], float(green[0]), rel_tol=5e-6, abs_tol=1e-25):
            raise ValueError('Proxy/Smeft C5 threshold mismatch: '+key)
    return {
        'format': EFT1_VISUAL_FORMAT,
        'model': centre['model'],
        'status': 'IllustrativeInterpolation_NotEFT1RGESolution',
        'mu_gev': x.tolist(), 'series': series,
        'c5_muS_source': 'centre.json numerical_details.c5_at_scalar_threshold',
        'formula': 'C5_proxy(mu)=C5_matched(muS)*ln(muF/mu)/ln(muF/muS)',
        'muF_gev': 1e5, 'muS_gev': 1e3,
        'chosen_proxy_C5_muF_gev_inv': 'zero (display convention, NOT matched Wilson boundary)',
        'physical_eft1_c5': False,
        'eft1_beta_function_solved': False,
        'matched_to_smeft_at_muS': True,
        'do_not_interpret_as_neutrino_mass_prediction': True,
    }


# Original template source (not the later shared-A2 custom-canvas renderer).
TEMPLATE = ROOT / 'Numerical' / 'plotting' / 'templates' / 'T3OptimisedReferenceRGE.html'
UV_EFT_FORMAT = 't3-shared-a2-general-uv-eft1-rge-trajectories-v1'
OUTPUT_FORMAT = 't3-original-optimised-reference-coloured-rge-general-v1'

# These quantities are present in the saved *renormalisable* UV/EFT1 state.
# Do not create nonexistent Wilson coefficient entries in EFT1.
RUNNING_KEYS = {
    'gY': ('Hypercharge coupling gY', ''),
    'g2': ('Weak coupling g2', ''),
    'g3': ('Strong coupling g3', ''),
    'lambdaH': ('Higgs quartic λH', ''),
    'lambdaT3_abs': ('|λT3|', ''),
    'lambdaT3_real': ('Re λT3', ''),
    'mS1_sq_gev2': ('mS1²', 'GeV²'),
    'mS2_sq_gev2': ('mS2²', 'GeV²'),
    'mS1_gev': ('√mS1²', 'GeV'),
    'mS2_gev': ('√mS2²', 'GeV'),
    'muS1_real_gev2': ('Re μS1²', 'GeV²'),
    'muS1_imag_gev2': ('Im μS1²', 'GeV²'),
    'muS2_real_gev2': ('Re μS2²', 'GeV²'),
    'muS2_imag_gev2': ('Im μS2²', 'GeV²'),
    'y1_norm': ('||Y1||F', ''),
    'y2_norm': ('||Y2||F', ''),
    'y1Conj_norm': ('||Y1conj||F', ''),
    'y2Conj_norm': ('||Y2conj||F', ''),
}

SMEFT_C5 = [key for key in SMEFT_KEYS if key.startswith('c5_')]
SMEFT_NEUTRINO = [
    'm1_ev', 'm2_ev', 'm3_ev',
    'delta_m21_sq_ev2', 'delta_m3l_sq_ev2', 'delta_m32_sq_ev2',
    'sin2_theta12', 'sin2_theta13', 'sin2_theta23',
]
SMEFT_ONLY = set(SMEFT_C5 + SMEFT_NEUTRINO)
# These are *only* the leading-log, LLSS -> LLHH source contributions.
# In particular they are not the complete C5(mu) or a neutrino mass prediction.
EFT1_SOURCE_KEYS = {
    'eft1_delta_c5_fnorm': ('EFT1 ||Delta C5 (one-loop log)||F', 'GeV⁻¹'),
}
for _name in _C5_INDEPENDENT:
    EFT1_SOURCE_KEYS['eft1_delta_c5_' + _name + '_abs'] = (
        'EFT1 |Delta C5_' + _name + ' (one-loop log)|', 'GeV⁻¹'
    )


def _certified_c5_matrix(raw: Any, name: str) -> np.ndarray:
    value = matrix(raw)
    if not np.allclose(value, value.T, rtol=1e-10, atol=1e-25):
        raise ValueError(name + ': coefficient is not symmetric in flavour')
    return 0.5 * (value + value.T)


def certified_eft1_source(model: str, boundary_path: Path, n_points: int) -> tuple[dict, dict]:
    """Run the already-implemented seven-term EFT1 mixing using a physical input.

    IMPORTANT: all nine external certifications are required by the installed
    RGE solver; this adapter neither generates nor asserts any certificates.
    A leading log is plotted alone, NEVER added to one-scale C5 at mu_S.
    """
    from RGE.running.intermediate.GeneralEFT1C5Mixing import (
        solve_boundary, CERTIFICATES,
    )
    raw = load(boundary_path)
    if raw.get('model') != model:
        raise ValueError('EFT1 source model does not match the selected benchmark')
    if raw.get('authority') != 'certified':
        raise ValueError('Synthetic or non-certified EFT1 inputs are forbidden in the thesis dashboard')
    cert = raw.get('certificates') or {}
    missing = [name for name in CERTIFICATES if cert.get(name) is not True]
    if missing:
        raise ValueError('Missing independent EFT1 master-Weyl certifications: ' + ', '.join(missing))
    result = solve_boundary(raw, allow_synthetic=False, trajectory_points=n_points)
    scope = result.get('physics_completeness') or {}
    if (result.get('format') != EFT1_SOURCE_FORMAT or result.get('status') != EFT1_PHYSICAL_STATUS
            or result.get('authority') != 'certified' or
            scope.get('physical_eft1_source_verified') is not True or
            scope.get('joint_two_scalar_finite_matching_included') is not False):
        raise ValueError('EFT1 master-RGE result failed its physical-source certification gate')
    scales = result.get('scales_gev') or {}
    if not math.isclose(float(scales.get('mu_F', float('nan'))), 1e5, rel_tol=1e-9) or not math.isclose(
            float(scales.get('mu_S', float('nan'))), 1e3, rel_tol=1e-9):
        raise ValueError('EFT1 source does not use this campaign fermion/scalar matching scales')
    samples = result.get('trajectory') or []
    if len(samples) != n_points:
        raise ValueError('Wrong number of leading-log samples')
    mu, entries = [], []
    for sample in samples:
        scale = float(sample['mu_gev'])
        if not math.isfinite(scale) or scale <= 0:
            raise ValueError('Invalid EFT1 source running scale')
        mu.append(scale)
        entries.append(_certified_c5_matrix(sample['delta_C5_log_gev_inverse'], 'EFT1 Delta C5'))
    if any(mu[i] <= mu[i+1] for i in range(len(mu)-1)) or not math.isclose(
            mu[0], 1e5, rel_tol=1e-9) or not math.isclose(mu[-1], 1e3, rel_tol=1e-9):
        raise ValueError('EFT1 leading-log scales must descend from 1e5 to 1e3 GeV')
    if np.linalg.norm(entries[0]) > 1e-13 * max(np.linalg.norm(entries[-1]), 1e-40):
        raise ValueError('EFT1 log source is not zero at muF')
    final = _certified_c5_matrix(result['total_delta_C5_log_gev_inverse'], 'EFT1 total Delta C5')
    rel = float(np.linalg.norm(entries[-1]-final) / max(np.linalg.norm(final), 1e-40))
    if rel > 1e-8:
        raise ValueError('EFT1 log endpoint does not match its integrated source')
    series = {'eft1_delta_c5_fnorm': [float(np.linalg.norm(z)) for z in entries]}
    for name,(i,j) in _C5_INDEPENDENT.items():
        series['eft1_delta_c5_'+name+'_abs'] = [float(abs(z[i,j])) for z in entries]
    return result, {'mu_gev':mu, 'series':series}


def remove_obsolete_dashboard_files(root: Path, output_dir: Path, *, clean_scripts: bool) -> list[str]:
    """Conservative, explicit cleanup after successful replacement is on disk."""
    removed = []
    for name in ('interactive_RGE_with_SMEFT.html', 'interactive_RGE_previous.html'):
        path = output_dir / name
        if path.is_file():
            path.unlink()
            removed.append(str(path))
    if clean_scripts:
        for name in OBSOLETE_LOCAL_FILES:
            path = root / name
            if path.is_file() and path.is_relative_to(root):
                path.unlink()
                removed.append(str(path))
    return removed


def display_name(model: str) -> str:
    import re
    m = re.fullmatch(r'([A-E])_([+-])(\d+)', model)
    if m is None:
        raise ValueError('Unrecognized numerical model key: ' + model)
    return m.group(1) + '(' + m.group(2) + m.group(3) + ')'


def _validated_stage(stage: dict, label: str) -> tuple[list[float], dict[str, list[float]]]:
    if not isinstance(stage, dict):
        raise ValueError(label + ': no stage payload')
    mu = stage.get('mu_gev', [])
    series = stage.get('series', {})
    if not isinstance(mu, list) or len(mu) < 3 or not isinstance(series, dict):
        raise ValueError(label + ': no dense ODE trajectory')
    scale = [float(x) for x in mu]
    if any(not math.isfinite(x) or x <= 0 for x in scale) or any(
        scale[i] <= scale[i+1] for i in range(len(scale)-1)
    ):
        raise ValueError(label + ': scales must be finite positive and descending')
    parsed = {}
    for k, seq in series.items():
        if not isinstance(seq, list) or len(seq) != len(scale):
            raise ValueError(label + ': series length mismatch ' + k)
        vals = [float(v) for v in seq]
        if any(not math.isfinite(v) for v in vals):
            raise ValueError(label + ': non-finite ' + k)
        parsed[k] = vals
    return scale, parsed


def _joint_trace(stages: dict, quantity: str) -> dict | None:
    """Only join genuine adjacent ODE samples with a compatible threshold boundary."""
    mu_all, y_all = [], []
    for stage_name in ('uv', 'eft1', 'eft1_c5', 'smeft'):
        stage = stages.get(stage_name)
        if not stage or quantity not in stage['series']:
            continue
        mu, vals = stage['mu_gev'], stage['series'][quantity]
        if mu_all:
            # If the parameter was absent in an intermediate stage, do NOT draw
            # a false connecting line across a missing physical RGE interval.
            if not math.isclose(mu_all[-1], mu[0], rel_tol=1e-9, abs_tol=0):
                raise ValueError(quantity + ': cannot join stages across an uncomputed interval')
            if not math.isclose(y_all[-1], vals[0], rel_tol=2e-5,
                                abs_tol=2e-9 * max(1., abs(y_all[-1]), abs(vals[0]))):
                raise ValueError(quantity + ': threshold endpoint mismatch; cannot join')
            mu, vals = mu[1:], vals[1:]
        mu_all.extend(mu)
        y_all.extend(vals)
    if not mu_all:
        return None
    return {'x': mu_all, 'y': y_all}


def _assemble_curves(trajectories: dict[str, dict]) -> tuple[dict, dict]:
    """Convert saved physical trajectories to exact original-template data shape."""
    if not trajectories:
        raise ValueError('No original UV/EFT1 results available')
    quantities = {}
    labels = {}
    sorted_models = sorted(trajectories)
    for raw_name in sorted_models:
        stages = trajectories[raw_name]
        model = display_name(raw_name)
        keys = set().union(*(set(v.get('series', {})) for v in stages.values()))
        for q in sorted(keys):
            if q not in RUNNING_KEYS and q not in SMEFT_KEYS and q not in EFT1_SOURCE_KEYS:
                continue
            try:
                curve = _joint_trace(stages, q)
            except ValueError as exc:
                # A matching discontinuity must not be drawn as a smooth curve.
                # Preserve the stage-resolved result, with an explicit label.
                for stage_name, stage in stages.items():
                    if q not in stage['series']:
                        continue
                    name = q + '__' + stage_name
                    quantities.setdefault(name, {})[model] = {
                        'x': stage['mu_gev'], 'y': stage['series'][q]}
                    title, unit = RUNNING_KEYS.get(q, EFT1_SOURCE_KEYS.get(q, (SMEFT_KEYS.get(q, q), '')))
                    labels[name] = {'label': title + ' (' + stage_name + '; unjoined threshold)',
                                    'unit': unit, 'block': 'saved_ode',
                                    'weinberg': q.startswith('c5_'),
                                    'smeft_only': stage_name == 'smeft',
                                    'eft1_only': stage_name == 'eft1_c5',
                                    'note': str(exc)}
                continue
            if curve:
                quantities.setdefault(q, {})[model] = curve
                title, unit = RUNNING_KEYS.get(q, EFT1_SOURCE_KEYS.get(q, (SMEFT_KEYS.get(q, q), '')))
                proxy = (q in SMEFT_C5 and q in stages.get('eft1', {}).get('series', {}))
                if proxy:
                    title += ' — yellow: DISPLAY PROXY; green: SMEFT ODE'
                labels[q] = {'label': title, 'unit': unit,
                             'block': 'eft1_visual_proxy_plus_smeft_ode' if proxy else 'saved_ode',
                             'weinberg': q.startswith('c5_'),
                             'visual_proxy': proxy,
                             'smeft_only': q in SMEFT_ONLY and not proxy,
                             'eft1_only': q in EFT1_SOURCE_KEYS}
    keys = set(quantities)
    if not keys.intersection(SMEFT_C5):
        raise RuntimeError('SMEFT C5 has no valid numerical samples; refusing to claim success')
    payload = {
        'format': OUTPUT_FORMAT,
        'models': [display_name(m) for m in sorted_models],
        'scales': {'mu_uv': 1e7, 'mu_fermion_threshold': 1e5,
                   'mu_scalar_threshold': 1e3, 'mu_low': 1e2},
        'quantities': quantities,
        'matched': {}, 'components': {},
        'continuation_methods': {}, 'continuation_issues': {},
        'eft1_wilson_c5_running_certified': False,
        'eft1_leading_log_source_only': bool(set(quantities) & set(EFT1_SOURCE_KEYS)),
        'eft1_c5_visual_proxy': any(v.get('visual_proxy') for v in labels.values()),
        'smeft_boundary_one_scale_conditional': True,
    }
    labels['llss_kernel'] = {
        'label': 'Universal LLSS log kernel η (NOT physical C5 running)',
        'unit': '', 'block': 'analytic_kernel', 'weinberg': False,
        'smeft_only': False,
    }
    return payload, labels


def _options(meta: dict) -> str:
    import html
    def group(name: str, keys: list[str]) -> str:
        options = ''.join(
            '<option value="' + html.escape(k, quote=True) + '">' +
            html.escape(meta[k]['label']) + '</option>'
            for k in keys if k in meta
        )
        return '<optgroup label="' + html.escape(name, quote=True) + '">' + options + '</optgroup>' if options else ''
    return (
        group('EFT1 C5 — certified one-loop source ONLY', list(EFT1_SOURCE_KEYS))
        + group('C5 — yellow DISPLAY PROXY / green SMEFT solver',
                [k for k in SMEFT_C5 if meta.get(k,{}).get('visual_proxy')])
        + group('SMEFT C5 — actual ODE samples',
                [k for k in SMEFT_C5 if not meta.get(k,{}).get('visual_proxy')])
        + group('SMEFT neutrino masses and mixing', SMEFT_NEUTRINO)
        + group('Renormalisable couplings — UV / EFT1 / SMEFT', list(RUNNING_KEYS))
        + group('Separate unjoined threshold stages', [k for k in meta if '__' in k])
        + group('EFT1 diagnostic only (NOT Wilson evolution)', ['llss_kernel'])
    )


def _replace_one(template: str, src: str, dest: str, label: str) -> str:
    count = template.count(src)
    if count != 1:
        raise ValueError('Installed original RGBeta dashboard template changed: ' + label +
                         '; count=' + str(count))
    return template.replace(src, dest, 1)


def actual_original_template_page(template: str, data: dict, meta: dict,
                                  title: str = 'Shared A(-2) — general T3 RGE running') -> str:
    """Use the original three-colour dashboard/plot JS, not the custom A2 canvas."""
    import html
    expected = [
        "const COLORS=['#2466ae','#c46c25'",
        "rectRegion(s.mu_uv,s.mu_fermion_threshold,regionColours.uv,'Full T3');",
        "rectRegion(s.mu_fermion_threshold,s.mu_scalar_threshold,regionColours.int,'Intermediate scalar EFT');",
        "rectRegion(s.mu_scalar_threshold,s.mu_low,regionColours.sm,'SM + Weinberg EFT');",
        'id="zoomUV"', 'id="zoomInt"', 'id="zoomFinal"',
        'id="yScaleMode"', 'const DATA=__DATA__;', 'const META=__META__;',
        '<select id="quantity">__OPTIONS__</select>',
    ]
    for token in expected:
        if token not in template:
            raise ValueError('This is not the original coloured T3 RGE template: ' + token)
    # Only correct labels / conditional matching remarks; never replace its CSS,
    # model palette, shaded EFT regions, zoom logic, or canvas algorithm.
    template = _replace_one(
        template,
        'Numerical integration direction: UV → IR. Saved benchmark trajectories; no matching or RGE recalculation.',
        'UV/EFT1 from saved general-RGBeta trajectories; SMEFT integrated using the existing SM + Weinberg ODE. '
        'Conditional one-scale C5 matching. Yellow C5 (when present) is an ILLUSTRATIVE endpoint-anchored log INTERPOLATION, NOT an EFT1 RGE. Physical EFT1 LLSS Wilson running is NOT available.',
        'header provenance')
    template = _replace_one(
        template,
        ' if(isLLSS(q)&&s.mu_fermion_threshold>s.mu_scalar_threshold)return [s.mu_fermion_threshold,s.mu_scalar_threshold];',
        ' if(META[q]&&META[q].eft1_only)return [s.mu_fermion_threshold,s.mu_scalar_threshold];\n'
        ' if(META[q]&&META[q].smeft_only)return [s.mu_scalar_threshold,s.mu_low];\n'
        ' if(isLLSS(q)&&s.mu_fermion_threshold>s.mu_scalar_threshold)return [s.mu_fermion_threshold,s.mu_scalar_threshold];',
        'default SMEFT zoom')
    template = _replace_one(
        template,
        'const showParts=weinberg&&showComponents.checked;',
        'const showParts=weinberg&&!(META[q]&&(META[q].smeft_only||META[q].visual_proxy))&&showComponents.checked;',
        'unavailable hard/direct components')
    template = _replace_one(
        template,
        "document.getElementById('weinbergLegend').hidden=!weinberg;",
        "document.getElementById('weinbergLegend').hidden=!weinberg||(META[q]&&(META[q].smeft_only||META[q].visual_proxy));",
        'legend for SMEFT-only')
    template = _replace_one(
        template,
        'componentControl.hidden=!weinberg;',
        'componentControl.hidden=!weinberg||(META[q]&&(META[q].smeft_only||META[q].visual_proxy));',
        'component selector for SMEFT-only')
    template = _replace_one(
        template,
        '}else if(weinberg){',
        "}else if(weinberg&&META[q]&&META[q].visual_proxy){\n"
        "  caption.textContent='ILLUSTRATIVE ONLY: yellow EFT1 C5 is logarithmic interpolation from a chosen zero at MF=10^5 GeV to the SAVED one-scale C5(MS) at MS=10^3 GeV. This is NOT an EFT1 Wilson RGE solution or physical matched C5 above MS. Green: numerically integrated SMEFT C5. Do NOT use the yellow part as a physics prediction.';\n"
        '}else if(weinberg&&META[q]&&META[q].smeft_only){\n'
        "  caption.textContent='Green region: genuinely integrated SMEFT C5(μ), starting from the saved conditional one-scale scalar-threshold boundary. Yellow region: physical EFT1 C5 / LLSS Wilson evolution has NOT been computed; no continuation is shown.';\n"
        '}else if(weinberg){',
        'honest C5 caption')
    template = _replace_one(
        template,
        "  caption.textContent=DATA.models.filter(m=>enabled[m]&&traces[m]).length+' model(s) shown. Curves are taken from saved diagnostics.';",
        "  if(META[q]&&META[q].eft1_only){caption.textContent='Yellow EFT1: externally certified one-loop LLSS-to-LLHH source ΔC5(log μ), NOT complete C5(μ). No LLSS self-running or finite joint-scalar matching. Do not join or add this curve to the conditional SMEFT C5 boundary.';}else{caption.textContent=DATA.models.filter(m=>enabled[m]&&traces[m]).length+' model(s) shown. Curves are taken from saved diagnostics.';}",
        'certified EFT1 source-only caption')
    text_blob = json.dumps(data, ensure_ascii=False, separators=(',', ':'), allow_nan=False).replace('</', '<\\/')
    meta_blob = json.dumps(meta, ensure_ascii=False, separators=(',', ':'), allow_nan=False).replace('</', '<\\/')
    for token, content in (
        ('__TITLE__', html.escape(title)),
        ('__OPTIONS__', _options(meta)),
        ('__DATA__', text_blob),
        ('__META__', meta_blob),
    ):
        if token not in template:
            raise ValueError('Original template missing: ' + token)
        template = template.replace(token, content)
    if any(k in template for k in ('__TITLE__','__OPTIONS__','__DATA__','__META__')):
        raise RuntimeError('Original template still has unresolved markers')
    return template


def _save_smeft_csv(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='', encoding='utf-8') as out:
        writer = csv.writer(out)
        writer.writerow(['mu_gev', *SMEFT_KEYS])
        for i, mu in enumerate(payload['mu_gev']):
            writer.writerow([mu] + [payload['series'][k][i] for k in SMEFT_KEYS])


def run(run_id: str = DEFAULT_RUN, points: int = 41, *, root: Path = STUDY_ROOT,
        template_path: Path = TEMPLATE, certified_root: Path = CERTIFIED_ROOT,
        clean_obsolete: bool = False, display_proxy: bool = False) -> dict:
    """Regenerate ORIGINAL styled interactive_RGE.html, discarding wrong separate view."""
    import hashlib
    import re
    import os
    if not re.fullmatch('[A-Za-z0-9_-]+', run_id):
        raise ValueError('Invalid --run-id')
    if not 3 <= points <= 501:
        raise ValueError('--points must lie in [3,501]')
    folder = root / run_id
    if not (folder / 'manifest.json').is_file():
        raise FileNotFoundError('Missing shared-A2 general run: ' + str(folder))
    manifest = load(folder / 'manifest.json')
    if manifest.get('interactive_rge_generated') is not True:
        raise ValueError('Saved general UV/EFT1 original running has not been generated')
    if not template_path.is_file():
        raise FileNotFoundError('Missing actual original coloured dashboard template: ' + str(template_path))
    template = template_path.read_text(encoding='utf-8')
    # Validate original style before doing any physics work or touching the report.
    actual_original_template_page(template, {'models':[]}, {'llss_kernel': {'label':'', 'unit':'', 'weinberg':False}})
    models = {}
    statuses = []
    eft1_readiness = []
    for model_dir in sorted(folder.iterdir()):
        if not model_dir.is_dir() or not re.fullmatch('[A-E]_[mp][0-9]+', model_dir.name):
            continue
        tag = model_dir.name
        source_path = model_dir / 'rge_trajectories.json'
        centre_path = model_dir / 'centre.json'
        certified_input = certified_root / tag / 'boundary.json'
        eft1_stage = None
        if not certified_input.is_file():
            # Never retain a previous derived physical source when the required
            # external calibration input has disappeared.
            (model_dir/'eft1_c5_source.json').unlink(missing_ok=True)
            eft1_readiness.append({'model':tag, 'status':'MissingCertifiedBoundary',
                                   'expected_path':str(certified_input),
                                   'note':'Synthetic preflight output is NOT a physical C5 boundary.'})
        elif not source_path.is_file() or not centre_path.is_file():
            eft1_readiness.append({'model':tag, 'status':'MissingCompatibleGeneralBenchmark',
                                   'boundary_path':str(certified_input)})
        else:
            try:
                raw_cert, eft1_stage = certified_eft1_source(
                    tag[0]+'_'+('-' if tag[2]=='m' else '+')+tag[3:], certified_input, points)
                save(model_dir/'eft1_c5_source.json', raw_cert)
                (model_dir/'eft1_c5_error.json').unlink(missing_ok=True)
                eft1_readiness.append({'model':tag, 'status':'CertifiedSourceComputed',
                                       'boundary_path':str(certified_input),
                                       'full_LLSS_self_running':False,
                                       'finite_scalar_matching':False})
            except Exception as exc:
                (model_dir/'eft1_c5_source.json').unlink(missing_ok=True)
                message=type(exc).__name__+': '+str(exc)
                eft1_readiness.append({'model':tag,'status':'RejectedCertifiedBoundary',
                                       'boundary_path':str(certified_input), 'error':message})
                save(model_dir/'eft1_c5_error.json', {'model':tag,'error':message})
        if not source_path.is_file():
            statuses.append({'model':tag,'status':'Unavailable',
                             'reason':'Original dense UV/EFT1 ODE trajectory not saved (B/C conditional)'})
            continue
        raw = load(source_path)
        if raw.get('format') != UV_EFT_FORMAT or raw.get('model') != (
                tag[0] + '_' + ('-' if tag[2]=='m' else '+') + tag[3:]) or raw.get('status')!='Success':
            statuses.append({'model':tag,'status':'InvalidUV_EFT1', 'reason':'Source format/model mismatch'})
            continue
        try:
            stages = {}
            for label, upper, lower in (
                ('uv',1e7,1e5),('eft1',1e5,1e3),
            ):
                mu,series = _validated_stage(raw[label], tag+'/'+label)
                if not math.isclose(mu[0],upper,rel_tol=1e-9) or not math.isclose(mu[-1],lower,rel_tol=1e-9):
                    raise ValueError(label + ': saved scale endpoints differ from original study')
                stages[label] = {'mu_gev':mu,'series':series}
            if not centre_path.is_file():
                raise FileNotFoundError('No gamma=1 centre to construct SMEFT C5 boundary')
            centre = load(centre_path)
            smeft = sample_smeft(centre, points)
            if smeft['model'] != raw['model']:
                raise ValueError('SMEFT and UV/EFT1 model tags disagree')
            mu, series = _validated_stage(smeft, tag+'/smeft')
            if not math.isclose(mu[0],1e3,rel_tol=1e-9) or not math.isclose(mu[-1],1e2,rel_tol=1e-9):
                raise ValueError('SMEFT sample endpoints not at 1e3 and 1e2 GeV')
            stages['smeft'] = {'mu_gev':mu,'series':series}
            if display_proxy:
                proxy = visual_eft1_c5_proxy(centre, stages['eft1']['mu_gev'], smeft)
                stages['eft1']['series'].update(proxy['series'])
                save(model_dir / 'eft1_c5_visual_proxy.json', proxy)
                for row in eft1_readiness:
                    if row['model']==tag:
                        row['visual_proxy_plotted'] = True
                        if row['status']=='MissingCertifiedBoundary':
                            row['status']='VisualProxyPlotted_NoCertifiedBoundary'
                        break
            if eft1_stage is not None:
                stages['eft1_c5'] = eft1_stage
            models[raw['model']] = stages
            save(model_dir / 'smeft_trajectories.json', {
                **smeft,
                'centre_sha256': hashlib.sha256(centre_path.read_bytes()).hexdigest(),
                'uv_eft1_sha256': hashlib.sha256(source_path.read_bytes()).hexdigest(),
            })
            _save_smeft_csv(model_dir / 'smeft_observables.csv', smeft)
            (model_dir / 'smeft_error.json').unlink(missing_ok=True)
            statuses.append({'model':tag,'status':'Success','smeft_points':points,
                             'c5_endpoint_relative_error':smeft['relative_endpoint_c5_residual']})
        except Exception as exc:
            # Retain actual renormalisable UV/EFT1 curves even when SMEFT fails.
            try:
                stages = {label: dict(zip(('mu_gev','series'),
                           _validated_stage(raw[label], tag+'/'+label)))
                          for label in ('uv','eft1')}
                models[raw['model']] = stages
            except Exception:
                pass
            message = type(exc).__name__+': '+str(exc)
            save(model_dir / 'smeft_error.json', {'model':tag,'status':'Failed','reason':message})
            statuses.append({'model':tag,'status':'FailedSMEFT','reason':message})
    payload, meta = _assemble_curves(models)
    html_body = actual_original_template_page(template, payload, meta,
                    'Shared A(-2) general T3 — original coloured RGE dashboard')
    destination = folder / 'interactive_RGE.html'
    tmp = destination.with_suffix('.html.tmp')
    tmp.write_text(html_body, encoding='utf-8')
    os.replace(tmp, destination)
    removed = remove_obsolete_dashboard_files(ROOT, folder, clean_scripts=clean_obsolete)
    save(folder/'eft1_c5_readiness.json', {
        'status':'CertifiedSourceFound' if any(r['status']=='CertifiedSourceComputed' for r in eft1_readiness)
                 else 'NoCompatibleCertifiedC5Source',
        'models':eft1_readiness,
        'source':'RGE/running/intermediate/GeneralEFT1C5Mixing.py',
        'physical_source_requires_nine_external_certificates':True,
        'required_external_certifications':list(EFT1_CERTIFICATE_GATES),
        'synthetic_or_legacy_source_accepted':False,
        'full_LLSS_wilson_running_computed':False,
        'finite_joint_scalar_matching_computed':False,
        'notes':('Visual-proxy mode draws an explicitly unphysical, endpoint-anchored C5 interpolation in yellow. '
                 'Only certified fixed-order Delta C5 is eligible for the physical-source plot.')})
    result = {
        'format': OUTPUT_FORMAT, 'run_id': run_id,
        'original_renderer_source': 'Numerical/plotting/templates/T3OptimisedReferenceRGE.html',
        'dashboard': str(destination),
        'real_ode_model_count': len(models),
        'smeft_model_count': sum(row['status']=='Success' for row in statuses),
        'model_status': statuses,
        'physical_eft1_wilson_running': False,
        'scalar_threshold_one_scale_matching_conditional': True,
        'wrong_parallel_dashboard_removed': not (folder / 'interactive_RGE_with_SMEFT.html').exists(),
        'eft1_certified_source_count': sum(r['status']=='CertifiedSourceComputed' for r in eft1_readiness),
        'eft1_c5_visual_proxy_model_count': sum(bool(r.get('visual_proxy_plotted')) for r in eft1_readiness),
        'eft1_c5_visual_proxy_not_physical': bool(display_proxy),
        'eft1_readiness_report': str(folder/'eft1_c5_readiness.json'),
        'removed_obsolete_files': removed,
    }
    save(folder / 'original_coloured_rge_status.json', result)
    return result


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run-id', default=DEFAULT_RUN)
    p.add_argument('--points', type=int, default=41)
    p.add_argument('--certified-boundary-root', type=Path, default=CERTIFIED_ROOT,
                   help='Only externally certified master-Weyl boundary.json files; synthetic rejected')
    p.add_argument('--eft1-display-proxy', action='store_true',
                   help='SHOW numerical but UNPHYSICAL endpoint-anchored C5 curves in original yellow EFT1 region')
    p.add_argument('--clean-obsolete', action='store_true',
                   help='Delete only the explicitly superseded additive dashboard scripts/tests/readmes')
    args = p.parse_args(argv)
    status = run(args.run_id,args.points,certified_root=args.certified_boundary_root,
                 clean_obsolete=args.clean_obsolete,display_proxy=args.eft1_display_proxy)
    print(json.dumps(status,ensure_ascii=False,indent=2))
    return 0 if status['smeft_model_count'] else 2


if __name__=='__main__':
    raise SystemExit(main())
