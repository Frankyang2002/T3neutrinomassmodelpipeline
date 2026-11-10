"""Strict source-gated completion path for A(-4)/D(-2) full-general T3.

This module never interprets Laurent-regular Matchete operators as renormalized
Wilson coefficients. A *numerical* complete scalar-threshold Weinberg boundary
must originate in a separate full-general calculation and carry matching-source
hashes. Once available, this module performs actual SMEFT RG, Takagi masses,
PMNS angles, and five NuFIT residuals with the repository's existing solvers.

Stage status is evidence-based. This is an execution interface, NOT a claim
that the presently committed matching artifacts meet its boundary contract.
"""
from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from Numerical.parallel_t3.background.general_batch1_api import (
    load_batch1, decode_native, encode_native,
)
from Numerical.parallel_t3.wilson.general_batch2_laurent import (
    check_saved_raw, check_audit_report, read_json,
)

MODELS = ('A_-4', 'D_-2')
C5_FORMAT = 't3-thesis-full-general-renormalized-c5-at-ms-v1'
REPORT_FORMAT = 't3-thesis-full-general-consolidated-candidate-v1'
STAGES = ('UV_general_RG', 'F_all_fermion_raw_one_loop',
          'F_Laurent', 'F_finite_one_loop_Wilson',
          'EFT1_full_Wilson_RGE', 'S_joint_full_one_loop',
          'SMEFT_RG', 'nu_observables', '18_real_Yukawa_fit')
OSCILLATION_KEYS = ('delta_m21_sq_ev2', 'delta_m3l_sq_ev2',
                    'sin2_theta12', 'sin2_theta13', 'sin2_theta23')


class CompletionBlocked(ValueError):
    """The requested full-general EFT calculation lacks a source-required stage."""


def insist(ok: bool, reason: str) -> None:
    if not ok:
        raise CompletionBlocked(reason)


def digest(path: Path) -> str:
    p = Path(path)
    insist(p.is_file(), 'Missing required file: ' + str(p))
    h = sha256()
    with p.open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def put_json(path: Path, obj: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + '.partial')
    tmp.write_text(json.dumps(obj, sort_keys=True, indent=2, allow_nan=False)
                   + '\n', encoding='utf-8')
    tmp.replace(path)


def from_complex_json(value: Any, *, label: str, shape=(3, 3)) -> np.ndarray:
    """Decode native complex list or {real,imag} matrix; refuse invalid shapes."""
    if isinstance(value, dict) and set(value) == {'real', 'imag'}:
        re = np.asarray(value['real'], dtype=float)
        im = np.asarray(value['imag'], dtype=float)
        insist(re.shape == im.shape, label + ': real/imag shapes differ')
        arr = re + 1j * im
    else:
        arr = np.asarray(decode_native(value), dtype=np.complex128)
    insist(arr.shape == tuple(shape) and np.isfinite(arr).all(),
           label + ': expected finite ' + str(shape) + ' complex matrix')
    return np.asarray(arr, dtype=np.complex128)


def require_oscillation_target(obj: Mapping[str, Any]) -> dict:
    insist(obj.get('ordering') == 'NO', 'Five-observable fit assumes NO')
    centers = obj.get('central_values', {})
    errors = obj.get('one_sigma_errors', {})
    insist(all(k in centers and k in errors for k in OSCILLATION_KEYS),
           'NuFIT five-observable target incomplete')
    for k in OSCILLATION_KEYS:
        insist(np.isfinite(float(centers[k])) and
               np.isfinite(float(errors[k])) and float(errors[k]) > 0,
               'Nonfinite or invalid target for ' + k)
    return {'central_values': centers, 'one_sigma_errors': errors,
            'dataset': obj.get('dataset')}


def score_nufit(prediction: Mapping[str, Any], target: Mapping[str, Any]) -> dict:
    """Five Gaussian residuals, not the official correlated NuFIT likelihood."""
    t = require_oscillation_target(target)
    req = dict(prediction)
    if 'delta_m3l_sq_ev2' not in req and 'dm31_sq_eV2' in req:
        req['delta_m3l_sq_ev2'] = req['dm31_sq_eV2']
    if 'delta_m21_sq_ev2' not in req and 'dm21_sq_eV2' in req:
        req['delta_m21_sq_ev2'] = req['dm21_sq_eV2']
    if req.get('mixing_defined') is False:
        raise CompletionBlocked('Neutrino PMNS mixing is numerically undefined')
    residuals = {}
    for k in OSCILLATION_KEYS:
        try:
            v = float(req[k])
        except (KeyError, ValueError, TypeError) as exc:
            raise CompletionBlocked('Missing predicted NuFIT observable ' + k) from exc
        insist(np.isfinite(v), 'Nonfinite observable ' + k)
        residuals[k] = float((v - t['central_values'][k]) / t['one_sigma_errors'][k])
    return {'order': list(OSCILLATION_KEYS), 'normalized_residuals': residuals,
            'chi2_gaussian_uncorrelated': float(sum(x*x for x in residuals.values())),
            'dataset': t['dataset']}


@dataclass(frozen=True)
class ExistingEvidence:
    model: str
    stage_statuses: dict[str, str]
    sources_sha256: dict[str, str]
    scalar_at_ms: dict[str, Any]
    rg_at_ms: dict[str, np.ndarray]
    ms_gev: float
    mf_gev: float


def existing_evidence(model: str, batch1: Path, raw: Path, laurent: Path) -> ExistingEvidence:
    """Validate saved outputs directly; do not infer a native/Matchete map."""
    insist(model in MODELS, 'Unknown model')
    b1 = load_batch1(batch1)
    insist(b1.model == model, 'Batch1 model mismatch')
    checked = check_saved_raw(model, batch1, raw)
    laurent_manifest = read_json(laurent / 'manifest.json')
    insist(laurent_manifest.get('model') == model and
           laurent_manifest.get('status') ==
           'FullOneLoopFormalLaurentAudited_TreeMassRemainderPreserved_NotFiniteWilsons',
           model + ': not a successful Laurent v2 run')
    insist(laurent_manifest.get('batch1_manifest_sha256') ==
           checked['batch1_manifest_sha256'] and
           laurent_manifest.get('raw_manifest_sha256') ==
           checked['raw_manifest_sha256'],
           model + ': Laurent and source provenance mismatch')
    audit = check_audit_report(read_json(laurent / 'laurent_audit.json'),
                               checked['one_loop_source_sha256'], laurent)
    src = b1.scalar_matching_input()
    insist(src.get('model') == model and src.get('no_heavy_fermions_active') is True
           and src.get('both_heavy_scalar_multiplets_to_integrate_out_together') is True,
           model + ': scalar matching handover violates required sequence')
    endpoints = b1.at_scalar_threshold()
    mapped = src.get('native_eft1_rgbeta_couplings_at_muS', {})
    insist(set(endpoints) == set(mapped), 'Lost scalar threshold RG inventory')
    for name, value in mapped.items():
        insist(np.allclose(endpoints[name], decode_native(value), rtol=1e-13, atol=1e-11),
               model + ': changed RG scalar-threshold state ' + name)
    spectrum = src.get('candidate_spherical_scalar_quadratic_spectrum', {})
    insist(spectrum.get('analytic_spectrum_verified') is True and
           spectrum.get('included_in_threshold_matching') is False,
           model + ': scalar masses must remain candidate until physical matching')
    evs = np.asarray(spectrum.get('mass_squared_eigenvalues_gev2', ()), float)
    insist(evs.shape == (8,) and np.isfinite(evs).all() and np.min(evs) > 0,
           model + ': invalid candidate heavy scalar squared masses')
    if model == 'D_-2':
        extras = [x for x in endpoints if x.startswith('lambdaT3GeneralQ')]
        insist(len(extras) == 5 and 'muT3GeneralS2Sq' in endpoints
               and src.get('quadratic_holomorphic_mass_symbol') == 'muT3GeneralS2Sq',
               'D(-2): fifth quartic or self-bilinear missing')
        insist(spectrum.get('additional_quadratic_structures') == 1,
               'D(-2): scalar self-bilinear omitted')
    else:
        insist(not any(x.startswith('lambdaT3GeneralQ') for x in endpoints)
               and spectrum.get('additional_quadratic_structures') == 0,
               'A(-4): unexpected additional structure')
    scale = b1.manifest['scales_gev']
    insist(scale == {'UV':1e7,'F':1e5,'S':1e3,'low':100},
           model + ': unexpected thesis scales')
    statuses = {
        'UV_general_RG': 'CompletedGeneralRenormalisableRG',
        'F_all_fermion_raw_one_loop': 'CompletedUnrenormalisedMatcheteFOnly',
        'F_Laurent': 'CompletedFormalLaurentNoSubtraction',
    }
    hashes = {
        'batch1_manifest': digest(batch1/'manifest.json'),
        'batch1_scalar_input': digest(batch1/'scalar_matching_input.json'),
        'raw_manifest': digest(raw/'manifest.json'),
        'raw_matching_ledger': checked['raw_ledger_sha256'],
        'laurent_manifest': digest(laurent/'manifest.json'),
        'laurent_audit': digest(laurent/'laurent_audit.json'),
        'laurent_unsubtracted_regular': audit['wolfram_output_sha256'].get(
            'epsilon_zero_regular_unsubtracted.wxf'),
    }
    return ExistingEvidence(model, statuses, hashes, src, endpoints,
                            float(scale['S']), float(scale['F']))


def verify_completed_c5_document(path: Path, evidence: ExistingEvidence) -> tuple[np.ndarray, dict]:
    """Require a truly independent finite physical matching input.

    Metadata-only declarations do not prove the physics of a source, so this
    interface does not call the calculation 'certified'. The input MUST
    provide an actual numerical MSbar coefficient at the scalar scale and
    all pertinent source hashes; the saved Laurent regular part alone fails.
    """
    document = read_json(path)
    insist(document.get('format') == C5_FORMAT, 'Wrong renormalized C5 input format')
    insist(document.get('model') == evidence.model, 'C5 model mismatch')
    insist(document.get('status') == 'SourceDerivedFiniteFullGeneralC5AtMS',
           'C5 is not a completed finite full-general matching result')
    insist(document.get('renormalization_scheme') == 'MSbar'
           and document.get('fermion_threshold_subtraction') == 'UV_hard_minus_EFT1_soft_plus_counterterms',
           'Finite UV/EFT1 subtraction and MSbar counterterms not supplied')
    insist(document.get('one_loop_full_EFT1_wilson_running_performed') is True
           and document.get('one_loop_joint_scalar_threshold_matching_performed') is True,
           'EFT1 Wilson RG and BOTH-scalar one-loop matching must be complete')
    insist(document.get('integrated_heavy_fields_in_order') ==
           [['F'], ['S1', 'S2']], 'EFT field-decoupling order wrong')
    insist(float(document.get('scale_GeV', 0)) == evidence.ms_gev,
           'C5 scale mismatch')
    hashes = document.get('input_sha256', {})
    for key in ('batch1_manifest','raw_manifest','laurent_manifest'):
        insist(hashes.get(key) == evidence.sources_sha256[key],
               'C5 source hash differs from completed RG/matching input: ' + key)
    # Hash an independently generated symbolic scalar-matching source, not
    # just its user-supplied file name.
    sources = document.get('symbolic_sources', {})
    insist(isinstance(sources, dict) and 'scalar_matching_eft_wxf' in sources,
           'Missing complete source-derived scalar threshold matching expression')
    source = sources['scalar_matching_eft_wxf']
    insist(isinstance(source, dict) and 'file' in source and 'sha256' in source,
           'Scalar matching source record incomplete')
    symbolic_path = Path(source['file'])
    if not symbolic_path.is_absolute():
        symbolic_path = path.parent / symbolic_path
    insist(digest(symbolic_path) == source['sha256'],
           'Scalar matching source WXF hash mismatch')
    c5 = from_complex_json(document.get('C5_MS_GeV_inverse'),label='MS C5')
    insist(np.allclose(c5,c5.T,rtol=1e-10,atol=1e-22),
           'Weinberg Wilson C5 must be complex symmetric, not Hermitian')
    return (c5+c5.T)/2, document


def finish_smeft(evidence: ExistingEvidence, c5: np.ndarray, *,
                 vev: float, target: Mapping[str, Any]) -> dict:
    """Actual one-loop SMEFT RGE and Takagi/PMNS with project-native solvers."""
    from Numerical.running.SMWeinbergEvolution import (
        SMWeinbergInitialConditions, evolve_sm_weinberg,
    )
    from Numerical.parallel_t3.threshold.rough_numerics import neutrino_observables
    from Numerical.parallel_t3.background.general_batch1_api import encode_native
    s = evidence.rg_at_ms
    initial = SMWeinbergInitialConditions(
        gY=float(np.real(s['gY'])),g2=float(np.real(s['g2'])),
        g3=float(np.real(s['g3'])),lambdaH=float(np.real(s['lambdaH'])),
        ye=np.asarray(s['ye'],complex),yu=np.asarray(s['yu'],complex),
        yd=np.asarray(s['yd'],complex),K=c5,
    )
    low = evolve_sm_weinberg(initial, evidence.ms_gev, 100.)
    insist(low.solver_success, 'SMEFT one-loop RG failed')
    mass=-.5*float(vev)**2*np.asarray(low.K,complex)
    observables, mixing_matrices = neutrino_observables(mass, low.ye)
    score = score_nufit(observables, target)
    return {
        'status': 'NumericalSMEFTEvolvedAndNeutrinoObservablesComputed',
        'model': evidence.model,
        'C5_MS_GeV_inverse':encode_native(c5),
        'C5_100GeV_GeV_inverse':encode_native(low.K),
        'mnu_100GeV_eV':encode_native(mass*1e9),
        'observables':observables,
        'nufit_gaussian_approximation': score,
        'PMNS_matrices': {name:encode_native(v) for name,v in mixing_matrices.items()},
        'smeft_solver_nfev':low.nfev,
        'smeft_scheme':'SM + Weinberg, one-loop existing repository solver',
        'physical_normalization_independently_certified':False,
        'source_contract': evidence.sources_sha256,
    }
