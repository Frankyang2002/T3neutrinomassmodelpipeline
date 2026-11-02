"""Fail-closed parser for *raw* Matchete simultaneous-F one-loop threshold evidence.

This module does NOT convert a raw Matchete result into an MSbar Wilson
coefficient, and does not infer any value for the GammaCC -> Weyl phase.
All outputs from its batch driver belong under the owned spinor output tree.
"""
from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
from typing import Any

from .boundary import MODELS, BoundaryError, require, tag

MIXED = {'B_-1', 'C_-1'}
ORDINARY_DRIVER = 'Lagrangian/RunT3GeneralFThreshold.wl'
MIXED_DRIVER = 'Lagrangian/RunGeneralT3FermionFirstEFT1OneLoop.wl'
MIXED_INPUT = 'output/general_t3_mass_eigen_matchete/{tag}_final_parameters.json'
FORMAT = 't3-f-only-raw-one-loop-native-evidence-v1'


def sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def model_command(root: Path, model: str, directory: Path,
                  executable: str = 'wolframscript') -> list[str]:
    """Native, *single-F-threshold*, one-process-per-model command.

    The drivers are existing repository code. No new matching conventions are
    introduced or guessed here. The two B/C models use their saved numeric
    mixed-doublet mass-eigenstate input; do not use the ordinary driver.
    """
    require(model in MODELS, f'Unknown model: {model}')
    directory = directory.resolve()
    root = root.resolve()
    if model in MIXED:
        driver = root / MIXED_DRIVER
        source = root / MIXED_INPUT.format(tag=tag(model))
        require(source.is_file(), f'{model}: missing mixed-mass input: {source}')
        obj = json.loads(source.read_text(encoding='utf-8-sig'))
        require(obj.get('model') == model and
                obj.get('format') == 't3-experimental-mass-eigen-matchete-parameters-v1' and
                obj.get('status') == 'NumericMassEigenstateBoundaryPreparedMatcheteNotRun',
                f'{model}: wrong saved mixed-mass input schema/status')
        require(obj.get('quartic_source_count') == 17 and
                obj.get('yukawa_source_count') == 2,
                f'{model}: incomplete mixed source interaction inventory')
        return [executable, '-file', str(driver), str(source),
                str(directory / f'{tag(model)}_native_one_loop.json')]
    return [executable, '-file', str(root / ORDINARY_DRIVER),
            str(directory), model, '1']


def preflight(root: Path, model: str) -> dict:
    """Return a hash-stamped read-only native driver/input plan."""
    cmd = model_command(root, model, root / 'output/general_observable/lo_certification_spinor/dryrun')
    driver = root / (MIXED_DRIVER if model in MIXED else ORDINARY_DRIVER)
    require(driver.is_file(), f'{model}: source driver missing: {driver}')
    return {
        'model': model, 'driver': driver.relative_to(root).as_posix(),
        'driver_sha256': sha(driver),
        'mass_basis_input': (str((root / MIXED_INPUT.format(tag=tag(model))).relative_to(root))
                             if model in MIXED else None),
        'mass_basis_input_sha256': (sha(root / MIXED_INPUT.format(tag=tag(model)))
                                    if model in MIXED else None),
        'evidence_kind': 'native_matchete_run_required',
        'integration': 'all BSM fermions simultaneously at mu_F; both scalar multiplets active',
    }


def _read_json(path: Path) -> dict:
    require(path.is_file(), f'Native Matchete result missing: {path}')
    data = json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(data, dict), f'Native result is not JSON object: {path}')
    return data


def _owned_filename(directory: Path, raw: str) -> Path:
    """Accept only basename output declarations, never absolute/parent paths."""
    require(isinstance(raw, str) and raw and Path(raw).name == raw and
            '/' not in raw and '\\' not in raw and raw not in ('.', '..'),
            'Unsafe native output filename')
    path = (directory / raw).resolve()
    require(path.is_relative_to(directory.resolve()), 'Native file escaped evidence directory')
    return path


def inspect_raw(model: str, directory: Path) -> dict:
    """Inspect actual native results; *never* certify a renormalized C5.

    Ordinary matching's native status says the match was performed at
    LoopOrder=1 but does NOT certify hard LLHH extraction. B/C additionally
    export a raw LLHH candidate sector. A 0 candidate does not prove zero
    physical hard matching: selection/poles/scheme still require checking.
    """
    require(model in MODELS, 'Unknown model')
    directory = directory.resolve()
    if model not in MIXED:
        path = directory / 'native_status.json'
        obj = _read_json(path)
        require(obj.get('model') == model and obj.get('loop_order') == 1 and
                obj.get('status') == 'FThresholdTreeWilsonExported' and
                obj.get('raw_one_loop_F_matching_included') is True and
                obj.get('hard_one_loop_MSbar_certified') is False and
                obj.get('integrated_fields') == ['F'] and
                obj.get('light_bsm_fields') == ['S1', 'S2'] and
                obj.get('heavy_flavour_count') == 3,
                f'{model}: native general F-only LoopOrder=1 evidence invalid')
        # This is the actual matching-generated one-loop *correction* export,
        # not the loop-order flag. Authenticate its presence and hash; do not
        # split Mathematica InputForm on '+' (nested expressions and signs).
        raw = directory / 'data' / 'eft1_after_F_one_loop_inputform.txt'
        require(raw.is_file() and raw.stat().st_size > 0,
                f'{model}: one-loop InputForm export missing; status flag alone insufficient')
        text = raw.read_text(encoding='utf-8-sig').strip()
        require(text and 'Field[NewFermion' not in text,
                f'{model}: exported loop sector empty or still contains heavy F')
        return {
            'status': 'MATCHeTE_F_ONLY_LOOP_ORDER_ONE_RETURNED_RAW_SECTOR_UNRESOLVED',
            'driver_type': 'general_14',
            'native_report_file': path.name,
            'native_report_sha256': sha(path),
            'one_loop_F_native_run_returned': True,
            'raw_loop_inputform_file': 'data/eft1_after_F_one_loop_inputform.txt',
            'raw_loop_inputform_sha256': sha(raw),
            'raw_loop_nonzero': text != '0',
            'raw_loop_contains_GammaCC': 'GammaCC' in text,
            'raw_loop_contains_Higgs_field': 'Field[H,' in text,
            'explicit_hard_LLHH_candidate_terms': None,
            'hard_sector_renormalized': False,
            'physical_wilson_coefficient': None,
        }
    path = directory / f'{tag(model)}_native_one_loop.json'
    obj = _read_json(path)
    require(obj.get('model') == model and
            obj.get('format') == 't3-mass-eigen-fermion-first-eft1-oneloop-v1' and
            obj.get('status') == 'FermionFirstOneLoopRawSectorsExtractedNotRenormalized' and
            obj.get('matching_EFT_order') == 5 and obj.get('matching_loop_order') == 1 and
            obj.get('matching_approach') == 'FermionFirstTreeAndOneLoopScalarsRemainLight' and
            obj.get('integrated_heavy_fields') == ['NewFermion'] and
            obj.get('scalar_light_fields') == ['NewScalarX', 'NewScalarY'] and
            obj.get('one_loop_matching_completed') is True and
            obj.get('loop_difference_hbar_present') is True and
            obj.get('loop_difference_fermion_absent') is True and
            obj.get('raw_hard_llhh_scheme_renormalized') is False and
            obj.get('raw_hard_llhh_pole_rge_certified') is False and
            obj.get('physical_C5_computed') is False,
            f'{model}: mixed native raw one-loop evidence invalid')
    n = obj.get('raw_hard_llhh_candidate_term_count')
    require(isinstance(n, int) and not isinstance(n, bool) and n >= 0,
            f'{model}: invalid raw hard candidate count')
    for key in ('one_loop_difference_file', 'raw_hard_llhh_file'):
        artifact = _owned_filename(directory, obj.get(key))
        require(artifact.is_file() and artifact.stat().st_size > 0,
                f'{model}: missing exported raw Matchete term: {artifact}')
    hard = _owned_filename(directory, obj['raw_hard_llhh_file'])
    hard_text = hard.read_text(encoding='utf-8-sig').strip()
    require('Field[NewFermion' not in hard_text,
            f'{model}: heavy fermion not eliminated from hard candidate')
    if n == 0:
        require(hard_text == '0', f'{model}: zero candidates but raw term is not 0')
    else:
        require(hard_text != '0' and 'GammaCC' in hard_text,
                f'{model}: candidate count but missing GammaCC LLHH terms')
    return {
        'status': 'MATCHeTE_RAW_MIXED_LOOP_LLHH_EXTRACTED_NOT_MSbar',
        'driver_type': 'mixed_BC',
        'native_report_file': path.name,
        'native_report_sha256': sha(path),
        'one_loop_F_native_run_returned': True,
        'explicit_hard_LLHH_candidate_terms': n,
        'raw_hard_llhh_file': hard.name,
        'raw_hard_llhh_sha256': sha(hard),
        'hard_sector_renormalized': False,
        'physical_wilson_coefficient': None,
    }


def blocked_record(model: str, reason: str, run_return_code: int | None = None) -> dict:
    return {
        'status': 'BLOCKED', 'model': model, 'reason': reason,
        'native_exit_code': run_return_code,
        'hard_sector_renormalized': False,
        'physical_wilson_coefficient': None,
    }


PHASE_DRIVER = 'Lagrangian/RunT3YukawaSpinorPhaseNative.wl'
PHASE_KEYS = (
    'gamma_cc_square_minus_one',
    'gamma_cc_transpose_antisymmetric',
    'cconj_field_definition_in_bilinear',
    'cconj_flips_left_projector',
    'hermitian_conjugate_matches_explicit_cconj',
    'distinct_native_y1_and_y1conj_present',
)


def inspect_native_phase(path: Path) -> dict:
    """Import ABSTRACT Matchete CConj identities; cannot fix absolute C phase."""
    obj = _read_json(path)
    require(obj.get('format') == 't3-workstream-a-a0-native-spinor-algebra-v1'
            and obj.get('model') == 'A_p0'
            and obj.get('status') == 'NativeMatcheteAbstractChargeConjugationAndHcPassed',
            'Native A_p0 abstract phase probe failed or wrong schema')
    checks = obj.get('checks')
    require(isinstance(checks, dict) and
            all(checks.get(name) is True for name in PHASE_KEYS),
            'Incomplete Matchete abstract charge conjugation identities')
    return {
        'status': 'ABSTRACT_MATCHeTE_CHARGE_CONJUGATION_VERIFIED_PHASE_OPEN',
        'native_probe_sha256': sha(path),
        'native_checks_passed': len(PHASE_KEYS),
        'GammaCC_to_master_Weyl_factor': None,
        'native_connected_UV_amplitude_certified': False,
        'numerical_C_matrix_phase_certified': False,
        'physical_master_wilson_boundary': None,
    }


FLAGS = (
    'T3_EXPERIMENTAL_GENERAL_FERMIONS',
    'T3_EXPERIMENTAL_GENERAL_QUARTICS',
    'T3_EXPERIMENTAL_GENERAL_QUADRATICS',
    'T3_EXPERIMENTAL_FLAVOR_MASSES',
)


def native_environment(model: str, parent: dict[str, str]) -> dict[str, str]:
    """Preserve each *existing* driver's independently checked UV build mode.

    The mixed eigenfield driver reconstructs its own seven extra scalar
    quartics and explicitly rejects GENERAL_QUARTICS=1. Keep all baseline
    experimental flags off for those B/C models; a fresh kernel prevents
    state cross-talk. The general nonmixed driver demands all four=1.
    """
    require(model in MODELS, 'Unknown T3 model')
    result = dict(parent)
    setting = '0' if model in MIXED else '1'
    result.update({key: setting for key in FLAGS})
    return result


def phase_environment(parent: dict[str, str]) -> dict[str, str]:
    """Pre-existing A_p0 abstract Matchete phase probe requires this mode."""
    result = dict(parent)
    result.update({
        'T3_EXPERIMENTAL_GENERAL_FERMIONS': '1',
        'T3_EXPERIMENTAL_GENERAL_QUARTICS': '1',
        'T3_EXPERIMENTAL_GENERAL_QUADRATICS': '0',
        'T3_EXPERIMENTAL_FLAVOR_MASSES': '1',
    })
    return result
