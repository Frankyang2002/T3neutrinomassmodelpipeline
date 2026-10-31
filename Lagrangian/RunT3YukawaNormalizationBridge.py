#!/usr/bin/env python3
"""Read-only Workstream A: exact 38-Yukawa CG norm and RGBeta prefactor gate.

Consumes live native Matchete CG exports and the successful component Ward
report. Optionally runs the isolated RGBeta companion in fresh Wolfram kernels
for eight representative model points (five representation classes; all six
exceptional Yukawas). Does not touch production RGEs, matching, or UV results.

What it verifies:
  * the exact native Matchete CG Hilbert--Schmidt norm for each of 38 Yukawas;
  * that the SU(2) singlet reference has squared norm 2 and the unit-normalised
    doublet x doublet -> triplet reference has squared norm 3;
  * candidate magnitude 1 (singlet) or sqrt(2)/3^(1/4) (triplet);
  * optionally, the *internal* RGBeta complex-scalar Yukawa tensor coefficient
    sqrt(2) multiplying Matrix[Bar[y]] for Chirality -> Right.

What it DOES NOT verify:
  The physical Matchete--RGBeta coupling dictionary (phase, complex
  conjugation, relative Lagrangian sign, spinor bilinear normalisation, global
  triplet component basis). A native RGBeta tensor coefficient is NOT the
  coefficient of a physical component interaction.

Examples from repository root:
  python Lagrangian/RunT3YukawaNormalizationBridge.py --self-test
  python Lagrangian/RunT3YukawaNormalizationBridge.py --all
  python Lagrangian/RunT3YukawaNormalizationBridge.py --all --no-wolfram
  python Lagrangian/RunT3YukawaNormalizationBridge.py --source-only
"""
from __future__ import annotations

import argparse
from collections import Counter
import json
import os
from pathlib import Path
import subprocess
import sys

import sympy as sp

from RunT3YukawaComponentWard import (
    EXTRAS, ROOT, model_labels, require, tensor_from_record,
    norm_squared,
)

REPORT_REL = Path('Reports/output/workstream_a/yukawa_normalization_bridge')
WARD_REL = Path('Reports/output/workstream_a/yukawa_component_ward/yukawa_component_ward.json')
NATIVE_REL = Path('Reports/output/workstream_a/native_cg_bridge/native')
EXTRAS_REL = Path('Reports/output/workstream_a/yukawa_component_ward/native_extras')
TRIPLET_REL = Path('output/general_t3_stage22/stage28_live_rgbeta_triplet_structure.json')
PROBE_NAMES = ('A_m2', 'A_p0', 'B_m1', 'C_m1', 'D_m2', 'D_p0', 'E_m2', 'E_p0')
PROBE_SUCCESS = 'NativeRGBetaYukawaInternalPrefactorsVerified'
CG_SCALE_TRIPLET = sp.sqrt(2) / sp.root(3, 4)


def read_json(path: Path) -> dict:
    require(path.is_file(), f'Missing required evidence: {path}')
    with path.open(encoding='utf-8') as stream:
        obj = json.load(stream)
    require(isinstance(obj, dict), f'Invalid object at {path}')
    return obj


def exact_equal(a: sp.Expr, b: sp.Expr) -> bool:
    return sp.simplify(a - b) == 0


def field_dims(model: str, which: int) -> tuple[int, int, int]:
    s1, s2, df, _ = model_labels()[model]
    return 2, df, (s1 if which == 1 else s2)


def native_record(root: Path, model: str, which: int, extra: bool) -> dict:
    if extra:
        p = root / EXTRAS_REL / f'{model}.json'
        data = read_json(p)
        require(data.get('status') == 'LiveExtraCGEqualityVerified'
                and data.get('model') == model
                and data.get('extra_equals_baseline_components_exactly') is True,
                f'{model}: extra Matchete CG not live-verified')
        return data['extra_tensor']
    p = root / NATIVE_REL / f'{model}.json'
    data = read_json(p)
    require(data.get('status') == 'LiveRegisteredBaselineCGsCaptured'
            and data.get('model') == model
            and data.get('full_matchete_check_passed') is True,
            f'{model}: baseline Matchete CG missing native proof')
    return data['native_tensors'][f'Yukawa{which}']


def cg_info(root: Path, model: str, which: int, extra: bool) -> dict:
    # Recompute independently from native CG vectors rather than trusting the
    # norm reported by the previous Ward stage.
    dims = field_dims(model, which)
    record = native_record(root, model, which, extra)
    name = f'{model}/y{which}' + ('Conj' if extra else '')
    values = tensor_from_record(record, dims, name)
    nsq = norm_squared(values)
    has_triplet = 3 in dims
    reference_norm_sq = sp.Integer(3 if has_triplet else 2)
    expected_norm_sq = 2 * sp.sqrt(3) if has_triplet else sp.Integer(2)
    require(exact_equal(nsq, expected_norm_sq),
            f'{name}: unexpected Matchete CG norm {nsq}; expected {expected_norm_sq}')
    ratio_abs = sp.simplify(sp.sqrt(nsq / reference_norm_sq))
    target = CG_SCALE_TRIPLET if has_triplet else sp.Integer(1)
    require(exact_equal(ratio_abs, target), f'{name}: exact CG magnitude mismatch')
    return {
        'model': model,
        'coupling': f'y{which}' + ('Conj' if extra else ''),
        'native_matchete_cg_norm_squared': str(nsq),
        'unit_reference_cg_norm_squared': str(reference_norm_sq),
        'conditional_absolute_cg_ratio': 'sqrt(2)/3**(1/4)' if has_triplet else '1',
        'representation_geometry': 'symmetric_triplet' if has_triplet else 'doublet_singlet',
        'reference_triplet_basis_is_independently_rgbeta_native': False,
        'physical_coupling_conversion_certified': False,
    }


def check_ward(root: Path) -> dict:
    ward = read_json(root / WARD_REL)
    require(ward.get('format') == 't3-workstream-a-exact-yukawa-component-ward-v1'
            and ward.get('status') == 'ExactSU2U1ComponentWardPassedPhysicalCrossPackageDictionaryPending'
            and ward.get('models_checked') == 16
            and ward.get('baseline_yukawas_checked') == 32
            and ward.get('new_native_extra_yukawas_checked') == 6
            and ward.get('su2_ward_component_equations_checked') == 1032
            and ward.get('physical_matchete_rgbeta_yukawa_normalization_certified') is False,
            'Previous component Ward audit is not the expected live 38-Yukawa result')
    require(len(ward.get('yukawas', [])) == 38
            and all(row.get('passed') is True for row in ward['yukawas']),
            'Incomplete Ward audit component rows')
    return ward


def check_stage28(root: Path) -> dict:
    data = read_json(root / TRIPLET_REL)
    require(data.get('models_checked') == 6
            and data.get('all_native_delS2_contractions_passed') is True
            and data.get('all_group_representation_checks_passed') is True
            and data.get('native_numerical_triplet_generators_available') is False
            and data.get('physical_coupling_mapping_certified') is False,
            'Stage28 native abstract symmetric-triplet evidence changed')
    return data


def probe_rgbeta(root: Path, timeout: int) -> None:
    script = root / 'Lagrangian/RunT3YukawaRGBetaNormalizationProbe.wl'
    require(script.is_file(), f'Missing native RGBeta probe: {script}')
    outdir = root / REPORT_REL / 'native_rgbeta'
    outdir.mkdir(parents=True, exist_ok=True)
    for label in PROBE_NAMES:
        family, tag = label.split('_')
        output = outdir / f'{label}.json'
        # A failed Wolfram process must never be diagnosed from stale JSON.
        if output.exists():
            output.unlink()
        process = subprocess.run(
            ['wolframscript', '-file', str(script), str(output), family, tag],
            cwd=root, env=os.environ.copy(), text=True,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            timeout=timeout, check=False,
        )
        (outdir / f'{label}.log').write_text(process.stdout, encoding='utf-8')
        if process.returncode != 0:
            details = ''
            if output.is_file():
                try:
                    failed = read_json(output)
                    details = (f"; step={failed.get('failed_step', 'unknown')}"
                               f"; reason={failed.get('reason', 'unknown')}")
                except (ValueError, OSError, json.JSONDecodeError):
                    details = '; native JSON was not readable'
            require(False,
                    f'{label} RGBeta native prefactor probe failed (exit {process.returncode})'
                    f'{details}; see {outdir / (label + ".log")}')
        print(f'Native RGBeta Yukawa prefactor: {label} PASS', flush=True)


def check_native_probes(root: Path) -> list[dict]:
    rows = []
    for label in PROBE_NAMES:
        data = read_json(root / REPORT_REL / 'native_rgbeta' / f'{label}.json')
        require(data.get('format') == 't3-rgbeta-native-yukawa-prefactor-v1'
                and data.get('status') == PROBE_SUCCESS
                and data.get('model') == label
                and data.get('internal_real_scalar_tensor_prefactor') == 'Sqrt[2]'
                and data.get('physical_action_equivalence_verified') is False
                and data.get('production_modified') is False
                and data.get('rgbeta_native_ward_invoked') is False
                and data.get('registration_scope') == 'IsolatedGaugeFieldsAndYukawas'
                and data.get('masses_and_quartics_registered') is False,
                f'{label}: unexpected RGBeta native probe status/scope')
        baseline = data.get('baseline_yukawas', [])
        require(len(baseline) == 2 and {x.get('coupling') for x in baseline} == {'y1', 'y2'},
                f'{label}: baseline couplings missing')
        extra = data.get('extra_yukawa')
        if label in EXTRAS:
            require(isinstance(extra, dict) and extra.get('coupling') ==
                    ('y1Conj' if label.endswith('_p0') else 'y2Conj'),
                    f'{label}: extra coupling missing')
        else:
            require(extra is None, f'{label}: unexpected extra Yukawa')
        all_records = baseline + ([extra] if extra else [])
        for item in all_records:
            require(item.get('chirality') == 'Right'
                    and item.get('coefficient_of_Matrix_Bar_y') in ('Sqrt[2]', '2^Rational[1, 2]')
                    and item.get('coefficient_of_Matrix_y_in_CouplingBar') in ('Sqrt[2]', '2^Rational[1, 2]')
                    and item.get('native_complex_scalar_yukawa_tensor_prefactor_verified') is True,
                    f'{label}: native RGBeta tensor-normalization probe failed')
        rows.append({'model': label, 'registered_yukawas_checked': len(all_records),
                     'internal_real_scalar_tensor_factor': 'sqrt(2)',
                     'rgbeta_native_ward_invoked': False,
                     'physical_action_equivalence_verified': False})
    return rows


def compile_report(root: Path, native: bool) -> dict:
    check_ward(root)
    check_stage28(root)
    all_rows = []
    for model in model_labels():
        for which in (1, 2):
            all_rows.append(cg_info(root, model, which, False))
        if model in EXTRAS:
            alpha = model_labels()[model][-1]
            all_rows.append(cg_info(root, model, 1 if alpha == 0 else 2, True))
    require(len(all_rows) == 38, 'Wrong total Yukawa count')
    counts = Counter(x['representation_geometry'] for x in all_rows)
    require(counts == {'doublet_singlet': 14, 'symmetric_triplet': 24},
            f'Unexpected representation inventory {counts}')
    probes = check_native_probes(root) if native else []
    return {
        'format': 't3-yukawa-normalization-bridge-stage1-v1',
        'status': ('ExactCGMagnitudesAndNativeRGBetaInternalPrefactorsVerifiedActionDictionaryPending'
                   if native else 'ExactCGMagnitudesVerifiedRGBetaInternalProbePending'),
        'models_checked': 16,
        'native_matchete_yukawas_checked': 38,
        'singlet_geometry_count': counts['doublet_singlet'],
        'triplet_geometry_count': counts['symmetric_triplet'],
        'su2_ward_equations_previously_verified': 1032,
        'exact_absolute_cg_ratio_singlet': '1',
        'conditional_exact_absolute_cg_ratio_triplet': 'sqrt(2)/3**(1/4)',
        'rgbeta_native_internal_prefactor_verified': native,
        'rgbeta_internal_real_scalar_prefactor': 'sqrt(2)' if native else None,
        'rgbeta_native_probe_models': len(probes),
        'native_rgbeta_yukawa_records_checked': sum(p['registered_yukawas_checked'] for p in probes),
        'matchete_rgbeta_physical_action_verified': False,
        'matchete_rgbeta_absolute_coupling_dictionary_certified': False,
        'native_rgbeta_ward_identity_certified': False,
        'native_rgbeta_triplet_component_basis_certified': False,
        'production_modified': False,
        'remaining': [
            'Fix one global gauge-field component map, including SU2 pseudoreal and triplet Bar orientation',
            'Match Matchete NCM/PlusHc bilinears to RGBeta left-handed Weyl bilinears',
            'Reconstruct full complex-to-real scalar action with signs, H.c. and identical-fermion factors',
            'Solve signed complex conjugation-aware physical coupling map and verify every component',
        ],
        'cg_rows': all_rows,
        'rgbeta_native_probe_rows': probes,
    }


def self_test() -> None:
    require(exact_equal(sp.sqrt((2 * sp.sqrt(3)) / 3), CG_SCALE_TRIPLET),
            'Triplet analytic normalization failed')
    require(exact_equal(sp.sqrt(sp.Integer(2) / 2), 1),
            'Singlet analytic normalization failed')
    # Negative controls: a unit off the singlet tensor and incorrectly
    # normalised triplet source tensor must fail.
    require(not exact_equal(sp.Integer(3), sp.Integer(2)), "Singlet negative control passed unexpectedly")
    require(not exact_equal(sp.Integer(3), 2 * sp.sqrt(3)), "Triplet negative control passed unexpectedly")
    print('PASS: exact singlet/triplet CG magnitudes and normalization-negative controls')


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root', type=Path, default=ROOT)
    parser.add_argument('--self-test', action='store_true')
    parser.add_argument('--all', action='store_true', help='Run 8 native RGBeta probes and build the full report')
    parser.add_argument('--no-wolfram', action='store_true', help='Use 8 previously saved native RGBeta probe files')
    parser.add_argument('--source-only', action='store_true', help='Exact CG census; do not claim RGBeta internal prefactor verification')
    parser.add_argument('--timeout', type=int, default=600, help='Seconds per native RGBeta probe')
    args = parser.parse_args()
    try:
        if args.self_test:
            self_test()
            return 0
        require(args.source_only != args.all, 'Choose exactly one of --all or --source-only')
        require(not args.no_wolfram or args.all, '--no-wolfram requires --all')
        root = args.repo_root.resolve()
        if args.all and not args.no_wolfram:
            probe_rgbeta(root, args.timeout)
        report = compile_report(root, native=args.all)
        output = root / REPORT_REL / 'yukawa_normalization_bridge.json'
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
        print(f"{report['status']}: {report['native_matchete_yukawas_checked']} CGs; "
              f"{report['singlet_geometry_count']} singlet, {report['triplet_geometry_count']} triplet")
        print('Absolute physical action/coupling dictionary: NOT CERTIFIED')
        print(f'Report: {output}')
        return 0
    except (ValueError, OSError, json.JSONDecodeError, subprocess.TimeoutExpired, KeyError) as exc:
        print('Yukawa normalization bridge BLOCKED: ' + str(exc), file=sys.stderr)
        return 2


if __name__ == '__main__':
    sys.exit(main())
