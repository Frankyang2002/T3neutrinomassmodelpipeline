"""Experimental ONE-STEP symbolic Matchete dimension-5 matching for B_-1/C_-1.

This is a one-scale all-heavy matching diagnostic, NOT the physical fermion-first
EFT1 threshold sequence and NOT a numerical or fully normalized C5 result.

Requires a successful previous native eigenfield HeavyMassBasis check.

python -m scripts.RunGeneralT3MassEigenOneLoop --model B_-1
python -m scripts.RunGeneralT3MassEigenOneLoop --model C_-1
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
import shutil
import subprocess
from typing import Any, Mapping

from Numerical.experimental.GeneralT3MixedScalarMassBasis import MIXED_MODELS
from scripts.RunGeneralT3MassEigenMatchete import (
    ROOT, INPUT, OUTPUT, _read_json, prepare_model,
    validate_wolfram_output,
)

WOLFRAM = ROOT / 'Lagrangian/RunGeneralT3MassEigenOneLoop.wl'
MATCHING_SOURCE = ROOT / 'Lagrangian/RunMatching.wl'
STATUS = 'OneLoopDirectSymbolicWeinbergProbePassedNumericalC5Pending'
FORMAT = 't3-mass-eigen-one-loop-matching-probe-v1'


def fail(message: str) -> None:
    raise ValueError('T3 exploratory one-loop matching fail-closed: ' + message)


def _stem(model: str, which: str) -> str:
    if model not in MIXED_MODELS or which not in ('initial', 'final'):
        fail(f'Unsupported model/endpoint: {model}, {which}')
    return model.replace('+', 'p').replace('-', 'm') + '_' + which


def _prior_certificate(model: str, *, which: str,
                       certificate_dir: Path, unitary: Any) -> Path:
    """Require the separate already-completed native eigenfield audit."""
    source = certificate_dir / (_stem(model, which) + '_matchete.json')
    report = _read_json(source)
    validate_wolfram_output(report, model, unitary)
    return source


def validate_probe_report(report: Mapping[str, Any], model: str,
                          outdir: Path, stem: str) -> tuple[Path, Path]:
    if not isinstance(report, Mapping) or report.get('status') != STATUS:
        fail(f'{model}: no successful symbolic Weinberg probe: '
             + str(report.get('status') if isinstance(report, Mapping) else report))
    if (report.get('model') != model or
        report.get('mass_basis_ordinary_check_passed') is not True or
        report.get('heavy_mass_basis_detailed_check_passed') is not True or
        report.get('one_loop_match_stage_completed') is not True or
        report.get('weinberg_extraction_status') != 'Success' or
        report.get('matching_EFT_order') != 5 or
        report.get('matching_loop_order') != 1 or
        report.get('matching_approach') != 'SimultaneousAllHeavyOneScaleExploratory' or
        report.get('numerical_c5_computed') is not False or
        report.get('production_sequential_matching_performed') is not False or
        report.get('physical_matching_performed') is not False or
        report.get('ordered_coefficient_is_canonical_flavour_symmetric_c5') is not False or
        report.get('matchete_rgbeta_quartic_normalization_certified') is not False):
        fail(model + ': symbolic matching result overclaims physical certification or lacks a required stage')
    n = report.get('weinberg_holomorphic_term_count')
    if not isinstance(n, int) or isinstance(n, bool) or n <= 0:
        fail(model + ': holomorphic Weinberg operator sector is empty')
    names = (
        ('matched_eft_inputform_file', '_matched_eft.inputform.txt'),
        ('weinberg_ordered_coefficient_inputform_file', '_weinberg_ordered_coefficient.inputform.txt'),
    )
    paths = []
    for key, extension in names:
        expected = stem + extension
        if report.get(key) != expected:
            fail(model + ': unexpected matching artifact filename: ' + key)
        path = outdir / expected
        if not path.is_file() or path.stat().st_size == 0:
            fail(model + ': missing/empty matching artifact: ' + str(path))
        paths.append(path)
    return paths[0], paths[1]


def run_model(model: str, *, which: str = 'final', folder: Path = INPUT,
              certificate_dir: Path = OUTPUT,
              outdir: Path | None = None,
              wolframscript: str = 'wolframscript',
              timeout_sec: int = 3600) -> dict:
    if timeout_sec <= 0:
        fail('timeout_sec must be positive')
    outdir = outdir if outdir is not None else OUTPUT / 'one_loop_probe'
    prepared = prepare_model(model, which=which, folder=folder, outdir=certificate_dir)
    original_check = _prior_certificate(model, which=which,
         certificate_dir=certificate_dir, unitary=prepared['payload']['unitary_U'])
    if not WOLFRAM.is_file() or not MATCHING_SOURCE.is_file():
        fail(f'Missing matching source {WOLFRAM} or {MATCHING_SOURCE}')
    if not shutil.which(wolframscript) and not Path(wolframscript).is_file():
        fail('wolframscript unavailable')
    outdir.mkdir(parents=True, exist_ok=True)
    stem = _stem(model, which)
    target = outdir / (stem + '_one_loop_probe.json')
    progress = outdir / (target.stem + '.progress.json')
    log = outdir / (stem + '_one_loop_probe.log')
    files = (target, progress, outdir / (stem + '_one_loop_probe_matched_eft.inputform.txt'),
             outdir / (stem + '_one_loop_probe_weinberg_ordered_coefficient.inputform.txt'))
    for path in files:
        if path.exists():
            path.unlink()
    env = dict(os.environ)
    env.pop('T3_EXPERIMENTAL_GENERAL_QUARTICS', None)
    command = [wolframscript, '-file', str(WOLFRAM), str(prepared['parameters']), str(target)]
    try:
        proc = subprocess.run(command, cwd=ROOT, env=env, capture_output=True,
                              text=True, timeout=timeout_sec, check=False)
    except subprocess.TimeoutExpired as exc:
        raise ValueError(
            f'T3 exploratory one-loop matching fail-closed: {model}: '
            f'Wolfram direct matching timed out after {timeout_sec}s; '
            f'check {progress} for the last completed stage') from exc
    log.write_text('command: ' + ' '.join(command) + '\nreturncode: ' + str(proc.returncode)
                   + '\nstdout:\n' + proc.stdout + '\nstderr:\n' + proc.stderr,
                   encoding='utf-8')
    last_stage = 'unknown'
    if progress.is_file():
        try:
            last_stage = _read_json(progress).get('last_stage', 'unknown')
        except (ValueError, OSError, json.JSONDecodeError):
            last_stage = 'invalid-progress-json'
    tail = (proc.stdout + '\n' + proc.stderr)[-2500:]
    if not target.is_file():
        fail(f'{model}: no matching JSON (exit {proc.returncode}), last stage {last_stage}; '
             f'see {log}\n{tail}')
    report = _read_json(target)
    if proc.returncode != 0:
        fail(f'{model}: matching not complete (exit {proc.returncode}, '
             f'last stage={last_stage}, status={report.get("status")}); see {log}\n{tail}')
    matched, coefficient = validate_probe_report(report, model, outdir, target.stem)
    return {
        'model': model, 'which': which, 'status': STATUS,
        'heavy_basis_certificate': str(original_check.relative_to(ROOT)),
        'probe_script_sha256': sha256(WOLFRAM.read_bytes()).hexdigest(),
        'matching_source_sha256': sha256(MATCHING_SOURCE.read_bytes()).hexdigest(),
        'matched_eft_file': str(matched.relative_to(ROOT)),
        'weinberg_ordered_coefficient_file': str(coefficient.relative_to(ROOT)),
        'numerical_c5_computed': False,
        'physical_matching_performed': False,
        'production_sequential_matching_performed': False,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    selector = parser.add_mutually_exclusive_group(required=True)
    selector.add_argument('--model', choices=MIXED_MODELS)
    selector.add_argument('--all', action='store_true')
    parser.add_argument('--which', choices=['initial', 'final'], default='final')
    parser.add_argument('--input-dir', type=Path, default=INPUT)
    parser.add_argument('--certificate-dir', type=Path, default=OUTPUT)
    parser.add_argument('--output-dir', type=Path, default=OUTPUT / 'one_loop_probe')
    parser.add_argument('--wolframscript', default='wolframscript')
    parser.add_argument('--timeout-sec', type=int, default=3600)
    args = parser.parse_args()
    if args.timeout_sec <= 0:
        parser.error('--timeout-sec must be positive')
    models = MIXED_MODELS if args.all else (args.model,)
    successes = []
    for model in models:
        item = run_model(model, which=args.which, folder=args.input_dir,
             certificate_dir=args.certificate_dir, outdir=args.output_dir,
             wolframscript=args.wolframscript, timeout_sec=args.timeout_sec)
        successes.append(item)
        print(f'{model}: exploratory one-loop symbolic Weinberg coefficient extracted; '
              'numerical C5 and sequential EFT1 matching NOT performed', flush=True)
    summary = {
        'format': FORMAT, 'status': 'ExperimentalSymbolicMatchingProbePassed',
        'models': successes, 'numerical_c5_computed': False,
        'physical_matching_performed': False,
        'production_sequential_matching_performed': False,
    }
    path = args.output_dir / ('summary_' + args.which + '.json')
    path.write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
    print('Summary:', path)


if __name__ == '__main__':
    main()
