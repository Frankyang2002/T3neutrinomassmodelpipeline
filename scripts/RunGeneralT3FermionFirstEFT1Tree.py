"""First actual Matchete fermion-first tree EFT1 matching for B_-1 and C_-1.

F is Heavy; mass-eigenstate X,Y are Light. A dimension-five LLSS tree
sector is extracted. This stage alone does NOT give physical C5 or a
complete two-threshold EFT1 prediction.

python -m scripts.RunGeneralT3FermionFirstEFT1Tree --all
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

WOLFRAM = ROOT / 'Lagrangian/RunGeneralT3FermionFirstEFT1Tree.wl'
MATCHING_SOURCE = ROOT / 'Lagrangian/RunMatching.wl'
OUTDIR = OUTPUT / 'eft1_fermion_first_tree'
FORMAT = 't3-mass-eigen-fermion-first-eft1-tree-v1'
STATUS = 'FermionFirstTreeLLSSExtractedNotSequentialEFT1'


def fail(message: str) -> None:
    raise ValueError('T3 F-first EFT1 tree fail-closed: ' + message)


def _stem(model: str, which: str) -> str:
    if model not in MIXED_MODELS or which not in {'initial', 'final'}:
        fail(f'Unsupported model/endpoint: {model}, {which}')
    return model.replace('+', 'p').replace('-', 'm') + '_' + which + '_eft1_fermion_first_tree'


def _prior_certificate(model: str, which: str, certificate_dir: Path, unitary: Any) -> Path:
    base = model.replace('+', 'p').replace('-', 'm') + '_' + which
    source = certificate_dir / (base + '_matchete.json')
    if not source.is_file():
        fail(f'{model}: missing prior live Matchete eigenfield check: {source}')
    validate_wolfram_output(_read_json(source), model, unitary)
    return source


def validate_tree_report(report: Mapping[str, Any], model: str,
                         outdir: Path, stem: str) -> tuple[Path, Path]:
    if not isinstance(report, Mapping):
        fail(f'{model}: EFT1 tree report is not a JSON object')
    required = {
        'format': FORMAT, 'status': STATUS, 'model': model,
        'matching_approach': 'FermionFirstTreeOnlyScalarsRemainLight',
        'matching_EFT_order': 5, 'matching_loop_order': 0,
        'mass_basis_ordinary_check_passed': True,
        'fermion_heavy_mass_basis_detailed_check_passed': True,
        'source_quartic_count': 17, 'source_yukawa_count': 2,
        'scalar_light_fields': ['NewScalarX', 'NewScalarY'],
        'integrated_heavy_fields': ['NewFermion'],
        'tree_matching_completed': True,
        'fermion_field_absent': True, 'x_scalar_present': True,
        'y_scalar_present': True, 'llss_contains_x': True,
        'llss_contains_y': True, 'eft1_rge_performed': False,
        'scalar_threshold_matching_performed': False,
        'one_loop_hard_fermion_matching_performed': False,
        'physical_C5_computed': False,
        'sequential_EFT1_completed': False,
        'physical_matching_performed': False,
        'matchete_rgbeta_quartic_normalization_certified': False,
    }
    for key, expected in required.items():
        if report.get(key) != expected or type(report.get(key)) is not type(expected):
            fail(f'{model}: missing/inconsistent stage certificate: {key}')
    n = report.get('llss_term_count')
    if not isinstance(n, int) or isinstance(n, bool) or n <= 0:
        fail(f'{model}: no dimension-five LLSS operators found')
    artefacts = []
    for key, suffix in (
        ('matched_eft_file', '_eft1_tree.inputform.txt'),
        ('llss_inputform_file', '_llss_tree.inputform.txt'),
    ):
        expected = stem + suffix
        if report.get(key) != expected:
            fail(f'{model}: unexpected artifact path for {key}')
        path = outdir / expected
        if not path.is_file() or path.stat().st_size == 0:
            fail(f'{model}: missing or empty artifact: {path}')
        artefacts.append(path)
    llss_text = artefacts[1].read_text(encoding='utf-8')
    if not all(x in llss_text for x in ('GammaCC', 'NewScalarX', 'NewScalarY')):
        fail(f'{model}: exported LLSS sector lacks certified field signatures')
    if 'Field[NewFermion' in llss_text:
        fail(f'{model}: EFT1 operator still contains integrated fermion field')
    return artefacts[0], artefacts[1]


def run_model(model: str, *, which: str = 'final', folder: Path = INPUT,
              certificate_dir: Path = OUTPUT, outdir: Path = OUTDIR,
              wolframscript: str = 'wolframscript', timeout_sec: int = 3600) -> dict:
    if timeout_sec <= 0:
        fail('timeout_sec must be positive')
    prepared = prepare_model(model, which=which, folder=folder, outdir=certificate_dir)
    prior = _prior_certificate(model, which, certificate_dir, prepared['payload']['unitary_U'])
    if not WOLFRAM.is_file() or not MATCHING_SOURCE.is_file():
        fail(f'Missing source {WOLFRAM} or {MATCHING_SOURCE}')
    if not shutil.which(wolframscript) and not Path(wolframscript).is_file():
        fail('wolframscript not available')
    outdir.mkdir(parents=True, exist_ok=True)
    stem = _stem(model, which)
    target = outdir / (stem + '.json')
    progress = outdir / (stem + '.progress.json')
    log = outdir / (stem + '.log')
    for item in (target, progress, outdir / (stem + '_eft1_tree.inputform.txt'),
                 outdir / (stem + '_llss_tree.inputform.txt')):
        if item.exists():
            item.unlink()
    env = dict(os.environ)
    env.pop('T3_EXPERIMENTAL_GENERAL_QUARTICS', None)
    command = [wolframscript, '-file', str(WOLFRAM), str(prepared['parameters']), str(target)]
    try:
        proc = subprocess.run(command, cwd=ROOT, env=env, capture_output=True,
                              text=True, timeout=timeout_sec, check=False)
    except subprocess.TimeoutExpired as e:
        fail(f'{model}: timed out after {timeout_sec}s; check progress {progress}')
    log.write_text('command: ' + ' '.join(command) + '\nreturncode: ' + str(proc.returncode)
                   + '\nstdout:\n' + proc.stdout + '\nstderr:\n' + proc.stderr,
                   encoding='utf-8')
    last = 'unknown'
    if progress.is_file():
        try:
            last = _read_json(progress).get('last_stage', 'unknown')
        except (ValueError, OSError, json.JSONDecodeError):
            last = 'invalid-progress-json'
    tail = (proc.stdout + '\n' + proc.stderr)[-3500:]
    if not target.is_file():
        fail(f'{model}: no EFT1 result JSON, exit={proc.returncode}, stage={last}; see {log}\n{tail}')
    report = _read_json(target)
    if proc.returncode != 0:
        fail(f'{model}: EFT1 incomplete (exit={proc.returncode}, stage={last}, '
             f'status={report.get("status")}); see {log}\n{tail}')
    eft, llss = validate_tree_report(report, model, outdir, stem)
    return {
        'model': model, 'status': STATUS,
        'heavy_basis_certificate': str(prior.relative_to(ROOT)),
        'first_threshold_wolfram_sha256': sha256(WOLFRAM.read_bytes()).hexdigest(),
        'existing_matching_source_sha256': sha256(MATCHING_SOURCE.read_bytes()).hexdigest(),
        'matched_eft_file': str(eft.relative_to(ROOT)),
        'llss_file': str(llss.relative_to(ROOT)),
        'llss_term_count': report['llss_term_count'],
        'physical_C5_computed': False, 'sequential_EFT1_completed': False,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    selector = parser.add_mutually_exclusive_group(required=True)
    selector.add_argument('--model', choices=MIXED_MODELS)
    selector.add_argument('--all', action='store_true')
    parser.add_argument('--which', choices=['initial', 'final'], default='final')
    parser.add_argument('--input-dir', type=Path, default=INPUT)
    parser.add_argument('--certificate-dir', type=Path, default=OUTPUT)
    parser.add_argument('--output-dir', type=Path, default=OUTDIR)
    parser.add_argument('--wolframscript', default='wolframscript')
    parser.add_argument('--timeout-sec', type=int, default=3600)
    args = parser.parse_args()
    if args.timeout_sec <= 0:
        parser.error('--timeout-sec must be positive')
    models = MIXED_MODELS if args.all else (args.model,)
    records = []
    for model in models:
        record = run_model(model, which=args.which, folder=args.input_dir,
                           certificate_dir=args.certificate_dir, outdir=args.output_dir,
                           wolframscript=args.wolframscript, timeout_sec=args.timeout_sec)
        records.append(record)
        print(f'{model}: Matchete F-first EFT1 TREE matching passed; '
              f'LLSS terms={record["llss_term_count"]}; '
              'one-loop hard matching, EFT1 RGEs, scalar threshold C5 NOT performed',
              flush=True)
    summary = {'format': FORMAT, 'status': 'FirstThresholdTreeLLSSObtained',
               'models': records, 'physical_C5_computed': False,
               'one_loop_hard_fermion_matching_performed': False,
               'sequential_EFT1_completed': False}
    path = args.output_dir / ('summary_' + args.which + '.json')
    path.write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
    print('Summary:', path)


if __name__ == '__main__':
    main()
