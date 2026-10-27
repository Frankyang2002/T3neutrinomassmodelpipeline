"""Experimental raw one-loop Matchete matching at the fermion-first EFT1 boundary.

The mass eigen-doublets X and Y stay Light; F is integrated out.
The output is the exact Matchete one-loop EFT *difference* and a raw
LLHH candidate sector, NOT a renormalized/physical Majorana coefficient.

python -m scripts.RunGeneralT3FermionFirstEFT1OneLoop --model B_-1
python -m scripts.RunGeneralT3FermionFirstEFT1OneLoop --model C_-1
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
)
from scripts.RunGeneralT3FermionFirstEFT1Tree import (
    OUTDIR as TREE_OUTDIR, _stem as tree_stem, _prior_certificate,
    validate_tree_report,
)

WOLFRAM = ROOT / 'Lagrangian/RunGeneralT3FermionFirstEFT1OneLoop.wl'
MATCHING_SOURCE = ROOT / 'Lagrangian/RunMatching.wl'
OUTDIR = OUTPUT / 'eft1_fermion_first_oneloop'
FORMAT = 't3-mass-eigen-fermion-first-eft1-oneloop-v1'
STATUS = 'FermionFirstOneLoopRawSectorsExtractedNotRenormalized'


def fail(message: str) -> None:
    raise ValueError('T3 F-first EFT1 one-loop fail-closed: ' + message)


def _stem(model: str, which: str) -> str:
    if model not in MIXED_MODELS or which not in ('initial', 'final'):
        fail(f'Unsupported model/endpoint {model} {which}')
    return model.replace('+', 'p').replace('-', 'm') + '_' + which + '_eft1_fermion_first_oneloop'


def _certify_prior_tree(model: str, which: str, tree_dir: Path) -> Path:
    stem = tree_stem(model, which)
    source = tree_dir / (stem + '.json')
    if not source.is_file():
        fail(f'{model}: run the successful tree EFT1 probe first: {source}')
    validate_tree_report(_read_json(source), model, tree_dir, stem)
    return source


def validate_one_loop_report(report: Mapping[str, Any], model: str,
                             outdir: Path, stem: str) -> dict[str, Path]:
    if not isinstance(report, Mapping):
        fail(f'{model}: one-loop report is not a JSON object')
    required = {
        'format': FORMAT, 'status': STATUS, 'model': model,
        'matching_approach': 'FermionFirstTreeAndOneLoopScalarsRemainLight',
        'matching_EFT_order': 5, 'matching_loop_order': 1,
        'mass_basis_ordinary_check_passed': True,
        'fermion_heavy_mass_basis_detailed_check_passed': True,
        'source_quartic_count': 17, 'source_yukawa_count': 2,
        'scalar_light_fields': ['NewScalarX', 'NewScalarY'],
        'integrated_heavy_fields': ['NewFermion'],
        'tree_matching_completed': True, 'fermion_field_absent': True,
        'x_scalar_present': True, 'y_scalar_present': True,
        'llss_contains_x': True, 'llss_contains_y': True,
        'one_loop_matching_completed': True,
        'loop_difference_hbar_present': True,
        'loop_difference_fermion_absent': True,
        'raw_hard_llhh_no_fermion_field': True,
        'raw_hard_llhh_scheme_renormalized': False,
        'raw_hard_llhh_pole_rge_certified': False,
        'eft1_rge_performed': False,
        'scalar_threshold_matching_performed': False,
        'one_loop_hard_fermion_matching_performed': False,
        'physical_C5_computed': False,
        'sequential_EFT1_completed': False,
        'physical_matching_performed': False,
        'matchete_rgbeta_quartic_normalization_certified': False,
    }
    for key, expected in required.items():
        observed = report.get(key)
        if observed != expected or type(observed) is not type(expected):
            fail(f'{model}: invalid result certificate {key}: {observed!r}')
    for key in ('llss_term_count', 'loop_difference_term_count',
                'raw_hard_llhh_candidate_term_count'):
        value = report.get(key)
        if not isinstance(value, int) or isinstance(value, bool):
            fail(f'{model}: invalid {key}')
        if value < (0 if key.startswith('raw_hard') else 1):
            fail(f'{model}: invalid term inventory {key}={value}')
    suffixes = {
        'matched_eft_file': '_eft1_tree.inputform.txt',
        'llss_inputform_file': '_llss_tree.inputform.txt',
        'one_loop_matched_eft_file': '_eft1_tree_plus_loop.inputform.txt',
        'one_loop_difference_file': '_eft1_one_loop_correction.inputform.txt',
        'raw_hard_llhh_file': '_raw_hard_llhh_sector.inputform.txt',
    }
    result: dict[str, Path] = {}
    for key, suffix in suffixes.items():
        expected = stem + suffix
        if report.get(key) != expected:
            fail(f'{model}: unexpected artifact name {key}')
        path = outdir / expected
        if not path.is_file() or path.stat().st_size == 0:
            fail(f'{model}: missing or empty artifact: {path}')
        result[key] = path
    tree_text = result['llss_inputform_file'].read_text(encoding='utf-8')
    if not all(t in tree_text for t in ('GammaCC', 'Field[NewScalarX', 'Field[NewScalarY')):
        fail(f'{model}: tree LLSS structure not found')
    difference = result['one_loop_difference_file'].read_text(encoding='utf-8')
    if 'hbar' not in difference or 'Field[NewFermion' in difference or not difference.strip():
        fail(f'{model}: one-loop difference missing hbar or includes F')
    hard = result['raw_hard_llhh_file'].read_text(encoding='utf-8').strip()
    if 'Field[NewFermion' in hard or 'Field[NewScalarX' in hard or 'Field[NewScalarY' in hard:
        fail(f'{model}: raw hard sector includes integrated F or unintegrated scalars')
    if report['raw_hard_llhh_candidate_term_count'] == 0:
        if hard != '0':
            fail(f'{model}: hard sector claimed empty but text nonzero')
    elif not ('GammaCC' in hard and 'Field[H,' in hard and 'hbar' in hard):
        fail(f'{model}: hard-sector candidate lacks expected operator signatures')
    return result


def run_model(model: str, *, which: str = 'final', folder: Path = INPUT,
              certificate_dir: Path = OUTPUT, tree_dir: Path = TREE_OUTDIR,
              outdir: Path = OUTDIR, wolframscript: str = 'wolframscript',
              timeout_sec: int = 7200) -> dict[str, Any]:
    if timeout_sec <= 0:
        fail('timeout_sec must be positive')
    prepared = prepare_model(model, which=which, folder=folder, outdir=certificate_dir)
    prior_mass = _prior_certificate(model, which, certificate_dir,
                                   prepared['payload']['unitary_U'])
    prior_tree = _certify_prior_tree(model, which, tree_dir)
    if not WOLFRAM.is_file() or not MATCHING_SOURCE.is_file():
        fail(f'Missing Wolfram source {WOLFRAM} or {MATCHING_SOURCE}')
    if not shutil.which(wolframscript) and not Path(wolframscript).is_file():
        fail('wolframscript not found')
    outdir.mkdir(parents=True, exist_ok=True)
    stem = _stem(model, which)
    target = outdir / (stem + '.json')
    progress = outdir / (stem + '.progress.json')
    log = outdir / (stem + '.log')
    for suffix in ('.json', '.progress.json', '_eft1_tree.inputform.txt',
                   '_llss_tree.inputform.txt', '_eft1_tree_plus_loop.inputform.txt',
                   '_eft1_one_loop_correction.inputform.txt',
                   '_raw_hard_llhh_sector.inputform.txt'):
        path = outdir / (stem + suffix)
        if path.exists():
            path.unlink()
    command = [wolframscript, '-file', str(WOLFRAM),
               str(prepared['parameters']), str(target)]
    env = dict(os.environ)
    env.pop('T3_EXPERIMENTAL_GENERAL_QUARTICS', None)
    try:
        completed = subprocess.run(command, cwd=ROOT, env=env, text=True,
                                   capture_output=True, timeout=timeout_sec,
                                   check=False)
    except subprocess.TimeoutExpired:
        fail(f'{model}: Matchete timed out after {timeout_sec}s; see {progress}')
    log.write_text('command: ' + ' '.join(command) + '\nexit: '
                   + str(completed.returncode) + '\nstdout:\n' + completed.stdout
                   + '\nstderr:\n' + completed.stderr, encoding='utf-8')
    last = 'unknown'
    if progress.is_file():
        try:
            last = _read_json(progress).get('last_stage', 'unknown')
        except (ValueError, OSError, json.JSONDecodeError):
            last = 'invalid-progress-json'
    tail = (completed.stdout + '\n' + completed.stderr)[-4500:]
    if not target.is_file():
        fail(f'{model}: no one-loop JSON (exit={completed.returncode}, last={last}); '
             f'see {log}\n{tail}')
    report = _read_json(target)
    if completed.returncode != 0:
        fail(f'{model}: one-loop Matchete incomplete (exit={completed.returncode}, '
             f'last={last}, status={report.get("status")}); see {log}\n{tail}')
    exports = validate_one_loop_report(report, model, outdir, stem)
    return {
        'model': model, 'status': STATUS,
        'heavy_mass_basis_certificate': str(prior_mass.relative_to(ROOT)),
        'fermion_first_tree_certificate': str(prior_tree.relative_to(ROOT)),
        'wolfram_source_sha256': sha256(WOLFRAM.read_bytes()).hexdigest(),
        'production_matching_source_sha256': sha256(MATCHING_SOURCE.read_bytes()).hexdigest(),
        'llss_tree_terms': report['llss_term_count'],
        'loop_difference_terms': report['loop_difference_term_count'],
        'raw_hard_llhh_candidate_terms': report['raw_hard_llhh_candidate_term_count'],
        'exports': {k: str(v.relative_to(ROOT)) for k, v in exports.items()},
        'raw_hard_llhh_scheme_renormalized': False,
        'physical_C5_computed': False, 'sequential_EFT1_completed': False,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--model', choices=MIXED_MODELS)
    group.add_argument('--all', action='store_true')
    parser.add_argument('--which', choices=('initial', 'final'), default='final')
    parser.add_argument('--input-dir', type=Path, default=INPUT)
    parser.add_argument('--certificate-dir', type=Path, default=OUTPUT)
    parser.add_argument('--tree-dir', type=Path, default=TREE_OUTDIR)
    parser.add_argument('--output-dir', type=Path, default=OUTDIR)
    parser.add_argument('--wolframscript', default='wolframscript')
    parser.add_argument('--timeout-sec', type=int, default=7200)
    args = parser.parse_args()
    selected = MIXED_MODELS if args.all else [args.model]
    results = []
    for model in selected:
        result = run_model(model, which=args.which, folder=args.input_dir,
                           certificate_dir=args.certificate_dir,
                           tree_dir=args.tree_dir, outdir=args.output_dir,
                           wolframscript=args.wolframscript, timeout_sec=args.timeout_sec)
        results.append(result)
        print(f'{model}: F-first tree LLSS terms={result["llss_tree_terms"]}; '
              f'one-loop correction terms={result["loop_difference_terms"]}; '
              f'raw LLHH candidates={result["raw_hard_llhh_candidate_terms"]}; '
              'UNRENORMALIZED / not final C5', flush=True)
    summary = {
        'format': FORMAT, 'status': 'FermionFirstOneLoopSectorsExtracted',
        'models': results, 'physical_C5_computed': False,
        'sequential_EFT1_completed': False,
    }
    dest = args.output_dir / ('summary_' + args.which + '.json')
    dest.write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
    print('Summary:', dest)


if __name__ == '__main__':
    main()
