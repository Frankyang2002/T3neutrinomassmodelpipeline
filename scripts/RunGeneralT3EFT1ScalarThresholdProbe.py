"""Experimental scalar-threshold matching on saved *tree-level* mass-eigenstate EFT1.

Fresh Matchete kernel: F absent, X/Y Heavy. Match the tree EFT1 once at
LoopOrder=0 and once at LoopOrder=1, then isolate the raw one-loop LLHH
sector. This is NOT a renormalized sequential EFT1 C5: no EFT1 running,
MSbar counterterms or hard-plus-soft matching have been combined.

    python -m scripts.RunGeneralT3EFT1ScalarThresholdProbe --model B_-1
    python -m scripts.RunGeneralT3EFT1ScalarThresholdProbe --model C_-1
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
from scripts.RunGeneralT3MassEigenMatchete import ROOT, INPUT, OUTPUT, _read_json, prepare_model
from scripts.RunGeneralT3FermionFirstEFT1Tree import (
    OUTDIR as TREE_OUTDIR, _prior_certificate, validate_tree_report, _stem as tree_stem,
)
from scripts.RunGeneralT3FermionFirstEFT1OneLoop import (
    OUTDIR as ONELOOP_OUTDIR, _stem as one_loop_stem, validate_one_loop_report,
)

WOLFRAM = ROOT / 'Lagrangian/RunGeneralT3EFT1ScalarThresholdProbe.wl'
MATCHING_SOURCE = ROOT / 'Lagrangian/RunMatching.wl'
OUTDIR = OUTPUT / 'eft1_scalar_threshold_probe'
FORMAT = 't3-mass-eigen-eft1-scalar-threshold-probe-v1'
STATUS = 'ScalarThresholdOneLoopRawSectorsExtractedNoEFT1Running'


def fail(message: str) -> None:
    raise ValueError('T3 EFT1 scalar-threshold probe fail-closed: ' + message)


def _stem(model: str, which: str, *, c5_relevant_only: bool = False) -> str:
    if model not in MIXED_MODELS or which not in ('final', 'initial'):
        fail(f'Unsupported model or endpoint: {model} {which}')
    return (model.replace('+', 'p').replace('-', 'm') + '_' + which +
            '_eft1_scalar_threshold_probe' + ('_c5_relevant' if c5_relevant_only else ''))


def _verified_tree_input(model: str, which: str, tree_dir: Path,
                         one_loop_dir: Path) -> tuple[Path, Path, Path]:
    tree_certificate = tree_dir / (tree_stem(model, which) + '.json')
    if not tree_certificate.is_file():
        fail(f'{model}: missing live tree EFT1 certificate: {tree_certificate}')
    tree_report = _read_json(tree_certificate)
    validate_tree_report(tree_report, model,
                         tree_dir, tree_stem(model, which))
    first_loop_certificate = one_loop_dir / (one_loop_stem(model, which) + '.json')
    if not first_loop_certificate.is_file():
        fail(f'{model}: missing one-loop F threshold certificate: {first_loop_certificate}')
    exports = validate_one_loop_report(
        _read_json(first_loop_certificate), model,
        one_loop_dir, one_loop_stem(model, which),
    )
    # EOMEFT is the production canonical continuation. The earlier raw
    # MatchedEFT is ONLY for diagnostics, not a second Match input.
    transition_name = tree_stem(model, which) + '_eft1_transition_eom.inputform.txt'
    if tree_report.get('transition_representation') != 'EOMEFT' or \
            tree_report.get('transition_eom_eft1_file') != transition_name:
        fail(f'{model}: first-stage certificate lacks canonical EOMEFT handoff; '
             'rerun the tree-only fermion-first probe')
    tree_path = tree_dir / transition_name
    if not tree_path.is_file():
        fail(f'{model}: missing canonical transition: {tree_path}')
    source = tree_path.read_text(encoding='utf-8')
    if any(x not in source for x in ('Field[NewScalarX', 'Field[NewScalarY', 'GammaCC')):
        fail('F-first tree EFT1 lacks its two scalar fields and Majorana operators')
    if 'Field[NewFermion' in source:
        fail('F-first tree EFT1 still contains the integrated fermion')
    if len(source) < 1000:
        fail('Tree-level EFT1 expression unexpectedly short')
    return tree_path, tree_certificate, first_loop_certificate


def validate_report(report: Mapping[str, Any], model: str,
                    outdir: Path, stem: str, tree_path: Path, *,
                    c5_relevant_only: bool = False) -> dict[str, Path]:
    if not isinstance(report, Mapping):
        fail('Non-object Matchete scalar-threshold report')
    expected: dict[str, Any] = {
        'format': FORMAT, 'status': STATUS, 'model': model,
        'input_eft1_tree': tree_path.name,
        'eigenfield_input_yukawas_and_quartics_reconstructed': True,
        'fermion_removed_before_scalar_matching': True,
        'heavy_scalar_fields': ['NewScalarX', 'NewScalarY'],
        'scalar_heavy_mass_basis_detailed_check_passed': True,
        'hermiticity_check_passed': True,
        'tree_matching_completed': True, 'one_loop_matching_completed': True,
        'scalar_one_loop_correction_hbar_present': True,
        'no_scalar_fields_in_matched_eft': True,
        'eft1_running_performed': False,
        'hard_pole_counterterms_combined': False,
        'hard_pole_rge_certified': False,
        'physical_sequential_C5_computed': False,
        'physical_matching_performed': False,
    }
    if c5_relevant_only:
        expected.update({
            'matching_scope': 'LeadingC5OneLoopSectorOnly',
            'full_eft_validated': False,
            'full_saved_eft_unchanged_on_disk': True,
            'hermiticity_applies_to_matching_subtheory_only': True,
            'excluded_four_heavy_scalar_quartic_terms': 122,
        })
        if report.get('full_saved_eft_hermiticity') is not False:
            fail(f'{model}: expected explicit full-EFT Hermiticity failure provenance')
    else:
        if report.get('matching_scope', 'FullEFT') != 'FullEFT' or \
                report.get('excluded_four_heavy_scalar_quartic_terms', 0) != 0 or \
                report.get('full_eft_validated', True) is not True:
            fail(f'{model}: cannot accept restricted C5 sector as full matching')
    for k, v in expected.items():
        if report.get(k) != v or type(report.get(k)) is not type(v):
            fail(f'{model}: mismatched result certificate {k}={report.get(k)!r}')
    for k in ('raw_scalar_llhh_term_count', 'scalar_loop_difference_term_count',
              'input_eft1_tree_character_count'):
        val = report.get(k)
        if not isinstance(val, int) or isinstance(val, bool) or val < (0 if k.startswith('raw_') else 1):
            fail(f'{model}: invalid integer {k}={val!r}')
    if report['input_eft1_tree_character_count'] != len(tree_path.read_text(encoding='utf-8')):
        fail('Saved tree EFT1 source length changed')
    method = report.get('canonical_EOMEFT_hermiticity_rewrite')
    if method not in ('None', 'OriginalRelabel', 'ExpandOnly', 'ExpandRelabel', 'HcSimplify', 'GreensSimplify', 'VerifiedRenormalizableReconstruction'):
        fail(f'{model}: unsupported equivalence-rewrite certificate: {method!r}')
    diagnostic = report.get('hermiticity_diagnostic_file')
    expected_diagnostic = stem + '_hermiticity_diagnostics.json'
    if diagnostic != expected_diagnostic:
        fail(f'{model}: missing fixed-name Hermiticity diagnostic certificate')
    diag_path = outdir / diagnostic
    if not diag_path.is_file():
        fail(f'{model}: Hermiticity diagnostic not exported')
    diag = _read_json(diag_path)
    if diag.get('equivalence_rewrite_accepted') != method or \
       diag.get('final_hermiticity') is not True or \
       diag.get('final_ordinary_CheckLagrangian') != 'True' or \
       diag.get('physical_sequential_C5_computed') is not False:
        fail(f'{model}: Hermiticity diagnostic is not a successful strict certificate')
    if c5_relevant_only and (diag.get('matching_scope') !=
          'LeadingC5OneLoopSectorOnly' or
          diag.get('full_saved_eft_hermiticity') is not False or
          diag.get('excluded_four_heavy_scalar_quartic_terms') != 122):
        fail(f'{model}: restricted C5 Hermiticity provenance not certified')
    if method == 'VerifiedRenormalizableReconstruction':
        name = stem + '_reconstruction_audit.json'
        if report.get('reconstruction_audit_file') != name:
            fail(f'{model}: reconstruction audit certificate missing')
        proof_path = outdir / name
        if not proof_path.is_file():
            fail(f'{model}: verified reconstruction audit missing')
        proof = _read_json(proof_path)
        if (proof.get('status') != 'VerifiedEquivalentReconstructionAccepted'
                or proof.get('algebraic_equality_proven') is not True
                or proof.get('independent_matchete_check_passed') is not True
                or proof.get('accepted') is not True
                or proof.get('physical_sequential_C5_computed') is not False):
            fail(f'{model}: unverified reconstruction cannot be accepted')
    else:
        if report.get('reconstruction_audit_file', '') != '':
            fail(f'{model}: unexpected reconstruction certificate')
    artifacts = {
        'sm_tree_file': '_sm_tree.inputform.txt',
        'sm_one_loop_file': '_sm_tree_plus_loop.inputform.txt',
        'scalar_correction_file': '_scalar_one_loop_correction.inputform.txt',
        'raw_scalar_llhh_file': '_raw_scalar_llhh.inputform.txt',
    }
    files: dict[str, Path] = {}
    for key, suffix in artifacts.items():
        expected_name = stem + suffix
        if report.get(key) != expected_name:
            fail(f'{model}: unexpected file name for {key}')
        path = outdir / expected_name
        if not path.is_file() or path.stat().st_size == 0:
            fail(f'{model}: empty or missing exported {key}: {path}')
        files[key] = path
    full = files['sm_one_loop_file'].read_text(encoding='utf-8')
    delta = files['scalar_correction_file'].read_text(encoding='utf-8')
    llhh = files['raw_scalar_llhh_file'].read_text(encoding='utf-8').strip()
    if any('Field[' + f in full for f in ('NewFermion', 'NewScalarX', 'NewScalarY')):
        fail(f'{model}: matching result retains a BSM field')
    if 'hbar' not in delta:
        fail(f'{model}: scalar correction lost its loop power counting')
    if report['raw_scalar_llhh_term_count'] == 0:
        if llhh != '0':
            fail(f'{model}: nonzero LLHH expression with zero reported count')
    elif any(x not in llhh for x in ('GammaCC', 'Field[H, Scalar', 'Field[l, Fermion', 'hbar')):
        fail(f'{model}: extracted LLHH is missing field/loop structure')
    elif any('Field[' + f in llhh for f in ('NewFermion', 'NewScalarX', 'NewScalarY')):
        fail(f'{model}: LLHH contains an integrated field')
    return files


def run_model(model: str, *, which: str = 'final', folder: Path = INPUT,
              certificate_dir: Path = OUTPUT, tree_dir: Path = TREE_OUTDIR,
              one_loop_dir: Path = ONELOOP_OUTDIR, outdir: Path = OUTDIR,
              wolframscript: str = 'wolframscript', timeout_sec: int = 7200,
              c5_relevant_only: bool = False) -> dict[str, Any]:
    if timeout_sec <= 0:
        fail('timeout_sec must be positive')
    prepared = prepare_model(model, which=which, folder=folder, outdir=certificate_dir)
    mass_basis_certificate = _prior_certificate(model, which, certificate_dir,
                                                prepared['payload']['unitary_U'])
    tree, tree_cert, loop_cert = _verified_tree_input(model, which, tree_dir, one_loop_dir)
    if not WOLFRAM.is_file() or not MATCHING_SOURCE.is_file():
        fail('Wolfram runner or production RunMatching.wl not found')
    if not shutil.which(wolframscript) and not Path(wolframscript).is_file():
        fail('wolframscript not found')
    outdir.mkdir(parents=True, exist_ok=True)
    stem = _stem(model, which, c5_relevant_only=c5_relevant_only)
    target = outdir / (stem + '.json')
    progress = outdir / (stem + '.progress.json')
    log = outdir / (stem + '.log')
    for suffix in ('.json', '.progress.json', '_hermiticity_diagnostics.json', '_hermiticity_sector_audit.json', '_reconstruction_audit.json', '_free_lagrangian_rebase.json', '_matching_stage_audit.json', '_sm_tree.inputform.txt',
                   '_sm_tree_plus_loop.inputform.txt',
                   '_scalar_one_loop_correction.inputform.txt',
                   '_raw_scalar_llhh.inputform.txt'):
        p = outdir / (stem + suffix)
        if p.is_file():
            p.unlink()
    cmd = [wolframscript, '-file', str(WOLFRAM),
           str(prepared['parameters']), str(tree), str(target)]
    env = dict(os.environ)
    env.pop('T3_EXPERIMENTAL_GENERAL_QUARTICS', None)
    env.pop('T3_EFT1_LEADING_C5_ONLY', None)
    if c5_relevant_only:
        env['T3_EFT1_LEADING_C5_ONLY'] = '1'
    try:
        completed = subprocess.run(cmd, cwd=ROOT, env=env, text=True,
                                   capture_output=True, timeout=timeout_sec, check=False)
    except subprocess.TimeoutExpired:
        fail(f'{model}: Matchete timed out; see {progress}')
    log.write_text('command: ' + ' '.join(cmd) + '\nexit: ' +
                   str(completed.returncode) + '\nstdout:\n' + completed.stdout +
                   '\nstderr:\n' + completed.stderr, encoding='utf-8')
    last = 'unknown'
    if progress.is_file():
        try:
            last = _read_json(progress).get('last_stage', 'unknown')
        except (ValueError, OSError):
            last = 'invalid-progress-file'
    tail = (completed.stdout + '\n' + completed.stderr)[-1800:]
    audit_path = outdir / (stem + '_matching_stage_audit.json')
    match_context = ''
    if audit_path.is_file():
        try:
            audit = _read_json(audit_path)
            stages = {}
            for stage in ('tree_stage', 'one_loop_stage'):
                obj = audit.get(stage, {})
                stages[stage] = {
                    k: {name: entry.get(name) for name in
                        ('head', 'characters', 'explicit_hbar', 'contains_X',
                         'contains_Y', 'contains_GammaCC')}
                    for k, entry in obj.items() if isinstance(entry, dict)
                }
            match_context = ('\nMatchete stage summaries: ' +
                             json.dumps(stages, separators=(',', ':')) +
                             '\nFull diagnostic: ' + str(audit_path))
        except (ValueError, OSError, TypeError):
            match_context = '\nMatchete stage audit was not readable: ' + str(audit_path)
    if not target.is_file():
        fail(f'{model}: no scalar-threshold report (exit={completed.returncode}, '
             f'last={last}), see {log}\n{tail}{match_context}')
    report = _read_json(target)
    if completed.returncode != 0:
        fail(f'{model}: failed scalar-threshold matching (exit={completed.returncode}, '
             f'last={last}, status={report.get("status")}), see {log}\n{tail}{match_context}')
    exports = validate_report(report, model, outdir, stem, tree,
                              c5_relevant_only=c5_relevant_only)
    return {
        'model': model, 'status': STATUS,
        'matching_scope': 'LeadingC5OneLoopSectorOnly' if c5_relevant_only else 'FullEFT',
        'full_eft_validated': not c5_relevant_only,
        'excluded_four_heavy_scalar_quartic_terms': 122 if c5_relevant_only else 0,
        'mass_basis_certificate': str(mass_basis_certificate.relative_to(ROOT)),
        'tree_eft1_certificate': str(tree_cert.relative_to(ROOT)),
        'fermion_one_loop_certificate': str(loop_cert.relative_to(ROOT)),
        'tree_eft1_sha256': sha256(tree.read_bytes()).hexdigest(),
        'wolfram_runner_sha256': sha256(WOLFRAM.read_bytes()).hexdigest(),
        'raw_scalar_llhh_term_count': report['raw_scalar_llhh_term_count'],
        'scalar_loop_difference_term_count': report['scalar_loop_difference_term_count'],
        'exports': {k: str(p.relative_to(ROOT)) for k, p in exports.items()},
        'eft1_running_performed': False,
        'physical_sequential_C5_computed': False,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--model', choices=MIXED_MODELS)
    group.add_argument('--all', action='store_true')
    parser.add_argument('--which', choices=('final', 'initial'), default='final')
    parser.add_argument('--output-dir', type=Path, default=OUTDIR)
    parser.add_argument('--wolframscript', default='wolframscript')
    parser.add_argument('--timeout-sec', type=int, default=7200)
    parser.add_argument('--c5-relevant-only', action='store_true',
        help='Explicit one-loop LLHH projection: retains all LLSS and Higgs-bearing quartics; excludes four-scalar spectator quartics. NOT full EFT matching.')
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    results = []
    for model in (MIXED_MODELS if args.all else (args.model,)):
        row = run_model(model, which=args.which, outdir=args.output_dir,
                        wolframscript=args.wolframscript, timeout_sec=args.timeout_sec,
                        c5_relevant_only=args.c5_relevant_only)
        results.append(row)
        print(f'{model}: scalar-threshold tree+one-loop matched; '
              f'raw scalar LLHH candidates={row["raw_scalar_llhh_term_count"]}; '
              f'scope={row["matching_scope"]}; '
              'EFT1 running and physical sequential C5 NOT performed', flush=True)
    summary = {
        'format': 't3-mass-eigen-eft1-scalar-threshold-probe-batch-v1',
        'status': 'ScalarThresholdRawMatchingCompletedNotFullSequentialC5',
        'matching_scope': 'LeadingC5OneLoopSectorOnly' if args.c5_relevant_only else 'FullEFT',
        'models': results, 'physical_sequential_C5_computed': False,
        'hard_pole_rge_certified': False,
    }
    p = args.output_dir / ('summary_' + args.which +
                           ('_c5_relevant' if args.c5_relevant_only else '') + '.json')
    p.write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
    print('Summary:', p)


if __name__ == '__main__':
    main()
