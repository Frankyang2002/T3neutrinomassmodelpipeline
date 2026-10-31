"""Recover completed native X-only raw Matchete outputs from the v3 reporting bug.

This is a STRICT forensic recovery, not a new Wolfram calculation. It does not
rewrite any prior Matchete results or certify sequential physical matching.
The v3 source had a five-argument AssociateTo at the very end of the X stage:
all Match and operator export stages ran, but the JSON was left 'in progress'.

python -u -X utf8 -m scripts.RecoverT3BCXRawMatching --all
"""
from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path
import zipfile

ROOT = Path(__file__).resolve().parents[1]
RELATIVE_V3 = Path('output/general_t3_mass_eigen_matchete/eft1_scalar_xy_fresh_kernel_v3')
RELATIVE_RECOVERED = Path('output/audits/t3_bc_x_raw_matching_recovery_v4')
MODELS = ('B_-1', 'C_-1')
FORMAT = 't3-bc-x-raw-native-matching-post-export-recovery-v4'
STATUS = 'XRawNativeMatchingOutputsRecoveredFromV3ReportBug'
V3_WL_SHA256 = 'ccf16e4c5aac78a8fc80c5f50fac766bec65b62279c40ae99db97cbe8e698ee0'
CHECK_FIELDS = {'ordinary_check', 'hermiticity', 'canonical', 'heavy_mass_basis',
                'gauge_invariant', 'all_objects_defined'}
ARTIFACTS = ('_x_heavy_rebased_eft2.inputform.txt',
             '_y_tree_eomeft.inputform.txt', '_y_full_eomeft.inputform.txt',
             '_x_tree_sm.inputform.txt', '_x_full_sm.inputform.txt',
             '_x_loop_difference.inputform.txt', '_x_raw_llhh.inputform.txt')
PROGRESS_STAGES = ('x-followup-tree-matching-started',
                   'x-followup-tree-matching-returned',
                   'x-followup-one-loop-matching-started',
                   'x-followup-one-loop-matching-returned',
                   'x-followup-operator-export-started',
                   'x-followup-matching-exports-returned',
                   'x-fresh-X-followup-returned')


def require(ok: bool, message: str) -> None:
    if not ok:
        raise RuntimeError('X native evidence recovery: ' + message)


def file_sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def read_json(path: Path) -> dict:
    require(path.is_file(), 'missing ' + str(path))
    value = json.loads(path.read_text(encoding='utf-8'))
    require(isinstance(value, dict), 'invalid JSON: ' + str(path))
    return value


def valid_six(check: object) -> bool:
    return (isinstance(check, dict) and set(check) == CHECK_FIELDS
            and all(value is True for value in check.values()))


def count_top_level_terms(expression: str) -> int:
    """Count additive Mathematica InputForm terms, preserving nested arguments.

    Fails on unbalanced delimiters rather than inferring a count from a keyword.
    This is not a full symbolic Mathematica expression parser.
    """
    require(bool(expression.strip()), 'empty Mathematica expression')
    openings = '([{'
    closings = ')]}'
    stack: list[str] = []
    quoted = False
    escaped = False
    count = 1
    for i, char in enumerate(expression):
        if quoted:
            if escaped:
                escaped = False
            elif char == '\\':
                escaped = True
            elif char == '"':
                quoted = False
            continue
        if char == '"':
            quoted = True
        elif char in openings:
            stack.append(char)
        elif char in closings:
            require(bool(stack) and openings.index(stack.pop()) == closings.index(char),
                    'unbalanced Mathematica delimiter')
        elif char in '+-' and not stack and i > 0 and expression[i-1] not in '+-*/^':
            count += 1
    require(not stack and not quoted, 'unterminated Mathematica expression')
    return count


def base_paths(model: str, root: Path) -> tuple[Path, str]:
    require(model in MODELS, 'unsupported model: ' + str(model))
    tag = model.replace('+', 'p').replace('-', 'm')
    return root / RELATIVE_V3, f'{tag}_final_x_from_saved_y_match'


def inspect_model(model: str, root: Path = ROOT) -> tuple[dict, dict[str, Path]]:
    """Read-only, fail-closed verification of the *existing* local v3 outputs."""
    directory, stem = base_paths(model, Path(root))
    files: dict[str, Path] = {
        'x_report': directory / (stem + '_x_followup.json'),
        'y_report': directory / (stem + '.json'),
        'provenance': directory / (stem + '.provenance.json'),
        'progress': directory / (stem + '.progress.json'),
        'log': directory / (stem + '.log'),
        'crosschecks': directory / (stem + '_x_followup_rebase_crosschecks.json'),
    }
    files.update({key: directory / (stem + '_x_followup' + key)
                  for key in ARTIFACTS})
    for name, path in files.items():
        require(path.is_file(), f'{model}: missing {name}: {path}')

    x = read_json(files['x_report'])
    y = read_json(files['y_report'])
    prov = read_json(files['provenance'])
    progress = read_json(files['progress'])
    crosschecks = read_json(files['crosschecks'])
    log = files['log'].read_text(encoding='utf-8', errors='replace')
    require(x.get('model') == y.get('model') == prov.get('model') ==
            progress.get('model') == model, f'{model}: mismatched model identifier')
    require(x.get('status') == 'XFollowupOneLoopMatchingInProgress' and
            x.get('x_tree_matching_completed') is True and
            x.get('physical_matching_certified') is False,
            f'{model}: not the known post-export v3 report failure')
    require(y.get('status') == 'SavedYTreeEOMConsumedForFreshXNoYRematch' and
            y.get('x_followup_status') == x.get('status') and
            y.get('physical_sequential_matching_certified') is False,
            f'{model}: unexpected outer fresh-kernel status')
    require(prov.get('format') == 't3-bc-fresh-x-only-native-attempt-v3' and
            prov.get('mode') == 'match' and prov.get('exit_code') == 0 and
            prov.get('fresh_wl_sha256') == V3_WL_SHA256 and
            prov.get('fresh_x_matchete_matching_completed') is False and
            prov.get('physical_matching_certified') is False,
            f'{model}: source fingerprint or provenance mismatch')
    require(x.get('matching_input_selection') == 'heavy_free_rebased_y_tree_eomeft'
            and x.get('rebased_check_passed') is True
            and x.get('selected_x_input_full_matchete_check_passed') is True
            and x.get('source_is_y_tree_EOMEFT') is True
            and x.get('y_loop_inherited_in_x_matching') is False
            and x.get('running_threshold_couplings_applied') is False
            and x.get('saved_fermion_EOMEFT_equivalence_certified') is False,
            f'{model}: unverified or changed X matching input')
    require(valid_six(x.get('rebased_x_heavy_eft2_detailed_check'))
            and valid_six(x.get('light_x_eft2_detailed_check')),
            f'{model}: X check did not pass all six conditions')
    require(crosschecks.get('status') == 'XHeavyRegistryVersusRebaseComparisonCompleted'
            and crosschecks.get('rebased_heavy_check_passed') is True
            and crosschecks.get('matched_eft_fallback_used') is False,
            f'{model}: rebase crosschecks mismatch')
    history = progress.get('stage_history', [])
    indices = [history.index(stage) if stage in history else -1
               for stage in PROGRESS_STAGES]
    require(all(i >= 0 for i in indices) and indices == sorted(indices)
            and progress.get('last_stage') == PROGRESS_STAGES[-1]
            and progress.get('physical_matching_performed') is False,
            f'{model}: Matchete matching/export checkpoints incomplete')
    require(log.count('Matching pipeline completed successfully.') == 2
            and log.count('AssociateTo::argrx') == 1
            and 'exit=0' in log,
            f'{model}: not the known one-loop-completed/AssociateTo-report bug')

    texts = {}
    for key in ARTIFACTS:
        text = files[key].read_text(encoding='utf-8')
        require(len(text) > 500 and 'Field[NewFermion' not in text,
                f'{model}: malformed or contaminated {key}')
        count_top_level_terms(text)
        texts[key] = text
    for key in ('_x_tree_sm.inputform.txt', '_x_full_sm.inputform.txt',
                '_x_loop_difference.inputform.txt', '_x_raw_llhh.inputform.txt'):
        require(all('Field[' + field not in texts[key]
                    for field in ('NewScalarX', 'NewScalarY', 'NewFermion')),
                f'{model}: integrated field remains in {key}')
    require('hbar' in texts['_x_loop_difference.inputform.txt']
            and 'hbar' in texts['_x_raw_llhh.inputform.txt']
            and 'GammaCC' in texts['_x_raw_llhh.inputform.txt']
            and 'Field[H, Scalar' in texts['_x_raw_llhh.inputform.txt']
            and 'Field[l, Fermion' in texts['_x_raw_llhh.inputform.txt'],
            f'{model}: missing one-loop LLHH operator structure')
    n_llhh = count_top_level_terms(texts['_x_raw_llhh.inputform.txt'])
    require(n_llhh == texts['_x_raw_llhh.inputform.txt'].count('GammaCC'),
            f'{model}: extracted LLHH term count and Majorana operator count differ')
    # The source archive shows 36 terms in each B/C benchmark; reject changes.
    require(n_llhh == 36, f'{model}: unexpected X LLHH term count {n_llhh}')
    ytree = texts['_y_tree_eomeft.inputform.txt']
    require('Field[NewScalarX' in ytree and 'Field[NewScalarY' not in ytree,
            f'{model}: saved Y tree EFT2 inventory invalid')
    result = {
        'format': FORMAT,
        'status': STATUS,
        'model': model,
        'original_v3_status': x['status'],
        'recovered_evidence_method': 'validated_completed_log_checkpoints_and_operator_exports',
        'source_wolfram_v3_sha256': V3_WL_SHA256,
        'x_tree_native_raw_match_completed': True,
        'x_one_loop_native_raw_match_completed': True,
        'x_raw_llhh_term_count': n_llhh,
        'x_raw_llhh_contains_dimreg_pole': '\u03f5' in texts['_x_raw_llhh.inputform.txt'],
        'x_one_loop_wilson_structure_present': True,
        'matchete_five_stages_succeeded_per_x_match': True,
        'verified_rebased_x_input': True,
        'source_is_y_tree_eomeft': True,
        'one_loop_y_contributions_inherited': False,
        'original_f_first_eomeft_equivalence_proven': False,
        'running_boundary_couplings_substituted': False,
        'renormalized_finite_scalar_matching_completed': False,
        'physical_sequential_matching_certified': False,
        'physical_neutrino_mass_prediction_certified': False,
        'input_sha256': {name: file_sha(path) for name, path in sorted(files.items())},
        'original_v3_run_provenance': {
            'source_sha256': prov.get('source_sha256'),
            'saved_y_tree_eomeft_sha256': prov.get('saved_y_tree_eomeft_sha256'),
            'saved_y_full_eomeft_sha256': prov.get('saved_y_full_eomeft_sha256'),
        },
    }
    return result, files


def recover(models: tuple[str, ...], *, root: Path = ROOT,
            output_dir: Path | None = None) -> dict:
    root = Path(root)
    output_dir = output_dir or (root / RELATIVE_RECOVERED)
    results = []
    evidence = {}
    for model in models:
        report, files = inspect_model(model, root)
        results.append(report)
        evidence[model] = files
    # No output is written until every requested model has passed validation.
    output_dir.mkdir(parents=True, exist_ok=True)
    for result in results:
        tag = result['model'].replace('+', 'p').replace('-', 'm')
        dest = output_dir / f'{tag}_x_raw_recovery.json'
        dest.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    summary = {
        'format': FORMAT, 'status': STATUS,
        'recovered_model_count': len(results), 'models': [r['model'] for r in results],
        'one_loop_x_raw_matching_completed_for_all': True,
        'physical_sequential_matching_certified': False,
        'results': results,
    }
    (output_dir / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
    with (output_dir / 'summary.csv').open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=[
            'model', 'status', 'x_raw_llhh_term_count',
            'x_one_loop_native_raw_match_completed',
            'renormalized_finite_scalar_matching_completed',
            'physical_sequential_matching_certified'])
        writer.writeheader()
        for row in results:
            writer.writerow({field: row[field] for field in writer.fieldnames})
    archive = output_dir / 'recovered_x_native_raw_evidence_v4.zip'
    with zipfile.ZipFile(archive, 'w', compression=zipfile.ZIP_DEFLATED) as zf:
        for result in results:
            tag = result['model'].replace('+', 'p').replace('-', 'm')
            zf.write(output_dir / f'{tag}_x_raw_recovery.json',
                     arcname=f'recovery/{tag}_x_raw_recovery.json')
            for name, path in sorted(evidence[result['model']].items()):
                zf.write(path, arcname=f'evidence/{tag}/{path.name}')
        zf.write(output_dir / 'summary.json', arcname='recovery/summary.json')
        zf.write(output_dir / 'summary.csv', arcname='recovery/summary.csv')
    print('Recovered raw X-only Matchete evidence:', archive)
    for row in results:
        print(f"{row['model']}: X tree=PASS; X one-loop=PASS; "
              f"raw LLHH terms={row['x_raw_llhh_term_count']}; "
              'physical=NOT CERTIFIED')
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    selection = parser.add_mutually_exclusive_group(required=True)
    selection.add_argument('--all', action='store_true')
    selection.add_argument('--model', choices=MODELS)
    parser.add_argument('--source-root', type=Path, default=ROOT)
    parser.add_argument('--output-dir', type=Path)
    args = parser.parse_args()
    recover(MODELS if args.all else (args.model,), root=args.source_root,
            output_dir=args.output_dir)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
