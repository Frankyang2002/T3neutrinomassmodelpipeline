"""Batch 19: first new general-scalar one-scale matching probe, B(+1).

Run the existing Matchete model builder once without and once with its single
extra quartic, then independently compare the extracted ordered C5 coefficients
at controlled points using Wolfram arithmetic. This is NOT sequential EFT
matching and does NOT perform an oscillation refit.

All outputs are new-only in output/studies/general_scalar_extended/matching_probes.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'output/studies/general_scalar_extended/matching_probes/B_p1'
REGISTRY = 'RGE/general/generated/general_rgbeta_quartic_registry.json'
WOLFRAM_MODEL = 'Lagrangian/RunModel.wl'
WOLFRAM_COMPARE = 'scripts/CompareT3BPlus1OneScaleC5.wl'
SOURCES = (
    REGISTRY, WOLFRAM_MODEL, WOLFRAM_COMPARE,
    'Lagrangian/model/LagrangianBuilder.wl',
    'Lagrangian/interactions/GeneralScalarOperators.wl',
    'Lagrangian/interactions/GeneralScalarVerifiedQuartics.wl',
    'Lagrangian/interactions/T3Topology.wl',
    'Lagrangian/RunMatching.wl',
)
MATCHING_SOURCES = tuple(name for name in SOURCES if name != WOLFRAM_COMPARE)
EXTRA = 'lambdaT3ExtraQ022000I1'
INTERACTION = 'T3ExtraQ_0_2_2_0_0_0_1'
FORMAT = 't3-b-plus1-one-scale-general-quartic-c5-probe-v1'


def require(ok: bool, message: str) -> None:
    if not ok:
        raise RuntimeError('Batch 19 B(+1) C5 probe: ' + message)


def read(path: Path) -> dict:
    obj = json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(obj, dict), 'expected JSON object at ' + str(path))
    return obj


def fingerprints(root: Path) -> dict[str, str]:
    result = {}
    for relative in SOURCES:
        path = root / relative
        require(path.is_file(), f'missing source {relative}')
        result[relative] = sha256(path.read_bytes()).hexdigest()
    return result


def check_registry(root: Path) -> dict:
    registry = read(root / REGISTRY)
    models = registry.get('models', {})
    require(isinstance(models, dict) and len(models) == 16, 'registry model inventory changed')
    candidate = models.get('B_+1', {})
    new = candidate.get('new_rgbeta_symbols', [])
    require(len(new) == 1 and new[0].get('name') == 'lambdaT3GeneralQ022000I1'
            and new[0].get('counts') == [0, 2, 2, 0, 0, 0],
            'B(+1) additional quartic differs from source registry')
    require(candidate.get('preexisting_rgbeta_symbols') == [] and
            candidate.get('quadratics_not_implemented') == 0,
            'B(+1) no-mixing/one-quartic basis has changed')
    return {'model': 'B_+1', 'rgbeta_symbol': new[0]['name'],
            'matchete_symbol': EXTRA, 'counts': [0, 2, 2, 0, 0, 0]}


def expected_destination(run_id: str, root: Path = ROOT) -> Path:
    require(re.fullmatch(r'[a-zA-Z0-9][a-zA-Z0-9_-]{0,70}', run_id) is not None,
            'invalid run ID')
    return root / 'output/studies/general_scalar_extended/matching_probes/B_p1' / run_id


def checked_artifacts(folder: Path, *, extended: bool) -> dict:
    summary = folder / 'comparison_summary.json'
    coefficient = folder / 'c5_coefficient.txt'
    require(summary.is_file() and coefficient.is_file(),
            f'incomplete Matchete output in {folder}')
    doc = read(summary)
    # The saved upstream B(+1) reference uses the representation label
    # "T3-d2-d2-F1". The user's live local RunModel.wl instead emitted
    # "T3-B" for the *same* representation. Neither string alone proves
    # the model: verify all dimensions, alpha and hypercharges independently.
    allowed_class_labels = ('T3-d2-d2-F1', 'T3-B')
    expected_model = {
        'Alpha': 1, 'Scalar1SU2': 2, 'Scalar2SU2': 2,
        'FermionSU2': 1, 'Scalar1Hypercharge': '1/2',
        'Scalar2Hypercharge': '3/2', 'FermionHypercharge': '1',
    }
    bad = {key: {'expected': value, 'actual': doc.get(key)}
           for key, value in expected_model.items() if doc.get(key) != value}
    if doc.get('ModelClass') not in allowed_class_labels:
        bad['ModelClass'] = {
            'allowed': list(allowed_class_labels),
            'actual': doc.get('ModelClass'),
        }
    require(not bad, f'{folder.name}: matching output is not B(+1): {bad}')
    for key in ('BuildStatus', 'MatchingStatus', 'WeinbergExtractionStatus'):
        require(doc.get(key) == 'Success', f'{folder.name}: {key} not successful')
    require(doc.get('T3IngredientsPresent') is True and
            doc.get('WeinbergOperatorPresent') is True and
            isinstance(doc.get('WeinbergHolomorphicTermCount'), int) and
            doc['WeinbergHolomorphicTermCount'] > 0,
            f'{folder.name}: incomplete Weinberg extraction')
    names = doc.get('AcceptedInteractions')
    require(isinstance(names, list) and all(isinstance(n, str) for n in names),
            'accepted interaction inventory missing')
    extras = sorted(n for n in names if n.startswith('T3ExtraQ_'))
    require(extras == ([INTERACTION] if extended else []),
            f'wrong extra quartic in {folder.name}: {extras}')
    coeff = coefficient.read_text(encoding='utf-8-sig').strip()
    require(len(coeff) > 30 and 'Coupling[lambdaT3' in coeff and
            'Coupling[y1' in coeff and 'Coupling[y2' in coeff,
            'C5 expression does not have expected topology/Yukawa symbols')
    # Symbol occurrence is an observation, not a proof of numerical contribution.
    seen_extra = EXTRA in coeff
    if not extended:
        require(not seen_extra, 'restricted C5 unexpectedly contains extra quartic')
    return {
        'reported_model_class': doc['ModelClass'],
        'accepted_extra_interactions': extras,
        'holomorphic_term_count': doc['WeinbergHolomorphicTermCount'],
        'coefficient_sha256': sha256(coefficient.read_bytes()).hexdigest(),
        'coefficient_bytes': coefficient.stat().st_size,
        'extra_quartic_symbol_occurs_in_ordered_c5': seen_extra,
    }


def checked_numeric_report(path: Path) -> dict:
    obj = read(path)
    require(obj.get('format') == 't3-b-plus1-c5-wolfram-numeric-audit-v1' and
            obj.get('status') == 'Success' and obj.get('sample_count') == 3,
            'native Wolfram numeric C5 check did not succeed')
    require(obj.get('baseline_agrees_with_general_at_zero_extra') is True,
            'general C5 does not reduce to restricted C5 at zero extra quartic')
    rows = obj.get('samples')
    require(isinstance(rows, list) and len(rows) == 3,
            'native Wolfram samples missing')
    for row in rows:
        require(row.get('zero_extra_agrees') is True and
                row.get('all_numeric_values_finite') is True,
                'native Wolfram sample failed a numerical check')
    require(isinstance(obj.get('nonzero_extra_changes_c5_at_any_sample'), bool),
            'numeric dependence flag missing')
    return obj


def summarize(probe: dict) -> str:
    m = probe['native_numeric_check']
    impact = ('Nonzero direct one-scale C5 sensitivity observed'
              if m['nonzero_extra_changes_c5_at_any_sample']
              else 'No nonzero direct one-scale C5 sensitivity in three samples')
    return '\n'.join([
        '# Batch 19 — B(+1) general-scalar one-scale C5 source probe', '',
        f"Status: `{probe['status']}`; model: `B_+1`.",
        f"Registered new Matchete quartic: `{EXTRA}` (Hdag Hdag S1 S1 + h.c.).",
        f"Baseline provenance: `{probe['baseline_source_provenance']}`.", '',
        f"- Restricted and general Matchete one-scale matching: **Success**.",
        f"- General C5 includes extra symbol: **{probe['general']['extra_quartic_symbol_occurs_in_ordered_c5']}**.",
        '- Three independent native numerical sample points reproduce the restricted limit at extra quartic zero.',
        f'- Nonzero quartic test: **{impact}**.', '',
        '**Not yet done:** new UV/EFT1 full 3-study B(+1) refits, model-specific general-parameter conversion, '
        'sequential Wilson matching, physical Matchete–RGBeta field map, or vacuum stability.',
        'Even an absent one-loop C5 response does not prove the new quartic has no effect on RG evolution.',
        '',
    ])


def run_process(command: list[str], env: dict[str, str], log: Path, timeout: int,
                *, cwd: Path) -> None:
    with log.open('w', encoding='utf-8') as handle:
        try:
            result = subprocess.run(command, cwd=cwd, env=env,
                                    stdout=handle, stderr=subprocess.STDOUT,
                                    timeout=timeout, check=False)
        except subprocess.TimeoutExpired as exc:
            raise RuntimeError('Wolfram timed out; examine ' + str(log)) from exc
    require(result.returncode == 0, f'Wolfram exit={result.returncode}; inspect {log}')


def _validate_resume(dest: Path, root: Path,
                     source_sha256: dict[str, str]) -> tuple[str, bool]:
    """Check matching-source provenance independently of the numeric audit code.

    The original v2 fingerprint covers the comparison script as well as all
    matching inputs.  A *comparison-only* correction may replay the already
    matched coefficients, but any change to RunModel/Matchete/registry sources
    rejects resume.  This is intentionally not full-manifest equality.
    """
    require(dest.is_dir(), 'cannot resume missing run ID: ' + str(dest))
    require(not (dest / 'report.json').exists() and
            not (dest / 'summary.md').exists(),
            'cannot resume an already finalised run')
    baseline = dest / 'baseline'
    require(baseline.is_dir() and (dest / 'baseline.log').is_file(),
            'resume requires the original completed baseline stage and log')
    checked_artifacts(baseline, extended=False)
    manifest = dest / 'source_fingerprints.json'
    if manifest.is_file():
        previous = read(manifest)
        require(set(previous) == set(source_sha256),
                'source SHA-256 inventory changed since interrupted run')
        require(all(previous[name] == source_sha256[name] for name in MATCHING_SOURCES),
                'matching source SHA-256 inputs changed since interrupted run')
        changed_compare = previous[WOLFRAM_COMPARE] != source_sha256[WOLFRAM_COMPARE]
        return ('MatchingSourcesSHA256Verified_ComparisonScriptUpdated' if changed_compare
                else 'SourceSHA256Verified', changed_compare)
    # The original v1 did not save a fingerprint manifest. Its timestamp-only
    # recovery is limited to the baseline-only incomplete layout.
    require(not (dest / 'general').exists() and
            not (dest / 'general.log').exists() and
            not (dest / 'numeric_audit.json').exists(),
            'legacy resume supported only immediately after baseline stage')
    baseline_time = min(
        (baseline / 'comparison_summary.json').stat().st_mtime_ns,
        (baseline / 'c5_coefficient.txt').stat().st_mtime_ns)
    require(all((root / name).stat().st_mtime_ns <= baseline_time
                for name in SOURCES),
            'source modified after original baseline; use new run ID')
    return 'LegacyBaselineReusedTimestampScreenOnlyNoSourceHashCertification', False


def _preserve_preexisting_audit(dest: Path) -> None:
    """Retain the original malformed native output and its Wolfram log."""
    for name, saved in (
        ('native_numeric_audit.json', 'native_numeric_audit.before_json_fix.json'),
        ('numeric_audit.log', 'numeric_audit.before_json_fix.log'),
    ):
        source = dest / name
        target = dest / saved
        if not source.is_file():
            continue
        # Preserve the FIRST failed audit across repeated recovery attempts;
        # never overwrite it with a later partial or regenerated file.
        if not target.exists():
            shutil.copy2(source, target)


def execute(run_id: str, *, root: Path = ROOT, wolframscript: str = 'wolframscript',
            timeout: int = 7200, runner=run_process, resume: bool = False) -> dict:
    require(timeout >= 10, 'timeout must be >=10 s')
    check_registry(root)
    source_sha256 = fingerprints(root)
    dest = expected_destination(run_id, root)
    if resume:
        baseline_provenance, comparison_script_updated = _validate_resume(
            dest, root, source_sha256)
    else:
        require(not dest.exists(), 'run ID already exists: ' + str(dest))
        dest.mkdir(parents=True, exist_ok=False)
        (dest / 'source_fingerprints.json').write_text(
            json.dumps(source_sha256, indent=2) + '\n', encoding='utf-8')
        baseline_provenance = 'SourceSHA256Verified'
        comparison_script_updated = False
    # The run directory is intentionally kept on a Matchete failure.
    for stage, flag in (('baseline', '0'), ('general', '1')):
        folder = dest / stage
        if resume and folder.is_dir():
            checked_artifacts(folder, extended=(stage == 'general'))
            continue
        require(not folder.exists(), 'incomplete stage directory: ' + str(folder))
        folder.mkdir()
        env = dict(os.environ)
        env['T3_EXPERIMENTAL_GENERAL_QUARTICS'] = flag
        command = [wolframscript, '-file', str(root / WOLFRAM_MODEL),
                   str(folder), '5', '1', 'B', '1']
        runner(command, env, dest / f'{stage}.log', timeout, cwd=root)
        checked_artifacts(folder, extended=(stage == 'general'))
    audit_path = dest / 'native_numeric_audit.json'
    audit_recomputed = False
    malformed_audit_recovered = False
    previously_valid = None
    if resume and audit_path.is_file() and not comparison_script_updated:
        try:
            previously_valid = checked_numeric_report(audit_path)
        except (json.JSONDecodeError, UnicodeDecodeError):
            # Never reinterpret an invalid native JSON as an accepted physics
            # result. Keep the exact original and rerun the corrected writer.
            malformed_audit_recovered = True
    if previously_valid is not None:
        native = previously_valid
    else:
        if audit_path.is_file():
            _preserve_preexisting_audit(dest)
        runner([wolframscript, '-file', str(root / WOLFRAM_COMPARE),
                str(dest / 'baseline' / 'c5_coefficient.txt'),
                str(dest / 'general' / 'c5_coefficient.txt'), str(audit_path)],
               dict(os.environ), dest / 'numeric_audit.log', timeout, cwd=root)
        audit_recomputed = True
        try:
            native = checked_numeric_report(audit_path)
        except (json.JSONDecodeError, UnicodeDecodeError) as exc:
            raise RuntimeError(
                'Batch 19 native JSON remains malformed after audit replay; '
                f'keep {audit_path} and inspect line/column: {exc}') from exc
    report = {
        'format': FORMAT, 'status': 'OneScaleSourceAndNumericAuditPassed',
        'model': 'B_+1', 'run_id': run_id,
        'registry': check_registry(root),
        'source_sha256': source_sha256,
        'baseline_source_provenance': baseline_provenance,
        'baseline_source_hash_certified': baseline_provenance in {
            'SourceSHA256Verified',
            'MatchingSourcesSHA256Verified_ComparisonScriptUpdated'},
        'full_original_source_manifest_unchanged': not comparison_script_updated and
            baseline_provenance == 'SourceSHA256Verified',
        'comparison_script_updated_since_match': comparison_script_updated,
        'comparison_script_sha256_at_matching': (
            read(dest / 'source_fingerprints.json').get(WOLFRAM_COMPARE)
            if (dest / 'source_fingerprints.json').is_file() else None),
        'numeric_audit_recomputed': audit_recomputed,
        'malformed_prior_numeric_audit_detected': malformed_audit_recovered,
        'resumed_interrupted_run': bool(resume),
        'baseline': checked_artifacts(dest / 'baseline', extended=False),
        'general': checked_artifacts(dest / 'general', extended=True),
        'native_numeric_check': native,
        'conditional_one_scale_only': True,
        'three_study_neutrino_refits_completed': False,
        'sequential_EFT_matching_completed': False,
        'physical_matchete_rgbeta_map_certified': False,
        'historical_outputs_overwritten': False,
    }
    (dest / 'report.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n',encoding='utf-8')
    (dest / 'summary.md').write_text(summarize(report),encoding='utf-8')
    return report

def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run-id', required=True)
    p.add_argument('--wolframscript', default='wolframscript')
    p.add_argument('--timeout-sec', type=int, default=7200)
    p.add_argument('--resume', action='store_true',
                   help='resume an incomplete run without redoing a valid Matchete stage')
    a = p.parse_args()
    r = execute(a.run_id, wolframscript=a.wolframscript,
                timeout=a.timeout_sec, resume=a.resume)
    print('B(+1):',r['status'])
    print('Nonzero quartic changes sampled C5:',
          r['native_numeric_check']['nonzero_extra_changes_c5_at_any_sample'])
    print('Saved:', expected_destination(a.run_id) / 'summary.md')


if __name__ == '__main__':
    main()
