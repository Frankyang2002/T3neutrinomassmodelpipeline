"""Batch 6: independently audit full B(-1)/C(-1) EFT1 scalar residuals.

This is a READ-ONLY, NO-MATCH component-equivalence diagnostic. It consumes
newly exported residuals from RunT3MixedScalarHardBatch6 v6 and reuses the
repository's existing exact T3BCSU2ComponentEquivalence.audit_model solver and
ResolveT3BCEFT1HermiticityAndYGate native-array verifier.

Passing component equivalence does NOT turn a failed Matchete CheckLagrangian
into success, and does not authorize physical one-loop matching.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import re
import sys
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / 'output/general_sequential_pilot/mixed_scalar_hard_batch6'
MODELS = {'B_-1': 'B_m1', 'C_-1': 'C_m1'}
FORMAT = 't3-batch6-full-eft1-native-component-audit-v1'


class AuditError(ValueError):
    pass


def load_json(path: Path) -> dict[str, Any]:
    try:
        val = json.loads(path.read_text(encoding='utf-8-sig'))
    except (OSError, UnicodeError, ValueError) as exc:
        raise AuditError(f'Missing or invalid JSON: {path}') from exc
    if not isinstance(val, dict):
        raise AuditError(f'JSON object required: {path}')
    return val


class DirectoryEvidence:
    """Existing exact solver's read-only evidence interface, without a ZIP."""

    def __init__(self, folder: Path):
        self.folder = folder

    def _safe_path(self, name: str) -> Path:
        if not isinstance(name, str) or not name or Path(name).name != name or \
                not re.fullmatch(r'[A-Za-z0-9_.-]+', name):
            raise AuditError(f'Invalid residual file name: {name!r}')
        return self.folder / name

    def text(self, name: str) -> str:
        try:
            return self._safe_path(name).read_text(encoding='utf-8-sig')
        except (OSError, UnicodeError) as exc:
            raise AuditError(f'Missing residual source file {name}') from exc

    def sha(self, name: str) -> str:
        try:
            return sha256(self._safe_path(name).read_bytes()).hexdigest()
        except OSError as exc:
            raise AuditError(f'Missing residual source file {name}') from exc


def inspect_inputs(folder: Path, model: str) -> dict[str, Any]:
    """Fail closed on wrong physics scope or incomplete residual partition."""
    if model not in MODELS:
        raise AuditError('Only B_-1 and C_-1 are supported')
    stem = MODELS[model]
    ev = DirectoryEvidence(folder)
    prefix = f'{stem}_full_eft1_residual'
    man = load_json(folder / f'{prefix}_full_residual_manifest.json')
    audit = load_json(folder / 'hard_scalar_reconstruction_audit.json')
    check = load_json(folder / 'hard_scalar_eft1_validation.json')
    progress = load_json(folder / 'native_result.progress.json')
    if man.get('format') != 't3-bc-complete-eft1-reconstruction-residual-v1' or \
       man.get('status') != 'CompleteSymbolicResidualExportedEquivalenceUnproven' or \
       man.get('model') != model or man.get('matching_scope') != 'FullEFT' or \
       man.get('coupling_group_partition_conserved') is not True or \
       man.get('contains_majorana_residual') is not False or \
       man.get('physical_matching_certified') is not False or \
       man.get('CG_and_index_tensor_equivalence_proven') is not False:
        raise AuditError(model + ': malformed/unverified residual manifest')
    if check.get('model') != model or check.get('full_EFT_validated') is not False:
        raise AuditError(model + ': unexpected EFT1 validation state')
    if check.get('initial_check', {}).get('checks', {}).get('Hermiticity') != 'False':
        raise AuditError(model + ': this audit requires a recorded Hermiticity-only failure')
    for key in ('CanonicallyNormalized', 'HeavyMassBasis', 'GaugeInvariant', 'AllObjectsDefined'):
        if check['initial_check']['checks'].get(key) != 'True':
            raise AuditError(model + ': another Matchete invariant failed: ' + key)
    if audit.get('algebraic_equality_proven') is not False or \
       audit.get('independent_matchete_check_passed') is not False or \
       audit.get('component_equivalence_proven') is not False or \
       audit.get('residual_manifest_file') != f'{prefix}_full_residual_manifest.json':
        raise AuditError(model + ': residual does not match the failed reconstruction')
    stages = progress.get('stage_history', [])
    if not isinstance(stages, list) or \
       ('batch6-full-component-residual-exported' not in stages or \
       'batch6-native-cross-CG-array-exported' not in stages) or \
       any('scalar-one-loop-match' in str(x) for x in stages):
        raise AuditError(model + ': residual not exported before scalar Match')
    groups = man.get('coupling_groups')
    if not isinstance(groups, list) or len(groups) != 25 or \
       man.get('residual_additive_terms') != 137 or \
       sum(g.get('additive_terms', 0) for g in groups if isinstance(g, dict)) != 137:
        raise AuditError(model + ': residual is not the repository-verified 137-term/25-sector problem')
    files = man.get('files')
    if not isinstance(files, dict) or set(files) != \
       {'saved_eomeft', 'reconstructed_candidate', 'symbolic_residual'}:
        raise AuditError(model + ': full source expressions missing')
    hashes = {key: ev.sha(name) for key, name in files.items()}
    for group in groups:
        if not isinstance(group, dict) or not isinstance(group.get('structural_sector'), str):
            raise AuditError(model + ': invalid residual group')
        ev.sha(group['file'])
    if audit.get('residual_additive_terms') != 137 or \
       audit.get('residual_structural_sectors') != 25:
        raise AuditError(model + ': audit counts disagree with residual manifest')
    return {'model': model, 'stem': stem, 'prefix': prefix, 'files_sha256': hashes,
            'sectors': len(groups), 'terms': man['residual_additive_terms']}


def interpret_solution(solutions: dict[str, Any]) -> dict[str, Any]:
    """Mirror the checked Workstream C exact SU(2) coefficient conditions."""
    from scripts.T3BCSU2ComponentEquivalence import SYM
    required = set(SYM) | {'T3MENativeCrossCGInv1', 'T3MENativeCrossCGInv2'}
    if set(solutions) != required:
        raise AuditError('Native CG family inventory changed: ' + repr(sorted(solutions)))
    noncross = [v for k, v in solutions.items() if k not in
                ('T3MENativeCrossCGInv1', 'T3MENativeCrossCGInv2')]
    if len(noncross) != 7 or any(v['matrix_rank'] != 0 or
        v['augmented_rank'] != 0 for v in noncross):
        raise AuditError('Scalar Bose-symmetric non-cross families did not cancel')
    a = solutions['T3MENativeCrossCGInv1']
    b = solutions['T3MENativeCrossCGInv2']
    if a['status'] != 'UNIQUE' or a['required_coefficients'] != ['1', '-1']:
        raise AuditError('CrossInv1 is not the required native A-B tensor')
    if b['status'] != 'NONUNIQUE' or b['matrix_rank'] != 1 or \
       not any('tau0' in str(x) for x in b['required_coefficients']):
        raise AuditError('CrossInv2 does not admit the required A+B tensor family')
    return {
        'exact_component_polynomial_conditions_solved': True,
        '7_noncross_families_cancelled': True,
        'CrossInv1_required_tensor': 'A-B',
        'CrossInv2_required_tensor': 'gamma*(A+B), gamma nonzero',
    }


def _enable_legacy_absolute_imports() -> None:
    """Make the repository's unchanged direct-script imports work with -m.

    ResolveT3BCEFT1HermiticityAndYGate.py imports
    ``T3BCSU2ComponentEquivalence`` as a *top-level* module. Python -m
    scripts.RunT3MixedComponentAuditBatch6 only puts the repository root on
    sys.path, so the legacy top-level name is otherwise unresolved. Expose
    the local scripts directory without editing or copying either dependency.
    """
    scripts_directory = str(ROOT / 'scripts')
    if scripts_directory not in sys.path:
        sys.path.insert(0, scripts_directory)


def audit_model_run(folder: Path, model: str) -> dict[str, Any]:
    inputs = inspect_inputs(folder, model)
    _enable_legacy_absolute_imports()
    from scripts.T3BCSU2ComponentEquivalence import audit_model as solve
    from scripts.ResolveT3BCEFT1HermiticityAndYGate import check_actual_native_arrays
    sectors, solutions = solve(DirectoryEvidence(folder), inputs['stem'])
    if len(sectors) != 25 or sum(item['raw_terms'] for item in sectors) != 137:
        raise AuditError(model + ': independent parser did not recover 25 sectors/137 terms')
    condition = interpret_solution(solutions)
    try:
        native = check_actual_native_arrays(folder, inputs['stem'])
        present = True
    except (ValueError, OSError) as exc:
        native = {'status': 'UnavailableOrInvalid', 'diagnostic': str(exc)}
        present = False
    # A previous component certificate is NOT a Matchete CheckLagrangian pass.
    return {
        'model': model, 'inputs': inputs,
        'status': ('NativeCGComponentConditionsVerifiedNoMatch' if present else
                   'ConditionalComponentIdentityNativeArraysUnavailable'),
        'component_conditions': condition,
        'native_array_evidence': native,
        'arrays_from_current_native_run': True,
        'native_arrays_independently_checked': present,
        'component_equivalence_under_native_CG': present,
        'matchete_Hermiticity_check_passed': False,
        'scalar_Match_performed': False,
        'physical_C5_certified': False,
        'full_pole_cancellation_certified': False,
    }


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--hard-run', required=True, help='Existing v6 native run ID')
    p.add_argument('--output', type=Path,
                   help='JSON report path (default: within hard-run directory)')
    args = p.parse_args(argv)
    if not re.fullmatch(r'[A-Za-z0-9_-]+', args.hard_run):
        p.error('Invalid run ID')
    root = OUTPUT / args.hard_run
    out = args.output or root / 'component_equivalence_audit.json'
    rows: list[dict[str, Any]] = []
    for model, stem in MODELS.items():
        try:
            row = audit_model_run(root / stem, model)
            rows.append(row)
            print(f'{model}: {row["status"]}; residual 137 terms/25 sectors; '
                  'scalar Match NOT RUN', flush=True)
        except Exception as exc:
            rows.append({'model': model, 'status': 'Blocked',
                         'reason': str(exc), 'physical_C5_certified': False})
            print(f'{model}: BLOCKED: {exc}', flush=True)
    complete = all(row.get('status') == 'NativeCGComponentConditionsVerifiedNoMatch'
                   for row in rows)
    report = {'format': FORMAT,
              'status': 'NativeCGComponentConditionsVerifiedNoMatch' if complete
                  else 'ConditionalOrBlocked',
              'models': rows, 'source_hard_run': args.hard_run,
              'matchete_Hermiticity_check_passed': False,
              'scalar_Match_performed': False,
              'fermion_Match_recomputed': False,
              'physical_C5_certified': False,
              'full_pole_cancellation_certified': False}
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print('Evidence report: ' + str(out), flush=True)
    # Nonzero if proof absent; even a 0 means ONLY component evidence verified.
    return 0 if complete else 2


if __name__ == '__main__':
    sys.exit(main())
