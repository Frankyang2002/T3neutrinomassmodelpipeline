#!/usr/bin/env python3
"""Read-only final-eigenbasis MF audit for T3 B_-1 and C_-1.

Run from the repository root, with the existing prepared mass-eigenstate
parameter JSON for each model. The original Wolfram builder is NEVER edited.

This runner makes a disposable copy of RunGeneralT3MassEigenMatchete.wl in an
OS temporary directory, inserts a strictly observational audit at the point
where that script has constructed its final `lag`, and runs the temporary copy.
The scriptDir assignment in the disposable copy points back at the repository
so that all model definitions are loaded from the local working tree.

The builder still performs both native Matchete checks. This wrapper accepts
no result unless the builder and the MF-index audit both succeed. Neither the
physically normalized component map nor EFT matching is tested here.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from typing import Any


SCRIPT_DIR_ANCHOR = 'scriptDir = DirectoryName @ ExpandFileName[$InputFileName];'
AUDIT_POINT_ANCHOR = ' T3MECheckpoint["ordinary-CheckLagrangian-started"];'
EXPECTED_BUILDER_OUTPUT = 'MatcheteMassEigenHeavyBasisPassedPhysicalMatchingNotPerformed'


class AuditError(RuntimeError):
    """The source, arguments or audit could not be validated."""


def wl_string(value: Path) -> str:
    """Encode an absolute Windows/Unix file path as a Wolfram string literal."""
    # Forward slashes are supported by Mathematica on Windows; JSON string
    # escaping is also compatible for these simple, absolute path literals.
    return json.dumps(value.resolve().as_posix(), ensure_ascii=True)


def build_probe(source: str, script_dir: Path, audit_json: Path) -> str:
    """Prepare an instrumented *copy*, failing closed on source changes."""
    for label, anchor in (
        ('script directory', SCRIPT_DIR_ANCHOR),
        ('final mass-eigenbasis validation', AUDIT_POINT_ANCHOR),
    ):
        if source.count(anchor) != 1:
            raise AuditError(f'Expected exactly one {label} anchor, got {source.count(anchor)}')
    if source.count('lag=(built["LSM"] + free + rotated) // Expand // RelabelIndices;') != 1:
        raise AuditError('Expected final lag assignment is missing or has changed')
    if source.index('lag=(built["LSM"] + free + rotated) // Expand // RelabelIndices;') > source.index(AUDIT_POINT_ANCHOR):
        raise AuditError('Final lag assignment unexpectedly follows audit point')

    # This code observes only the already-assembled `lag`. It neither replaces
    # any Lagrangian term nor changes Matchete registries, fields or couplings.
    # The MF detector matches Matchete's *internal* Coupling[MF,...] objects,
    # following the existing RunGeneralFermionMassBasisAudit.wl convention.
    audit = r'''
 (* T3MFAUDIT: diagnostic ONLY; final `lag` is unchanged. *)
 Module[{mfIndexLists, mfTerms, mfTermContainsFermion, indexedQ,
         unindexedQ, noncanonicalQ, rawMFPresentQ, expandedLag, terms,
         mfFields, fermionIndicesQ, mfLabelAtoms, report, mfIndexStrings},
   mfIndexLists = Cases[
     lag,
     HoldPattern[Coupling[name_, indices_List, ___]] /; name === MF :> indices,
     {0, Infinity}
   ];
   indexedQ = mfIndexLists =!= {} &&
     AllTrue[mfIndexLists, MatchQ[#, {Index[_, NFlavor]}] &];
   unindexedQ = AnyTrue[mfIndexLists, # === {} &];
   noncanonicalQ = AnyTrue[mfIndexLists,
     !MatchQ[#, {Index[_, NFlavor]}] &];
   rawMFPresentQ = !FreeQ[lag, Verbatim[MF], {0, Infinity}];
   mfLabelAtoms = Count[lag, Verbatim[MF], {0, Infinity}];
   expandedLag = Expand[lag];
   terms = If[Head[expandedLag] === Plus, List @@ expandedLag, {expandedLag}];
   mfTerms = Select[terms,
     !FreeQ[#, Verbatim[MF], {0, Infinity}] &];
   mfTermContainsFermion = mfTerms =!= {} &&
     AllTrue[mfTerms,
       !FreeQ[#, Verbatim[NewFermion], {0, Infinity}] &];
   mfFields = Cases[mfTerms,
     HoldPattern[Field[NewFermion, Fermion, indices_List, ___]] :> indices,
     {0, Infinity}];
   fermionIndicesQ = mfFields =!= {} &&
     AllTrue[mfFields, !FreeQ[#, Index[_, NFlavor], {0, Infinity}] &];
   mfIndexStrings = ToString[#, InputForm] & /@
     DeleteDuplicates[mfIndexLists];
   report = <|
     "format" -> "t3-bc-final-eigenbasis-mf-index-audit-v1",
     "model" -> modelName,
     "status" -> If[
       indexedQ && !unindexedQ && !noncanonicalQ &&
       rawMFPresentQ && mfTermContainsFermion && fermionIndicesQ &&
       mfLabelAtoms === Length[mfIndexLists],
       "PassedFinalEigenbasisIndexedMF",
       "FinalEigenbasisMFIndexAssertionFailed"
     ],
     "final_eigenbasis_mass_label_present" -> rawMFPresentQ,
     "final_eigenbasis_indexed_mf" -> indexedQ,
     "final_eigenbasis_unindexed_mf" -> unindexedQ,
     "final_eigenbasis_noncanonical_mf_indices" -> noncanonicalQ,
     "mf_coupling_node_count" -> Length[mfIndexLists],
     "mf_raw_label_atom_count" -> mfLabelAtoms,
     "mf_unique_index_lists_inputform" -> mfIndexStrings,
     "mf_lagrangian_term_count" -> Length[mfTerms],
     "all_mf_terms_contain_new_fermion" -> mfTermContainsFermion,
     "new_fermion_field_instances_in_mf_terms" -> Length[mfFields],
     "new_fermion_flavour_indices_present" -> fermionIndicesQ,
     "three_flavour_mass_option_expected" -> True,
     "generic_off_diagonal_mass_tested" -> False,
     "physical_matching_performed" -> False,
     "component_basis_map_certified" -> False,
     "matchete_rgbeta_quartic_normalization_certified" -> False
   |>;
   If[Quiet @ Check[Export[T3MFAUDIT_OUTPUT_PATH, report, "RawJSON"], $Failed] === $Failed,
      Print["T3MFAUDIT: failed to export result"]];
   Print["T3MFAUDIT: ", report["status"],
     "; MF indices=", InputForm[mfIndexStrings]];
 ];
'''
    audit = audit.replace('T3MFAUDIT_OUTPUT_PATH', wl_string(audit_json))
    patched = source.replace(
        SCRIPT_DIR_ANCHOR,
        f'scriptDir = {wl_string(script_dir)};'
    )
    patched = patched.replace(AUDIT_POINT_ANCHOR, audit + '\n' + AUDIT_POINT_ANCHOR)
    return patched


def write_json_exclusive(path: Path, data: dict[str, Any], overwrite: bool) -> None:
    """Avoid accidentally replacing a previous experiment's output."""
    path.parent.mkdir(parents=True, exist_ok=True)
    mode = 'w' if overwrite else 'x'
    with path.open(mode, encoding='utf-8') as f:
        json.dump(data, f, indent=2, sort_keys=True, allow_nan=False)
        f.write('\n')


def run(args: argparse.Namespace) -> int:
    script_dir = Path(__file__).resolve().parent
    builder = script_dir / 'RunGeneralT3MassEigenMatchete.wl'
    parameters = Path(args.input).resolve(strict=True)
    output = Path(args.output).resolve()
    log_path = output.with_suffix('.log')

    if not builder.is_file():
        raise AuditError(f'Local source not found: {builder}')
    if output == parameters or log_path == parameters or output == log_path:
        raise AuditError('Output and log paths must differ from the input')
    if (output.exists() or log_path.exists()) and not args.overwrite:
        raise AuditError('Output or log already exists; choose a new path, or use --overwrite')
    if not output.parent.is_dir() and args.create_output_dir:
        output.parent.mkdir(parents=True, exist_ok=True)
    elif not output.parent.is_dir():
        raise AuditError('Output directory does not exist; use --create-output-dir')

    with parameters.open(encoding='utf-8') as f:
        param_json = json.load(f)
    model = param_json.get('model')
    if model not in {'B_-1', 'C_-1'}:
        raise AuditError('This audit is only for B_-1/C_-1 prepared inputs')
    if param_json.get('format') != 't3-experimental-mass-eigen-matchete-parameters-v1':
        raise AuditError('Wrong prepared parameter JSON format')
    if os.environ.get('T3_EXPERIMENTAL_FLAVOR_MASSES') != '1':
        raise AuditError('T3_EXPERIMENTAL_FLAVOR_MASSES must be set to 1')
    if os.environ.get('T3_EXPERIMENTAL_GENERAL_FERMIONS') != '1':
        raise AuditError('T3_EXPERIMENTAL_GENERAL_FERMIONS must be set to 1')
    if os.environ.get('T3_EXPERIMENTAL_GENERAL_QUARTICS') == '1':
        raise AuditError('Disable T3_EXPERIMENTAL_GENERAL_QUARTICS for the eigenbasis builder')
    if os.environ.get('T3_EXPERIMENTAL_GENERAL_QUADRATICS') == '1':
        raise AuditError('Disable T3_EXPERIMENTAL_GENERAL_QUADRATICS for the eigenbasis builder')

    source = builder.read_text(encoding='utf-8-sig')
    builder_hash = hashlib.sha256(source.encode('utf-8')).hexdigest()
    wolframscript = args.wolframscript
    if shutil.which(wolframscript) is None and not Path(wolframscript).is_file():
        raise AuditError(f'wolframscript executable not found: {wolframscript}')

    with tempfile.TemporaryDirectory(prefix='t3_bc_mf_audit_') as td:
        tmp = Path(td)
        probe = tmp / 'RunGeneralT3MassEigenMatchete_MFProbe.wl'
        mf_result_path = tmp / 'mf_observation.json'
        builder_result_path = tmp / 'builder_result.json'
        probe.write_text(build_probe(source, script_dir, mf_result_path), encoding='utf-8')
        cmd = [wolframscript, '-file', str(probe), str(parameters), str(builder_result_path)]
        with (tmp / 'matchete.log').open('w', encoding='utf-8') as log_file:
            completed = subprocess.run(cmd, stdout=log_file, stderr=subprocess.STDOUT, check=False)

        def load_or_none(path: Path) -> dict[str, Any] | None:
            if not path.is_file():
                return None
            try:
                value = json.loads(path.read_text(encoding='utf-8-sig'))
            except (UnicodeError, json.JSONDecodeError):
                return None
            return value if isinstance(value, dict) else None

        observation = load_or_none(mf_result_path)
        original_result = load_or_none(builder_result_path)
        original_ok = bool(original_result) and (
            original_result.get('status') == EXPECTED_BUILDER_OUTPUT
            and original_result.get('mass_basis_ordinary_check_passed') is True
            and original_result.get('heavy_mass_basis_detailed_check_passed') is True
        )
        mass_ok = bool(observation) and observation.get('status') == 'PassedFinalEigenbasisIndexedMF'
        valid = completed.returncode == 0 and original_ok and mass_ok
        report = {
            'format': 't3-bc-final-eigenbasis-mf-index-audit-wrapper-v1',
            'model': model,
            'status': 'PassedFinalEigenbasisIndexedMFAndMatchete' if valid else 'MFOrMatcheteAuditFailed',
            'builder_source': str(builder),
            'builder_sha256': builder_hash,
            'parameter_input': str(parameters),
            'matchete_exit_code': completed.returncode,
            'builder_result_status': original_result.get('status') if original_result else None,
            'builder_ordinary_check_passed': (original_result or {}).get('mass_basis_ordinary_check_passed') is True,
            'builder_heavy_basis_check_passed': (original_result or {}).get('heavy_mass_basis_detailed_check_passed') is True,
            'mf_observation': observation,
            'physical_matching_performed': False,
            'component_basis_map_certified': False,
            'matchete_rgbeta_quartic_normalization_certified': False,
            'log_path': str(log_path),
        }
        # Preserve the full Matchete log, even when the run fails.
        copy_mode = 'wb' if args.overwrite else 'xb'
        with log_path.open(copy_mode) as destination, (tmp / 'matchete.log').open('rb') as source_log:
            shutil.copyfileobj(source_log, destination)
        write_json_exclusive(output, report, args.overwrite)

    print(f'{model}: {report["status"]}')
    print(f'JSON: {output}')
    print(f'LOG:  {log_path}')
    return 0 if valid else 1


def self_test() -> int:
    example = (SCRIPT_DIR_ANCHOR + '\n'
               + 'lag=(built["LSM"] + free + rotated) // Expand // RelabelIndices;\n'
               + AUDIT_POINT_ANCHOR + '\n')
    result = build_probe(example, Path('/tmp/local/Lagrangian'), Path('/tmp/audit.json'))
    assert 'mfIndexLists = Cases[' in result
    assert 'lag=(built["LSM"] + free + rotated) // Expand // RelabelIndices;' in result
    assert result.index('mfIndexLists = Cases[') < result.index(AUDIT_POINT_ANCHOR)
    assert result.count('T3MFAUDIT: diagnostic ONLY') == 1
    print('PASS: source-anchor and non-mutating diagnostic instrumentation self-test')
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', help='existing B/C prepared mass-eigenstate parameter JSON')
    parser.add_argument('--output', help='new read-only MF audit JSON output')
    parser.add_argument('--wolframscript', default='wolframscript')
    parser.add_argument('--create-output-dir', action='store_true')
    parser.add_argument('--overwrite', action='store_true', help='explicitly allow replacing audit JSON and log')
    parser.add_argument('--self-test', action='store_true', help='local instrumentation check; no Matchete needed')
    args = parser.parse_args()
    if args.self_test:
        return self_test()
    if not args.input or not args.output:
        parser.error('--input and --output are required unless --self-test is used')
    try:
        return run(args)
    except (AuditError, FileNotFoundError, OSError, ValueError, json.JSONDecodeError) as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        return 2


if __name__ == '__main__':
    sys.exit(main())
