"""Locally rename stage-numbered research source files, not scientific artifacts.

This is an intentionally conservative, transactional, Git-aware migration.
Default is a read-only plan. Use --apply --confirm RENAME_SOURCE_FILES to apply.
No git operations that write to Git, no Wolfram runs, and no changes under output/.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
from datetime import datetime, timezone

EXPECTED_HEAD = 'fe0f511104cce8d7d914d114e3fd884eb47c6a95'
CONFIRM = 'RENAME_SOURCE_FILES'

# These are concise descriptions of the source *work*, not claims about
# whether the corresponding historical physics calculation is certified.
TOPICS = {
    6: 'ScalarQuarticResponseAudit',
    7: 'C5LinearResponseCorrection',
    8: 'WolframSelfTestArgumentFix',
    9: 'MatcheteTensorInitializationFix',
    10: 'MatchingRGSourcePilots',
    11: 'C5FiniteMassLogarithms',
    12: 'ExtraC5PhaseResponse',
    13: 'CoupledUVRGE',
    14: 'ConditionalUVToC5',
    15: 'QuarticConventionC5',
    16: 'FermionFirstThresholds',
    17: 'EFT1ScalarPairMixing',
    18: 'EFT1WilsonSourceTransport',
    19: 'WeinbergSourceTransport',
    20: 'GeneralQuarticRGBeta',
    21: 'QuarticExportAudit',
    22: 'ScalarTensorCertification',
    23: 'QuarticPrefactorAudit',
    24: 'BCQuarticCovariance',
    25: 'NativeQuarticAction',
    26: 'QuarticActionIdentity',
    27: 'TripletIntertwiner',
    28: 'RGBetaTripletStructure',
    29: 'UVConventionHandoff',
    30: 'CoupledA0UVRunning',
    31: 'CoupledUVRunning',
    33: 'FermionThresholdHandoff',
    34: 'IntermediateEFTRunning',
    35: 'ScalarThresholdHandoff',
    36: 'FinalWeinbergRunning',
}

# Some titles contain only STATUS or contain an opaque numeric SU(2) key.
# These special titles get an explicit descriptive filename.
DOC_OVERRIDES = {
    'GENERAL_T3_STAGE6_STATUS.md': 'GENERAL_T3_SCALAR_QUARTIC_VALIDATION.md',
    'GENERAL_T3_STAGE7_STATUS.md': 'GENERAL_T3_C5_LINEAR_RESPONSE_CORRECTION.md',
    'GENERAL_T3_STAGE8_STATUS.md': 'GENERAL_T3_WOLFRAMSCRIPT_SELFTEST_FIX.md',
    'GENERAL_T3_STAGE9_STATUS.md': 'GENERAL_T3_MATCHETE_TENSOR_INITIALIZATION.md',
    'GENERAL_T3_STAGE10_STATUS.md': 'GENERAL_T3_MATCHING_RG_SOURCE_PILOTS.md',
    'GENERAL_T3_STAGE11_STATUS.md': 'GENERAL_T3_C5_FINITE_MASS_LOGARITHMS.md',
    'GENERAL_T3_STAGE12_STATUS.md': 'GENERAL_T3_EXTRA_C5_PHASE_RESPONSE.md',
    'GENERAL_T3_STAGE13_STATUS.md': 'GENERAL_T3_COUPLED_UV_RUNNING.md',
    'GENERAL_T3_STAGE14_STATUS.md': 'GENERAL_T3_CONDITIONAL_UV_TO_C5.md',
    'GENERAL_T3_STAGE15_STATUS.md': 'GENERAL_T3_QUARTIC_CONVENTION_C5.md',
    'GENERAL_T3_STAGE16_STATUS.md': 'GENERAL_T3_FERMION_FIRST_THRESHOLDS.md',
    'GENERAL_T3_STAGE17_STATUS.md': 'GENERAL_T3_EFT1_SCALAR_PAIR_MIXING.md',
    'GENERAL_T3_STAGE18_STATUS.md': 'GENERAL_T3_EFT1_WILSON_SOURCE_TRANSPORT.md',
    'GENERAL_T3_STAGE19_STATUS.md': 'GENERAL_T3_WEINBERG_SOURCE_TRANSPORT.md',
    'GENERAL_T3_STAGE20_STATUS.md': 'GENERAL_T3_GENERAL_QUARTIC_RGBETA_PROGRAMME.md',
    'GENERAL_T3_STAGE21_EXPORT_AUDIT.md': 'GENERAL_T3_RGBETA_EXPORT_COVERAGE_AUDIT.md',
    'GENERAL_T3_STAGE22_TENSOR_CERTIFICATION.md': 'GENERAL_T3_SCALAR_TENSOR_CERTIFICATION.md',
    'GENERAL_T3_STAGE22_110101_LOCAL_CORRECTION.md': 'GENERAL_T3_TENSOR_110101_LOCAL_CORRECTION.md',
    'GENERAL_T3_STAGE22_110101_SU2_COVARIANCE.md': 'GENERAL_T3_TENSOR_110101_SU2_COVARIANCE.md',
    'GENERAL_T3_STAGE22_ALIGNED_GEOMETRY.md': 'GENERAL_T3_TENSOR_ALIGNED_GEOMETRY.md',
    'GENERAL_T3_STAGE22_EXPLICIT_GEOMETRY.md': 'GENERAL_T3_TENSOR_EXPLICIT_GEOMETRY.md',
    'GENERAL_T3_STAGE31_16_MODEL_UV_RUNNING.md': 'GENERAL_T3_16_MODEL_UV_RUNNING.md',
    'GENERAL_T3_STAGE36_EXISTING_C5_AND_WEINBERG.md': 'GENERAL_T3_FINAL_C5_AND_WEINBERG_RUNNING.md',
}

# Preferred full names for Python modules that are widely imported in the
# current numerical study. Only names change; the solver code is untouched.
SCRIPT_OVERRIDES = {
    'ExportGeneralT3Stage22Matchete': 'ExportGeneralT3ScalarTensorMatchete',
    'ExportGeneralT3Stage22RGBetaDeclarations': 'ExportGeneralT3ScalarTensorRGBetaDeclarations',
    'ExportGeneralT3Stage22RGBetaInvariants': 'ExportGeneralT3ScalarTensorRGBetaInvariants',
    'ExportGeneralT3Stage22SU2NativeCovariance': 'ExportGeneralT3SU2NativeCovariance',
    'ExportGeneralT3Stage24BCNativeCovariance': 'ExportGeneralT3BCQuarticNativeCovariance',
    'ExportGeneralT3Stage24BCPostFix': 'ExportGeneralT3BCQuarticCovarianceCorrection',
    'ExportGeneralT3Stage25NativeQuarticAssembly': 'ExportGeneralT3NativeQuarticAssembly',
    'ExportGeneralT3Stage28NativeTripletStructure': 'ExportGeneralT3NativeTripletStructure',
    'ExportGeneralT3Stage30A0OneLoop': 'ExportGeneralT3A0UVOneLoop',
    'ExportGeneralT3Stage31OneLoop': 'ExportGeneralT3AllModelsUVOneLoop',
    'ExportGeneralT3Stage34EFT1OneLoop': 'ExportGeneralT3EFT1OneLoop',
    'VerifyGeneralT3Stage22Live110101Fix': 'VerifyGeneralT3SU2Tensor110101Correction',
    'VerifyGeneralT3Stage24BCPostFix': 'VerifyGeneralT3BCQuarticCovarianceCorrection',
    'VerifyGeneralT3Stage24LiveBC': 'VerifyGeneralT3BCQuarticNativeInvariants',
    'Stage22RGBeta110101Covariance': 'RGBeta110101Covariance',
    'Stage22RGBetaAlignedGeometry': 'RGBetaAlignedGeometry',
    'Stage22RGBetaExplicitBasis': 'RGBetaExplicitBasis',
    'Stage22RGBetaNativeAccess': 'RGBetaNativeAccess',
    'Stage22TensorAlgebra': 'RGBetaTensorAlgebra',
    'Stage23QuarticPrefactorAudit': 'QuarticPrefactorAudit',
    'Stage24BCRemainderCovariance': 'BCQuarticRemainderCovariance',
    'Stage25QuarticActionAudit': 'NativeQuarticActionAudit',
    'Stage26QuarticActionIdentity': 'QuarticActionIdentity',
    'Stage27TripletIntertwiner': 'TripletIntertwiner',
    'Stage28TripletNativeStructureAudit': 'RGBetaTripletNativeStructureAudit',
    'Stage29ConventionUVHandoff': 'UVConventionHandoff',
    'Stage30CoupledA0UVRunning': 'CoupledA0UVRunning',
    'Stage31CoupledUVRunning': 'CoupledUVRunning',
    'Stage33FermionThresholdHandoff': 'FermionThresholdHandoff',
    'Stage34IntermediateEFTRunning': 'IntermediateEFTRunning',
    'Stage35ScalarThresholdHandoff': 'ScalarThresholdHandoff',
    'Stage36FinalWeinbergRunning': 'FinalWeinbergRunning',
}

ELIGIBLE_SUFFIXES = {
    '.py', '.ps1', '.wl', '.md', '.rst', '.txt', '.sh', '.bat', '.ini',
    '.toml', '.yml', '.yaml',
}
EXCLUDED_PREFIXES = (
    'output/', 'Reports/output/', 'docs/archive/', 'tests/fixtures/',
    '.github/workflows/archive/', 'RGE/general/generated/',
    'Numerical/benchmarks/',
)


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def snake(name: str) -> str:
    return re.sub(r'(?<!^)(?=[A-Z])', '_', name).lower()


def rename_name(rel: str) -> str | None:
    """Return renamed relative path or None for paths left untouched."""
    p = Path(rel)
    if p.suffix.lower() not in ('.py', '.ps1', '.wl', '.md'):
        return None
    if rel.startswith('scripts/') and p.parent == Path('scripts'):
        stem = p.stem
        if stem in SCRIPT_OVERRIDES:
            updated = SCRIPT_OVERRIDES[stem]
        else:
            m = re.fullmatch(r'RunGeneralT3Stage(\d+)', stem)
            if m and int(m.group(1)) in TOPICS:
                updated = 'RunGeneralT3' + TOPICS[int(m.group(1))]
            else:
                m = re.search(r'Stage(\d+)', stem)
                if not m or int(m.group(1)) not in TOPICS:
                    return None
                # Preserve already descriptive words after the stage number.
                # For export/verification scripts, a short family label is useful.
                prefix = stem[:m.start()]
                suffix = stem[m.end():]
                if prefix.startswith(('ExportGeneralT3', 'VerifyGeneralT3')):
                    family = ('Tensor' if int(m.group(1)) == 22 else
                              'BCQuartic' if int(m.group(1)) == 24 else '')
                    updated = prefix + family + suffix
                else:
                    updated = prefix + suffix
        return str(p.with_name(updated + p.suffix)).replace('\\', '/')

    if rel.startswith('docs/') and p.parent == Path('docs'):
        if p.name in DOC_OVERRIDES:
            return str(p.with_name(DOC_OVERRIDES[p.name])).replace('\\', '/')
        m = re.match(r'GENERAL_T3_STAGE(\d+)_(.+)\.md$', p.name)
        if m and int(m.group(1)) in TOPICS:
            stage, suffix = int(m.group(1)), m.group(2)
            family = 'TENSOR_' if stage == 22 else 'GENERAL_QUARTIC_' if stage == 20 else ''
            updated = f'GENERAL_T3_{family}{suffix}.md'
            return str(p.with_name(updated)).replace('\\', '/')
        return None

    if rel.startswith('tests/python/') and p.parent == Path('tests/python'):
        m = re.match(r'(test_general_t3_)stage(\d+)(.*)', p.stem)
        if m and int(m.group(2)) in TOPICS:
            stage = int(m.group(2))
            tail = m.group(3).lstrip('_')
            # Keep detailed test names; give ambiguous numeric tests a family.
            if tail:
                # Stages 25 and 26 both have a test named 'quartic_action'.
                # Removing the stage number verbatim would collapse two
                # distinct scientific checks onto one file.
                if stage == 25 and tail == 'quartic_action':
                    tail = 'native_quartic_action'
                elif stage == 26 and tail == 'quartic_action':
                    tail = 'quartic_action_identity'
                elif stage == 22:
                    tail = 'tensor_' + tail
                elif stage == 20:
                    tail = 'general_quartic_' + tail
                updated = m.group(1) + tail
            else:
                updated = m.group(1) + snake(TOPICS[stage])
            return str(p.with_name(updated + p.suffix)).replace('\\', '/')
    return None


def tracked_paths(root: Path) -> list[str]:
    p = subprocess.run(['git', 'ls-files', '-z'], cwd=root,
                       capture_output=True, check=True)
    return sorted(s.decode('utf-8', errors='surrogateescape').replace('\\', '/')
                  for s in p.stdout.split(b'\0') if s)


def head(root: Path) -> str:
    return subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=root, text=True).strip()


def replacement_pairs(rename: dict[str, str]) -> list[tuple[str, str]]:
    """Replace qualified filename references and Python/Wolfram stem tokens.

    Names for history data output/general_t3_stageNN remain unchanged. We do
    NOT globally replace 'Stage 22', 'stage22', scientific statuses or run IDs.
    """
    replace = {}
    for old, new in rename.items():
        op, np = Path(old), Path(new)
        replace[op.name] = np.name
        # Both dotted/qualified module invocation and Python import identifier.
        if old.startswith('scripts/'):
            replace[op.stem] = np.stem
        # Test file stem in explicit pytest module/file path invocations.
        if old.startswith('tests/python/'):
            replace[op.stem] = np.stem
        # Docs path references may omit .md, but avoid broad substrings.
    return sorted(replace.items(), key=lambda kv: (-len(kv[0]), kv[0]))


def replace_text(src: bytes, pairs: list[tuple[str, str]], rel: str) -> bytes:
    try:
        value = src.decode('utf-8')
    except UnicodeDecodeError:
        # BOM/legacy encodings are not silently rewritten.
        return src
    result = value
    for old, new in pairs:
        # only whole identifiers or filenames; prevent accidental replacement
        # inside longer Python names and source-generated identifiers.
        pattern = r'(?<![A-Za-z0-9_])' + re.escape(old) + r'(?![A-Za-z0-9_])'
        # Filenames and identifiers contain underscores, so exact matching is
        # chosen: Python references must be full qualified names.
        result = re.sub(pattern, lambda _: new, result)
    if rel.startswith('docs/GENERAL_T3_') and re.search(r'STAGE\d+', rel):
        # Titles are reading-facing. Internal stage numbering remains history.
        result = re.sub(r'^(#\s+)(?:Experimental\s+)?Stage\s+\d+\s*[:—-]\s*',
                        r'\1General T3: ', result, count=1, flags=re.MULTILINE)
    return result.encode('utf-8') if result != value else src


def eligible(rel: str) -> bool:
    return not rel.startswith(EXCLUDED_PREFIXES) and Path(rel).suffix.lower() in ELIGIBLE_SUFFIXES


def make_plan(root: Path) -> dict:
    tracked = tracked_paths(root)
    valid = [rel for rel in tracked if (root/rel).is_file() and eligible(rel)]
    rename = {rel: new for rel in valid if (new := rename_name(rel)) and new != rel}
    destinations: dict[str, str] = {}
    for old, new in sorted(rename.items()):
        previous = destinations.get(new)
        if previous is not None:
            raise RuntimeError(
                f'Renaming collision: {previous} and {old} -> {new}'
            )
        destinations[new] = old
    existing = set(tracked)
    for old, new in rename.items():
        if new in existing or (root/new).exists():
            raise RuntimeError(f'Destination already exists: {old} -> {new}')
    pairs = replacement_pairs(rename)
    updates = {}
    for rel in valid:
        source = (root/rel).read_bytes()
        result = replace_text(source, pairs, rel)
        if source != result or rel in rename:
            updates[rel] = (source, result)
    return {'rename': rename, 'updates': updates,
            'tracked_files_scanned': len(valid), 'pairs': pairs}


def dirty_tracked_paths(root: Path) -> set[str]:
    raw = subprocess.check_output(['git', 'status', '--porcelain=v1', '-z',
                                   '--untracked-files=no'], cwd=root)
    # Porcelain -z entries contain two status bytes, one space, and a path;
    # rename entries have an additional origin path.
    entries = raw.split(b'\0')
    changed = set()
    i = 0
    while i < len(entries):
        entry = entries[i]
        i += 1
        if not entry:
            continue
        status = entry[:2].decode('ascii', errors='replace')
        path = entry[3:].decode('utf-8', errors='surrogateescape').replace('\\', '/')
        changed.add(path)
        if 'R' in status or 'C' in status:
            if i < len(entries) and entries[i]:
                changed.add(entries[i].decode('utf-8', errors='surrogateescape').replace('\\', '/'))
                i += 1
    return changed


def apply_plan(root: Path, plan: dict, backup: Path) -> dict:
    if head(root) != EXPECTED_HEAD:
        raise RuntimeError('HEAD changed; rerun against inspected revision')
    changed = dirty_tracked_paths(root)
    conflict = sorted(changed.intersection(plan['updates']))
    if conflict:
        raise RuntimeError('Affected tracked files have local edits: '+', '.join(conflict[:15]))
    if backup.exists():
        raise RuntimeError(f'Backup already exists; inspect it before rerunning: {backup}')
    changed_files = sorted(plan['updates'])
    manifest = {'baseline_commit': EXPECTED_HEAD,
                'created_utc': datetime.now(timezone.utc).isoformat(),
                'backup_root': str(backup),
                'source_files': [
                    {'path': rel, 'sha256': sha(plan['updates'][rel][0]),
                     'renamed_to': plan['rename'].get(rel),
                     'post_sha256': sha(plan['updates'][rel][1])}
                    for rel in changed_files],
                'renames': plan['rename'],
                'note': 'Restore original relative paths from originals/ after removing renamed destinations.'}
    # A complete, SHA-verified external backup is created BEFORE any writes.
    for rel in changed_files:
        original = backup/'originals'/rel
        original.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(root/rel, original)
        if sha(original.read_bytes()) != sha(plan['updates'][rel][0]):
            raise IOError(f'Backup verification failed for {rel}')
    (backup/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n', encoding='utf-8')
    # Stage all target bytes to temporary files; only then replace/move.
    temps = {}
    try:
        for index, rel in enumerate(changed_files):
            dest_rel = plan['rename'].get(rel, rel)
            dest = root/dest_rel
            dest.parent.mkdir(parents=True, exist_ok=True)
            tmp = dest.with_name(dest.name + f'.rename-prepared-{index}.tmp')
            if tmp.exists():
                raise RuntimeError(f'Preparation file already exists: {tmp}')
            tmp.write_bytes(plan['updates'][rel][1])
            temps[rel] = tmp
        for rel in changed_files:
            dest_rel = plan['rename'].get(rel, rel)
            os.replace(temps[rel], root/dest_rel)
            if rel in plan['rename']:
                (root/rel).unlink()
        for rel in changed_files:
            dest = root/plan['rename'].get(rel, rel)
            if sha(dest.read_bytes()) != sha(plan['updates'][rel][1]):
                raise IOError(f'Post-migration SHA mismatch: {dest}')
    except BaseException:
        # The original copy of every affected file is already verified.
        # Restore on interrupted/failed application; keep the backup too.
        for rel in changed_files:
            dest_rel = plan['rename'].get(rel, rel)
            if dest_rel != rel and (root/dest_rel).is_file():
                (root/dest_rel).unlink()
            original = backup/'originals'/rel
            if original.is_file():
                (root/rel).parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(original, root/rel)
        raise
    finally:
        for tmp in temps.values():
            if tmp.exists():
                tmp.unlink()
    return manifest


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    p.add_argument('--apply', action='store_true', help='Execute after backing up originals outside the repository')
    p.add_argument('--confirm', default='')
    p.add_argument('--report', type=Path, help='Optional JSON manifest of the read-only plan')
    args = p.parse_args(argv)
    root = args.root.resolve()
    if head(root) != EXPECTED_HEAD:
        raise SystemExit('STOP: checkout HEAD differs from inspected '+EXPECTED_HEAD)
    plan = make_plan(root)
    rename = plan['rename']
    updates = plan['updates']
    report = {'baseline_commit': EXPECTED_HEAD,
              'renamed_count': len(rename),
              'edited_or_renamed_files': len(updates),
              'tracked_text_scanned': plan['tracked_files_scanned'],
              'renames': rename,
              'modified_without_rename': sorted(set(updates)-set(rename)),
              'excluded': ['output/**', 'Reports/output/**', 'docs/archive/**',
                           'tests/fixtures/**'],
              'dry_run': not args.apply}
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print('SOURCE FILE RENAMES:', len(rename))
    print('OTHER FILES WITH UPDATED REFERENCES:', len(updates)-len(rename))
    print('ELIGIBLE TRACKED FILES SCANNED:', plan['tracked_files_scanned'])
    if not args.apply:
        print('DRY RUN: no source or numerical files changed')
        for old, new in sorted(rename.items()):
            print(old, '->', new)
        return 0
    if args.confirm != CONFIRM:
        raise SystemExit('STOP: --apply requires --confirm '+CONFIRM)
    backup = root.parent / ('T3_stage_filename_backup_'+EXPECTED_HEAD[:7])
    result = apply_plan(root, plan, backup)
    print('APPLIED:', len(result['renames']), 'renames')
    print('ORIGINAL FILE BACKUP:', backup)
    print('BACKUP MANIFEST:', backup/'manifest.json')
    print('Run Python tests now; external Matchete/RGBeta tools not exercised.')
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
