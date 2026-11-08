"""Read-only A(-4) EFT1 RGBeta group-contraction / quartic-projector audit.

Uses the freshly generated native Mathematica diagnostic, the earlier verified
9-operator witness, and the exact local model source.  Optionally scans the
*installed* RGBeta .m/.wl package files to locate contraction and projector
rules; its excerpts are evidence, never assumed executable identities.

NO physics file is modified.  This does not certify a Matchete coupling map.
"""
from __future__ import annotations

import argparse
import ast
from hashlib import sha256
import json
from pathlib import Path
import re
from typing import Any, Mapping

ROOT = Path(__file__).resolve().parents[1]
SOURCE_REL = Path('RGE/running/rgbeta/T3RGBetaModel.wl')
FORMAT = 'a-minus4-rgbeta-symbolic-contractions-v1'
STATUS = 'NATIVE_SYMBOLIC_CONTRACTION_DIAGNOSTICS_CAPTURED'
OLD_FORMAT = 'a-minus4-rgbeta-operator-witness-v1'
GROUP_IDS = (
    'fund_casimir_3_over_4',
    'triplet_casimir_2',
    'triplet_fully_contracted_generator_norm_6',
    'fund_fully_contracted_generator_norm_3_over_2',
)
QUARTIC_NAMES = (
    'lambdaS2', 'lambdaS2Adj', 'lambdaH2', 'lambdaH2Adj', 'lambdaT3',
)
GROUP_HEADS = {
    'lambdaS2': ('del[SU2L[S2]',),
    'lambdaS2Adj': ('tGen[SU2L[S2]',),
    'lambdaH2': ('del[SU2L[fund]', 'del[SU2L[S2]'),
    'lambdaH2Adj': ('tGen[SU2L[fund]', 'tGen[SU2L[S2]'),
    'lambdaT3': ('delS2[SU2L',),
}
TARGET_STATUSES = {
    'EVALUATED_MATCH', 'SYMBOLICALLY_UNRESOLVED_OR_DIFFERENT',
    'EVALUATION_BLOCKED', 'TIMED_OUT',
}
# Exact installed RGBeta source signatures to search for. This is a source
# locator, not an assertion about which internal algorithms are implemented.
KEYWORDS = re.compile(
    r'tGen\s*\[|delS2\s*\[|CheckPotentialCouplingMixing|BetaTerm\s*\[|'
    r'Quartic|[Pp]roject|[Ff]ierz|[Cc]ontract|[Cc]asimir|'
    r'Lam\s*\[|[Gg]roup[Ii]nv',
)


class ContractionAuditError(ValueError):
    """Missing/changed native evidence or unsafe certification assertion."""


def require(test: Any, message: str) -> None:
    if not test:
        raise ContractionAuditError(message)


def read_json(path: Path) -> dict:
    data = json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(data, dict), f'{path}: expected JSON object')
    return data


def sha_file(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def normalize(value: str) -> str:
    require(isinstance(value, str), 'Missing native InputForm')
    return re.sub(r'\s+', '', value)


def parse_installed_path(value: str) -> Path:
    require(isinstance(value, str) and 'RGBeta' in value,
            'Missing installed RGBeta package provenance')
    # OPString[FindFile["RGBeta`"]] returns an InputForm-quoted string;
    # literal_eval decodes its escaping without executing arbitrary code.
    if '"' in value:
        require(value.startswith('"') and value.endswith('"'),
                'Malformed quoted package path InputForm')
        try:
            value = ast.literal_eval(value)
        except (ValueError, SyntaxError, UnicodeError) as exc:
            raise ContractionAuditError('Malformed package path InputForm') from exc
    require(isinstance(value, str) and 'RGBeta' in value and
            '"' not in value and '\x00' not in value,
            'Invalid installed package path')
    return Path(value)


def scan_installed_package(package_init: Path) -> dict[str, Any]:
    """Capture bounded snippets from the installed package's exact sources.

    Do not recursively traverse outside the installed RGBeta package root.
    Only inspect .m/.wl files smaller than 2 MB and skip binary/large files.
    """
    if not package_init.is_file():
        return {'status': 'INSTALLED_PACKAGE_NOT_ACCESSIBLE',
                'package_init': str(package_init), 'source_files': [],
                'relevant_source_snippets': []}
    root = package_init.parent.parent if package_init.parent.name.lower() == 'kernel' \
        else package_init.parent
    sources = sorted(set(root.rglob('*.m')) | set(root.rglob('*.wl')))
    sources = [p for p in sources if p.is_file() and 0 < p.stat().st_size <= 2_000_000]
    require(len(sources) <= 300,
            'Unexpectedly large installed RGBeta source inventory')
    manifest = []
    snippets = []
    for path in sources:
        content = path.read_text(encoding='utf-8-sig', errors='replace').splitlines()
        relative = path.relative_to(root).as_posix()
        matches = [(i + 1, line.strip()) for i, line in enumerate(content)
                   if KEYWORDS.search(line)]
        manifest.append({'path': relative, 'sha256': sha_file(path),
                         'line_count': len(content), 'matching_lines': len(matches)})
        # Keep lines especially useful for exact generator/projector diagnosis.
        priority = sorted(matches, key=lambda pair: (
            0 if re.search(r'CheckPotentialCouplingMixing|[Pp]roject|[Ff]ierz|'
                           r'tGen\s*\[|delS2\s*\[', pair[1]) else 1,
            pair[0]))
        for number, line in priority[:12]:
            snippets.append({'source': relative, 'line': number,
                             'code_excerpt': line[:400]})
    snippets = sorted(snippets, key=lambda x: (
        0 if re.search(r'[Pp]roject|[Ff]ierz|CheckPotentialCouplingMixing',
                       x['code_excerpt']) else 1,
        x['source'], x['line']))[:120]
    return {'status': 'INSTALLED_PACKAGE_SOURCE_SCANNED',
            'package_init': str(package_init),
            'package_source_root': str(root),
            'source_files': manifest,
            'relevant_source_snippets': snippets}


def analyze(data: Mapping[str, Any], previous: Mapping[str, Any],
            root: Path, *, source_scan: bool = True) -> dict[str, Any]:
    require(isinstance(data, Mapping) and isinstance(previous, Mapping),
            'Native/previous records must be JSON mappings')
    require(data.get('format') == FORMAT and data.get('status') == STATUS,
            'Wrong/blocked native contraction witness')
    require(previous.get('format') == OLD_FORMAT and
            previous.get('status') ==
            'NATIVE_SYMBOLIC_OPERATOR_TENSOR_CAPTURED_COMPONENT_MAP_PENDING',
            'Previous nine-operator witness missing/blocked')
    for record in (data, previous):
        require(record.get('model') == 'A_-4' and record.get('stage') == 'EFT1'
                and record.get('dimensions_alpha') == [1, 3, 2, -4],
                'Wrong model or EFT stage')
        require(record.get('source_model_path') == SOURCE_REL.as_posix(),
                'Incorrect native model-source path')
    source = root / SOURCE_REL
    require(source.is_file(), f'Missing model source {source}')
    source_sha = sha_file(source)
    require(data.get('source_model_sha256') == source_sha and
            previous.get('source_model_sha256') == source_sha,
            'Native evidence is stale relative to the current model source')
    pkg = data.get('rgbeta_package_file')
    require(pkg == previous.get('rgbeta_package_file'),
            'Native probes used different installed RGBeta package paths')
    init_path = parse_installed_path(pkg)
    require(previous.get('native_operator_registry_complete') is True and
            previous.get('native_Lam_formal_tensor_available') is True,
            'Previous native quartic registry or Lam tensor incomplete')
    for flag in ('native_group_reductions_certified',
                 'quartic_projection_certified',
                 'physical_matchete_coupling_map_certified',
                 'native_to_reference_axes_certified', 'source_modified'):
        require(data.get(flag) is False,
                f'Unsafe or missing native witness flag: {flag}')
    require(previous.get('matchete_physical_coupling_map_certified') is False
            and previous.get('native_to_reference_component_axes_certified') is False,
            'Prior witness claims unsupported physical certification')
    group = data.get('generator_contraction_probes')
    require(isinstance(group, list) and len(group) == 4,
            'Expected four native generator probes')
    names = [p.get('id') for p in group if isinstance(p, Mapping)]
    require(len(names) == 4 and set(names) == set(GROUP_IDS),
            'Incorrect/duplicate generator probe inventory')
    checked = []
    for probe in group:
        identity = probe['id']
        status = probe.get('status')
        require(status in TARGET_STATUSES, f'{identity}: invalid status')
        truth = probe.get('native_exact_target_match')
        require(type(truth) is bool and (status == 'EVALUATED_MATCH') == truth,
                f'{identity}: inconsistent native equality claim')
        if truth:
            require(normalize(probe.get('native_difference_inputform')) == '0',
                    f'{identity}: exact match without a zero symbolic difference')
        for key in ('source_expression_inputform',
                    'native_evaluated_inputform', 'target_inputform',
                    'native_difference_inputform'):
            val = probe.get(key)
            require(isinstance(val, str) and len(val) > 0,
                    f'{identity}: missing native expression {key}')
        src = normalize(probe['source_expression_inputform'])
        require('tGen[SU2L[' in src and 'HoldComplete[' in src,
                f'{identity}: not a native held generator contraction')
        checked.append({'id': identity, 'status': status,
                        'native_exact_target_match': truth,
                        'native_difference_inputform':
                        probe['native_difference_inputform']})
    quartic = data.get('quartic_invariant_samples')
    require(isinstance(quartic, list) and len(quartic) == 5,
            'Missing A(-4) quartic native invariant samples')
    require({r.get('name') for r in quartic if isinstance(r, Mapping)} ==
            set(QUARTIC_NAMES), 'Quartic sample inventory mismatch')
    for row in quartic:
        name = row['name']
        require(row.get('captured') is True,
                f'{name}: native quartic invariant could not be read')
        expr = normalize(row.get('native_invariant_inputform'))
        require('$Failed' not in expr and '$TimedOut' not in expr,
                f'{name}: uncomputed invariant')
        for head in GROUP_HEADS[name]:
            require(head in expr, f'{name}: missing {head} in native invariant')
    double_calls = data.get('native_dummy_index_double_call_probes')
    required_pairs = {
        'S2Adj_times_S2Adj', 'H2Adj_times_H2Adj',
        'S2Adj_times_H2Adj', 'S2Singlet_times_S2Singlet_control',
    }
    require(isinstance(double_calls, list) and len(double_calls) == 4 and
            {r.get('id') for r in double_calls if isinstance(r, Mapping)} ==
            required_pairs, 'Missing/duplicate native dummy-index probes')
    dummy_count = 0
    dummy_rows = []
    for row in double_calls:
        sharing = row.get('same_literal_A_present_in_both_raw_calls')
        require(type(sharing) is bool,
                'Native dummy-index probe did not report a boolean')
        a = normalize(row.get('first_call_inputform'))
        b = normalize(row.get('second_call_inputform'))
        require('$Failed' not in a+b and len(a) > 0 and len(b) > 0,
                'Native double-call invariant unavailable')
        # A shared raw dummy is evidence of an indexing risk, not proof that
        # RGBeta's subsequent quartic projector fails to alpha-rename it.
        if sharing:
            require('tGen[' in a and 'tGen[' in b and
                    ',A,' in a and ',A,' in b,
                    'Raw A-sharing claim lacks the literal shared generator A')
        dummy_count += int(sharing)
        dummy_rows.append({'id': row['id'],
                           'same_literal_A_present_in_both_raw_calls': sharing})
    sym_names = data.get('native_projector_related_symbol_names')
    sym_rows = data.get('native_projector_symbol_definitions')
    require(isinstance(sym_names, list) and len(sym_names) <= 80 and
            len(sym_names) == len(set(sym_names)) and
            all(isinstance(s, str) and s.startswith('RGBeta`') for s in sym_names),
            'Invalid package projector symbol inventory')
    require(isinstance(sym_rows, list) and len(sym_rows) <= 25 and
            [p.get('symbol') for p in sym_rows] == sym_names[:len(sym_rows)],
            'Inconsistent package projector introspection')
    for row in sym_rows:
        require(isinstance(row.get('definition_prefix_inputform'), str) and
                isinstance(row.get('full_definition_length'), int) and
                type(row.get('truncated')) is bool,
                'Malformed native projector definition record')
    scan = scan_installed_package(init_path) if source_scan else {
        'status': 'INSTALLED_PACKAGE_SCAN_NOT_REQUESTED', 'source_files': [],
        'relevant_source_snippets': []}
    successful = sum(p['native_exact_target_match'] for p in checked)
    return {
        'format': 'a-minus4-rgbeta-symbolic-contraction-check-v1',
        'status': 'NATIVE_GROUP_AND_PROJECTOR_EVIDENCE_RECORDED_UNCERTIFIED',
        'model': 'A_-4', 'stage': 'EFT1',
        'model_source_sha256': source_sha,
        'native_package_init': str(init_path),
        'native_operator_registry_verified': 9,
        'native_Lam_formal_tensor_recorded': True,
        'generator_probes_evaluated_to_reference_targets': successful,
        'generator_probes_total': 4,
        'generator_probe_details': checked,
        'quartic_native_invariants_recorded': [r['name'] for r in quartic],
        'native_raw_dummy_index_reuse_probes': dummy_rows,
        'native_raw_dummy_index_reuse_count': dummy_count,
        'native_projector_alpha_renaming_verified': False,
        'projector_native_symbol_count': len(sym_names),
        'projector_native_definitions_recorded': len(sym_rows),
        'installed_package_source': scan,
        'native_to_Matchete_component_map_certified': False,
        'scalar_beta_polynomial_certified': False,
        'physical_quartic_normalization_certified': False,
        'production_modified': False,
        'next_gate': 'Inspect installed RGBeta generator contraction and quartic projection definitions; compare raw Lam-loop contraction in native axes to globally projected beta coordinates',
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--native', type=Path, required=True)
    parser.add_argument('--previous-operator', type=Path, required=True)
    parser.add_argument('--report', type=Path, required=True)
    parser.add_argument('--no-package-source-scan', action='store_true')
    args = parser.parse_args()
    try:
        result = analyze(read_json(args.native), read_json(args.previous_operator),
                         ROOT, source_scan=not args.no_package_source_scan)
        require(args.report.resolve() not in {
            args.native.resolve(), args.previous_operator.resolve()},
            'Refusing to overwrite native evidence')
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(result, indent=2) + '\n',
                               encoding='utf-8')
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        print(f'BLOCKED: {exc}')
        return 2
    print('A(-4) native quartic registry: 9/9 retained')
    print('Native group contractions matching test targets:',
          result['generator_probes_evaluated_to_reference_targets'], '/ 4')
    for row in result['generator_probe_details']:
        print(' ', row['id'], ':', row['status'])
    scan = result['installed_package_source']
    print('Installed RGBeta source:', scan['status'],
          'files=', len(scan['source_files']))
    print('Raw literal adjoint A index reused:',
          result['native_raw_dummy_index_reuse_count'], '/ 4 probes')
    print('Native quartic projector symbols:',
          result['projector_native_symbol_count'])
    print('Physical scalar beta / Matchete map certified: NO')
    print('Report:', args.report)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
