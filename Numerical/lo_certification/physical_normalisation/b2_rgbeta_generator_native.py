"""B2 G2 stage 2: validate live RGBeta *symbolic* gauge-generator/field index evidence.

No RGBeta numerical eps/delS2 components, SU(2) axis map, action coefficient,
Matchete-to-RGBeta physical dictionary or physical certificates follow.
Source is guarded against all six reviewed production physics blobs. This file
never modifies the existing native G1/G2 reports or physics sources.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import re
import subprocess

from b2_gauge_coupling import AuditBlocked, CLASSES, SOURCE_BLOBS, require
from b2_rgbeta_ward_witness import source_provenance
from b2_rgbeta_internal_action import ALLOWED_B2_PATHS

FORMAT = 'b2-rgbeta-native-symbolic-generator-orientation-v1'
STATUS = 'NATIVE_SYMBOLIC_GENERATOR_RULES_CAPTURED_NUMERIC_AXES_PENDING'
SCRIPT = Path('scripts/RunPhysicalB2RGBetaGeneratorProbe.wl')
OUTPUT_DEFAULT = Path('output/general_observable/physical_normalisation/b2_rgbeta_generator_native_02/orientation_witness.json')
ENV_NAMES = {
    'RGE/running/rgbeta/T3RGBetaModel.wl': 'T3_B2_PROVENANCE_RGBETA',
    'Lagrangian/interactions/T3Topology.wl': 'T3_B2_PROVENANCE_TOPOLOGY',
    'Lagrangian/model/SU2Invariants.wl': 'T3_B2_PROVENANCE_SU2',
}
FLAGS = (
    'T3_EXPERIMENTAL_GENERAL_FERMIONS',
    'T3_EXPERIMENTAL_GENERAL_QUADRATICS',
    'T3_EXPERIMENTAL_GENERAL_QUARTICS',
    'T3_EXPERIMENTAL_FLAVOR_MASSES',
)


def child_environment(provenance: dict, parent: dict[str, str]) -> dict[str, str]:
    """Native Wolfram gets only Python-verified Git metadata, not inherited values."""
    require(provenance.get('all_six_source_blobs_unchanged') is True,
            'six-source checkout guard did not complete')
    head = provenance.get('repository_head_measured')
    require(isinstance(head, str) and re.fullmatch('[0-9a-f]{40}', head),
            'source guard returned an invalid Git HEAD')
    source_rows = provenance.get('validated_git_blobs')
    require(isinstance(source_rows, list) and len(source_rows) == len(SOURCE_BLOBS),
            'source guard did not cover all six reviewed paths')
    by_path = {}
    for row in source_rows:
        require(isinstance(row, dict), 'malformed source-guard row')
        path = row.get('path')
        require(path in SOURCE_BLOBS and path not in by_path
                and row.get('git_blob_sha') == SOURCE_BLOBS[path]
                and row.get('validated_with_git_clean_filter') is True,
                'unverified, changed or duplicate source-guard row')
        by_path[path] = row
    require(set(by_path) == set(SOURCE_BLOBS), 'incomplete six-source inventory')
    require(set(ENV_NAMES) == set(ALLOWED_B2_PATHS), 'G1/G2 source paths changed')
    require(all(parent.get(f, '0') != '1' for f in FLAGS),
            'unrestricted experimental T3 flag enabled')
    env = dict(parent)
    for name in ('T3_B2_PROVENANCE_HEAD', *ENV_NAMES.values()):
        env.pop(name, None)
    env['T3_B2_PROVENANCE_HEAD'] = head
    for path, key in ENV_NAMES.items():
        env[key] = SOURCE_BLOBS[path]
    return env


def validate_data(data: dict, provenance: dict) -> dict:
    """Pure JSON validation; never equate abstract generator rules to numeric axes."""
    require(isinstance(data, dict), 'native probe is not a JSON object')
    require(data.get('format') == FORMAT and data.get('status') == STATUS,
            'wrong native probe status / incomplete native capture')
    require(data.get('git_head_at_capture') == provenance.get('repository_head_measured'),
            'native probe captured at a different Git HEAD')
    require(data.get('git_source_blobs_at_capture') ==
            {path: SOURCE_BLOBS[path] for path in ALLOWED_B2_PATHS},
            'native source inventory does not match pinned physics blobs')
    require(data.get('git_provenance_capture_method') == 'Python six-source clean-filter guard',
            'native output lacks Python source-provenance handoff')
    require(isinstance(data.get('rgbeta_package_file'), str)
            and bool(data['rgbeta_package_file'].strip()),
            'installed RGBeta package path not captured')
    index_names = data.get('field_index_map_symbol_names')
    require(isinstance(index_names, list) and len(index_names) == 1
            and isinstance(index_names[0], str)
            and index_names[0].startswith('RGBeta`')
            and index_names[0].endswith('`$fieldIndexMap')
            and data.get('field_index_map_unique_package_scope') is True,
            'native $fieldIndexMap was not uniquely resolved in RGBeta package scope')
    sym = data.get('Tscal_symbol_names')
    require(isinstance(sym, list) and len(sym) == 1
            and isinstance(sym[0], str) and sym[0].endswith('`Tscal'),
            'could not uniquely find native Tscal symbol')
    down = data.get('Tscal_downvalues_inputform')
    require(isinstance(down, str) and len(down) > 50
            and 'tGen[' in down and 'Bar[' in down
            and '$fieldIndexMap' in down,
            'native Tscal definitions were not captured')
    require(data.get('Tscal_rule_has_conjugate_scalar_blocks') is True,
            'native scalar generator did not show both barred and unbarred blocks')
    require(data.get('tGen_bar_transpose_identity_native') is True,
            'native tGen[Bar[rep]] transpose relation did not hold')
    require(data.get('eps_sample_remains_symbolic') is True
            and data.get('delS2_sample_remains_symbolic') is True,
            'native sample behaviour changed: inspect before using convention')
    rows = data.get('classes')
    require(isinstance(rows, list) and len(rows) == 5,
            'native field registry did not cover exactly five classes')
    seen = set()
    for row in rows:
        require(isinstance(row, dict), 'non-object native class')
        family = row.get('family')
        require(family in CLASSES and family not in seen,
                'unknown or duplicate native class')
        seen.add(family)
        d1, d2, df, alphas = CLASSES[family]
        rep = row.get('dimensions_alpha')
        require(isinstance(rep, list) and rep[:3] == [d1,d2,df]
                and len(rep) == 4 and rep[3] in alphas,
                'incorrect class representation or alpha')
        require(row.get('registered') is True
                and row.get('scalar_representations_match_registration') is True
                and row.get('fermion_representations_match_registration') is True
                and row.get('all_scalar_field_and_conjugate_indices_present') is True
                and row.get('native_package_field_index_map_available') is True
                and row.get('yukawa_field_orders_match') is True,
                f'{family}: native field/field-index registration failed')
        maps = row.get('field_index_ids')
        require(isinstance(maps, dict) and set(maps) ==
                {'S1', 'Bar[S1]', 'S2', 'Bar[S2]'}
                and all(isinstance(v, int) and not isinstance(v, bool) and v > 0
                        for v in maps.values())
                and len(set(maps.values())) == 4,
                f'{family}: scalar/antiscalar registry index blocks ambiguous')
        for field, dimension in (('S1',d1), ('S2',d2), ('F',df)):
            expression = row.get(f'{field}_gauge_reps')
            require(isinstance(expression, str) and 'U1Y[' in expression,
                    f'{family}: missing registered {field} gauge representations')
            expected_su2 = {1: None, 2: 'SU2L[fund]', 3: 'SU2L[S2]'}[dimension]
            require((('SU2L[' not in expression) if expected_su2 is None
                    else expected_su2 in expression),
                    f'{family}: registered {field} SU2 dimensionality mismatch')
        require(row.get('physical_action_certified') is False
                and row.get('numeric_SU2_axes_certified') is False,
                f'{family}: conditional evidence falsely promoted')
    require(seen == set(CLASSES), 'native class registry coverage incomplete')
    require(data.get('RGBeta_numeric_eps_components_certified') is False
            and data.get('RGBeta_numeric_delS2_components_certified') is False
            and data.get('RGBeta_to_Matchete_basis_identified') is False
            and data.get('physical_coupling_dictionary_certified') is False
            and data.get('physical_source_certified_models') == 0
            and data.get('safe_to_integrate_into_production') is False
            and data.get('EFT_performed') is False,
            'native symbolic evidence overclaimed a physical result')
    return {
        'format': 'b2-g2-native-symbolic-generator-validation-v1',
        'status': STATUS,
        'local_git_head': provenance['repository_head_measured'],
        'families_registered': 5,
        'native_symbolic_scalar_generator_rule_captured': True,
        'native_barred_generator_transpose_verified': True,
        'numeric_SU2_component_axes_identified': False,
        'B_C_index_orientation_resolved_from_native_numeric_axes': False,
        'physical_source_certified_models': 0,
        'next_gate': 'Construct native field-index/action contraction in explicit common SU2 basis; do not infer it from tGen/eps symbolic names.',
    }


def read_json(path: Path) -> dict:
    require(path.is_file() and path.stat().st_size > 16,
            f'fresh native RGBeta JSON missing or empty: {path}')
    try:
        data = json.loads(path.read_text(encoding='utf-8-sig'))
    except (UnicodeError, json.JSONDecodeError) as exc:
        raise AuditBlocked(f'bad native JSON: {exc}') from exc
    require(isinstance(data, dict), 'native JSON is not an object')
    return data


def run_native(root: Path, output: Path, *, wolfram: str = 'wolframscript',
               parent_env: dict[str, str] | None = None) -> dict:
    root, output = root.resolve(), output.resolve()
    allowed = root/'output/general_observable/physical_normalisation'
    require(output.is_relative_to(allowed) and output.suffix.lower() == '.json',
            'native output must be a new JSON under physical_normalisation')
    require(not output.exists() and not Path(str(output)+'.partial').exists(),
            'refusing existing native output or leftover partial')
    script = root/SCRIPT
    require(script.is_file(), f'missing Wolfram source probe: {script}')
    parent = dict(os.environ if parent_env is None else parent_env)
    require(all(parent.get(f, '0') != '1' for f in FLAGS),
            'non-restricted T3 environment')
    proven = source_provenance(root)
    env = child_environment(proven, parent)
    try:
        result = subprocess.run([wolfram, '-file', str(script), str(output)],
                                cwd=root, env=env, check=False)
    except OSError as exc:
        raise AuditBlocked(f'cannot launch wolframscript: {exc}') from exc
    require(result.returncode == 0,
            f'native generator probe exited {result.returncode}; inspect {output}')
    native = read_json(output)
    validated = validate_data(native, proven)
    validated['native_output'] = str(output)
    validated['native_output_sha256'] = hashlib.sha256(output.read_bytes()).hexdigest()
    return validated
