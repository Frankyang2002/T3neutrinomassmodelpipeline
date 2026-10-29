"""Synthetic, filesystem-isolated regression tests for the final docs consolidation."""
from __future__ import annotations

import importlib.util
from pathlib import Path
import subprocess

import pytest

SOURCE = Path(__file__).resolve().parents[2] / 'tools' / 'FinishHistoricalDocsConsolidation.py'
spec = importlib.util.spec_from_file_location('FinishHistoricalDocsConsolidation', SOURCE)
assert spec and spec.loader
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def init_repo(tmp_path: Path) -> Path:
    root = tmp_path / 'repo'
    (root / 'docs').mkdir(parents=True)
    for v in m.VOLUMES.values():
        (root / 'docs' / v).write_text('# Existing volume\n\n## Contents\n\n- [Existing](#existing)\n\n---\n\n<a id="existing"></a>\n\nPreserved existing record.\n', encoding='utf-8')
    for n in m.CANDIDATES:
        (root / 'docs' / n).write_text('# Research note: '+n+'\n\n'
            'Reproducibility detail, with preserved equation $x=1$.\n', encoding='utf-8')
    index = '# Historical index\n\n## Original-file lookup\n\n| Original path | New location | SHA-256 |\n|---|---|---|\n| `docs/old.md` | [old](old.md) | `abc` |\n\n## Retained individual documents\n\n| File | Reason |\n|---|---|\n'
    index += ''.join(f'| [{n}]({n}) | `RETAIN_REFERENCED` |\n' for n in m.CANDIDATES)
    (root / 'docs' / m.INDEX).write_text(index, encoding='utf-8')
    (root / 'docs' / 'PHYSICS.md').write_text('# Core Physics\n\n[Research note](GENERAL_T3_TENSOR_ALIGNED_GEOMETRY.md)\n', encoding='utf-8')
    (root / 'docs' / 'general_t3_threshold_benchmark.json').write_text('{"physics": 1}\n', encoding='utf-8')
    subprocess.run(['git', 'init', '-q', str(root)], check=True)
    subprocess.run(['git', '-C', str(root), 'config', 'user.email', 'test@example.org'], check=True)
    subprocess.run(['git', '-C', str(root), 'config', 'user.name', 'Synthetic Test'], check=True)
    subprocess.run(['git', '-C', str(root), 'add', '.'], check=True)
    subprocess.run(['git', '-C', str(root), 'commit', '-qm', 'baseline'], check=True)
    return root


def test_candidate_inventory_exact():
    assert len(m.RENAME_ONLY) == 23
    assert len(m.AUDIT_ONLY) == 4
    assert len(m.CANDIDATES) == 27
    assert len(set(m.CANDIDATES)) == 27


def test_protected_live_dependencies_excluded():
    assert 'HIGH_SCALAR_HALF_FERMION_START.md' not in m.CANDIDATES
    assert 'T3_REDUCED_SCALAR_POTENTIAL_CONCLUSION.md' not in m.CANDIDATES
    assert 'TEST_BASELINE_FIX.md' not in m.CANDIDATES


def test_volume_assignment_consistent():
    assert m.volume_for('GENERAL_T3_TENSOR_ALIGNED_GEOMETRY.md') == m.VOLUMES['tensors']
    assert m.volume_for('GENERAL_T3_EFT1_SCALAR_PAIR_MIXING.md') == m.VOLUMES['eft1']
    assert m.volume_for('ONE_TIME_ARTIFACT_CLEANUP.md') == m.VOLUMES['maintenance']
    assert m.volume_for('T3_ZERO_BOUNDARY_LOW_ENERGY_BATCH8.md') == m.VOLUMES['fits']


def test_historical_mapping_literal_is_accepted_but_active_python_reference_rejected():
    line = "'GENERAL_T3_STAGE31_16_MODEL_UV_RUNNING.md': 'GENERAL_T3_16_MODEL_UV_RUNNING.md',"
    assert m.static_reference_is_historical('tools/RenameNumberedResearchFiles.py',
        'GENERAL_T3_16_MODEL_UV_RUNNING.md', line)
    assert not m.static_reference_is_historical('Numerical/example.py',
        'GENERAL_T3_16_MODEL_UV_RUNNING.md', 'DOC = "docs/GENERAL_T3_16_MODEL_UV_RUNNING.md"')
    assert not m.static_reference_is_historical('tools/RenameNumberedResearchFiles.py',
        'GENERAL_T3_16_MODEL_UV_RUNNING.md', 'open("docs/GENERAL_T3_16_MODEL_UV_RUNNING.md")')


def test_rewrite_markdown_links_preserves_fenced_historical_commands():
    name = 'GENERAL_T3_TENSOR_ALIGNED_GEOMETRY.md'
    text = f'[Active](docs/{name})\nPlain historical mention: {name}\n```powershell\ncat docs/{name}\n```\n'
    changed = m.rewrite_links(text, set(m.CANDIDATES))
    assert '[Active](' + m.target_for(name) + ')' in changed
    assert 'Plain historical mention: ' + name in changed
    assert 'cat docs/' + name in changed


def test_rewrite_external_links_unchanged():
    name = 'GENERAL_T3_TENSOR_ALIGNED_GEOMETRY.md'
    text = f'[External](https://example.org/{name})\n'
    assert m.rewrite_links(text, set(m.CANDIDATES)) == text


def test_index_reclassification_removes_old_retained_row():
    name = 'GENERAL_T3_TENSOR_ALIGNED_GEOMETRY.md'
    text = '# Index\n\n## Original-file lookup\n| A | B | C |\n|---|---|---|\n\n## Retained individual documents\n\n| File | Reason |\n|---|---|\n| [' + name + '](' + name + ') | `RETAIN_REFERENCED` |\n'
    updated = m.index_update(text, {name: {'volume':m.volume_for(name),'sha256':'ab'*32}})
    assert f'| [{name}]({name}) |' not in updated
    assert f'| `docs/{name}` |' in updated
    assert '## Retained individual documents' in updated


def test_index_missing_expected_row_fails_closed():
    with pytest.raises(RuntimeError, match='Retained index entry'):
        m.index_update('## Original-file lookup\n\n## Retained individual documents\n',
                       {'GENERAL_T3_TENSOR_ALIGNED_GEOMETRY.md':{'volume':'v','sha256':'0'}})


def test_append_maintains_existing_document_text():
    n = 'GENERAL_T3_TENSOR_ALIGNED_GEOMETRY.md'
    old = '# Volume\n\n## Contents\n\n- [Old](#old)\n\n---\n\nOld scientific data.\n'
    new = m.append_records(old, {n:b'# Title\n\nOriginal scientific equation $x=1$.\n'})
    assert 'Old scientific data.' in new
    assert 'Original scientific equation $x=1$.' in new
    assert '#'+m.anchor(n) in new


def test_make_plan_and_archive_apply(tmp_path):
    root = init_repo(tmp_path)
    plan = m.make_plan(root, expected_head=None)
    assert plan['ready']
    assert plan['candidate_count'] == 27
    assert plan['root_markdown_after'] == plan['root_markdown_before'] - 27
    archived = tmp_path / 'verified_originals'
    out = m.apply(root, plan, archived)
    assert out['originals_archived'] == 27
    assert not (root/'docs'/m.CANDIDATES[0]).exists()
    assert (archived/'originals/docs'/m.CANDIDATES[0]).is_file()
    assert (archived/'before/docs'/m.INDEX).is_file()
    assert (root/'docs'/'PHYSICS.md').read_text().count('RESEARCH_GENERAL_T3_TENSORS_AND_QUARTICS.md') == 1
    assert (root/'docs'/'general_t3_threshold_benchmark.json').read_text() == '{"physics": 1}\n'
    assert (root/'docs'/m.INDEX).read_text().count('docs/GENERAL_T3_TENSOR_ALIGNED_GEOMETRY.md') >= 1
    assert (root/'docs'/m.volume_for(m.CANDIDATES[0])).read_text().find(m.anchor(m.CANDIDATES[0])) != -1


def test_modified_target_is_blocked(tmp_path):
    root = init_repo(tmp_path)
    (root / 'docs' / m.CANDIDATES[0]).write_text('# Locally modified\n')
    with pytest.raises(RuntimeError, match='Dirty required path'):
        m.make_plan(root, expected_head=None)


def test_unknown_runtime_reference_blocks_apply(tmp_path):
    root = init_repo(tmp_path)
    p = root/'Numerical'/'consumer.py'
    p.parent.mkdir(parents=True)
    p.write_text('DOC = "docs/GENERAL_T3_TENSOR_ALIGNED_GEOMETRY.md"\n')
    plan = m.make_plan(root, expected_head=None)
    assert not plan['ready']
    assert any(x['source']=='Numerical/consumer.py' for x in plan['reference_audit']['blocking_external_references'])
    with pytest.raises(RuntimeError, match='Unknown references'):
        m.apply(root, plan, tmp_path/'archive')


def test_head_mismatch_is_rejected(tmp_path):
    root = init_repo(tmp_path)
    with pytest.raises(RuntimeError, match='Unexpected HEAD'):
        m.make_plan(root, expected_head='0000000000000000')


def test_prior_consolidation_fixture_literals_are_nonblocking(tmp_path):
    """The seven known references are only synthetic test literals."""
    root = init_repo(tmp_path)
    previous_test = root / m.PREVIOUS_CONSOLIDATION_FIXTURE_PATH
    previous_test.parent.mkdir(parents=True, exist_ok=True)
    previous_test.write_text('\n'.join(sorted(m.PREVIOUS_CONSOLIDATION_FIXTURE_LINES)) + '\n', encoding='utf-8')
    plan = m.make_plan(root, expected_head=None)
    assert plan['ready'] is True
    assert plan['reference_audit']['blocking_external_references'] == []
    seen = [r for r in plan['reference_audit']['historical_literal_references']
            if r['source'] == m.PREVIOUS_CONSOLIDATION_FIXTURE_PATH]
    assert len(seen) == 7
    assert {r['name'] for r in seen} == {
        'GENERAL_T3_TENSOR_ALIGNED_GEOMETRY.md',
        'GENERAL_T3_FINAL_C5_AND_WEINBERG_RUNNING.md',
        'ONE_TIME_ARTIFACT_CLEANUP.md',
    }


def test_previous_fixture_path_is_not_globally_exempt(tmp_path):
    root = init_repo(tmp_path)
    previous_test = root / m.PREVIOUS_CONSOLIDATION_FIXTURE_PATH
    previous_test.parent.mkdir(parents=True, exist_ok=True)
    previous_test.write_text('\n'.join(sorted(m.PREVIOUS_CONSOLIDATION_FIXTURE_LINES)) +
        '\nPath("docs/GENERAL_T3_TENSOR_ALIGNED_GEOMETRY.md").read_text()\n', encoding='utf-8')
    plan = m.make_plan(root, expected_head=None)
    assert plan['ready'] is False
    assert len(plan['reference_audit']['blocking_external_references']) == 1
    assert plan['reference_audit']['blocking_external_references'][0]['source'] == \
        m.PREVIOUS_CONSOLIDATION_FIXTURE_PATH


def test_changed_historical_fixture_line_is_blocked():
    path = m.PREVIOUS_CONSOLIDATION_FIXTURE_PATH
    name = 'GENERAL_T3_TENSOR_ALIGNED_GEOMETRY.md'
    assert m.static_reference_is_historical(path, name,
        "(root/'docs/GENERAL_T3_TENSOR_ALIGNED_GEOMETRY.md').write_text(")
    assert not m.static_reference_is_historical(path, name,
        "(root/'docs/GENERAL_T3_TENSOR_ALIGNED_GEOMETRY.md').read_text(")
    assert not m.static_reference_is_historical('Numerical/some_module.py', name,
        "(root/'docs/GENERAL_T3_TENSOR_ALIGNED_GEOMETRY.md').write_text(")
