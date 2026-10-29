"""Unit/regression checks for fail-closed, local-only docs consolidation."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from tools.ConsolidateHistoricalDocs import (
    INDEX, archive_and_apply, group_for, make_plan, rewrite_archived_doc_links,
    render_volume, write_report,
)


def git(root, *args):
    return subprocess.run(['git', '-C', str(root), *args],
                          check=True, capture_output=True)


@pytest.fixture
def repo(tmp_path):
    root=tmp_path/'repo'
    root.mkdir()
    (root/'docs').mkdir()
    (root/'scripts').mkdir()
    (root/'tests').mkdir()
    (root/'output').mkdir()
    (root/'docs/PHYSICS.md').write_text('# Physics\n',encoding='utf8')
    (root/'docs/RESULTS.md').write_text('# Results\n',encoding='utf8')
    (root/'docs/GENERAL_T3_TENSOR_ALIGNED_GEOMETRY.md').write_text(
        '# Aligned geometry\n\nExact finite geometry statement.\n',encoding='utf8')
    (root/'docs/GENERAL_T3_ACTION_IDENTITY.md').write_text(
        '# Action identity\n\n[Geometry](GENERAL_T3_TENSOR_ALIGNED_GEOMETRY.md)\n\n'
        '```python\n# [keep](GENERAL_T3_TENSOR_ALIGNED_GEOMETRY.md)\n```\n',encoding='utf8')
    (root/'docs/RENAMED_QUARTIC_PROVENANCE_COMPATIBILITY.md').write_text(
        '# Rename compatibility\n\nHistorical correction.\n',encoding='utf8')
    (root/'docs/general_c5_toy_benchmark.json').write_text('{"a": 1}\n',encoding='utf8')
    (root/'output/science.json').write_text('{"c5": 1}\n',encoding='utf8')
    git(root,'init','-q')
    git(root,'config','user.email','test@example.org')
    git(root,'config','user.name','test')
    git(root,'add','-A')
    git(root,'commit','-qm','synthetic fixtures')
    return root


def rows_by_name(plan):
    return {Path(r['path']).name:r for r in plan['rows']}


def test_topic_group_assignments():
    assert group_for('GENERAL_T3_TENSOR_ALIGNED_GEOMETRY.md')=='tensors'
    assert group_for('GENERAL_T3_FINAL_C5_AND_WEINBERG_RUNNING.md')=='running'
    assert group_for('T3_MATCHETE_NATIVE_H2_PORTAL_CG_FIX.md')=='matchete'
    assert group_for('T3_GENERAL_EFT1_HERMITICITY_SECOND_THRESHOLD.md')=='eft1'
    assert group_for('T3_BC_UV_RUNNING_MASS_HIERARCHY_BATCH14.md')=='fits'
    assert group_for('ONE_TIME_ARTIFACT_CLEANUP.md')=='maintenance'


def test_safe_fixture_plan_and_no_implicit_filesystem_mutation(repo):
    plan=make_plan(repo)
    assert plan['totals']['root_markdown']==5
    assert plan['totals']['eligible']==3
    assert plan['totals']['new_volumes']==2
    assert plan['totals']['projected_root_markdown']==5
    assert (repo/'docs/GENERAL_T3_ACTION_IDENTITY.md').exists()
    assert not (repo/'docs'/INDEX).exists()


def test_rewrite_links_not_source_code_or_external():
    s=('[local](GENERAL_T3_ACTION_IDENTITY.md) '
       '[external](https://example.org/GENERAL_T3_ACTION_IDENTITY.md)\n'
       '```sh\n[no](GENERAL_T3_ACTION_IDENTITY.md)\n```\n'
       '[target]: ./GENERAL_T3_ACTION_IDENTITY.md\n')
    t=rewrite_archived_doc_links(s,{'GENERAL_T3_ACTION_IDENTITY.md':'RESEARCH.md#action'})
    assert '[local](RESEARCH.md#action)' in t
    assert '[external](https://example.org/GENERAL_T3_ACTION_IDENTITY.md)' in t
    assert '```sh\n[no](GENERAL_T3_ACTION_IDENTITY.md)\n```' in t
    assert '[target]: RESEARCH.md#action' in t


def test_reference_in_non_moving_test_blocks(repo):
    (repo/'tests/test_marker.py').write_text(
        "DOC = 'docs/GENERAL_T3_ACTION_IDENTITY.md'\n",encoding='utf8')
    git(repo,'add','-A');git(repo,'commit','-qm','add reference')
    plan=make_plan(repo)
    assert rows_by_name(plan)['GENERAL_T3_ACTION_IDENTITY.md']['classification']=='RETAIN_REFERENCED'
    assert rows_by_name(plan)['GENERAL_T3_TENSOR_ALIGNED_GEOMETRY.md']['classification']=='RETAIN_REFERENCED'


def test_non_link_prose_reference_nonblocking(repo):
    (repo/'docs/PHYSICS.md').write_text(
        '# Physics\nThe old GENERAL_T3_ACTION_IDENTITY.md was a research note.\n')
    git(repo,'add','-A');git(repo,'commit','-qm','historical prose')
    plan=make_plan(repo)
    assert rows_by_name(plan)['GENERAL_T3_ACTION_IDENTITY.md']['classification']=='ARCHIVE_AND_CONSOLIDATE'


def test_active_markdown_link_blocks(repo):
    (repo/'docs/PHYSICS.md').write_text(
        '# Physics\n[Active](GENERAL_T3_ACTION_IDENTITY.md)\n')
    git(repo,'add','-A');git(repo,'commit','-qm','active link')
    plan=make_plan(repo)
    assert rows_by_name(plan)['GENERAL_T3_ACTION_IDENTITY.md']['classification']=='RETAIN_REFERENCED'


def test_dirty_historical_doc_retained(repo):
    (repo/'docs/GENERAL_T3_ACTION_IDENTITY.md').write_text('Local work not committed\n')
    plan=make_plan(repo)
    assert rows_by_name(plan)['GENERAL_T3_ACTION_IDENTITY.md']['classification']=='RETAIN_LOCALLY_MODIFIED'


def test_untracked_doc_remains_untouched(repo):
    (repo/'docs/TEMP_RESEARCH.md').write_text('Local scratch\n')
    plan=make_plan(repo)
    assert rows_by_name(plan)['TEMP_RESEARCH.md']['classification']=='RETAIN_UNTRACKED'


def test_external_archive_and_exact_original_sha(repo,tmp_path):
    plan=make_plan(repo)
    original=(repo/'docs/GENERAL_T3_ACTION_IDENTITY.md').read_bytes()
    output=tmp_path/'archive'
    result=archive_and_apply(repo,plan,output)
    assert result['archived']==3
    assert result['volumes']==2
    assert not (repo/'docs/GENERAL_T3_ACTION_IDENTITY.md').exists()
    assert (repo/'docs/PHYSICS.md').exists()
    assert (repo/'docs/general_c5_toy_benchmark.json').read_text()=='{"a": 1}\n'
    assert (repo/'output/science.json').read_text()=='{"c5": 1}\n'
    assert (output/'docs/GENERAL_T3_ACTION_IDENTITY.md').read_bytes()==original
    manifest=json.loads((output/'manifest.json').read_text())
    expected=hashlib.sha256(original).hexdigest()
    assert any(x['sha256']==expected for x in manifest['files'])
    volume=(repo/'docs/RESEARCH_GENERAL_T3_TENSORS_AND_QUARTICS.md').read_text()
    assert '[Geometry](RESEARCH_GENERAL_T3_TENSORS_AND_QUARTICS.md#source-general-t3-tensor-aligned-geometry-md)' in volume
    assert 'Exact finite geometry statement.' in volume
    assert (repo/'docs'/INDEX).is_file()


def test_archive_target_inside_repository_refused(repo):
    plan=make_plan(repo)
    with pytest.raises(RuntimeError,match='outside'):
        archive_and_apply(repo,plan,repo/'internal_archive')
    assert (repo/'docs/GENERAL_T3_ACTION_IDENTITY.md').exists()


def test_stale_plan_refused_before_mutation(repo,tmp_path):
    plan=make_plan(repo)
    (repo/'docs/PHYSICS.md').write_text('# changed\n')
    with pytest.raises(RuntimeError,match='state changed'):
        archive_and_apply(repo,plan,tmp_path/'archive')
    assert not (tmp_path/'archive').exists()
    assert (repo/'docs/GENERAL_T3_ACTION_IDENTITY.md').exists()


def test_report_is_readonly_for_docs(repo,tmp_path):
    plan=make_plan(repo)
    write_report(repo,plan,tmp_path/'report')
    assert (tmp_path/'report/summary.md').is_file()
    assert (tmp_path/'report/inventory.csv').is_file()
    assert (repo/'docs/GENERAL_T3_ACTION_IDENTITY.md').is_file()


def test_report_scan_includes_reports_python_code(repo):
    d=repo/'Reports';d.mkdir()
    (d/'BuildReports.py').write_text(
      'INFO="docs/GENERAL_T3_ACTION_IDENTITY.md"\n',encoding='utf8')
    git(repo,'add','-A');git(repo,'commit','-qm','report script reference')
    plan=make_plan(repo)
    assert rows_by_name(plan)['GENERAL_T3_ACTION_IDENTITY.md']['classification']=='RETAIN_REFERENCED'
