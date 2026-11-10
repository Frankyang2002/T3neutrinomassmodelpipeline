"""Isolated regression checks for the local-only documentation updater."""
from __future__ import annotations

import importlib.util
from pathlib import Path
import tempfile

SOURCE = Path(__file__).resolve().parents[1] / 'scripts' / 'UpdateT3DocumentationStatus.py'
spec = importlib.util.spec_from_file_location('t3_docs_updater', SOURCE)
updater = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(updater)


def test_idempotent_and_history_preserved():
    old = '# Numerical pipeline, fitting, diagnostics and reproducibility\n\n**Scope:** ordinary T3 neutrino-mass studies, with emphasis on saved results.\n\n## Results\n\n32/32 accepted historically.\n'
    a = updater.update_text(old, 'docs/NUMERICAL.md')
    b = updater.update_text(a, 'docs/NUMERICAL.md')
    assert a == b
    assert '32/32 accepted historically.' in b
    assert '**Historical scope:**' in b
    assert b.count(updater.BEGIN) == 1
    assert 'full-general, complete-one-loop' in b


def test_old_result_wording_not_mistaken_for_current():
    old = '# T3 numerical results: accepted benchmarks, comparisons and limits\n\n**Authoritative results-level summary for the honours project.**\n'
    a = updater.update_text(old, 'docs/RESULTS.md')
    assert a.startswith('# Historical T3 results: reduced and exploratory benchmarks')
    assert 'Historical results-level summary' in a
    assert '32/32' not in a


def test_bad_marker_rejected():
    try:
        updater.update_text('# Documentation\n\n<!-- T3_CURRENT_SCOPE_BEGIN: old -->', 'docs/README.md')
    except ValueError as exc:
        assert 'marker' in str(exc)
    else:
        raise AssertionError('Incomplete marker must be rejected')


def test_apply_then_check_all_targets_and_links():
    with tempfile.TemporaryDirectory() as folder:
        root = Path(folder)
        for rel in updater.TARGETS:
            path = root / rel
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(f'# {rel}\n\nHistorical content.\n', encoding='utf-8')
        assert updater.run(root, 'check') == 1
        assert updater.run(root, 'apply') == 0
        assert updater.run(root, 'check') == 0
        assert updater.run(root, 'apply') == 0
        root_readme = (root / 'README.md').read_text(encoding='utf-8')
        num = (root / 'Numerical/parallel_t3/wilson/README.md').read_text(encoding='utf-8')
        assert '(docs/IMPLEMENTATION_STATUS.md)' in root_readme
        assert '(../../../docs/IMPLEMENTATION_STATUS.md)' in num
