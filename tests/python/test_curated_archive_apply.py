"""Synthetic regression tests for selective, reversible local archive."""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

_SCRIPT_DIR = Path(__file__).resolve().parents[2] / 'scripts' / 'maintenance'
sys.path.insert(0, str(_SCRIPT_DIR))
import ApplyCuratedArchive as curated  # noqa: E402
import ExtremeRepositoryCleanup as core  # noqa: E402


def git(root: Path, *args: str) -> str:
    command = subprocess.run(['git', '-C', str(root), *args], check=True,
                             text=True, capture_output=True)
    return command.stdout.strip()


class CuratedArchiveTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name) / 'repo'
        (self.root / 'Lagrangian').mkdir(parents=True)
        (self.root / 'docs').mkdir()
        (self.root / 'pipeline.py').write_text('pass\n', encoding='utf-8')
        (self.root / 'README_OLD_EXAMPLE.md').write_text('Outdated installation note\n', encoding='utf-8')
        (self.root / 'docs' / 'TECHNICAL_EXAMPLE.md').write_text('Important derivation\n', encoding='utf-8')
        git(self.root, 'init', '-q')
        git(self.root, 'config', 'user.email', 'synthetic@example.invalid')
        git(self.root, 'config', 'user.name', 'Fixture')
        git(self.root, 'add', '.')
        git(self.root, 'commit', '-qm', 'Fixture')
        self.head = git(self.root, 'rev-parse', 'HEAD')
        inventory = core.build_plan(self.root)
        candidate = next(x for x in inventory['archive_candidates']
                         if x['path'] == 'README_OLD_EXAMPLE.md')
        self.curated = {
            'schema': curated.SCHEMA,
            'source_git_head': self.head,
            'options': {'include_nested_notes': False, 'include_presentations': False},
            'archive_now': [candidate],
            'deferred_technical_notes': [x for x in inventory['archive_candidates']
                                         if x['category'] == 'historical-document'],
            'held_unchanged': [],
        }
        self.selection = Path(self.tmp.name) / 'curated.json'
        self._write()

    def _write(self):
        self.selection.write_text(json.dumps(self.curated), encoding='utf-8')

    def test_selective_archive_and_restore(self):
        plan = curated.prepare(self.root, self.selection)
        self.assertEqual(plan['counts']['approved'], 1)
        self.assertEqual(plan['counts']['deferred'], 1)
        archive = core.apply_plan(self.root, plan)
        self.assertIsNotNone(archive)
        self.assertFalse((self.root / 'README_OLD_EXAMPLE.md').exists())
        self.assertTrue((self.root / 'docs' / 'TECHNICAL_EXAMPLE.md').exists())
        self.assertEqual(core.restore(self.root, archive), 1)
        self.assertTrue((self.root / 'README_OLD_EXAMPLE.md').exists())

    def test_changed_content_is_not_archived(self):
        (self.root / 'README_OLD_EXAMPLE.md').write_text('Modified locally\n', encoding='utf-8')
        with self.assertRaises(RuntimeError):
            curated.prepare(self.root, self.selection)

    def test_git_head_mismatch_fails_closed(self):
        self.curated['source_git_head'] = '0' * 40
        self._write()
        with self.assertRaisesRegex(RuntimeError, 'HEAD'):
            curated.prepare(self.root, self.selection)

    def test_no_document_under_docs_can_be_approved(self):
        live = core.build_plan(self.root)
        self.curated['archive_now'] = [next(x for x in live['archive_candidates']
                                             if x['path'] == 'docs/TECHNICAL_EXAMPLE.md')]
        self._write()
        with self.assertRaisesRegex(ValueError, 'Only root artifacts'):
            curated.prepare(self.root, self.selection)

    def test_path_traversal_rejected(self):
        item = dict(self.curated['archive_now'][0])
        item['path'] = '../outside.txt'
        self.curated['archive_now'] = [item]
        self._write()
        with self.assertRaisesRegex(ValueError, 'Invalid project-relative'):
            curated.prepare(self.root, self.selection)

    def test_apply_requires_explicit_confirmation(self):
        with self.assertRaises(SystemExit) as caught:
            curated.main(['--root', str(self.root), '--plan', str(self.selection), '--apply'])
        self.assertEqual(caught.exception.code, 2)
        self.assertTrue((self.root / 'README_OLD_EXAMPLE.md').exists())


if __name__ == '__main__':
    unittest.main()
