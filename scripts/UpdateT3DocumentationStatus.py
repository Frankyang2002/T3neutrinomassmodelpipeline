#!/usr/bin/env python3
"""Update current-scope notices in old T3 project notes without erasing history.

Complete standalone file. Run from the project root:
  python scripts/UpdateT3DocumentationStatus.py --apply
  python scripts/UpdateT3DocumentationStatus.py --check

The updater preserves the contents and historical tables of existing documents,
changing only named stale scope statements and one versioned status block near
its title. It does not modify physics code, existing numerical results, or GitHub.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import re
import sys

BEGIN = '<!-- T3_CURRENT_SCOPE_BEGIN: 2026-11-10 -->'
END = '<!-- T3_CURRENT_SCOPE_END -->'
BASE_DIR = Path(__file__).resolve().parents[1]

COMMON = (
    '**Current implementation status (10 November 2026):** '
    'The requested A(-4)/D(-2) **full-general, complete-one-loop, '
    'fermion-first -> joint-scalar -> SMEFT Yukawa fit is not yet an '
    'executable end-to-end numerical workflow**. The repository already '
    'implements general Matchete/RGBeta components and an actual symbolic '
    '`RunSequentialT3Matching` stage engine, but its general Wilson matching '
    'outputs are not fully connected to the two-model numerical fit. '
    'The current thesis Sobol pilot is **restricted/conditional** '
    '(analytic fermion triangle and selected scalar bubble), even if it reaches '
    'a low chi-square. See [implementation status]({status}), '
    '[A(-4)/D(-2) feasibility]({feas}) and '
    '[implementation batches]({batches}). '
    'Earlier numerical results and verification claims below refer to '
    'their **dated, historical configurations** unless an individual run '
    'manifest establishes otherwise.'
)

# Documents kept for their old numerical data, technical conventions, or history.
# No complete file is replaced from the default GitHub version by this script.
TARGETS = {
    'README.md': 'Project overview; current A(-4)/D(-2) task is not finished.',
    'README_T3_FINAL_REPRO_BATCH8.md': 'Batch 8 is a historical conditional report, not current full-general output.',
    'docs/README.md': 'Read the current implementation-status and batch plan before old studies.',
    'docs/NUMERICAL.md': 'Historical reduced/frozen numerical methods; not the new full-general evaluator.',
    'docs/RESULTS.md': 'Historical reduced/conditional results, not newly fitted general-model values.',
    'docs/SCALAR_SECTOR.md': 'Detailed reduced-potential analysis, not a statement that all general couplings are numerically matched.',
    'docs/MODEL_AND_PIPELINE.md': 'The required EFT sequence is a physics contract; its full-general implementation remains in progress.',
    'docs/PROJECT_GUIDE.md': 'The cited production smoke test is historical and does not demonstrate full-general one-loop matching.',
    'docs/T3_THESIS_PLAN.md': 'Earlier thesis planning does not supersede the current code-audit status.',
    'docs/T3_COMPLETE_GENERAL_SCALAR_SYMBOLIC.md': 'General symbolic exports exist; general numerical two-threshold Wilson matching is not implied.',
    'docs/T3_FULL_SCALAR_CONTINUOUS_BATCH.md': 'Full scalar renormalisable RGE running does not imply full Wilson matching and SMEFT.',
    'Numerical/parallel_t3/wilson/README.md': 'Current pilot uses the candidate Wilson ADM; no complete general one-loop matching input is connected.',
}

# Safe, exact text corrections; if a locally modified document differs, its
# unchanged contents remain intact and the visible scope banner still applies.
REPLACEMENTS = {
    'docs/NUMERICAL.md': {
        '**Scope:** ordinary T3 neutrino-mass studies, with emphasis':
            '**Historical scope:** ordinary T3 neutrino-mass studies, with emphasis',
        '**Source status (29 October 2026):**':
            '**Historical source status (29 October 2026):**',
    },
    'docs/RESULTS.md': {
        '# T3 numerical results: accepted benchmarks, comparisons and limits':
            '# Historical T3 results: reduced and exploratory benchmarks',
        '**Authoritative results-level summary for the honours project.**':
            '**Historical results-level summary for the reduced-scalar and exploratory campaigns.**',
        '## 3. Final reduced-scalar frozen comparison — primary thesis result':
            '## 3. Historical reduced-scalar frozen comparison',
        'The primary honours result is the **restricted leading-one-loop, zero-additional-scalar-boundary comparison**.':
            'The historical reduced-scalar study result is the **restricted leading-one-loop, zero-additional-scalar-boundary comparison**.',
    },
    'docs/SCALAR_SECTOR.md': {
        '**Authoritative scalar-sector explanation for the honours thesis.**':
            '**Historical scalar-sector analysis of the reduced-potential approximation.**',
        '### 6.2 Final reduced-potential frozen-UV fits (primary thesis results)':
            '### 6.2 Historical reduced-potential frozen-UV fits',
    },
    'docs/PROJECT_GUIDE.md': {
        '## Verified production workflow (one-point smoke test, 2026-11-01)':
            '## Historical operational smoke test (one point, 2026-11-01)',
    },
}


def status_notice(path: str, description: str) -> str:
    prefix = '../../../docs/' if path.startswith('Numerical/') else ''
    if path in ('README.md', 'README_T3_FINAL_REPRO_BATCH8.md'):
        prefix = 'docs/'
    # docs/*.md uses same-directory paths.
    return (
        f'{BEGIN}\n\n'
        f'**Document scope:** {description}\n\n'
        + COMMON.format(
            status=prefix + 'IMPLEMENTATION_STATUS.md',
            feas=prefix + 'T3_A4_D2_GENERAL_PIPELINE_FEASIBILITY_2026_11_10.md',
            batches=prefix + 'T3_A4_D2_FULL_GENERAL_IMPLEMENTATION_BATCHES.md',
        )
        + f'\n\n{END}\n'
    )


def update_text(text: str, path: str) -> str:
    notice = status_notice(path, TARGETS[path])
    previous = re.compile(re.escape(BEGIN) + r'.*?' + re.escape(END) + r'\n?', re.S)
    if previous.search(text):
        text = previous.sub('', text)
    else:
        # Old-version banner marker? Fail rather than risk duplicate ambiguous status.
        if 'T3_CURRENT_SCOPE_BEGIN:' in text or END in text:
            raise ValueError(f'{path}: unknown/partial current-scope marker')
    for source, replacement in REPLACEMENTS.get(path, {}).items():
        if source in text:
            text = text.replace(source, replacement)
    # Preserve original Markdown beyond this first-title insertion.
    title = re.search(r'^# [^\n]+\n', text, re.M)
    if not title:
        raise ValueError(f'{path}: missing Markdown H1 title')
    return text[:title.end()] + '\n' + notice + text[title.end():].lstrip('\n')


def run(root: Path, mode: str) -> int:
    failed = False
    changed = 0
    for rel in TARGETS:
        path = root / rel
        if not path.is_file():
            print(f'MISSING {rel}', flush=True)
            failed = True
            continue
        raw = path.read_bytes()
        newline = '\r\n' if b'\r\n' in raw else '\n'
        try:
            current = raw.decode('utf-8-sig').replace('\r\n', '\n')
            updated = update_text(current, rel)
        except (UnicodeError, ValueError) as exc:
            print(f'ERROR {rel}: {exc}', flush=True)
            failed = True
            continue
        expected = updated.replace('\n', newline).encode('utf-8')
        if raw.startswith(b'\xef\xbb\xbf'):
            expected = b'\xef\xbb\xbf' + expected
        if expected != raw:
            changed += 1
            if mode == 'apply':
                path.write_bytes(expected)
                print(f'UPDATED {rel}', flush=True)
            else:
                print(f'OUTDATED {rel}', flush=True)
        elif mode == 'apply':
            print(f'UNCHANGED {rel}', flush=True)
    if failed:
        return 2
    if mode == 'check':
        if changed:
            print(f'DOCUMENTATION OUTDATED: {changed} files; run --apply')
            return 1
        print(f'DOCUMENTATION CURRENT: {len(TARGETS)} page-level notices')
    else:
        print(f'DOCUMENTATION UPDATED: {changed} / {len(TARGETS)} documents')
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    choice = ap.add_mutually_exclusive_group(required=True)
    choice.add_argument('--apply', action='store_true')
    choice.add_argument('--check', action='store_true')
    ap.add_argument('--root', type=Path, default=BASE_DIR,
                    help='Project root (defaults to parent of installed scripts directory)')
    args = ap.parse_args()
    return run(args.root.resolve(), 'apply' if args.apply else 'check')


if __name__ == '__main__':
    sys.exit(main())
