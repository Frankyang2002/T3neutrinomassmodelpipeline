#!/usr/bin/env python3
"""Guarded local full-file replacement of RGBeta's heavy-doublet Dirac mass.

This utility reads the user's CURRENT T3RGBetaModel.wl, verifies the exact
known-bad mass definition, and writes a COMPLETE corrected replacement file.
It preserves all unrelated source text and makes a timestamped backup.

    python scripts/InstallT3RGBetaDoubletMassFix.py --self-test
    python scripts/InstallT3RGBetaDoubletMassFix.py --all
    python scripts/InstallT3RGBetaDoubletMassFix.py --verify

No GitHub writes or production beta-function execution are performed here.
"""
from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path
import re
import shutil
import sys

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'RGE' / 'running' / 'rgbeta' / 'T3RGBetaModel.wl'

# Anchor includes the distinct charged Dirac mass {F,Fc}; it cannot match
# the neutral Majorana branch {F,F}, which must remain unchanged.
OLD = '''        fermionInvariant = If[
            dF === 1,
            (1 &),
            (del[T3RGBetaRep[dF], #1, #2] &)
        ];

        AddFermionMass[
            MF,
            {F, Fc},'''
NEW = '''        (* F and Fc are both registered as unbarred SU(2) fundamentals
           when dF == 2. Their invariant Dirac mass needs epsilon, not
           delta. The singlet and real triplet cases retain their forms. *)
        fermionInvariant = Which[
            dF === 1, (1 &),
            dF === 2, (eps[T3RGBetaRep[dF], #1, #2] &),
            dF === 3, (del[T3RGBetaRep[dF], #1, #2] &),
            True, Return[$Failed]
        ];

        AddFermionMass[
            MF,
            {F, Fc},'''

NEUTRAL = '''    If[yF === 0,
        fermionInvariant = If[
            dF === 1,
            (1 &),
            (del[T3RGBetaRep[dF], #1, #2] &)
        ];

        AddFermionMass[
            MF,
            {F, F},'''


def require(ok: bool, detail: str) -> None:
    if not ok:
        raise ValueError(detail)


def normalize(text: str) -> str:
    return text.replace('\r\n', '\n')


def mass_section(text: str) -> tuple[int, int]:
    """Locate one complete function, so no other If/dF mass can be altered."""
    key = 'T3RGBetaAddMasses[\n'
    start = text.find(key)
    require(start >= 0 and text.find(key, start + 1) < 0,
            'Could not uniquely locate T3RGBetaAddMasses definition')
    end_marker = '\n(* ---------------------------------------------------------------------- *)'
    end = text.find(end_marker, start)
    require(end >= 0, 'Could not locate end of mass definition')
    return start, end


def check_structure(text: str, patched: bool) -> None:
    start, end = mass_section(text)
    mass = text[start:end]
    require(mass.count(NEUTRAL) == 1, 'Neutral Majorana branch changed')
    require(mass.count('{F, F}') == 1 and mass.count('{F, Fc}') == 1,
            'Expected one Majorana and one Dirac mass registration')
    require(mass.count('T3RGBetaAddScalarMasses[dS1, dS2]') == 1,
            'Scalar mass setup changed')
    require(mass.count(NEW) == (1 if patched else 0) and
            mass.count(OLD) == (0 if patched else 1),
            'Charged mass invariant differs from expected source')
    require(text.count('T3RGBetaAddMasses[dS1, dS2, dF, yF];') == 1,
            'T3RGBetaBuild no longer invokes mass builder as expected')


def transform(source: str) -> str:
    text = normalize(source)
    check_structure(text, patched=False)
    start, end = mass_section(text)
    before, mass, after = text[:start], text[start:end], text[end:]
    require(mass.count(OLD) == 1, 'Source changed: refusal to modify')
    patched = before + mass.replace(OLD, NEW, 1) + after
    check_structure(patched, patched=True)
    require(patched.replace(NEW, OLD, 1) == text,
            'Round-trip validation failed: other source bytes would change')
    if '\r\n' in source:
        return patched.replace('\n', '\r\n')
    return patched


def self_test() -> None:
    sample = '''T3RGBetaAddMasses[\n    dS1_Integer, dS2_Integer, dF_Integer, yF_\n] := Module[\n    {fermionInvariant},\n\n    T3RGBetaAddScalarMasses[dS1, dS2];\n\n''' + NEUTRAL + '''\n            GroupInvariant -> fermionInvariant\n        ],\n\n''' + OLD + '''\n            GroupInvariant -> fermionInvariant\n        ];\n    ];\n];\n(* ---------------------------------------------------------------------- *)\nT3RGBetaBuild[] := T3RGBetaAddMasses[dS1, dS2, dF, yF];\n'''
    fixed = transform(sample)
    check_structure(normalize(fixed), patched=True)
    require(fixed.count('eps[T3RGBetaRep[dF], #1, #2]') == 1,
            'Doublet epsilon absent')
    for bad in (sample.replace('{F, Fc}', '{F, F}'),
                sample.replace(OLD, OLD.replace('del[', 'eps[')),
                sample.replace(NEUTRAL, NEUTRAL.replace('del[', 'eps['))):
        try:
            transform(bad)
        except ValueError:
            pass
        else:
            raise AssertionError('Malformed source unexpectedly accepted')
    require(transform(sample.replace('\n', '\r\n')).count('\r\n') > 0,
            'CRLF not preserved')
    print('PASS: guarded single-branch replacement, neutral branch preserved, '
          'round trip, negative controls, CRLF preservation')


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    op = parser.add_mutually_exclusive_group(required=True)
    op.add_argument('--self-test', action='store_true')
    op.add_argument('--all', action='store_true', help='Write a complete replacement')
    op.add_argument('--verify', action='store_true', help='Check installed replacement')
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return 0
    require(SOURCE.is_file(), 'Missing local source: ' + str(SOURCE))
    contents = SOURCE.read_bytes()
    require(not contents.startswith(b'\xef\xbb\xbf'), 'Unexpected UTF-8 BOM')
    original = contents.decode('utf-8')
    if args.verify:
        check_structure(normalize(original), patched=True)
        print('PASS: installed doublet mass epsilon; neutral mass branch unchanged')
        return 0
    if NEW in original:
        check_structure(normalize(original), patched=True)
        print('Already fixed: ' + str(SOURCE))
        return 0
    result = transform(original)
    stamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    backup = SOURCE.with_name('T3RGBetaModel.before_doublet_mass_fix.' + stamp + '.wl')
    require(not backup.exists(), 'Backup exists: ' + str(backup))
    shutil.copy2(SOURCE, backup)
    try:
        SOURCE.write_bytes(result.encode('utf-8'))
        check_structure(normalize(SOURCE.read_text(encoding='utf-8')),
                        patched=True)
    except Exception:
        shutil.copy2(backup, SOURCE)
        raise
    print('PASS: created complete replacement: ' + str(SOURCE))
    print('Backup: ' + str(backup))
    print('Changed only charged dF=2 mass invariant: del -> eps')
    print('Requires local live RGBeta Ward test; UV beta outputs NOT recomputed')
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except (ValueError, OSError, UnicodeError) as exc:
        print('BLOCKED (source left unchanged): ' + str(exc), file=sys.stderr)
        sys.exit(2)
