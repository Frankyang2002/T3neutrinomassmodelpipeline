#!/usr/bin/env python3
"""Static regression for the Matchete load-before-parse issue.

The preceding RunT3ExtendedC5Native.wl erroneously placed Needs["Matchete`"]
inside a Module containing Bar[Coupling[...]]. Wolfram parsed the entire
Module first, creating Global`Bar and Global`Coupling and failing Needs with
Bar::shdw and Coupling::shdw. The fixed driver loads Matchete in its OWN
completed top-level expression. This is structural, not a native Matchete test.

Project path: scripts/TestT3ExtendedMatcheteLoadBoundary.py
Run: python scripts/TestT3ExtendedMatcheteLoadBoundary.py
"""
from __future__ import annotations

from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
DRIVER = ROOT / 'Lagrangian/RunT3ExtendedC5Native.wl'


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def check_load_boundary(source: str) -> None:
    top = 'res = CheckAbort[Catch[Module['
    require(source.count(top) == 1, 'Native matching body missing/duplicated')
    prefix, body = source.split(top, 1)
    # Explicit declaration of the real load operation, before the matching body
    # is parsed. Do not merely search for the name in a comment.
    load = re.search(
        r'nativeMatcheteLoad\s*=\s*CheckAbort\[\s*'
        r'Check\[\s*UsingFrontEnd\[Needs\["Matchete`"\];\s*True\]\s*,\s*'
        r'\$Failed\s*\]\s*,\s*\$Aborted\s*\]\s*;',
        prefix,
        re.S,
    )
    require(load is not None, 'Matchete load must be separate and complete')
    require('NativeStep["load-Matchete"];' in prefix, 'Missing load checkpoint')
    require('If[nativeMatcheteLoad =!= True,' in prefix,
            'Failed load must block before matching starts')
    require('Exit[6]' in prefix, 'Failed load must exit nonzero')
    require('Needs["Matchete`"]' not in body,
            'Needs reintroduced inside the parsed matching body')
    require('HoldPattern[Bar[Coupling[y1Conj,___]]]' in body and
            'HoldPattern[Bar[Coupling[y2Conj,___]]]' in body,
            'Native extra-coupling rules are missing')
    require('RunT3Matching[built["LUV"],5,1]' in body and
            'ExtractWeinbergCoefficient[matched["MatchedEFT"]]' in body,
            'Original Matchete matching/extraction calls changed')
    require('"absolute_matchete_rgbeta_dictionary_certified"->False' in source,
            'Uncertified physical dictionary must remain explicit')


def check_delimiters(source: str) -> None:
    """Check Wolfram brackets while ignoring nested comments and strings."""
    stack: list[tuple[str, int]] = []
    pairs = {']': '[', '}': '{', ')': '('}
    pos, depth, string = 0, 0, False
    while pos < len(source):
        two = source[pos:pos + 2]
        if string:
            if source[pos] == '\\':
                pos += 2
                continue
            if source[pos] == '"':
                string = False
            pos += 1
            continue
        if two == '(*':
            depth += 1
            pos += 2
            continue
        if depth:
            if two == '*)':
                depth -= 1
                pos += 2
                continue
            pos += 1
            continue
        if source[pos] == '"':
            string = True
        elif source[pos] in '([{':
            stack.append((source[pos], pos))
        elif source[pos] in ')]}':
            require(bool(stack) and stack[-1][0] == pairs[source[pos]],
                    f'Unbalanced Wolfram delimiter at offset {pos}')
            stack.pop()
        pos += 1
    require(not stack and not string and depth == 0,
            'Unclosed Wolfram bracket, comment or string')


def main() -> int:
    text = DRIVER.read_text(encoding='utf-8-sig')
    check_load_boundary(text)
    check_delimiters(text)
    # Deliberately reintroduce the previous bug; the guard must reject it.
    bad = text.replace('res = CheckAbort[Catch[Module[',
                       'res = CheckAbort[Catch[Module[\n  '
                       'Needs["Matchete`"];', 1)
    try:
        check_load_boundary(bad)
    except AssertionError:
        pass
    else:
        raise AssertionError('Regression accepted in-Module Matchete load')
    print('PASS: separate top-level Matchete load precedes matching body')
    print('PASS: blocked native-load failure; unchanged matching calls; bracket syntax')
    print('PASS: negative control detects Needs reintroduced inside Module')
    print('STATIC ONLY: native Matchete execution must be run on Windows')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
