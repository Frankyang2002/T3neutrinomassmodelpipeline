#!/usr/bin/env python3
"""Regression for Matchete's unsupported NFlavor=1 native matching pilot.

Project path: scripts/TestT3ExtendedFlavorIndexFix.py
This is an executable Python/mock test; live Matchete matching is not simulated
as physical evidence.
"""
from __future__ import annotations

import importlib.util
import json
import re
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
WL = ROOT / 'Lagrangian/RunT3ExtendedC5Native.wl'
RUNNER = ROOT / 'scripts/RunT3ExtendedMatchingBatch.py'
SYNTH = ROOT / 'scripts/TestT3ExtendedPhysicalStackSynthetic.py'


def require(ok: bool, why: str) -> None:
    if not ok:
        raise AssertionError(why)


def syntax_balance(source: str) -> None:
    """WL bracket check with nested comments and literal strings skipped."""
    stack: list[str] = []
    pairs = {')': '(', ']': '[', '}': '{'}
    i, comment, in_string = 0, 0, False
    while i < len(source):
        pair = source[i:i + 2]
        if in_string:
            if source[i] == '\\':
                i += 2
                continue
            if source[i] == '"':
                in_string = False
        elif pair == '(*':
            comment += 1
            i += 2
            continue
        elif comment:
            if pair == '*)':
                comment -= 1
                i += 2
                continue
        elif source[i] == '"':
            in_string = True
        elif source[i] in '([{':
            stack.append(source[i])
        elif source[i] in ')]}':
            require(bool(stack) and stack.pop() == pairs[source[i]],
                    'Unbalanced Wolfram bracket')
        i += 1
    require(not stack and not comment and not in_string,
            'Unterminated WL bracket, comment or string')


def test_source_contracts() -> None:
    wl = WL.read_text(encoding='utf-8-sig')
    python = RUNNER.read_text(encoding='utf-8-sig')
    syntax_balance(wl)
    require(re.search(r'nHeavy\s*=\s*If\[\s*MemberQ\[\s*\{"2",\s*"3"\}', wl) is not None,
            'Wolfram does not restrict native NFlavor to {2,3}')
    require('T3ModelFromClass[cl,allowed[key],nHeavy]' in wl,
            'Heavy flavour count not passed to model')
    require('Needs["Matchete`"]' in wl and
            wl.index('nativeMatcheteLoad = CheckAbort[') < wl.index('res = CheckAbort[Catch[Module['),
            'Matchete load boundary regressed')
    require('RunT3Matching[built["LUV"],5,1]' in wl and
            'ExtractWeinbergCoefficient[matched["MatchedEFT"]]' in wl,
            'Matching/C5 extraction changed unexpectedly')
    require('choices=(2,3),default=3' in python,
            'CLI should default to physical 3 heavy flavors')
    require('need(heavy_flavours in (2,3)' in python,
            'Python must reject invalid imported/native invocation')
    require('absolute_matchete_rgbeta_dictionary_certified"->False' in wl,
            'Physical dictionary must not be falsely certified')
    # Confirm upstream numerical state and Matchete input source contract
    fields = ROOT / 'Lagrangian/model/T3Fields.wl'
    if fields.is_file():
        fs = fields.read_text(encoding='utf-8-sig')
        require('DefineFlavorIndex[' in fs and 'Lookup[f, "Multiplicity", 3]' in fs,
                'T3 native heavy-flavor registration has changed')
    require('heavy_flavours=1' not in SYNTH.read_text(encoding='utf-8-sig'),
            'Synthetic test still exercises invalid NFlavor=1')


def test_mocked_native_invocations() -> None:
    spec = importlib.util.spec_from_file_location('nativebatchfix', RUNNER)
    require(spec is not None and spec.loader is not None, 'Could not import batch runner')
    batch = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(batch)
    batch.self_test()
    with TemporaryDirectory() as folder:
        root = Path(folder)
        # A true native process cannot be run here; verify invocation arguments
        # and reporting schema without pretending a result is Matchete evidence.
        def fake_run(cmd, **kwargs):
            heavy = int(cmd[-1])
            require(heavy in (2, 3), 'Invalid requested heavy-flavour dimension')
            require(cmd[1] == '-file', 'Changed command-line interface')
            require(kwargs['env']['T3_EXPERIMENTAL_FLAVOR_MASSES'] == '1',
                    'Flavour-indexed mass switch lost')
            out = Path(cmd[3])
            out.parent.mkdir(parents=True, exist_ok=True)
            key = f'{cmd[-3]}_{cmd[-2]}'
            datum = {
                'format': 't3-extended-c5-native-matchete-one-step-v1',
                'model': key, 'heavy_flavours': heavy,
                'status': batch.SUCCESS,
                'additional_yukawa': batch.SPECS[key][2],
                'one_loop_matching_performed': True,
                'native_matchete_lagrangian_check_passed': True,
                'matched_c5_in_matchete_native_coordinates': True,
                'absolute_matchete_rgbeta_dictionary_certified': False,
                'hierarchical_matching_certified': False,
                'cross_package_numeric_substitution_authorized': False,
                'physical_neutrino_mass_prediction': False,
            }
            for k in ('native_coefficient_full_path',
                      'native_coefficient_zero_extra_path',
                      'native_coefficient_extra_difference_path'):
                path = out.parent / f'{out.stem}.{k}.txt'
                path.write_text('synthetic C5', encoding='utf-8')
                datum[k] = str(path)
            out.write_text(json.dumps(datum), encoding='utf-8')
            return SimpleNamespace(returncode=0, stdout='MOCK RESULT ONLY\n')

        with patch.object(batch.subprocess, 'run', side_effect=fake_run):
            for flavors in (2, 3):
                for model in batch.SPECS:
                    result = batch.run_one(model, heavy_flavours=flavors,
                                           wl='fake-wolfram', timeout=1,
                                           root=root, output=root / 'out')
                    require(result['status'] == batch.SUCCESS,
                            'Mocked native schema rejected')
                    require(result['heavy_flavours'] == flavors,
                            'Heavy-flavor dimension lost')
            try:
                batch.run_one('A_p0', heavy_flavours=1,
                              wl='fake-wolfram', timeout=1,
                              root=root, output=root / 'out')
            except ValueError:
                pass
            else:
                raise AssertionError('Unsupported Matchete NFlavor=1 accepted')


def main() -> None:
    test_source_contracts()
    test_mocked_native_invocations()
    print('PASS: Wolfram syntax, native load boundary and physical matching source guards')
    print('PASS: NFlavor=1 rejected; all six models mocked at NFlavor=2 and 3')
    print('PASS: matching reports remain explicitly unnormalised and nonhierarchical')
    print('STATIC/MOCK ONLY: Matchete must be run locally for physical evidence')


if __name__ == '__main__':
    main()
