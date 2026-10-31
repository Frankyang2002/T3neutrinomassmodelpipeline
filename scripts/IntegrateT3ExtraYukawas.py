#!/usr/bin/env python3
"""Guarded local full-file replacement for six opt-in T3 RGBeta Yukawas.

The original local T3RGBetaModel.wl is read, structurally audited, and written
back as ONE COMPLETE replacement (never a patch/diff) with a timestamped backup.
It MUST already contain the independently Ward-checked charged doublet MF fix.
The extension is opt-in via T3_EXPERIMENTAL_GENERAL_FERMIONS=1 in a process.
No GitHub access, UV RGEs, or numerical results are modified by this script.

python scripts/IntegrateT3ExtraYukawas.py --self-test
python scripts/IntegrateT3ExtraYukawas.py --install
python scripts/IntegrateT3ExtraYukawas.py --verify
"""
from __future__ import annotations

import argparse
from datetime import datetime
import os
from pathlib import Path
import shutil
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
MODEL = ROOT / 'RGE/running/rgbeta/T3RGBetaModel.wl'

MASS_CORRECT = '''        fermionInvariant = Which[
            dF === 1, (1 &),
            dF === 2, (eps[T3RGBetaRep[dF], #1, #2] &),
            dF === 3, (del[T3RGBetaRep[dF], #1, #2] &),
            True, Return[$Failed]
        ];

        AddFermionMass[
            MF,
            {F, Fc},'''
MASS_OLD = '''        fermionInvariant = If[
            dF === 1,
            (1 &),
            (del[T3RGBetaRep[dF], #1, #2] &)
        ];

        AddFermionMass[
            MF,
            {F, Fc},'''

BASE_FUNCTION_END = '''    ];
];

(* We add scalar and fermion mass term *)'''
# Exact old production build anchor.
BUILD_ANCHOR = '''    T3RGBetaAddT3Yukawas[dS1, dS2, dF, yF];
    T3RGBetaAddMasses[dS1, dS2, dF, yF];'''
BUILD_NEW = '''    T3RGBetaAddT3Yukawas[dS1, dS2, dF, yF];
    (* Extension is opt-in; disabled mode produces the original five Yukawas. *)
    If[T3RGBetaAddGeneralFermionYukawas[dS1, dS2, dF, yF] === $Failed,
        Return[$Failed]
    ];
    T3RGBetaAddMasses[dS1, dS2, dF, yF];'''

# Place this immediately before the mass definitions; its GaugeRep and Fields
# choices are identical to those verified in the isolated RGBeta registration.
EXTRA_DEFINITION = '''
(* Opt-in Workstream A general-fermion Yukawas.  The original restricted
   topology is preserved unless T3_EXPERIMENTAL_GENERAL_FERMIONS=1.
   CheckInvariance is deliberately not invoked: RGBeta v1.2.0's symbolic
   diagnostic is not a certified Ward check for these component tensors.
   Independent exact component Ward evidence lives in Workstream A reports.
   This DOES NOT establish a physical Matchete/RGBeta coupling dictionary. *)
T3RGBetaAddGeneralFermionYukawas[
    dS1_Integer, dS2_Integer, dF_Integer, yF_
] := Module[{which, extra, fields, invariant, scalarY, fermionY,
             expectedName, countBefore, registered},
    If[Environment["T3_EXPERIMENTAL_GENERAL_FERMIONS"] =!= "1",
        Return[{}]
    ];
    (* Hypercharges of S1 and S2 are yF-1/2 and yF+1/2. *)
    which = Which[
        TrueQ[yF === 1/2], 1,
        TrueQ[yF === -1/2], 2,
        True, 0
    ];
    If[which === 0, Return[{}]];
    If[dF =!= 2 || !MemberQ[{1, 3}, If[which === 1, dS1, dS2]],
        Print["Invalid general-fermion Yukawa representation: ",
            {dS1, dS2, dF, yF}];
        Return[$Failed]
    ];
    invariant = T3RGBetaYukawaInvariant[If[which === 1, dS1, dS2], dF];
    If[which === 1,
        extra = y1Conj;
        fields = {S1, l, F};
        scalarY = yF - 1/2;
        fermionY = yF,
        extra = y2Conj;
        fields = {Bar[S2], l, Fc};
        scalarY = -(yF + 1/2);
        fermionY = -yF
    ];
    If[!TrueQ[FullSimplify[-1/2 + scalarY + fermionY == 0]],
        Print["General Yukawa U1Y test failed"];
        Return[$Failed]
    ];
    expectedName = If[which === 1, y1, y2];
    If[!AssociationQ[$yukawas] || Length[$yukawas] =!= 5 ||
       !KeyExistsQ[$yukawas, expectedName] || KeyExistsQ[$yukawas, extra],
        Print["Unexpected RGBeta Yukawa registry before optional extension"];
        Return[$Failed]
    ];
    countBefore = Length[$yukawas];
    registered = CheckAbort[Check[
        AddYukawa[extra, fields,
            GroupInvariant -> invariant,
            CouplingIndices -> ({gen[#2], heavy[#3]} &),
            Chirality -> Right]; True,
        $Failed], $Aborted];
    If[registered =!= True || Length[$yukawas] =!= countBefore + 1 ||
       !KeyExistsQ[$yukawas, extra],
        Print["Extra RGBeta Yukawa registration failed"];
        Return[$Failed]
    ];
    {extra}
];
'''
BETAS_END = '''    result
];

(* ---------------------------------------------------------------------- *)
(* If shared scalars we use the following instead'''
BETAS_NEW = '''    (* These one-loop beta functions must be included in the same UV system
       if their independent flavour matrices have been registered. *)
    If[KeyExistsQ[$yukawas, y1Conj],
        AssociateTo[result, "y1Conj" -> Quiet[BetaTerm[y1Conj, 1]]]
    ];
    If[KeyExistsQ[$yukawas, y2Conj],
        AssociateTo[result, "y2Conj" -> Quiet[BetaTerm[y2Conj, 1]]]
    ];
    result
];

(* ---------------------------------------------------------------------- *)
(* If shared scalars we use the following instead'''


def require(ok: bool, detail: str) -> None:
    if not ok:
        raise ValueError(detail)


def normalize(src: str) -> str:
    return src.replace('\r\n', '\n')


def function_window(text: str, name: str, next_marker: str) -> str:
    start = text.find(name)
    require(start >= 0 and text.find(name, start + 1) == -1,
            f'Source no longer has one unique {name!r} definition')
    end = text.find(next_marker, start)
    require(end >= 0, f'Expected end marker missing after {name!r}')
    return text[start:end]


def validate_base(text: str) -> None:
    require(text.count(MASS_CORRECT) == 1, 'The previously verified doublet epsilon mass correction is missing')
    require(text.count(MASS_OLD) == 0, 'Uncorrected charged doublet MF delta is present')
    require(text.count(BUILD_ANCHOR) == 1, 'Unexpected UV build call site')
    require(text.count(BASE_FUNCTION_END) == 1, 'Cannot identify end of baseline T3 Yukawa function')
    require(text.count(BETAS_END) == 1, 'Cannot identify end of UV beta function')
    subset = function_window(text, 'T3RGBetaAddT3Yukawas[\n', '(* We add scalar and fermion mass term *)')
    require(subset.count('AddYukawa[') == 2,
            'Original two T3 Yukawa registrations changed')
    neutral = function_window(text, 'T3RGBetaAddMasses[\n', '(* ---------------------------------------------------------------------- *)')
    require(neutral.count('{F, F}') == 1 and neutral.count('{F, Fc}') == 1,
            'Neutral Majorana/charged Dirac mass branches changed')
    require(text.count('T3RGBetaOneLoopBetas[] := Module[') == 1,
            'Expected one UV beta function definition')
    require(text.count('T3RGBetaAddGeneralFermionYukawas[') == 0,
            'Extension is already installed or source partly changed')


def validate_installed(text: str) -> None:
    require(text.count(MASS_CORRECT) == 1 and MASS_OLD not in text,
            'Corrected charged-doublet mass invariant was lost')
    require(text.count(BUILD_NEW) == 1,
            'Opt-in extension must be invoked exactly once in UV build')
    require(text.count(EXTRA_DEFINITION) == 1,
            'Source must contain exactly one complete optional Yukawa definition')
    require(text.count(BETAS_NEW) == 1,
            'Full UV beta system does not include optional Yukawa beta keys')
    require(text.count('AddFermionMass[\n            MF,\n            {F, F},') == 1,
            'Neutral MF branch no longer present')
    require(text.count('T3RGBetaBuildShared[') >= 1,
            'Shared-model definitions unexpectedly lost')


def transform(source: str) -> str:
    text = normalize(source)
    validate_base(text)
    result = text.replace(BASE_FUNCTION_END, '    ];\n];\n' + EXTRA_DEFINITION + '\n(* We add scalar and fermion mass term *)', 1)
    result = result.replace(BUILD_ANCHOR, BUILD_NEW, 1)
    result = result.replace(BETAS_END, BETAS_NEW, 1)
    validate_installed(result)
    recovered = result.replace(BETAS_NEW, BETAS_END, 1).replace(BUILD_NEW, BUILD_ANCHOR, 1).replace(
        '    ];\n];\n' + EXTRA_DEFINITION + '\n(* We add scalar and fermion mass term *)',
        BASE_FUNCTION_END, 1)
    require(recovered == text, 'Unrelated RGBeta source text would change')
    return result.replace('\n', '\r\n') if '\r\n' in source else result


def fixture() -> str:
    # Minimal but contract-equivalent structure for self-tests. Real installed
    # source must pass the same guards before any bytes are modified.
    return (
        'T3RGBetaAddT3Yukawas[\n' +
        'AddYukawa[y1]; AddYukawa[y2];\n' + BASE_FUNCTION_END + '\n' +
        'T3RGBetaAddMasses[\n' +
        '    AddFermionMass[\n            MF,\n            {F, F},\n        ];\n' +
        MASS_CORRECT + '\n        ];\n];\n(* ---------------------------------------------------------------------- *)\n' +
        'T3RGBetaBuild[\n' + BUILD_ANCHOR + '\n];\n' +
        'T3RGBetaOneLoopBetas[] := Module[\n' + BETAS_END + '\n' +
        'T3RGBetaBuildShared[foo_] := foo;\n'
    )


def self_test() -> None:
    src = fixture()
    dst = transform(src)
    validate_installed(dst)
    require(dst.replace('\r\n', '\n').count('AddYukawa[extra, fields,') == 1,
            'Extra Yukawa registration missing')
    require(transform(src.replace('\n', '\r\n')).count('\r\n') > 5,
            'CRLF line endings not retained')
    for invalid in (
        src.replace(MASS_CORRECT, MASS_OLD),
        src.replace('AddYukawa[y2]', 'AddYukawa[y3] AddYukawa[y2]'),
        src.replace(BUILD_ANCHOR, BUILD_ANCHOR + '\n' + BUILD_ANCHOR),
    ):
        try:
            transform(invalid)
        except ValueError:
            continue
        raise AssertionError('Malformed source unexpectedly accepted')
    print('PASS: guarded complete-file replacement; 6 opt-in cases, neutral MF preserved, CRLF and negative controls')


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument('--self-test', action='store_true')
    g.add_argument('--install', action='store_true')
    g.add_argument('--verify', action='store_true')
    opts = p.parse_args()
    if opts.self_test:
        self_test()
        return 0
    require(MODEL.is_file(), 'Missing current local model source: '+str(MODEL))
    original = MODEL.read_bytes()
    require(not original.startswith(b'\xef\xbb\xbf'), 'Unexpected UTF-8 BOM')
    src = original.decode('utf-8')
    if opts.verify:
        validate_installed(normalize(src))
        print('PASS: complete opt-in RGBeta general-Yukawa model installed, mass epsilon retained')
        return 0
    if 'T3RGBetaAddGeneralFermionYukawas[' in src:
        validate_installed(normalize(src))
        print('Already integrated, no modifications: '+str(MODEL))
        return 0
    rendered = transform(src).encode('utf-8')
    stamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    backup = MODEL.with_name('T3RGBetaModel.before_general_yukawa_'+stamp+'.wl')
    require(not backup.exists(), 'Backup name collision; no changes made')
    shutil.copy2(MODEL, backup)
    try:
        with tempfile.NamedTemporaryFile(dir=MODEL.parent, suffix='.wl', delete=False) as tmp:
            tmp.write(rendered)
            temp = Path(tmp.name)
        os.replace(temp, MODEL)
        validate_installed(normalize(MODEL.read_text(encoding='utf-8')))
    except Exception:
        shutil.copy2(backup, MODEL)
        raise
    print('PASS: local full replacement installed: '+str(MODEL))
    print('Backup: '+str(backup))
    print('Opt-in switch: T3_EXPERIMENTAL_GENERAL_FERMIONS=1')
    print('All generated RGE outputs must be regenerated; physical coupling dictionary NOT certified')
    return 0


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except (ValueError, OSError, UnicodeError) as exc:
        print('BLOCKED (no new source committed): '+str(exc), file=sys.stderr)
        raise SystemExit(2)
