#!/usr/bin/env python3
"""Read-only B/C scalar-threshold coupling-support closure (not a CT).

Match the *coupling monomials* of each actual native Y/LLHH pole term against
products of the EFT1 LLSS Wilson and HHSS quartic terms.  Restrict internal
scalar pairs to YY or XY, require scalar-complex-conjugation matching, the
lepton-bilinear chirality and the external Higgs field orientation.

No SU(2) index contraction, numerical Wilson coefficient, beta-function
normalization, UV/IR pole assignment, or physical C5 is computed here.

The program is intentionally independent of the earlier toy beta calculator.
It never overwrites input expressions.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from hashlib import sha1, sha256
import json
from pathlib import Path
import re
import sys

SOURCE_BASE = Path('output/general_t3_mass_eigen_matchete')
KNOWN_GIT_BLOBS = {
    'B_m1': {
        'eft1_transition_eom': 'f6a3fd1fb2da9468909fdd89d0eeb0bf11635297',
        'llss_tree': '76067102a585b93969c063158fcdd3a9ebce85f7',
        'y_raw_llhh': '1e0eb2c51e20604b1c4923919911aaff6c4cbb47',
    },
    'C_m1': {
        'eft1_transition_eom': 'f33da5d1033dc3ddc94ad74da83c0bc7f40a4e66',
        'llss_tree': 'aae6a68ce1ef9a39a1c08350f0017e2aa4e18ff5',
        'y_raw_llhh': '3f47f81e04c4a9d1189b55f6d42bc698184e5add',
    },
}
COUPLING = re.compile(r'Coupling\[([A-Za-z][A-Za-z0-9_]*),')
FIELDS = ('H', 'NewScalarX', 'NewScalarY')
MASS_COUPLINGS = {'MF', 'MSX', 'MSY'}
U12 = 'Coupling[T3MEU12, {}, 0]'
U21 = 'Coupling[T3MEU21, {}, 0]'
U11 = 'Coupling[T3MEU11, {}, 0]'
U22 = 'Coupling[T3MEU22, {}, 0]'


class AuditError(ValueError):
    pass


def require(test: bool, message: str) -> None:
    if not test:
        raise AuditError(message)


def git_blob(data: bytes) -> str:
    return sha1(b'blob ' + str(len(data)).encode('ascii') + b'\0' + data).hexdigest()


def split_additive(expression: str) -> list[str]:
    """Split top-level InputForm +/-, never inside parentheses, lists or calls."""
    expression = expression.strip()
    require(bool(expression), 'Empty source')
    stack: list[str] = []
    pairs = {')': '(', ']': '[', '}': '{'}
    out: list[str] = []
    start = 0
    in_string = False
    escaped = False
    for k, c in enumerate(expression):
        if in_string:
            if escaped:
                escaped = False
            elif c == '\\':
                escaped = True
            elif c == '"':
                in_string = False
            continue
        if c == '"':
            in_string = True
            continue
        if c in '([{':
            stack.append(c)
        elif c in ')]}':
            require(bool(stack) and stack.pop() == pairs[c],
                    f'Unbalanced InputForm at position {k}')
        elif not stack and k > start and c in '+-' and expression[k-1] == ' ':
            out.append(expression[start:k].strip())
            start = k
    require(not stack and not in_string, 'Unbalanced InputForm at end')
    out.append(expression[start:].strip())
    require(all(out), 'Empty additive term')
    return out


def field_count(expr: str, label: str) -> int:
    return expr.count('Field[' + label + ', Scalar')


def barred_field_count(expr: str, label: str) -> int:
    return expr.count('Bar[Field[' + label + ', Scalar')


def chirality(expr: str) -> str:
    pr = 'Proj[1]' in expr
    pl = 'Proj[-1]' in expr
    require(not (pr and pl), 'Two chirality projectors in a term')
    return 'PR' if pr else 'PL' if pl else 'NONE'


def pairing(expr: str) -> str:
    x, y = field_count(expr, 'NewScalarX'), field_count(expr, 'NewScalarY')
    require(x + y == 2, 'Not a two-scalar operator')
    return ('X' * x) + ('Y' * y)


def expand_known_hhxy_quartic(expr: str) -> list[str]:
    """Expand exactly the two B/C HHXY lambdaT3 mixing products.

    This is NOT general symbolic expansion.  Reject any new/unknown additive
    expression rather than silently parse its factors as one monomial.
    """
    for barred in (True, False):
        def cv(s: str) -> str:
            return 'Bar[' + s + ']' if barred else s
        a = cv(U12) + '*' + cv(U21)
        b = cv(U11) + '*' + cv(U22)
        target = '(' + a + ' + ' + b + ')'
        if target in expr:
            require(expr.count(target) == 1, 'Repeated mixed-quartic expansion')
            return [expr.replace(target, a), expr.replace(target, b)]
    require(' + ' not in expr, 'Unexpected unexpanded additive quartic expression')
    return [expr]


def coupling_signature(expr: str) -> tuple[str, ...]:
    """Extract factors, exact conjugation and integer powers from InputForm.

    Deliberately ignores all CGs, SU(2), rational factors, flavour indices and
    denominator masses. It is therefore ONLY a monomial-support projection.
    """
    out: list[str] = []
    for match in COUPLING.finditer(expr):
        name = match.group(1)
        is_barred = expr[max(0, match.start()-4):match.start()] == 'Bar['
        depth = 0
        end = None
        for j in range(match.start(), len(expr)):
            if expr[j] == '[':
                depth += 1
            elif expr[j] == ']':
                depth -= 1
                if depth == 0:
                    end = j
                    break
        require(end is not None, 'Unterminated Coupling call')
        after = end + 1
        if is_barred:
            require(after < len(expr) and expr[after] == ']',
                    'Unclosed Bar[Coupling] expression')
            after += 1
        power = 1
        if after < len(expr) and expr[after] == '^':
            digits = re.match(r'[0-9]+', expr[after+1:])
            require(digits is not None, 'Unsupported coupling power')
            power = int(digits.group(0))
            require(1 <= power <= 8, 'Unsupported coupling power size')
        if name in MASS_COUPLINGS:
            continue
        out.extend([('bar_' if is_barred else '') + name] * power)
    require(bool(out), 'No coupling factors in expected operator')
    return tuple(sorted(out))


def matching_key(expr: str, *, projected_higgs_from: str | None = None,
                 bilinear_from: str | None = None) -> tuple:
    return (
        coupling_signature(expr),
        chirality(bilinear_from if bilinear_from is not None else expr),
        barred_field_count(projected_higgs_from if projected_higgs_from is not None else expr, 'H'),
    )


def source_paths(root: Path, stem: str) -> dict[str, Path]:
    return {
        'eft1_transition_eom': root / 'eft1_fermion_first_tree' /
            (stem + '_final_eft1_fermion_first_tree_eft1_transition_eom.inputform.txt'),
        'llss_tree': root / 'eft1_fermion_first_tree' /
            (stem + '_final_eft1_fermion_first_tree_llss_tree.inputform.txt'),
        'y_raw_llhh': root / 'eft1_scalar_Y_only_canonical_probe' /
            (stem + '_final_y_only_canonical_match_raw_scalar_llhh.inputform.txt'),
    }


def analyze_model(root: Path, stem: str, strict_snapshot: bool = False) -> dict:
    files = source_paths(root, stem)
    data = {}
    provenance = {}
    for key, path in files.items():
        require(path.is_file(), 'Missing native source: ' + str(path))
        blob = path.read_bytes()
        got_blob = git_blob(blob)
        expected = KNOWN_GIT_BLOBS[stem][key]
        if strict_snapshot:
            require(got_blob == expected,
                    f'{stem}/{key}: different source snapshot {got_blob} != {expected}')
        data[key] = blob.decode('utf-8').strip()
        provenance[key] = {
            'project_path': str(path).replace('\\', '/'),
            'sha256': sha256(blob).hexdigest(),
            'git_blob_sha1': got_blob,
            'matches_inspected_git_snapshot': got_blob == expected,
        }
    full_terms = split_additive(data['eft1_transition_eom'])
    llss_terms = split_additive(data['llss_tree'])
    raw_terms = split_additive(data['y_raw_llhh'])
    wilson = [(n, t) for n, t in enumerate(llss_terms, 1)
              if field_count(t, 'NewScalarX')+field_count(t, 'NewScalarY') == 2
              and t.count('Field[l, Fermion') == 2 and 'GammaCC' in t]
    quartics = [(n, t) for n, t in enumerate(full_terms, 1)
                if field_count(t, 'H') == 2
                and field_count(t, 'NewScalarX')+field_count(t, 'NewScalarY') == 2
                and 'Field[l, Fermion' not in t]
    poles = [(n, t) for n, t in enumerate(raw_terms, 1) if 'ϵ' in t]
    logs = [(n, t) for n, t in enumerate(raw_terms, 1)
            if 'Log[' in t and 'ϵ' not in t]
    constants = [(n, t) for n, t in enumerate(raw_terms, 1)
                 if 'Log[' not in t and 'ϵ' not in t]
    expected_llss = 20 if stem == 'B_m1' else 24
    require(len(wilson) == expected_llss, f'{stem}: unexpected LLSS count {len(wilson)}')
    require(len(quartics) == 42, f'{stem}: unexpected HHSS count {len(quartics)}')
    require((len(raw_terms), len(poles), len(logs), len(constants)) == (132,50,50,32),
            f'{stem}: unexpected Y/LLHH Laurent inventory')
    require(all('hbar' in t and 'GammaCC' in t and t.count('Field[l, Fermion') == 2
                and field_count(t, 'H') == 2 for _,t in poles),
            f'{stem}: unsupported raw pole operator')

    # Project EXACT supported coupling products, retaining source coordinates.
    supports: dict[tuple, list[dict]] = defaultdict(list)
    n_expanded = 0
    for wi, w in wilson:
        pair = pairing(w)
        if pair == 'XX':
            continue  # purely light X loops, not hard Y-specific support
        require(chirality(w) in ('PR', 'PL'), 'LLSS chirality missing')
        for qi, q in quartics:
            if pairing(q) != pair:
                continue
            if any(barred_field_count(w, sc) + barred_field_count(q, sc)
                   != field_count(w, sc) for sc in ('NewScalarX','NewScalarY')):
                continue  # complex-field Wick-contraction orientation
            variations = expand_known_hhxy_quartic(q)
            n_expanded += len(variations) - 1
            for variant_no, expanded in enumerate(variations, 1):
                key = matching_key(w + ' ' + expanded,
                                   projected_higgs_from=q, bilinear_from=w)
                supports[key].append({
                    'pair': pair,
                    'LLSS_term': wi,
                    'HHSS_term': qi,
                    'HHSS_mixing_variant': variant_no,
                })

    by_pair = Counter()
    records = []
    unmatched = []
    for ri, raw in poles:
        key = matching_key(raw)
        origins = supports.get(key, [])
        classes = sorted({a['pair'] for a in origins})
        for c in classes:
            by_pair[c] += 1
        entry = {
            'native_pole_term': ri,
            'chirality': key[1],
            'external_higgs_bar_count': key[2],
            'coupling_factors': list(key[0]),
            'supported_pair_types': classes,
            'source_candidates': origins,
        }
        records.append(entry)
        if not origins:
            unmatched.append(entry)
    # Fail closed for this independently inspected B/C source snapshot.
    require(not unmatched, f'{stem}: {len(unmatched)} native pole monomials unaccounted')
    require(by_pair == Counter({'YY':18,'XY':32}),
            f'{stem}: unexpected YY/XY decomposition: {dict(by_pair)}')
    require(not any('XX' in e['supported_pair_types'] for e in records),
            'Pure-XX contribution incorrectly accepted')
    require(n_expanded >= 2, 'HHXY lambdaT3 sum not expanded')
    return {
        'model': {'B_m1':'B_-1','C_m1':'C_-1'}[stem],
        'source_provenance': provenance,
        'native_counts': {
            'eft1_tree_llss_terms': len(wilson),
            'eft1_tree_hhss_terms': len(quartics),
            'y_llhh_raw_terms': len(raw_terms),
            'y_llhh_pole_terms': len(poles),
            'y_llhh_log_terms': len(logs),
            'y_llhh_constant_terms': len(constants),
            'native_poles_supported_by_tree_coupling_products': len(records),
            'native_poles_without_support': len(unmatched),
        },
        'pole_support_by_internal_pair': dict(by_pair),
        'HHXY_quartic_mixing_sum_expansions': n_expanded,
        'all_native_pole_monomials_have_source_coupling_support': True,
        'unambiguously_identifies_all_diagrams': False,
        'SU2_component_coefficient_comparison_performed': False,
        'majorana_GammaCC_normalization_certified': False,
        'EFT1_EFT2_counterterms_computed': False,
        'UV_IR_poles_assigned': False,
        'physical_C5': None,
        'pole_records': records,
    }


def self_test() -> None:
    s = '-(2*Foo[1]+3*Foo[2]) + X[]^2 - Y[]'
    assert len(split_additive(s)) == 3
    expr = '-Bar[Coupling[lambdaT3, {}, 0]]*Coupling[T3MEU12, {}, 0]^2*' \
           'Coupling[y1, {}, 0]*Coupling[y1, {}, 0]/Coupling[MF, {}, 0]'
    assert coupling_signature(expr) == ('T3MEU12','T3MEU12','bar_lambdaT3','y1','y1')
    q = '-' + 'Bar[Coupling[lambdaT3, {}, 0]]*(' + \
        'Bar['+U12+']*Bar['+U21+'] + Bar['+U11+']*Bar['+U22+'])'
    variants = expand_known_hhxy_quartic(q)
    assert len(variants) == 2 and all(' + ' not in t for t in variants)
    assert len({coupling_signature(t) for t in variants}) == 2
    # A fabricated or conjugation-inverted coupling cannot match unchanged.
    ref = matching_key(expr)
    assert matching_key(expr.replace('y1','y3')) != ref
    assert matching_key(expr.replace('Bar[Coupling[lambdaT3, {}, 0]]',
                                     'Coupling[lambdaT3, {}, 0]')) != ref
    assert git_blob(b'abc') == 'f2ba8f84ab5c1bce84a7b441cb1959cfc7093b7f'
    print('PASS: InputForm splitting and conjugated/exponent coupling extraction')
    print('PASS: exact HHXY lambdaT3 sum alternatives and tamper rejection')
    print('PASS: no beta/CT, SU2 normalization, or physical C5 asserted')


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--self-test', action='store_true')
    p.add_argument('--root', type=Path, default=SOURCE_BASE,
                   help='directory containing eft1_fermion_first_tree and EFT1 Y-only output')
    p.add_argument('--output', type=Path)
    p.add_argument('--strict-known-snapshot', action='store_true',
                   help='require exact known Git blob hashes for all source files')
    args = p.parse_args()
    try:
        if args.self_test:
            self_test()
        if args.output:
            report = {
                'format':'t3-bc-y-mixed-scalar-native-pole-coupling-support-v1',
                'scope':'coupling monomial and scalar-orientation support only',
                'reference':'actual native B/C EFT1 tree and Y/LLHH hard-pole expressions',
                'models':{},
                'SU2_coefficient_level_pole_cancellation_verified':False,
                'MSbar_counterterm_computed':False,
                'physical_C5':None,
            }
            for stem in ('B_m1','C_m1'):
                entry = analyze_model(args.root, stem,
                                      strict_snapshot=args.strict_known_snapshot)
                report['models'][entry['model']] = entry
                counts = entry['pole_support_by_internal_pair']
                print(f"{entry['model']}: 50/50 native pole monomials supported; "
                      f"YY={counts.get('YY',0)}, XY={counts.get('XY',0)}, "
                      'XX=0; counterterms=NOT CALCULATED')
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n',
                                   encoding='utf-8')
            print('PASS: source-level support audit -> ' + str(args.output))
        require(args.self_test or args.output is not None,
                'Specify --self-test and/or --output')
        return 0
    except (AuditError, OSError, UnicodeError, AssertionError) as exc:
        print('FAIL (no physical matching certificate): ' + str(exc), file=sys.stderr)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
