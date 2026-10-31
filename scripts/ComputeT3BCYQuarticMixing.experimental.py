"""Independent Y-related (YY and mixed XY) scalar-quartic contribution to LLHH beta tensor.

Reference: Aebischer, Bresciani, Selimovic, arXiv:2512.16890v1,
Eq. (4.85), scalar-quartic term, and Eq. (3.3).

At fixed fermion indices ij in a REAL scalar basis, with fully symmetric
lambda_abcd (V4=lambda_abcd phi_a phi_b phi_c phi_d / 4!) and C_ijab
(L5=C_ijab psi_i psi_j phi_a phi_b + h.c.),

  (16 pi^2) [beta(C_ijab)]_(quartic, Y-related)
       = sum_{c,d in X+Y; c or d in Y} lambda_abcd C_ijcd.

When no X components are supplied, this reduces to the old Y-only contraction.

The result is an independently computed *contribution* to the EFT1 anomalous
dimension. It is NOT the complete beta function, not a counterterm, and not a
physical Y-threshold matching coefficient. The input must already be mapped
to the paper's real-scalar and two-component Weyl normalization. Native
Matchete GammaCC objects cannot be inserted without a verified conversion.

All coefficients are represented as products of exact rational numbers and
formal symbols, not as parsed Python or Mathematica expressions.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction
import json
from pathlib import Path
import re
from typing import Any

import sympy as sp

FORMAT = 't3bc-y-scalar-quartic-beta-v1'
TOKEN = re.compile(r'^[A-Za-z][A-Za-z0-9_]*$')


class PhysicsInputError(ValueError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise PhysicsInputError(message)


def exact_sum(term_list: Any) -> sp.Expr:
    """Parse a restricted exact multiplicative-polynomial JSON encoding."""
    require(isinstance(term_list, list), 'Each coefficient must be a list of exact terms')
    expression = sp.Integer(0)
    for item in term_list:
        require(isinstance(item, dict), 'Coefficient term must be an object')
        require(set(item) <= {'rational', 'factors'}, 'Unknown coefficient term key')
        raw_rational = item.get('rational', '1')
        require(isinstance(raw_rational, str) and
                bool(re.fullmatch(r'[+-]?[0-9]+(?:/[0-9]+)?', raw_rational)),
                'Rational must be an exact integer or p/q string')
        try:
            rational = Fraction(raw_rational)
        except (ValueError, TypeError, ZeroDivisionError) as exc:
            raise PhysicsInputError('Bad rational coefficient') from exc
        value = sp.Rational(rational.numerator, rational.denominator)
        factors = item.get('factors', {})
        require(isinstance(factors, dict), 'factors must be an object')
        for name, power in factors.items():
            require(isinstance(name, str) and bool(TOKEN.fullmatch(name)),
                    'Invalid factor label: ' + str(name))
            require(type(power) is int and -12 <= power <= 12,
                    'Factor power must be an integer between -12 and 12')
            # 'ImaginaryUnit' is explicit; arbitrary factors are free symbols.
            atom = sp.I if name == 'ImaginaryUnit' else (
                sp.sqrt(3) if name == 'Sqrt3' else sp.Symbol(name))
            value *= atom**power
        expression += value
    return sp.expand(expression)


def idx_tuple(raw: Any, n: int, arity: int, label: str) -> tuple[int, ...]:
    require(isinstance(raw, list) and len(raw) == arity,
            label + ': incorrect number of scalar indices')
    require(all(type(x) is int and 0 <= x < n for x in raw),
            label + ': scalar index must be an integer in range')
    return tuple(sorted(raw))


def flavor_pair(raw: Any, allowed: set[str]) -> tuple[str, str]:
    require(isinstance(raw, list) and len(raw) == 2 and
            all(isinstance(x, str) and x in allowed for x in raw),
            'Invalid fermion/flavour pair')
    return tuple(sorted(raw))


def build_tensor(entries: Any, n: int, arity: int,
                 allowed_flavours: set[str] | None, name: str) -> dict:
    require(isinstance(entries, list), name + ' must be an array')
    tensor: dict[tuple, sp.Expr] = defaultdict(lambda: sp.Integer(0))
    seen: set[tuple] = set()
    for entry in entries:
        require(isinstance(entry, dict), name + ' entry must be an object')
        require(set(entry) == ({'scalars', 'terms'} if allowed_flavours is None else
                               {'flavours', 'scalars', 'terms'}),
                name + ' entry has missing or extra keys')
        idx = idx_tuple(entry['scalars'], n, arity, name)
        key = idx if allowed_flavours is None else (
            flavor_pair(entry['flavours'], allowed_flavours), idx)
        require(key not in seen, name + ': repeated canonical tensor entry ' + str(key))
        seen.add(key)
        tensor[key] = exact_sum(entry['terms'])
    return dict(tensor)


def check_input(data: dict) -> tuple:
    require(data.get('format') == FORMAT, 'Wrong input format')
    labels = data.get('scalar_labels')
    flavours = data.get('fermion_labels')
    require(isinstance(labels, list) and len(labels) >= 2 and
            all(isinstance(x, str) and bool(x) for x in labels) and
            len(set(labels)) == len(labels), 'scalar_labels must be distinct strings')
    require(isinstance(flavours, list) and bool(flavours) and
            all(isinstance(x, str) and bool(x) for x in flavours) and
            len(set(flavours)) == len(flavours),
            'fermion_labels must be distinct strings')
    n = len(labels)
    h = data.get('higgs_scalar_indices')
    y = data.get('heavy_y_scalar_indices')
    x = data.get('light_x_scalar_indices', [])
    for name, value in [('higgs_scalar_indices', h), ('heavy_y_scalar_indices', y)]:
        require(isinstance(value, list) and bool(value) and
                len(set(value)) == len(value) and
                all(type(z) is int and 0 <= z < n for z in value),
                name + ' must contain distinct valid scalar indices')
    require(isinstance(x, list) and len(set(x)) == len(x) and
            all(type(z) is int and 0 <= z < n for z in x),
            'light_x_scalar_indices must contain distinct valid scalar indices')
    require(not(set(h) & set(y)) and not(set(x) & set(h)) and
            not(set(x) & set(y)), 'H/X/Y scalar subsets overlap')
    allowed = set(flavours)
    quartic = build_tensor(data.get('quartic_tensor'), n, 4, None, 'quartic_tensor')
    wilson = build_tensor(data.get('wilson_tensor'), n, 2, allowed, 'wilson_tensor')
    require(bool(quartic) and bool(wilson), 'Non-empty quartic and Wilson tensors required')
    require(all((i in y and j in y) or (i in x and j in y) or
                (i in y and j in x) for _, (i, j) in wilson),
            'Wilson terms must contain YY or XY, not pure XX or H')
    native = None
    if 'native_LLHH_pole_without_hbar' in data:
        native = build_tensor(data['native_LLHH_pole_without_hbar'], n, 2, allowed,
                              'native_LLHH_pole_without_hbar')
        require(all(a in h and b in h for _, (a, b) in native),
                'Native poles must already be in the chosen HH real-scalar basis')
    return (labels, flavours, tuple(sorted(h)), tuple(sorted(y)),
            tuple(sorted(x)), quartic, wilson, native)


def calculate(data: dict) -> dict:
    labels, flavours, h, y, x, quartic, wilson, native = check_input(data)
    # Eq. (4.85) sums over ORDERED real-scalar indices.
    # Here keep YY and XY terms (at least one index in Y), not XX.
    # For a symmetric pair (c,d) with c != d, multiplicity is two.
    out: dict[tuple, sp.Expr] = defaultdict(lambda: sp.Integer(0))
    nonzero_input_pair_counts = {'YY': 0, 'XY': 0}
    for (ij, cd), value in wilson.items():
        if value == 0:
            continue
        nonzero_input_pair_counts['YY' if cd[0] in y and cd[1] in y else 'XY'] += 1
        multiplicity = 1 if cd[0] == cd[1] else 2
        for h_position, a in enumerate(h):
            for b in h[h_position:]:
                key = (ij, (a, b))
                q = quartic.get(tuple(sorted((a, b, *cd))), sp.Integer(0))
                out[key] += multiplicity * q * value
    out = {key: sp.expand(value) for key, value in out.items()
           if sp.expand(value) != 0}

    result: dict[str, Any] = {
        'format': 't3bc-y-quartic-beta-contribution-result-v1',
        'reference': 'arXiv:2512.16890v1, Eqs. (3.3) and (4.85)',
        'normalization_required':
            'L5 = C_ijab psi_i psi_j phi_a phi_b + h.c.; '
            'V4 = lambda_abcd phi_a phi_b phi_c phi_d/4!; all phi real',
        'independent_beta_calculation':
            '(16*pi**2)*beta(C_ijab)|quartic,Y-related = '
            'sum_(c,d in X+Y; c or d in Y) lambda_abcd*C_ijcd',
        'higgs_scalar_labels': [labels[i] for i in h],
        'heavy_y_scalar_labels': [labels[i] for i in y],
        'light_x_scalar_labels': [labels[i] for i in x],
        'included_internal_pair_types': ('YY, XY' if x else 'YY only'),
        'nonzero_input_wilson_pair_counts': nonzero_input_pair_counts,
        'beta_y_terms': [
            {'flavours': list(ij), 'higgs_scalars': [labels[x] for x in ab],
             'higgs_scalar_indices': list(ab),
             '16pi2_beta_y': str(val)}
            for (ij, ab), val in sorted(out.items())
        ],
        'requires_mapped_paper_tensor_normalization': True,
        'gamma_cc_to_paper_normalization_checked_here': False,
        'includes_all_EFT1_beta_contributions': False,
        'independent_y_MSbar_counterterm_certified': False,
        'raw_Matchete_pole_UV_IR_origin_assigned': False,
        'net_pole_cancellation_certified': False,
        'physical_C5_Y': None,
    }
    if native is not None:
        keys = sorted(set(native) | set(out))
        comparison = []
        all_same = True
        all_opposite = True
        for key in keys:
            b = out.get(key, sp.Integer(0))
            p = native.get(key, sp.Integer(0))
            plus = sp.expand(2 * p - b)
            minus = sp.expand(2 * p + b)
            all_same &= (plus == 0)
            all_opposite &= (minus == 0)
            comparison.append({
                'flavours': list(key[0]),
                'higgs_scalar_indices': list(key[1]),
                '2P_minus_16pi2_beta_y': str(plus),
                '2P_plus_16pi2_beta_y': str(minus),
            })
        result['optional_pole_beta_comparison'] = {
            'same_sign_all_components': all_same,
            'opposite_sign_all_components': all_opposite,
            'differences': comparison,
            'note': 'P must have the explicit Matchete hbar=1/(16pi^2) removed. '
                    'Comparison assumes an independently established identical '
                    'real-scalar, flavour, chirality and operator normalization. '
                    'Even an exact equality is not generally required by this subset '
                    'and cannot establish UV/IR assignment, any counterterm, or C5.',
            'physical_matching_certified': False,
        }
    else:
        result['optional_pole_beta_comparison'] = None
    return result


def toy() -> dict:
    # This is an analytic toy with one H and two Y real components, NOT T3 B/C.
    return {
        'format': FORMAT,
        'scalar_labels': ['H', 'Y1', 'Y2'],
        'fermion_labels': ['e', 'mu'],
        'higgs_scalar_indices': [0],
        'heavy_y_scalar_indices': [1, 2],
        'quartic_tensor': [
            {'scalars': [0, 0, 1, 1],
             'terms': [{'rational': '1', 'factors': {'lam1': 1}}]},
            {'scalars': [0, 0, 1, 2],
             'terms': [{'rational': '1/3', 'factors': {'lam2': 1}}]},
        ],
        'wilson_tensor': [
            {'flavours': ['e', 'mu'], 'scalars': [1, 1],
             'terms': [{'rational': '1/2', 'factors': {'c1': 1}}]},
            {'flavours': ['e', 'mu'], 'scalars': [1, 2],
             'terms': [{'rational': '3/4', 'factors': {'c2': 1}}]},
        ],
    }


def self_test() -> None:
    a = toy()
    _, _, _, _, _, q, w, _ = check_input(a)
    require(q[(0, 0, 1, 1)] == sp.Symbol('lam1'), 'Quartic tensor encoding failed')
    require(w[(('e', 'mu'), (1, 2))] == sp.Rational(3, 4)*sp.Symbol('c2'),
            'Wilson tensor encoding failed')
    result = calculate(a)
    actual = result['beta_y_terms'][0]['16pi2_beta_y']
    wanted = (sp.Symbol('c1')*sp.Symbol('lam1')/2 +
              sp.Symbol('c2')*sp.Symbol('lam2')/2)
    require(sp.expand(sp.sympify(actual) - wanted) == 0,
            'Eq. (4.85) quartic contraction or ordered-pair multiplicity failed')
    a['native_LLHH_pole_without_hbar'] = [
        {'flavours': ['e', 'mu'], 'scalars': [0, 0],
         'terms': [
             {'rational': '1/4', 'factors': {'c1': 1, 'lam1': 1}},
             {'rational': '1/4', 'factors': {'c2': 1, 'lam2': 1}},
         ]}
    ]
    r = calculate(a)
    require(r['optional_pole_beta_comparison']['same_sign_all_components'],
            'Synthetic exact comparison failed')
    require(not r['net_pole_cancellation_certified'],
            'Synthetic comparison must not claim physical certification')
    a['native_LLHH_pole_without_hbar'][0]['terms'][0]['rational'] = '1/3'
    r = calculate(a)
    require(not r['optional_pole_beta_comparison']['same_sign_all_components'],
            'Mismatch guard did not reject a changed pole')
    bad = toy()
    bad['wilson_tensor'][0]['scalars'] = [0, 1]
    try:
        calculate(bad)
    except PhysicsInputError:
        pass
    else:
        raise AssertionError('Non-YY Wilson operator was accepted')
    # Actual new physics-case regression: the mixed X/Y propagator pair contributes
    # twice when its real scalar indices differ, while a pure XX pair is excluded.
    mixed = {
        'format': FORMAT,
        'scalar_labels': ['H', 'X', 'Y'],
        'fermion_labels': ['e', 'mu'],
        'higgs_scalar_indices': [0],
        'light_x_scalar_indices': [1],
        'heavy_y_scalar_indices': [2],
        'quartic_tensor': [{'scalars':[0,0,1,2],
                            'terms':[{'rational':'1', 'factors':{'lamXY':1}}]}],
        'wilson_tensor': [{'flavours':['e','mu'], 'scalars':[1,2],
                           'terms':[{'rational':'1/4', 'factors':{'cXY':1}}]}],
    }
    xy_result = calculate(mixed)
    require(sp.expand(sp.sympify(xy_result['beta_y_terms'][0]['16pi2_beta_y']) -
                      sp.Symbol('cXY')*sp.Symbol('lamXY')/2) == 0,
            'XY real-scalar ordered-pair multiplicity failed')
    mixed['wilson_tensor'][0]['scalars'] = [1,1]
    try:
        calculate(mixed)
    except PhysicsInputError:
        pass
    else:
        raise AssertionError('Pure XX Wilson pair was accepted')
    print('PASS: Eq. (4.85) exact YY/XY real-scalar mixing and pair multiplicity')
    print('PASS: exact synthetic pole comparison, XX exclusion and tamper rejection')
    print('PASS: no counterterm, UV/IR assignment, or physical C5 certification')


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--self-test', action='store_true')
    ap.add_argument('--write-toy-input', type=Path,
                    help='Writes a toy JSON template; not B/C benchmark data')
    ap.add_argument('--input', type=Path, help='Real-scalar tensor JSON')
    ap.add_argument('--output', type=Path, help='JSON report location')
    args = ap.parse_args()
    try:
        if args.self_test:
            self_test()
        if args.write_toy_input:
            args.write_toy_input.parent.mkdir(parents=True, exist_ok=True)
            args.write_toy_input.write_text(json.dumps(toy(), indent=2)+'\n',encoding='utf-8')
            print('Toy template written: ' + str(args.write_toy_input))
        if args.input or args.output:
            require(bool(args.input and args.output),
                    '--input and --output must both be provided')
            data = json.loads(args.input.read_text(encoding='utf-8'))
            result = calculate(data)
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
            print('Independent Y quartic beta diagnostic: ' + str(args.output))
            print('Physical C5: NOT CERTIFIED')
        require(bool(args.self_test or args.write_toy_input or args.input),
                'Specify --self-test, --write-toy-input, or --input/--output')
        return 0
    except (PhysicsInputError, ValueError, KeyError, OSError,
            json.JSONDecodeError) as exc:
        print('FAIL (physical matching not certified): '+str(exc))
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
