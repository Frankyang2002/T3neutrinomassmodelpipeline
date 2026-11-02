"""B1: structural audit of Matchete's already-exported native InputForm action.

This is deliberately a small *parser*, not a loose token/regex grep, and it
never evaluates Mathematica. It accepts a tightly restricted multiplicative
language of symbolic action terms; unfamiliar syntax fails closed.

Native source factors / chirality / conjugation are verified separately from
an unmeasured concrete GammaCC phase and the B2 tensor basis.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import re
from typing import Any


class ActionContractError(ValueError):
    """Native input disagrees with the asserted B1 operator contract."""


def require(cond: bool, reason: str) -> None:
    if not cond:
        raise ActionContractError(reason)


@dataclass(frozen=True)
class Atom:
    name: str


@dataclass(frozen=True)
class Call:
    name: str
    args: tuple[Any, ...]


@dataclass(frozen=True)
class Sequence:
    items: tuple[Any, ...]


@dataclass(frozen=True)
class Operator:
    name: str
    args: tuple[Any, ...]


_TOKEN = re.compile(r"\s*([A-Za-z$`_][A-Za-z0-9$`_]*|\d+|[\[\]{}(),+*/^-])")


class Parser:
    def __init__(self, source: str):
        require(isinstance(source, str) and len(source) < 200_000, 'Invalid InputForm text')
        self.tokens = []
        offset = 0
        while offset < len(source):
            match = _TOKEN.match(source, offset)
            require(match is not None, 'Unrecognised InputForm token at position ' + str(offset))
            self.tokens.append(match.group(1))
            offset = match.end()
        require(len(self.tokens) < 30_000, 'InputForm token limit exceeded')
        self.i = 0

    def peek(self) -> str | None:
        return self.tokens[self.i] if self.i < len(self.tokens) else None

    def take(self, expected: str | None = None) -> str:
        current = self.peek()
        require(current is not None and (expected is None or current == expected),
                f'Expected {expected!r}, got {current!r}')
        self.i += 1
        return current

    def parse(self) -> Any:
        value = self.sum()
        require(self.peek() is None, 'Trailing InputForm tokens')
        return value

    def sum(self) -> Any:
        value = self.multiply()
        while self.peek() in ('+', '-'):
            op = self.take()
            value = Operator('add' if op == '+' else 'subtract', (value, self.multiply()))
        return value

    def multiply(self) -> Any:
        value = self.unary()
        while self.peek() in ('*', '/'):
            op = self.take()
            value = Operator('multiply' if op == '*' else 'divide', (value, self.unary()))
        return value

    def unary(self) -> Any:
        if self.peek() == '-':
            self.take('-')
            value = self.unary()
            return -value if isinstance(value, Fraction) else Operator('negate', (value,))
        if self.peek() == '+':
            self.take('+')
            return self.unary()
        value = self.primary()
        if self.peek() == '^':
            self.take('^')
            value = Operator('power', (value, self.unary()))
        return value

    def arguments(self, closer: str) -> tuple[Any, ...]:
        out = []
        if self.peek() != closer:
            while True:
                out.append(self.sum())
                if self.peek() != ',':
                    break
                self.take(',')
        self.take(closer)
        return tuple(out)

    def primary(self) -> Any:
        tok = self.take()
        if tok == '(':
            value = self.sum()
            self.take(')')
            return value
        if tok == '{':
            return Sequence(self.arguments('}'))
        if tok.isdigit():
            return Fraction(int(tok))
        require(bool(re.fullmatch(r'[A-Za-z$`_][A-Za-z0-9$`_]*', tok)),
                'Unexpected InputForm primary ' + tok)
        if self.peek() == '[':
            self.take('[')
            return Call(tok, self.arguments(']'))
        return Atom(tok)


def parse(source: str) -> Any:
    return Parser(source).parse()


def multiplicative(value: Any) -> tuple[Fraction, list[Call]]:
    """Exact rational numeric factor and unordered symbolic action factors.

    Deliberately reject every symbolic denominator/addition/power, which would
    require a different interpretation of the native action coefficient.
    """
    if isinstance(value, Fraction):
        return value, []
    if isinstance(value, Call):
        return Fraction(1), [value]
    if isinstance(value, Operator):
        if value.name == 'negate':
            coeff, factors = multiplicative(value.args[0])
            return -coeff, factors
        if value.name == 'multiply':
            a, af = multiplicative(value.args[0]); b, bf = multiplicative(value.args[1])
            return a*b, af+bf
        if value.name == 'divide':
            a, af = multiplicative(value.args[0]); b, bf = multiplicative(value.args[1])
            require(not bf and b, 'Symbolic or zero denominator in action coefficient')
            return a/b, af
    raise ActionContractError('Unsupported additive/power/atomic native operator')


def iscall(node: Any, name: str, count: int | None = None) -> bool:
    return isinstance(node, Call) and node.name == name and (count is None or len(node.args) == count)


def named(node: Any, name: str) -> bool:
    return isinstance(node, Atom) and node.name == name


def wrapped(node: Any, wrapper: str, inner: Any = None) -> bool:
    return iscall(node, wrapper, 1) and (inner is None or node.args[0] == inner)


def field(node: Any, name: str) -> bool:
    return iscall(node, 'Field', 4) and named(node.args[0], name) and named(node.args[1],
             'Fermion' if name in ('NewFermion', 'l') else 'Scalar')


def field_form(node: Any, name: str, modifiers: tuple[str, ...] = ()) -> bool:
    while modifiers:
        require(modifiers[0] in ('Bar', 'Transp'), 'Invalid expected wrapper')
        if not iscall(node, modifiers[0], 1):
            return False
        node = node.args[0]
        modifiers = modifiers[1:]
    return field(node, name)


def coupling(node: Any, name: str, conjugated: bool = False) -> bool:
    if conjugated:
        if not iscall(node, 'Bar', 1):
            return False
        node = node.args[0]
    return (iscall(node, 'Coupling', 3) and named(node.args[0], name)
            and isinstance(node.args[1], Sequence) and node.args[2] == 0)


def cg(node: Any, name: str, conjugated: bool = False) -> bool:
    if not iscall(node, 'CG', 2) or not isinstance(node.args[1], Sequence):
        return False
    first = node.args[0]
    return wrapped(first, 'Bar', Atom(name)) if conjugated else named(first, name)


def projectors(node: Any, gamma_cc: bool, proj: int | None) -> bool:
    if not iscall(node, 'DiracProduct'):
        return False
    actual = node.args
    if proj is None:
        return gamma_cc and actual == (Atom('GammaCC'),)
    required = ((Atom('GammaCC'),) if gamma_cc else ()) + (Call('Proj', (Fraction(proj),)),)
    return actual == required


def ncm(factors: list[Call]) -> Call:
    items = [f for f in factors if f.name == 'NCM']
    require(len(items) == 1, 'Expected exactly one NCM bilinear')
    return items[0]


def verify_spinor(node: Call, model: str, vertex: str, hc: bool) -> dict:
    require(len(node.args) in (2, 3), 'Unexpected NCM arity')
    majorana = model in ('B_m1','C_m1')
    expected = None
    # For an anticommuting Weyl field, the HC is not interchangeable with
    # the original: charge conjugation and projector must be transported.
    if vertex == 'mass':
        if majorana:
            expected = (('NewFermion', ('Transp',)), (True, None), ('NewFermion', ()))
        else:
            expected = (('NewFermion', ('Bar',)), None, ('NewFermion', ()))
    elif vertex == 'y1' and not majorana:
        expected = ((('l', ('Transp',)), (True, -1), ('NewFermion', ())) if hc
                    else (('l', ('Bar',)), (True, 1), ('NewFermion', ('Transp', 'Bar'))))
    elif vertex == 'y2' and not majorana:
        expected = ((('NewFermion', ('Bar',)), (False, -1), ('l', ())) if hc
                    else (('l', ('Bar',)), (False, 1), ('NewFermion', ())))
    else:
        expected = ((('NewFermion', ('Transp',)), (True, -1), ('l', ())) if hc
                    else (('l', ('Bar',)), (False, 1), ('NewFermion', ())))
    left, middle, right = expected
    if middle is None:
        require(len(node.args) == 2, f'{model} {vertex}: Dirac mass NCM arity')
        actual_left, actual_right = node.args
    else:
        require(len(node.args) == 3 and projectors(node.args[1], *middle),
                f'{model} {vertex}: wrong GammaCC/projector')
        actual_left, actual_right = node.args[0], node.args[2]
    require(field_form(actual_left, *left) and field_form(actual_right, *right),
            f'{model} {vertex}: wrong NCM field ordering/conjugation')
    return {'kind': vertex, 'hc': hc, 'gamma_cc': bool(middle and middle[0]),
            'projector': middle[1] if middle else None, 'field_order_verified': True}


def verify_term(source: str, model: str, vertex: str, hc: bool = False) -> dict:
    coeff, factors = multiplicative(parse(source))
    majorana = model in ('B_m1','C_m1')
    if vertex == 'mass':
        require(not hc, 'Mass HC is not separately enumerated by this native witness')
        require(coeff == (Fraction(-1, 2) if majorana else Fraction(-1)),
                f'{model}: wrong fermion mass prefactor {coeff}')
        require(len(factors) == 2 and sum(coupling(f, 'MF') for f in factors) == 1,
                f'{model}: mass operator factors are not exactly MF*NCM')
    else:
        require(vertex in ('y1', 'y2'), 'Unsupported Yukawa')
        require(coeff == Fraction(1), f'{model} {vertex}: wrong Yukawa prefactor {coeff}')
        expected_scalar = 'NewScalar1' if vertex == 'y1' else 'NewScalar2'
        scalar_bar = (hc == (vertex == 'y1'))
        wanted = [
            lambda f: coupling(f, vertex, hc),
            lambda f: cg(f, 'T3Y1CG' if vertex == 'y1' else 'T3Y2CG', hc),
            lambda f: field_form(f, expected_scalar, ('Bar',) if scalar_bar else ())]
        require(len(factors) == 4 and all(sum(bool(test(f)) for f in factors) == 1 for test in wanted),
                f'{model} {vertex}: incorrect coupling/CG/scalar orientation or extra factors')
    spin = verify_spinor(ncm(factors), model, vertex, hc)
    return {'coefficient': f'{coeff.numerator}/{coeff.denominator}',
            **spin}


def verify_record(record: dict, model: str) -> dict:
    require(model in ('A_p0','B_m1','C_m1'), 'Unknown representative model')
    require(isinstance(record, dict) and record.get('model') == model,
            'Wrong model in native evidence')
    require(record.get('format') == 't3-b1-native-restricted-action-v1' and
            record.get('status') == 'NativeActionTermsIsolated_PhysicalDictionaryPending',
            'Wrong native format or status')
    require(record.get('physical_source_certified') is False and
            record.get('source_gammaCC_concrete_matrix_measured') is False and
            record.get('cg_component_coefficients_measured_here') is False,
            'Incorrect source certificate or scope claim')
    tests = record.get('tests')
    require(isinstance(tests, dict) and tests and all(v is True for v in tests.values()),
            'Native Matchete checks not successful')
    flag_names = {'T3_EXPERIMENTAL_GENERAL_FERMIONS',
                  'T3_EXPERIMENTAL_GENERAL_QUADRATICS',
                  'T3_EXPERIMENTAL_GENERAL_QUARTICS',
                  'T3_EXPERIMENTAL_FLAVOR_MASSES'}
    flags = record.get('effective_flags')
    require(isinstance(flags, dict) and set(flags) == flag_names and
            all(v == '0' for v in flags.values()),
            'Native run not restricted to original model flags')
    masses = record.get('mass_terms_inputform')
    require(isinstance(masses, list) and len(masses) == 1 and record.get('mass_term_count') == 1,
            'Expected exactly one mass operator')
    result = {'mass': verify_term(masses[0], model, 'mass')}
    for vertex in ('y1', 'y2'):
        row = record.get('yukawas', {}).get(vertex)
        require(isinstance(row, dict) and row.get('original_term_count') == 1 and
                row.get('hc_term_count') == 1, 'Native Yukawa original/HC count changed')
        orig, hc = row.get('original_terms_inputform'), row.get('hc_terms_inputform')
        require(isinstance(orig, list) and len(orig) == 1 and isinstance(hc, list) and len(hc) == 1,
                'Native Yukawa term lists changed')
        result[vertex] = {'original': verify_term(orig[0], model, vertex),
                          'hc': verify_term(hc[0], model, vertex, True)}
    require(record.get('source_action_LInt_inputform') and record.get('source_action_LFree_inputform'),
            'Full native action not recorded')
    # The native probe already performed symbolic term extraction from LInt/LFree.
    # Here ensure exact terms are present in the separately exported full actions.
    lint, lfree = record['source_action_LInt_inputform'], record['source_action_LFree_inputform']
    require(all(term in lint for row in record['yukawas'].values()
                for key in ('original_terms_inputform', 'hc_terms_inputform')
                for term in row[key]), 'Extracted Yukawa term missing from full LInt text')
    # LFree's top-level mass may be printed without the same outer parentheses;
    # verify its MF symbol, while term-extraction provenance is native-verified.
    require('Coupling[MF,' in lfree, 'Extracted MF missing from full LFree')
    return {'model': model, 'native_structure_verified': True,
            'couplings_and_chiral_projectors': result,
            'absolute_GammaCC_matrix_orientation_measured': False,
            'B2_relative_CG_phases_measured_here': False,
            'physical_source_certified': False}
