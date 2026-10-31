"""Lossless additive InputForm pole/log/constant splitting for T3 B/C scalar probes.

This is text-level bookkeeping of ALREADY computed Matchete local operator
expressions, not loop renormalization or a map to numerical Wilson tensors.
In particular an omitted 1/eps term is NOT a certified MS-bar counterterm.

Supported audited syntax has top-level additive monomials, epsilon exactly
once in a multiplicative denominator, and exactly one specified mass log.
Reject any other epsilon/log structure rather than guess its finite part.
"""
from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from typing import Iterable
import re

EPSILON = '\u03f5'  # Matchete InputForm: ϵ (not the eps[SU2L] Clebsch)
MASS_LOGS = {
    'Y': 'Log[\u03bcbar2/Coupling[MSY, {}, 0]^2]',
    'X': 'Log[\u03bcbar2/Coupling[MSX, {}, 0]^2]',
}
VALID_PAIRS = {'Y': ('LLHH', 'LLHX', 'LLXX'), 'X': ('LLHH',)}


class AuditError(ValueError):
    """Fail closed rather than silently erase an unrecognized term."""


def _balanced(s: str) -> bool:
    stack: list[str] = []
    pairs = {')': '(', ']': '[', '}': '{'}
    in_string = False
    escaped = False
    for c in s:
        if in_string:
            if escaped:
                escaped = False
            elif c == '\\':
                escaped = True
            elif c == '"':
                in_string = False
        elif c == '"':
            in_string = True
        elif c in '([{':
            stack.append(c)
        elif c in ')]}':
            if not stack or stack.pop() != pairs[c]:
                return False
    return not stack and not in_string


def additive_terms(expr: str) -> tuple[str, ...]:
    """Split only at top-level + and -, preserving Matchete term strings."""
    s = expr.strip()
    if not s:
        raise AuditError('Empty expression')
    if s == '0':
        return ()
    stack: list[str] = []
    partners = {')': '(', ']': '[', '}': '{'}
    in_string = False
    escaped = False
    start = 0
    terms = []
    for idx, c in enumerate(s):
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
            if not stack or stack.pop() != partners[c]:
                raise AuditError(f'Unbalanced expression near offset {idx}')
        elif c in '+-' and not stack and idx > start:
            term = s[start:idx].strip()
            if not term:
                raise AuditError('Empty additive term')
            terms.append(term)
            start = idx
    if stack or in_string:
        raise AuditError('Unbalanced expression')
    final = s[start:].strip()
    if not final:
        raise AuditError('Missing final term')
    terms.append(final)
    if any(t[0] not in '+-' and i > 0 for i, t in enumerate(terms)):
        raise AuditError('Malformed additive sign')
    if any('ϵ' in t and 'ϵ*' not in t for t in terms):
        raise AuditError('Unrecognized epsilon structure')
    return tuple(terms)


def _sum(terms: Iterable[str]) -> str:
    return ' '.join(terms) if (terms := list(terms)) else '0'


@dataclass(frozen=True)
class Split:
    sector: str
    stage: str
    source_sha256: str
    terms: tuple[str, ...]
    pole_terms: tuple[str, ...]
    log_terms: tuple[str, ...]
    constant_terms: tuple[str, ...]
    pole_residue: str
    pole_free: str
    scale_equal_mass: str

    def summary(self) -> dict:
        return {
            'sector': self.sector, 'threshold': self.stage,
            'source_sha256': self.source_sha256,
            'terms_total': len(self.terms),
            'pole_terms': len(self.pole_terms),
            'log_terms': len(self.log_terms),
            'constant_terms': len(self.constant_terms),
            'pole_residue_sha256': sha256(self.pole_residue.encode('utf8')).hexdigest(),
            'pole_free_sha256': sha256(self.pole_free.encode('utf8')).hexdigest(),
            'scale_equal_mass_sha256': sha256(self.scale_equal_mass.encode('utf8')).hexdigest(),
            'pole_free_is_renormalized_wilson': False,
            'physical_matching_certified': False,
        }


def split_inputform(source: str, *, stage: str, sector: str) -> Split:
    if stage not in VALID_PAIRS or sector not in VALID_PAIRS[stage]:
        raise AuditError(f'Unsupported stage/sector: {stage}/{sector}')
    raw = source.strip()
    if not _balanced(raw):
        raise AuditError('Unbalanced InputForm')
    terms = additive_terms(raw)
    if not terms and sector != 'LLHX':
        raise AuditError(f'Zero input for required sector {stage}/{sector}')
    mass_log = MASS_LOGS[stage]
    poles: list[str] = []
    logs: list[str] = []
    constants: list[str] = []
    for idx, t in enumerate(terms):
        n_eps = t.count(EPSILON)
        n_log = t.count('Log[')
        if 'hbar' not in t:
            raise AuditError(f'Unexpected tree term at {stage}/{sector} term {idx}')
        if 'GammaCC' not in t or 'Field[l, Fermion' not in t:
            raise AuditError(f'Unexpected non-Majorana operator at term {idx}')
        if stage == 'X' and 'Field[NewScalarX' in t:
            raise AuditError('X-heavy field retained in LLHH')
        if 'Field[NewScalarY' in t:
            raise AuditError('Y-heavy field retained after matching')
        if sector == 'LLHH' and (
                t.count('Field[H, Scalar') != 2 or
                'Field[NewScalarX' in t):
            raise AuditError(f'Invalid LLHH scalar field inventory at term {idx}')
        if sector == 'LLXX' and (
                t.count('Field[NewScalarX, Scalar') != 2 or
                'Field[H, Scalar' in t):
            raise AuditError(f'Invalid LLXX scalar field inventory at term {idx}')
        if sector == 'LLHX' and (
                t.count('Field[H, Scalar') != 1 or
                t.count('Field[NewScalarX, Scalar') != 1):
            raise AuditError(f'Invalid LLHX field inventory at term {idx}')
        if n_eps:
            if (n_eps != 1 or t.count('ϵ*') != 1 or n_log or
                    not re.search(r'/\([^()]*ϵ\*', t)):
                raise AuditError(f'Unexpected mixed / higher epsilon term at {idx}')
            poles.append(t)
        elif n_log:
            if n_log != 1 or t.count(mass_log) != 1:
                raise AuditError(f'Unexpected log argument at term {idx}')
            logs.append(t)
        else:
            constants.append(t)
    if len(terms) != len(poles) + len(logs) + len(constants):
        raise AuditError('Partition terms not conserved')
    pole_residue = _sum(t.replace('ϵ*', '', 1) for t in poles)
    pole_free = _sum(t for t in terms if 'ϵ' not in t)
    at_mass = _sum(constants)
    if 'ϵ' in pole_residue or 'ϵ' in pole_free or 'Log[' in at_mass:
        raise AuditError('Residual pole or log in split output')
    for candidate in [pole_residue, pole_free, at_mass]:
        if not _balanced(candidate):
            raise AuditError('Split output is not bracket-balanced')
    return Split(sector, stage, sha256(source.encode('utf8')).hexdigest(),
                 terms, tuple(poles), tuple(logs), tuple(constants),
                 pole_residue, pole_free, at_mass)
