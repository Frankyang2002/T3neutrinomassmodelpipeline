"""Numerical coupling dictionary from *explicit* native component polynomials.

This module never identifies Matchete and RGBeta coupling symbols by name.
It solves V_M(x;m)=V_R(x;r) at fixed scalar field coordinates for random
complex fields x, with r = map @ m.  Both inputs must provide the FULL
native polynomial tensor, prefactor, conjugation and field ordering.  A
native RGBeta *declaration string* is not a component tensor.

The result is a numerical field-basis dictionary only if the two inputs have
already been exported in the same complex scalar coordinates.  Spinor factors
and masses require separate certificates, and no physical claim is implied.
"""
from __future__ import annotations

from dataclasses import dataclass
import ast
import itertools
import math
import re
from typing import Any, Mapping, Sequence

import numpy as np


class TensorMapError(ValueError):
    pass


def numeric_component(value: Any) -> complex:
    """Parse restricted Mathematica InputForm numbers; never eval arbitrary code."""
    if isinstance(value, (int, float, complex)) and not isinstance(value, bool):
        return complex(value)
    if not isinstance(value, str):
        raise TensorMapError(f"Non-numeric component: {value!r}")
    s = value.strip().replace(' ', '').replace('*^', 'e')
    s = re.sub(r'Sqrt\[([^\[\]]+)\]', r'sqrt(\1)', s)
    s = re.sub(r'Rational\[(-?\d+),(-?\d+)\]', r'(\1/\2)', s)
    s = s.replace('^', '**')
    # InputForm can write imaginary literals without an explicit multiply.
    s = re.sub(r'(?<=[0-9)])I\b', '*I', s)
    s = re.sub(r'(?<=\))I\b', '*I', s)
    try:
        node = ast.parse(s, mode='eval')
    except SyntaxError as exc:
        raise TensorMapError(f"Unsupported tensor component {value!r}") from exc

    def run(n: ast.AST) -> complex:
        if isinstance(n, ast.Expression): return run(n.body)
        if isinstance(n, ast.Constant) and type(n.value) in (int, float, complex):
            return complex(n.value)
        if isinstance(n, ast.Name) and n.id == 'I': return 1j
        if isinstance(n, ast.Name) and n.id == 'Pi': return complex(math.pi)
        if isinstance(n, ast.UnaryOp) and isinstance(n.op, ast.USub): return -run(n.operand)
        if isinstance(n, ast.UnaryOp) and isinstance(n.op, ast.UAdd): return run(n.operand)
        if isinstance(n, ast.BinOp):
            a, b = run(n.left), run(n.right)
            if isinstance(n.op, ast.Add): return a + b
            if isinstance(n.op, ast.Sub): return a - b
            if isinstance(n.op, ast.Mult): return a * b
            if isinstance(n.op, ast.Div): return a / b
            if isinstance(n.op, ast.Pow):
                if abs(b.imag) > 1e-15 or abs(b.real) > 20:
                    raise TensorMapError('Unsupported exponent')
                return a ** b.real
        if (isinstance(n, ast.Call) and isinstance(n.func, ast.Name)
                and n.func.id == 'sqrt' and len(n.args) == 1 and not n.keywords):
            return np.sqrt(run(n.args[0]))
        raise TensorMapError(f"Unsupported numeric expression {value!r}")

    z = complex(run(node))
    if not np.isfinite(z): raise TensorMapError('Non-finite component')
    return z


def canonical_species(name: str) -> tuple[str, bool]:
    """Canonical species names: H/S1/S2 and Bar[H]/Bar[S1]/Bar[S2]."""
    t = name.replace(' ', '')
    conjugated = False
    # Actual T3ExtendedSpeciesList in GeneralScalarOperators.wl uses Hdag,
    # S1dag, S2dag (not Bar[...] strings).
    if t in ('Hdag', 'S1dag', 'S2dag'):
        return t[:-3], True
    for left, right in [('Bar[', ']'), ('Conjugate[', ']')]:
        if t.startswith(left) and t.endswith(right):
            t, conjugated = t[len(left):-len(right)], True
            break
    if t not in ('H', 'S1', 'S2'):
        raise TensorMapError(f'Unknown species {name!r}')
    return t, conjugated


@dataclass(frozen=True)
class Polynomial:
    symbol: str
    species: tuple[tuple[str, bool], ...]
    tensor: np.ndarray
    prefactor: complex
    add_hc: bool
    # If True, mathematical potential is the negative of source Lagrangian.
    source_is_lagrangian: bool

    @staticmethod
    def from_record(record: Mapping[str, Any]) -> 'Polynomial':
        if not isinstance(record.get('symbol'), str):
            raise TensorMapError('Missing tensor symbol')
        raw_species = record.get('species')
        dims = record.get('dimensions')
        if not isinstance(raw_species, list) or len(raw_species) != 4:
            raise TensorMapError('Exactly four scalar species required')
        if not isinstance(dims, list) or len(dims) != 4 or any(
                type(x) is not int or x < 1 for x in dims):
            raise TensorMapError('Invalid four scalar tensor dimensions')
        values = record.get('components')
        if not isinstance(values, list) or len(values) != math.prod(dims):
            raise TensorMapError('Incomplete component tensor')
        tensor = np.array([numeric_component(v) for v in values], dtype=complex).reshape(dims)
        prefactor = numeric_component(record.get('potential_prefactor', '1'))
        phase = numeric_component(record.get('phase', '1'))
        sign = -1 if record.get('source_is_lagrangian', False) else 1
        if record.get('add_hc') not in (True, False):
            raise TensorMapError('Explicit add_hc true/false required')
        if not np.isfinite(tensor).all() or np.max(np.abs(tensor)) == 0:
            raise TensorMapError(f"Zero or invalid tensor {record['symbol']}")
        return Polynomial(record['symbol'], tuple(map(canonical_species, raw_species)),
                          tensor, prefactor * phase * sign, record['add_hc'],
                          bool(record.get('source_is_lagrangian', False)))

    def evaluate(self, fields: Mapping[str, np.ndarray]) -> float:
        factors = []
        for (species, conjugated), dim in zip(self.species, self.tensor.shape):
            if species not in fields: raise TensorMapError(f'Missing field {species}')
            x = np.asarray(fields[species], dtype=complex)
            if x.shape != (dim,): raise TensorMapError(f'Bad field dimension {species}')
            factors.append(x.conj() if conjugated else x)
        val = np.einsum('abcd,a,b,c,d->', self.tensor, *factors, optimize=True)
        val *= self.prefactor
        if self.add_hc: val += val.conjugate()
        if abs(val.imag) > 1e-9 * max(1.0, abs(val)):
            raise TensorMapError(f'{self.symbol} is not a real potential polynomial')
        return float(val.real)


def _random_fields(dims: Mapping[str, int], rng: np.random.Generator) -> dict[str, np.ndarray]:
    return {s: (rng.normal(size=d) + 1j*rng.normal(size=d))/math.sqrt(2)
            for s, d in dims.items()}


def solve_dictionary(
    matchete_records: Sequence[Mapping[str, Any]],
    rgbeta_records: Sequence[Mapping[str, Any]],
    *, samples: int = 256, seed: int = 1907, tolerance: float = 1e-8,
) -> dict[str, Any]:
    """Solve map r=A@m using field-polynomial equality for all original symbols.

    Assumes same canonical complex-component basis, not just same SU(2) irrep.
    Rank deficiency is a BLOCKER rather than a guessed gauge-invariant mapping.
    All Matchete/RGBeta coefficients are treated as real parameters here.
    Complex parameter families need separate Re/Im component records.
    """
    ms = [Polynomial.from_record(q) for q in matchete_records]
    rs = [Polynomial.from_record(q) for q in rgbeta_records]
    if not ms or not rs: raise TensorMapError('No component polynomials')
    if len(set(x.symbol for x in ms)) != len(ms) or len(set(x.symbol for x in rs)) != len(rs):
        raise TensorMapError('Repeated coupling symbols')
    dims: dict[str, int] = {}
    for p in itertools.chain(ms, rs):
        for (species, _), d in zip(p.species, p.tensor.shape):
            if species in dims and dims[species] != d:
                raise TensorMapError(f'{species} component-basis mismatch')
            dims[species] = d
    rng = np.random.default_rng(seed)
    xs = [_random_fields(dims, rng) for _ in range(samples)]
    vm = np.asarray([[p.evaluate(x) for p in ms] for x in xs], dtype=float)
    vr = np.asarray([[p.evaluate(x) for p in rs] for x in xs], dtype=float)
    # vr @ A = vm; each column of A maps one source real coupling to RGBeta.
    a, _, rank_r, singular = np.linalg.lstsq(vr, vm, rcond=1e-12)
    rank_m = int(np.linalg.matrix_rank(vm, tol=1e-10 * max(1.,np.linalg.norm(vm,ord=2))))
    residual = float(np.linalg.norm(vr @ a - vm) / max(np.linalg.norm(vm), 1e-300))
    holdout = [_random_fields(dims,rng) for _ in range(max(32,samples//4))]
    mh = np.asarray([[p.evaluate(x) for p in ms] for x in holdout], dtype=float)
    rh = np.asarray([[p.evaluate(x) for p in rs] for x in holdout], dtype=float)
    hold_resid = float(np.linalg.norm(rh @ a - mh) / max(np.linalg.norm(mh),1e-300))
    independent = rank_r == len(rs) and rank_m == len(ms)
    exact_space = (rank_r == rank_m == len(ms) == len(rs) and
                   residual < tolerance and hold_resid < tolerance)
    return {
        'status': 'NumericalPolynomialMapComplete' if exact_space else 'Blocked_PolynomialMapIncomplete',
        'matchete_symbols': [p.symbol for p in ms],
        'rgbeta_symbols': [p.symbol for p in rs],
        'map_rgbeta_from_matchete': a.tolist(),
        'rank_matchete': rank_m, 'rank_rgbeta': int(rank_r),
        'independent_bases': bool(independent),
        'fit_relative_residual': residual,
        'holdout_relative_residual': hold_resid,
        'samples': samples, 'seed': seed,
        'scalar_coordinates_certified_by_this_step': False,
        'spinor_normalisation_certified_by_this_step': False,
        'note': 'Use only if both native tensor exports explicitly share the same component field basis.'
    }


def scan_hermitian_phases(records: Sequence[Mapping[str, Any]], *,
                          samples: int = 48, seed: int = 1403,
                          tolerance: float = 1e-9) -> list[dict[str, Any]]:
    """Detect quartics killed by PlusHc with an unphased real coupling.

    Re(O) identically zero but Im(O) nonzero signals an anti-Hermitian
    contraction requiring a factor I inside PlusHc. This is a *numerical*
    potential-level diagnostic and never silently rewrites source couplings.
    """
    rng = np.random.default_rng(seed)
    output = []
    for record in records:
        # Use the native unphased complex invariant; inspect it before + h.c.
        p = Polynomial.from_record({**record, 'add_hc': False,
                                    'phase': '1', 'source_is_lagrangian':False})
        dims = {}
        for (species, _), dim in zip(p.species, p.tensor.shape):
            if species in dims and dims[species] != dim:
                raise TensorMapError(f'Inconsistent scalar dimension {species}')
            dims[species] = dim
        vals = []
        for _ in range(samples):
            fields = _random_fields(dims,rng)
            legs = [np.conjugate(fields[name]) if barred else fields[name]
                    for name,barred in p.species]
            z = p.prefactor * np.einsum('abcd,a,b,c,d->',p.tensor,*legs,optimize=True)
            vals.append(complex(z))
        scale = max(max(abs(z) for z in vals), 1e-300)
        rr = max(abs(z.real) for z in vals) / scale
        ii = max(abs(z.imag) for z in vals) / scale
        if rr < tolerance and ii > tolerance: classification = 'AntiHermitianRequiresI'
        elif ii < tolerance and rr > tolerance: classification = 'HermitianPhase1'
        elif rr < tolerance and ii < tolerance: classification = 'ZeroPolynomialBlocked'
        else: classification = 'ComplexInvariantRequiresExplicitHc'
        output.append({'symbol':p.symbol, 'classification':classification,
                       'unphased_real_part_ratio':float(rr),
                       'unphased_imag_part_ratio':float(ii),
                       'samples':samples})
    return output
