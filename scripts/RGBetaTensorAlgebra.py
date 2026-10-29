"""Exact sparse tensor change-of-basis algebra for Stage 22.

No model-specific conventions, numerical tolerances, or identity-map fallback.
The caller must supply tensors in the SAME ordered field/component basis.
"""
from __future__ import annotations

import ast
from collections import Counter
from itertools import product
from typing import Any, Mapping, Sequence

import sympy as sp


class TensorCertificateError(ValueError):
    pass


def exact(value: Any) -> sp.Expr:
    """Parse a small, explicit subset of exact Wolfram InputForm; no eval."""
    if isinstance(value, bool):
        raise TensorCertificateError('Boolean tensor entry')
    if isinstance(value, int):
        return sp.Integer(value)
    if isinstance(value, str):
        text = value.replace('^', '**').replace('[', '(').replace(']', ')').strip()
        try:
            node = ast.parse(text, mode='eval').body
        except (SyntaxError, ValueError) as exc:
            raise TensorCertificateError(f'Invalid exact scalar {value!r}') from exc

        def walk(n: ast.AST) -> sp.Expr:
            if isinstance(n, ast.Constant) and isinstance(n.value, int) and not isinstance(n.value, bool):
                return sp.Integer(n.value)
            if isinstance(n, ast.Name) and n.id == 'I':
                return sp.I
            if isinstance(n, ast.UnaryOp) and isinstance(n.op, (ast.USub, ast.UAdd)):
                x = walk(n.operand)
                return -x if isinstance(n.op, ast.USub) else x
            if isinstance(n, ast.BinOp) and isinstance(n.op, (ast.Add, ast.Sub, ast.Mult, ast.Div, ast.Pow)):
                a, b = walk(n.left), walk(n.right)
                if isinstance(n.op, ast.Add): return a+b
                if isinstance(n.op, ast.Sub): return a-b
                if isinstance(n.op, ast.Mult): return a*b
                if isinstance(n.op, ast.Div):
                    if b == 0: raise TensorCertificateError('Division by zero')
                    return a/b
                if isinstance(n.op, ast.Pow):
                    if not b.is_Rational: raise TensorCertificateError('Nonrational exponent')
                    return a**b
            if isinstance(n, ast.Call) and isinstance(n.func, ast.Name) and not n.keywords:
                a = [walk(x) for x in n.args]
                if n.func.id in ('Sqrt','sqrt') and len(a) == 1: return sp.sqrt(a[0])
                if n.func.id == 'Rational' and len(a) == 2 and all(x.is_Integer for x in a):
                    return sp.Rational(a[0], a[1])
                if n.func.id == 'Complex' and len(a) == 2: return a[0]+sp.I*a[1]
            raise TensorCertificateError(f'Unsupported exact scalar {value!r}')

        ans = sp.simplify(walk(node))
        if ans.free_symbols or ans.has(sp.zoo, sp.nan, sp.oo, -sp.oo):
            raise TensorCertificateError(f'Nonfinite/nonconstant tensor entry {value!r}')
        return ans
    raise TensorCertificateError(f'Tensor scalar not exact string/integer: {value!r}')


def encoding(x: sp.Expr) -> str:
    return str(sp.simplify(x))


def normalized_record(record: Mapping[str, Any], label: str) -> dict[str, Any]:
    """Validate one fully explicit native basis tensor record."""
    if not isinstance(record, Mapping):
        raise TensorCertificateError(f'{label}: record not an object')
    counts = record.get('counts')
    if not (isinstance(counts, list) and len(counts) == 6 and
            all(isinstance(x, int) and x >= 0 for x in counts) and sum(counts) == 4):
        raise TensorCertificateError(f'{label}: invalid scalar field counts')
    dims = record.get('dimensions')
    if not (isinstance(dims, list) and len(dims) == 4 and
            all(isinstance(x, int) and x in (1, 2, 3) for x in dims)):
        raise TensorCertificateError(f'{label}: invalid tensor dimensions')
    vals = record.get('components')
    size = 1
    for d in dims: size *= d
    if not isinstance(vals, list) or len(vals) != size:
        raise TensorCertificateError(f'{label}: expected {size} components, got {len(vals) if isinstance(vals,list) else type(vals)}')
    comp = sp.Matrix([exact(x) for x in vals])
    if not any(comp):
        raise TensorCertificateError(f'{label}: identically zero tensor')
    if not isinstance(record.get('symbol'), str) or not record['symbol']:
        raise TensorCertificateError(f'{label}: absent coupling name')
    return {'counts': tuple(counts), 'dimensions': tuple(dims),
            'symbol': record['symbol'], 'components': comp}


def tensor_rank(columns: Sequence[sp.Matrix]) -> int:
    if not columns: return 0
    return int(sp.Matrix.hstack(*columns).rank())


def same_species_symmetry(record: Mapping[str, Any]) -> None:
    """Check Bose exchange invariance for every pair of identical complex species."""
    n = normalized_record(record, str(record.get('symbol', 'unnamed')))
    species = [a for a, count in enumerate(n['counts']) for _ in range(count)]
    dims = n['dimensions']
    for i in range(4):
        for j in range(i+1, 4):
            if species[i] != species[j]: continue
            if dims[i] != dims[j]: raise TensorCertificateError('identical species, unequal dimensions')
            for idx in product(*(range(d) for d in dims)):
                swapped = list(idx)
                swapped[i], swapped[j] = swapped[j], swapped[i]
                def offset(t: Sequence[int]) -> int:
                    k = 0
                    for x, d in zip(t, dims): k = k*d + x
                    return k
                if sp.simplify(n['components'][offset(idx)]-n['components'][offset(swapped)]) != 0:
                    raise TensorCertificateError(f"{n['symbol']}: nonsymmetric bosonic components ({i},{j})")


def solve_basis(mat: Sequence[Mapping[str, Any]], rgb: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    """Solve R = M A exactly. Coupling vectors satisfy m = A r.

    No prefactors or complex conjugation are silently applied here. Both input
    lists must already be in one common field/component order.
    """
    if not mat or not rgb:
        raise TensorCertificateError('Empty tensor block')
    mm = [normalized_record(x,'Matchete') for x in mat]
    rr = [normalized_record(x,'RGBeta') for x in rgb]
    c0, d0 = mm[0]['counts'], mm[0]['dimensions']
    for x in (*mm, *rr):
        if x['counts'] != c0 or x['dimensions'] != d0:
            raise TensorCertificateError('Incompatible field counts or tensor component order')
    if len(mm) != len(rr):
        raise TensorCertificateError(f'Incompatible basis sizes ({len(mm)} vs {len(rr)})')
    if len({x['symbol'] for x in mm}) != len(mm) or len({x['symbol'] for x in rr}) != len(rr):
        raise TensorCertificateError('Duplicate coupling names')
    for r in (*mat, *rgb): same_species_symmetry(r)
    M = sp.Matrix.hstack(*(x['components'] for x in mm))
    R = sp.Matrix.hstack(*(x['components'] for x in rr))
    n = len(mm)
    if M.rank() != n or R.rank() != n:
        raise TensorCertificateError(f'Rank deficient tensor bases (Matchete {M.rank()}, RGBeta {R.rank()}, needed {n})')
    sols = []
    for j in range(n):
        try:
            answer, params = M.gauss_jordan_solve(R[:,j])
        except ValueError as exc:
            raise TensorCertificateError('RGBeta tensor outside Matchete span') from exc
        if params.rows:
            raise TensorCertificateError('Underdetermined tensor transformation')
        sols.append(answer)
    A = sp.Matrix.hstack(*sols)
    if A.det() == 0 or any(sp.simplify(x) != 0 for x in (M*A-R)):
        raise TensorCertificateError('Noninvertible or nonzero exact tensor residual')
    return {'matchete_symbols':[x['symbol'] for x in mm],
            'rgbeta_symbols':[x['symbol'] for x in rr],
            'counts':list(c0), 'dimensions':list(d0),
            'r_in_m':[[encoding(A[i,j]) for j in range(n)] for i in range(n)],
            'm_in_r':[[encoding((A.inv())[i,j]) for j in range(n)] for i in range(n)],
            'matchete_rank':n,'rgbeta_rank':n,'maximum_exact_component_residual':'0',
            'zero_coupling_recovery': True,
            'tensor_equality':True}


def potential_map(tensor_map: dict[str, Any], mat_prefactors: Sequence[Any], rgb_prefactors: Sequence[Any],
                  mat_to_rgb_orientation: str = 'direct') -> dict[str, Any]:
    """Convert tensor map into coefficient map, with signed Lagrangian factors.

    L_M = sum_i pM_i m_i M_i + h.c.
    L_R = sum_j pR_j r_j R_j + h.c.
    Therefore r = diag(1/pR) A^{-1} diag(pM) m.
    If different source conventions conjugate whole holomorphic block, the
    coupling map is antilinear; this is EXPLICIT, never guessed.
    """
    if mat_to_rgb_orientation not in ('direct','conjugate'):
        raise TensorCertificateError('Orientation must be direct or conjugate')
    n = len(tensor_map['matchete_symbols'])
    if len(mat_prefactors) != n or len(rgb_prefactors) != n:
        raise TensorCertificateError('Missing any potential prefactor')
    pm = [exact(x) for x in mat_prefactors]
    pr = [exact(x) for x in rgb_prefactors]
    if any(sp.simplify(x)==0 for x in pm+pr):
        raise TensorCertificateError('Potential prefactors must be nonzero')
    A = sp.Matrix([[exact(x) for x in row] for row in tensor_map['r_in_m']])
    B = sp.diag(*[1/x for x in pr]) * A.inv() * sp.diag(*pm)
    if B.det() == 0:
        raise TensorCertificateError('Singular physical coupling transformation')
    return {'rgbeta_from_matchete':[[encoding(B[i,j]) for j in range(n)] for i in range(n)],
            'matchete_from_rgbeta':[[encoding((B.inv())[i,j]) for j in range(n)] for i in range(n)],
            'orientation':mat_to_rgb_orientation,
            'coupling_reconstruction': 'r = B*m' if mat_to_rgb_orientation == 'direct' else 'r = B*conjugate(m)',
            'zero_coupling_recovery':True}
