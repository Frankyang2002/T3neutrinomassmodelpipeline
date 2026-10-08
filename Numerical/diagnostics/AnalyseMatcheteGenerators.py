"""Diagnose Matchete's SU(2) generator ordering and diagonalize Cartan matrices.

Run: python -m Numerical.diagnostics.AnalyseMatcheteGenerators

This does not certify that the generator basis equals the exported CG tensor basis.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as sp

DEFAULT_INPUT = Path('Reports/output/full/multiplicity/matchete_generator_probe.json')
DEFAULT_OUTPUT = Path('Reports/output/full/multiplicity/matchete_generator_analysis.json')


def _matrix(rows):
    return sp.Matrix([[sp.sympify(entry, locals={'I': sp.I}) for entry in row] for row in rows])


def _string_matrix(matrix):
    return [[str(sp.simplify(value)) for value in matrix.row(i)] for i in range(matrix.rows)]


def analyse_representation(record):
    """Use the directly observed Matchete ordering (G1,G3,G2)."""
    dimension = record['dimension']
    raw = [_matrix(m) for m in record['generators_input_form']]
    if len(raw) != 3 or any(m.shape != (dimension, dimension) for m in raw):
        raise ValueError(f'Invalid generator dimensions for d={dimension}')
    t1, t2, t3 = raw[0], raw[2], raw[1]
    commutators = [
        sp.simplify(t1*t2 - t2*t1 - sp.I*t3),
        sp.simplify(t2*t3 - t3*t2 - sp.I*t1),
        sp.simplify(t3*t1 - t1*t3 - sp.I*t2),
    ]
    algebra_ok = all(c == sp.zeros(dimension) for c in commutators)
    hermitian_ok = all(sp.simplify(t-t.conjugate().T) == sp.zeros(dimension) for t in (t1,t2,t3))
    casimir = sp.simplify(t1*t1+t2*t2+t3*t3)
    j = sp.Rational(dimension-1,2)
    casimir_ok = casimir == j*(j+1)*sp.eye(dimension)
    # Eigenvectors of a Hermitian Cartan matrix form the physical weight basis.
    eigenvectors = []
    for eigenvalue, multiplicity, vectors in t3.eigenvects():
        if multiplicity != len(vectors):
            raise ValueError('Incomplete Cartan eigensystem')
        for v in vectors:
            eigenvectors.append((sp.simplify(eigenvalue), sp.simplify(v/sp.sqrt((v.conjugate().T*v)[0]))))
    eigenvectors.sort(key=lambda pair: float(sp.re(pair[0])), reverse=True)
    u = sp.Matrix.hstack(*(vector for _, vector in eigenvectors))
    diagonal = sp.simplify(u.conjugate().T*t3*u)
    unitary_ok = sp.simplify(u.conjugate().T*u-sp.eye(dimension)) == sp.zeros(dimension)
    spectrum = [str(value) for value, _ in eigenvectors]
    expected = [str(j-i) for i in range(dimension)]
    verified = algebra_ok and hermitian_ok and casimir_ok and unitary_ok and spectrum == expected
    return {
        'dimension': dimension,
        'status': 'generator_algebra_verified' if verified else 'validation_failed',
        'matchete_to_standard_order_zero_based': [0,2,1],
        'hermitian': hermitian_ok,
        'positive_su2_commutators': algebra_ok,
        'casimir': casimir_ok,
        'diagonalizer_unitary': unitary_ok,
        'cartan_spectrum_descending': spectrum,
        'T3_in_original_matchete_basis': _string_matrix(t3),
        'U_columns_are_weight_eigenvectors': _string_matrix(u),
        'Udag_T3_U': _string_matrix(diagonal),
        'U_phase_is_not_fixed_by_diagonalization': True,
        'CG_tensor_covariance_verified': False,
    }


def analyse(payload):
    entries = [analyse_representation(rec) for rec in payload['representations']]
    return {
        'status': 'generator_algebra_verified_CG_basis_unverified' if all(e['status']=='generator_algebra_verified' for e in entries) else 'validation_failed',
        'representations': entries,
        'warning': 'The generator matrices have not been checked against Matchete CG tensor invariance. Eigenvector phases are arbitrary. Do not assign physical loop-channel weights or coupling compensation yet.',
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=DEFAULT_INPUT)
    parser.add_argument('--output', type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    payload = json.loads(args.input.read_text(encoding='utf-8'))
    result = analyse(payload)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print('Status:',result['status'])
    for rec in result['representations']:
        print('dimension=',rec['dimension'],'T3 spectrum=',rec['cartan_spectrum_descending'],'algebra=',rec['status'])
    print('Output:',args.output)
    if result['status'] == 'validation_failed':
        raise SystemExit(1)


if __name__ == '__main__':
    main()
