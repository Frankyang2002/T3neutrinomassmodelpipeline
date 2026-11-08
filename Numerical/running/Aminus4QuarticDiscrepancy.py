"""Pure quadratic-polynomial diagnostics for the A(-4) RGBeta/reference mismatch.

These identities concern *only* the scalar-loop contribution in the declared
12-real-scalar reference convention. They are not physical Matchete coupling
conversions and must not be applied to the production beta source.

Coordinates are ordered as in Aminus4EFT1Audit.SCALAR_NAMES.
"""
from __future__ import annotations

from fractions import Fraction
from itertools import combinations
from typing import Callable, Mapping, Any

import numpy as np

COORDINATES = (
    'lambdaH', 'lambdaS1', 'lambdaS2', 'lambdaS2Adj',
    'lambdaH1', 'lambdaH2', 'lambdaH2Adj', 'lambda12',
    'lambdaT3.real', 'lambdaT3.imag',
)

# Native RGBeta minus the fourth-derivative scalar-loop reference, derived
# from the archived 55-probe source audit, not a physical counterterm.
EXPECTED_MONOMIALS = {
    ('lambdaS2Adj', 'lambdaS2Adj'): {'lambdaS2': Fraction(8),
                                     'lambdaS2Adj': Fraction(14)},
    ('lambdaS2Adj', 'lambdaH2Adj'): {'lambdaH2Adj': Fraction(10)},
    ('lambdaH2Adj', 'lambdaH2Adj'): {'lambdaH': Fraction(2),
                                     'lambdaS2Adj': Fraction(2),
                                     'lambdaH2': Fraction(4)},
    ('lambdaT3.real', 'lambdaT3.real'): {'lambdaH2Adj': Fraction(-16, 3)},
    ('lambdaT3.imag', 'lambdaT3.imag'): {'lambdaH2Adj': Fraction(-16, 3)},
}

# Equivalent compact polynomial, using the installed Python RGBeta InputForm
# evaluator's Bar[complex_scalar] syntax. The sign is NATIVE MINUS REFERENCE.
DELTA_INPUTFORM = {
    'lambdaH': '2*lambdaH2Adj^2',
    'lambdaS2': '8*lambdaS2Adj^2',
    'lambdaS2Adj': '14*lambdaS2Adj^2 + 2*lambdaH2Adj^2',
    'lambdaH2': '4*lambdaH2Adj^2',
    'lambdaH2Adj': '10*lambdaS2Adj*lambdaH2Adj - (16/3)*Bar[lambdaT3]*lambdaT3',
}


def coordinates_of(point: Mapping[str, Any]) -> np.ndarray:
    """Ten *real* coordinates; accepts a complex or JSON complex lambdaT3."""
    value = point['lambdaT3']
    if isinstance(value, Mapping):
        value = complex(value['real'], value['imag'])
    value = complex(value)
    return np.asarray([float(point[n]) for n in COORDINATES[:8]] +
                      [float(value.real), float(value.imag)], dtype=float)


def eval_polynomial(values: np.ndarray, monomials: Mapping) -> np.ndarray:
    """Return ten beta-component discrepancies in COORDINATES order."""
    x = np.asarray(values, dtype=float)
    if x.shape != (10,) or not np.all(np.isfinite(x)):
        raise ValueError('Expected ten finite real quartic coordinates')
    lookup = {v:i for i,v in enumerate(COORDINATES)}
    beta = np.zeros(10, dtype=float)
    for (a,b), row in monomials.items():
        for dst, coeff in row.items():
            beta[lookup[dst]] += float(coeff) * x[lookup[a]] * x[lookup[b]]
    return beta


def discrepancy_vector(comparison: Mapping) -> np.ndarray:
    rows = comparison['coupling_beta_comparison']
    if set(rows) != set(COORDINATES):
        raise ValueError('Comparison lacks ten scalar reference coordinates')
    return np.asarray([
        float(rows[n]['native_rgbeta_scalar_loop']) -
        float(rows[n]['reference_scalar_loop']) for n in COORDINATES], float)


def extract_monomials(
    run_probe: Callable[[tuple[str, ...]], Mapping], *, tolerance: float = 2e-8,
) -> dict:
    """Polarise all 55 probes into 10 squares + 45 cross coefficients.

    run_probe(names) evaluates the *actual* native/reference scalar-only beta
    at the point where every named coordinate equals one and others are zero.
    No source status flags or fitted coupling redefinitions are used.
    """
    singles = {}
    pair_rows = {}
    comparisons = {}
    for name in COORDINATES:
        row = run_probe((name,))
        comparisons[(name,)] = row
        singles[name] = discrepancy_vector(row)
    for a, b in combinations(COORDINATES, 2):
        row = run_probe((a, b))
        comparisons[(a,b)] = row
        pair_rows[(a,b)] = discrepancy_vector(row) - singles[a] - singles[b]
    monomials = {}
    for a, d in singles.items():
        coeffs = {dst: float(v) for dst,v in zip(COORDINATES,d) if abs(v)>tolerance}
        if coeffs: monomials[(a,a)] = coeffs
    for k,d in pair_rows.items():
        coeffs = {dst: float(v) for dst,v in zip(COORDINATES,d) if abs(v)>tolerance}
        if coeffs: monomials[k] = coeffs
    reconstructed_errors = []
    for probe, row in comparisons.items():
        vals = np.asarray([1.0 if n in probe else 0.0 for n in COORDINATES])
        observed = discrepancy_vector(row)
        predicted = eval_polynomial(vals,monomials)
        reconstructed_errors.append(float(np.max(np.abs(observed-predicted))))
    failed = sum(not c.get('reference_and_native_scalar_loop_match', False)
                 for c in comparisons.values())
    closure = max(float(c['scalar_loop_projection_closure_relative'])
                  for c in comparisons.values())
    return {'monomials':monomials, 'probe_count':len(comparisons),
            'failing_probes':failed,
            'max_scalar_reference_closure_relative':closure,
            'max_probed_reconstruction_abs':max(reconstructed_errors),
            'raw_comparisons':comparisons}


def assert_expected_five(monomials: Mapping, *, tolerance: float = 2e-8) -> float:
    if set(monomials) != set(EXPECTED_MONOMIALS):
        raise ValueError(f'Five-monomial discrepancy support changed: '
                         f'found={set(monomials)}, expected={set(EXPECTED_MONOMIALS)}')
    max_abs = 0.0
    for key, target in EXPECTED_MONOMIALS.items():
        row = monomials[key]
        if set(row) != set(target):
            raise ValueError(f'Discrepancy beta channel support changed for {key}')
        for n, value in target.items():
            err = abs(float(row[n]) - float(value))
            max_abs = max(max_abs, err)
            if err > tolerance:
                raise ValueError(f'Discrepancy coefficient changed: {key} -> {n}: '
                                 f'{row[n]} != {value}')
    return max_abs


def monomials_to_json(monomials: Mapping) -> list[dict]:
    out = []
    for pair, terms in monomials.items():
        row = {'couplings':list(pair), 'beta_differences':{}}
        for name, val in terms.items():
            rational=Fraction(float(val)).limit_denominator(4096)
            row['beta_differences'][name] = {
                'numerical':float(val),
                'rational_if_within_tolerance':str(rational)
                    if abs(float(rational)-float(val))<2e-8 else None,
            }
        out.append(row)
    return out
