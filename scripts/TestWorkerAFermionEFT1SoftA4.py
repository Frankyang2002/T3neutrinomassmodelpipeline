"""Tests A(-4) candidate scalar bubble and direct source provenance.

Uses the two original evidence ZIPs; no writes to repository or production data.
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from WorkerAFermionEFT1SoftA4 import evaluate, scalar_bubble_finite


def check(sources_zip: Path, soft_zip: Path) -> None:
    a, b, q = 1.4e6, 1.2e6, 1e10
    direct = 1 + (a * math.log(q / a) - b * math.log(q / b)) / (a - b)
    np.testing.assert_allclose(scalar_bubble_finite(a, b, q), direct, rtol=2e-14)
    np.testing.assert_allclose(scalar_bubble_finite(a, a, q), math.log(q / a), rtol=2e-14)
    for t in (1e-8, -1e-8, 1e-6, -1e-6):
        approx = scalar_bubble_finite(a * (1 + t), a, q)
        exact = 1 + ((1 + t) * math.log(q / (a * (1 + t))) - math.log(q / a)) / t
        np.testing.assert_allclose(approx, exact, rtol=2e-7, atol=2e-7)
    data = evaluate(sources_zip, soft_zip, 1e5)
    assert data['soft_symbolic_source']['symbolic_reconstruction_exact'] is True
    assert len(data['direct_hard_source_indices']) == 2
    assert data['fermion_threshold_finite_wilson_array_available'] is False
    assert data['direct_ordered_pole_sum_frobenius_GeV_inv'] == 0.0
    np.testing.assert_allclose(data['hard_direct_pole_frobenius_GeV_inv'],
                               2.445231389e-16, rtol=1e-9)
    # Explicit log(mu^2) terms cancel between source hard and EFT1 bubble,
    # with masses, flavour couplings and lambda held fixed (not an RGE test).
    def complex_decode(value):
        return np.asarray(value['real']) + 1j * np.asarray(value['imag'])
    high = data
    # Re-use the fixed initial values, applying the scalar-loop scale shift
    # algebraically because the benchmark is anchored at mu_F = 1e5 GeV.
    for factor in (0.5, 2.0, 10.0):
        shift = 2 * math.log(factor)
        pole = complex_decode(high['hard_direct_pole_GeV_inv'])
        hard = complex_decode(high['hard_direct_finite_GeV_inv']) + pole * shift
        soft = complex_decode(high['soft_candidate_finite_GeV_inv']) - pole * shift
        paired = complex_decode(high['paired_direct_finite_diagnostic_GeV_inv'])
        np.testing.assert_allclose(hard + soft, paired, atol=1e-27, rtol=1e-12)
    print('PASS source reconstruction, A4 numerical benchmark, mass limit, direct pole pairing, fixed-coupling scale cancellation')


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sources-zip', type=Path, required=True)
    p.add_argument('--soft-zip', type=Path, required=True)
    a = p.parse_args()
    check(a.sources_zip, a.soft_zip)


if __name__ == '__main__':
    main()
