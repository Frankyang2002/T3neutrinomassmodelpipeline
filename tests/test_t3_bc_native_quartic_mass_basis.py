"""Source-tensor rotation checks; no live Matchete installation needed.

Requires the repository's saved Stage22 native quartics and B/C mass-basis
reports. Missing upstream reports are a failed preflight, not a skipped test.
"""
from copy import deepcopy

import numpy as np
import pytest

from Numerical.experimental.GeneralT3MixedScalarMassBasis import MixedDoubletMassBasis
from scripts.AuditT3BCNativeQuarticMassBasis import (
    EXPECTED, MASS_BASIS, NATIVE, load_basis, native_tensors, read_json,
    rotate_tensor, shared_independent_block, verify_polynomial,
)


@pytest.mark.parametrize('model', ('B_-1', 'C_-1'))
def test_saved_native_quartics_rotate_and_preserve_independent_contractions(model):
    source = native_tensors(read_json(NATIVE), model)
    assert set(source) == set(EXPECTED)
    basis = load_basis(MASS_BASIS / f'{model.replace("-", "m")}_final_mass_basis.json', model)
    rotated = {}
    for name, (fields, tensor) in source.items():
        rotated[name] = rotate_tensor(tensor, fields, basis.U)
        assert len(rotated[name]) >= 1
        assert verify_polynomial(fields, tensor, rotated[name], basis, trials=5) < 3e-12
    evidence = shared_independent_block(rotated['lambdaT3ExtraQ110101I1'],
                                        rotated['lambdaT3ExtraQ110101I2'])
    assert evidence['rank'] == 2
    assert evidence['small_to_large_singular_value_ratio'] > 1e-9


def test_complex_bilinear_phase_native_quartic_rotation():
    source = native_tensors(read_json(NATIVE), 'B_-1')
    basis = MixedDoubletMassBasis.build('B_-1', 1.21e6, 1.43e6, 22000+15000j)
    assert not basis.tachyonic
    for name in ('lambdaT3ExtraQ110101I1', 'lambdaT3ExtraQ110101I2',
                 'lambdaT3ExtraQ020020I1'):
        fields, tensor = source[name]
        rotated = rotate_tensor(tensor, fields, basis.U)
        assert verify_polynomial(fields, tensor, rotated, basis, trials=8, seed=1001) < 3e-12


def test_native_inventory_fails_closed_when_one_tensor_is_missing():
    bad = deepcopy(read_json(NATIVE))
    bad['models']['B_-1']['tensors'].pop()
    with pytest.raises(ValueError, match='exactly seven'):
        native_tensors(bad, 'B_-1')


def test_independent_contraction_rank_test_detects_duplicate():
    basis = MixedDoubletMassBasis.build('C_-1', 1.2e6, 1.5e6, 5000+9000j)
    fields, tensor = native_tensors(read_json(NATIVE), 'C_-1')[
        'lambdaT3ExtraQ110101I1']
    a = rotate_tensor(tensor, fields, basis.U)
    with pytest.raises(ValueError, match='collapsed'):
        shared_independent_block(a, a)
