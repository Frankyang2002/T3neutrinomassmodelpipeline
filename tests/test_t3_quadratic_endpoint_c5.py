"""Fast isolated guards for conditional EFT1-endpoint C5 bridge."""
from __future__ import annotations

import numpy as np
import pytest

from scripts.RunT3QuadraticEndpointC5 import pair, numerical_c5, require
from Numerical.experimental.GeneralT3MixedScalarMassBasis import MixedDoubletMassBasis


def test_complex_pair_and_fail_closed():
    assert pair({'real':2,'imag':-3}) == 2-3j
    with pytest.raises(ValueError):
        pair({'real':1})
    with pytest.raises(ValueError):
        pair({'real':float('nan'),'imag':0})


def test_mass_basis_is_used_in_numeric_expression():
    b=MixedDoubletMassBasis.build('B_-1',1.8e6,3.2e6,10000+5000j)
    z={'real':0.01,'imag':0.0}
    y=[[z] for _ in range(3)]
    config={'y1':y,'y2':y,'fermion_masses_gev':[1e5],'vev_gev':246.0}
    # Simple checked source expression tests the interface, not Matchete physics.
    c5,mnu,s=numerical_c5('lambdaT3/(MSX^2+MSY^2)',config,b,
                         {'lambdaT3':0.2+0j})
    expected=2*0.2/sum(b.eigenvalues_sq)
    np.testing.assert_allclose(c5,np.full((3,3),expected,dtype=complex),rtol=1e-13)
    np.testing.assert_allclose(mnu,-246**2/2*1e9*c5)
    assert np.isfinite(s).all()


def test_disallowed_mass_basis():
    with pytest.raises(ValueError):
        MixedDoubletMassBasis.build('A_+0',1e6,2e6,500+0j)
