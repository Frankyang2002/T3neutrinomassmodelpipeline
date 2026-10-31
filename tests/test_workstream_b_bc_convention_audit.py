"""Independent algebraic controls; does not promote a physical map certificate."""
from __future__ import annotations

import numpy as np
import pytest
import sympy as sp

from scripts import WorkstreamBBCConventionAudit as audit


def test_two_invariant_110101_matrix_is_exact_and_not_diagonal():
    # R_j = sum_i A_ij M_i.  L_M=-V_RGB => g=-A^{-1} h.
    A = audit.A110101
    C = audit.EXPECTED_110101
    assert audit.same_matrix(A*C, -sp.eye(2))
    assert audit.same_matrix(C, sp.Matrix([[0,2*sp.sqrt(3)/3],
                                           [-1,-sp.sqrt(3)/3]]))
    assert not audit.same_matrix(C,-sp.eye(2))
    h=sp.Matrix([2+3*sp.I,-4+sp.I])
    assert audit.same_matrix(A*(C*h),-h)
    assert not audit.same_matrix(A*(sp.conjugate(C*h)),-h)


def test_cartesian_hessian_matches_independent_finite_difference():
    for mu in (0.3+0.2j, -0.17+0.41j, -0.5+0j):
        H= audit.expected_hessian(3.2,5.1,mu)
        fd=audit.finite_difference_hessian(3.2,5.1,mu)
        assert np.allclose(H,fd,atol=1e-9,rtol=0)
        assert np.allclose(H,H.T,atol=1e-14,rtol=0)
        delta=np.sqrt(((3.2-5.1)/2)**2+abs(mu)**2)
        pred=np.array([(3.2+5.1)/2-delta]*4+[(3.2+5.1)/2+delta]*4)
        assert np.allclose(np.linalg.eigvalsh(H),pred,atol=1e-12,rtol=0)


def test_wrong_mass_prefactor_and_conjugation_are_detected():
    mu=0.3+0.2j
    correct=audit.expected_hessian(3.2,5.1,mu)
    assert not np.allclose(correct,audit.expected_hessian(3.2,5.1,mu/2))
    assert not np.allclose(correct,audit.expected_hessian(3.2,5.1,mu.conjugate()))
    assert np.allclose(np.linalg.eigvalsh(correct),
                       np.linalg.eigvalsh(audit.expected_hessian(3.2,5.1,mu.conjugate())))
    # Spectra alone cannot resolve the phase/conjugation convention.


def test_eft1_mixed_mass_zero_boundary_source_and_reality():
    bmu=("6*lambdaS1S1bar2S2bar*mS1Sq + 6*lambdaS1barS2S2bar2*mS2Sq"
         " - (9*g2^2*muT3GeneralS1S2)/2 - (3*gY^2*muT3GeneralS1S2)/2"
         " + 2*lambda12*muT3GeneralS1S2"
         " - (3*lambda12Adj*muT3GeneralS1S2)/2"
         " + 12*lambdaS1bar2S2bar2*Bar[muT3GeneralS1S2]")
    beta=audit.eft1_poly(bmu)
    zero=sp.simplify(beta.subs({audit.symbol(audit.MASS):0,
                                sp.conjugate(audit.symbol(audit.MASS)):0}))
    expected=6*(audit.symbol(audit.L1)*audit.symbol(audit.M1)
                +audit.symbol(audit.L2)*audit.symbol(audit.M2))
    assert sp.simplify(zero-expected)==0
    assert zero != 0
    with pytest.raises(ValueError,match='Non-polynomial'):
        audit.eft1_poly('Tr[nonScalar]')


def test_full_standalone_mass_algebra():
    result=audit.mathematical_tests()
    assert result['hessian'] < 5e-9
    assert result['rephasing'] < 5e-12
