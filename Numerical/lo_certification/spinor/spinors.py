"""Exact matrix and Grassmann checks in a DECLARED Weyl convention.

These checks are conditional mathematics, not evidence that installed Matchete
uses a particular representation/sign for GammaCC.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction
from typing import Mapping
import numpy as np


def convention_probe() -> dict:
    z = np.zeros((2, 2), dtype=complex)
    eye = np.eye(2, dtype=complex)
    s1 = np.array([[0, 1], [1, 0]], dtype=complex)
    s2 = np.array([[0, -1j], [1j, 0]], dtype=complex)
    s3 = np.diag([1, -1]).astype(complex)
    gamma0 = np.block([[z, eye], [eye, z]])
    gamma2 = np.block([[z, s2], [-s2, z]])
    gamma5 = np.block([[-eye, z], [z, eye]])
    pl = (np.eye(4) - gamma5) / 2
    pr = (np.eye(4) + gamma5) / 2
    C = 1j * gamma2 @ gamma0
    eps = 1j * s2
    assert np.array_equal(C[:2, :2], eps)
    assert np.array_equal(C[2:, 2:], -eps)
    assert np.array_equal(C @ C, -np.eye(4))
    assert np.array_equal(C.T, -C)
    assert np.array_equal(C @ pl, np.block([[eps, z], [z, z]]))
    assert np.array_equal(C @ pr, np.block([[z, z], [z, -eps]]))
    assert np.array_equal(np.linalg.inv(C) @ pl, np.block([[-eps, z], [z, z]]))
    # Two arbitrary complex spinors, testing the bilinear including its phase.
    a = np.array([1+2j, 3-1j])
    b = np.array([2-1j, -1+4j])
    A = np.r_[a, [0j, 0j]]
    B = np.r_[b, [0j, 0j]]
    assert A.T @ C @ pl @ B == a.T @ eps @ b
    assert A.T @ np.linalg.inv(C) @ pl @ B == -a.T @ eps @ b
    return {
        'metric': '+---', 'basis': 'Weyl', 'epsilon_12': 1,
        'charge_conjugation_choice': 'C = i gamma^2 gamma^0',
        'PL': '(1-gamma5)/2', 'PR': '(1+gamma5)/2',
        'transpose_C_PL_over_weyl_epsilon': '+1',
        'transpose_Cinv_PL_over_weyl_epsilon': '-1',
        'exact_matrices_and_complex_spinor_checks_passed': True,
        'matchete_GammaCC_identified_with_this_C': False,
        'matchete_gamma_cc_to_master_weyl_factor': None,
    }


@dataclass(frozen=True)
class G:
    """Minimal sparse exterior (Grassmann) algebra with exact coefficients."""
    terms: Mapping[tuple[str, ...], Fraction]

    @staticmethod
    def var(name: str) -> 'G':
        return G({(name,): Fraction(1)})

    @staticmethod
    def scalar(value) -> 'G':
        v = Fraction(value)
        return G({(): v} if v else {})

    def __add__(self, other) -> 'G':
        if not isinstance(other, G):
            other = G.scalar(other)
        out = dict(self.terms)
        for k, v in other.terms.items():
            out[k] = out.get(k, Fraction()) + v
            if not out[k]:
                del out[k]
        return G(out)

    __radd__ = __add__

    def __neg__(self) -> 'G':
        return G({k: -v for k, v in self.terms.items()})

    def __sub__(self, other) -> 'G':
        return self + -other if isinstance(other, G) else self + -Fraction(other)

    def __mul__(self, other) -> 'G':
        if not isinstance(other, G):
            other = G.scalar(other)
        out: dict[tuple[str, ...], Fraction] = {}
        for a, av in self.terms.items():
            for b, bv in other.terms.items():
                if len(set(a+b)) < len(a+b):
                    continue
                inversions = sum(x > y for x in a for y in b)
                key = tuple(sorted(a+b))
                out[key] = out.get(key, Fraction()) + (-1 if inversions % 2 else 1)*av*bv
                if not out[key]:
                    del out[key]
        return G(out)

    __rmul__ = __mul__

    def __eq__(self, other) -> bool:
        if not isinstance(other, G):
            other = G.scalar(other)
        return dict(self.terms) == dict(other.terms)


def wedge(a: tuple[G, G], b: tuple[G, G]) -> G:
    return a[0]*b[1] - a[1]*b[0]


def exchange_probe() -> dict:
    """Substitute exact heavy-field EOM into Dirac/Majorana UV actions.

    Majorana: -M/2 h eps h + sum y_i h eps ell_i
    Dirac:    -M u eps v + u eps A + v eps B
    """
    M = Fraction(11)
    fields = [(G.var(f'e{i}a'), G.var(f'e{i}b')) for i in range(4)]
    weights = [Fraction(2), Fraction(-3), Fraction(5), Fraction(-7)]
    j = tuple(sum((weights[i]*fields[i][k] for i in range(4)), G.scalar(0)) for k in range(2))
    # For the majorana EOM, h = j/M, with epsilon12=+1.
    h = (j[0]*(1/M), j[1]*(1/M))
    l_uv = -(M/2)*wedge(h, h) + wedge(h, j)
    l_contact = wedge(j, j)*(1/(2*M))
    assert l_uv == l_contact
    for a in range(4):
        for b in range(4):
            if a != b:
                assert wedge(fields[a], fields[b]) == wedge(fields[b], fields[a])
    A, B = fields[0], fields[1]
    # Heavy Dirac EOM: u=-B/M, v=A/M.
    u = tuple(-x*(1/M) for x in B)
    v = tuple(x*(1/M) for x in A)
    dirac_uv = -M*wedge(u, v) + wedge(u, A) + wedge(v, B)
    dirac_contact = wedge(A, B)*(1/M)
    assert dirac_uv == dirac_contact
    assert wedge(A, B) == wedge(B, A)
    assert G.var('x')*G.var('y') == -(G.var('y')*G.var('x'))
    return {
        'external_independent_currents': 4,
        'majorana_identical_channel_factor': '1/(2M)',
        'majorana_mixed_channel_factor': '1/M',
        'majorana_distinct_cross_channels': 6,
        'dirac_mixed_channel_factor': '1/M',
        'exact_grassmann_saddle_signs_verified': True,
        'external_crossing_and_anticommutation_checked': True,
        'scope': 'declared two-component p=0 UV action; not native Matchete four-point amplitude',
        'gamma_cc_to_master_weyl_factor': None,
    }
