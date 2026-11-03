"""B1 conditional internal-field phase transport, independent of loop engine.

A symbolic source-action check, NOT a physical Matchete/RGBeta dictionary.

Chosen left-handed source monomials (external Higgs convention held fixed):
 y1 * nu chi S1, y2 * nu xi S2^*, M[xi,chi] xi chi,
 lambda * H H S1 S2^*; Majorana imposes xi=chi.

All basis changes use OLD field = U * NEW field. Coefficients therefore
transform with U (not U^dagger). The reverse active field convention reverses
all phase exponents; the invariant statements are unchanged.
"""
from __future__ import annotations
import math
from typing import Any
import numpy as np
import sympy as sp

class PhaseAuditError(ValueError):
    pass


def require(condition: bool, msg: str) -> None:
    if not condition:
        raise PhaseAuditError(msg)


OPERATORS = ('y1', 'y2', 'MF', 'lambdaT3')
FIELDS_DIRAC = ('nu', 'chi', 'xi', 'S1', 'S2')
FIELDS_MAJORANA = ('nu', 'chi', 'S1', 'S2')


def exponent_matrix(majorana: bool) -> sp.Matrix:
    """Rows are monomial field charges; conjugate scalar has charge -1."""
    if majorana:
        return sp.Matrix([[1, 1, 1, 0],
                          [1, 1, 0, -1],
                          [0, 2, 0, 0],
                          [0, 0, 1, -1]])
    return sp.Matrix([[1, 1, 0, 1, 0],
                      [1, 0, 1, 0, -1],
                      [0, 1, 1, 0, 0],
                      [0, 0, 0, 1, -1]])


def exact_invariant_witness(majorana: bool) -> dict[str, Any]:
    q = exponent_matrix(majorana)
    internal = q[:, 1:]
    invariant = sp.Matrix([1, 1, -1, -1])
    require(internal.T * invariant == sp.zeros(internal.cols, 1),
            'Source spurion not invariant under internal fields')
    require(q.T * invariant == sp.Matrix([2] + [0] * internal.cols),
            'External nu-nu covariance is incorrect')
    ns = internal.T.nullspace()
    require(len(ns) == 1 and ns[0] / ns[0][0] == invariant,
            'Expected unique invariant exponent up to scalar multiplication')
    for wrong in (sp.Matrix([1, 1, -1, 1]),
                  sp.Matrix([1, 1, 1, -1]),
                  sp.Matrix([1, -1, -1, -1]),
                  sp.Matrix([1, 0, 0, 0])):
        require(internal.T * wrong != sp.zeros(internal.cols, 1),
                'Incorrect alternative spurion accepted')
    return {
        'fermion_case': 'Majorana' if majorana else 'Dirac',
        'field_order': list(FIELDS_MAJORANA if majorana else FIELDS_DIRAC),
        'operator_order': list(OPERATORS),
        'field_exponent_matrix': [list(map(int, row)) for row in q.tolist()],
        'unique_internal_invariant_exponents': [1, 1, -1, -1],
        'invariant_source_spurion': 'y1*y2*conjugate(MF)*conjugate(lambdaT3)',
        'external_neutrino_charge': 2,
        'internal_invariant_witness_exact': True,
        'external_nu_pair_covariance_exact': True,
    }


def complex_matrix(raw: Any, name: str, shape=(3, 3)) -> np.ndarray:
    if isinstance(raw, dict):
        require(set(raw) == {'real', 'imag'}, name + ' unsupported complex encoding')
        z = np.asarray(raw['real'], dtype=float) + 1j * np.asarray(raw['imag'], dtype=float)
    else:
        z = np.asarray(raw, dtype=complex)
    require(z.shape == shape and np.all(np.isfinite(z)), name + ' must be finite '+str(shape))
    return z


def complex_scalar(raw: Any, name: str) -> complex:
    if isinstance(raw, dict):
        require(set(raw) == {'real', 'imag'}, name + ' unsupported complex encoding')
        z = complex(float(raw['real']), float(raw['imag']))
    else:
        require(not isinstance(raw, (list, tuple, np.ndarray)), name+' is not scalar')
        z = complex(raw)
    require(math.isfinite(z.real) and math.isfinite(z.imag), name+' nonfinite')
    return z


def unitary(z: Any, name: str, dimension=3) -> np.ndarray:
    z = complex_matrix(z, name, (dimension, dimension))
    require(np.linalg.norm(z.conj().T @ z - np.eye(dimension)) < 1e-11,
            name+' must be unitary')
    return z


def source_flavour_numerator(y1: Any, y2: Any, mf: Any, lam: complex) -> np.ndarray:
    """Source-coupling phase numerator, NOT the complete one-loop self-energy.

    MF is indexed [xi,chi] in the Dirac case, so the correct contraction is
    y1[nu,chi] @ MF.conj().T[chi,xi] @ y2[nu,xi].T.
    For Majorana, MF is symmetric and the same contraction applies.
    """
    a = complex_matrix(y1, 'y1')
    b = complex_matrix(y2, 'y2')
    m = complex_matrix(mf, 'MF')
    l = complex_scalar(lam, 'lambdaT3')
    return a @ m.conj().T @ b.T * l.conjugate()


def change_basis(y1: Any, y2: Any, mf: Any, lam: complex,
                 u_chi: Any, u_xi: Any | None, u_nu: Any,
                 s1_phase: float, s2_phase: float,
                 majorana: bool) -> tuple[np.ndarray, np.ndarray, np.ndarray, complex]:
    """Field basis: old chi=Uchi new chi; old xi=Uxi new xi.

    old nu=Unu new nu, old S1=e^(i c) new S1, old S2=e^(i d) new S2.
    For a Majorana field, Uxi is absent and identically Uchi.
    """
    a = complex_matrix(y1, 'y1'); b = complex_matrix(y2, 'y2')
    m = complex_matrix(mf, 'MF'); l = complex_scalar(lam, 'lambdaT3')
    uc = unitary(u_chi, 'Uchi'); un = unitary(u_nu, 'Unu')
    if majorana:
        require(u_xi is None, 'Majorana must use the same heavy flavour basis')
        require(np.linalg.norm(m-m.T) < 1e-11*max(1.0,np.linalg.norm(m)),
                'Majorana mass tensor must be symmetric')
        ux = uc
    else:
        require(u_xi is not None, 'Dirac xi rotation missing')
        ux = unitary(u_xi, 'Uxi')
    require(math.isfinite(s1_phase) and math.isfinite(s2_phase), 'Invalid scalar rephase')
    a2 = np.exp(1j*s1_phase) * un.T @ a @ uc
    b2 = np.exp(-1j*s2_phase) * un.T @ b @ ux
    m2 = ux.T @ m @ uc
    l2 = np.exp(1j*(s1_phase-s2_phase))*l
    if majorana:
        require(np.linalg.norm(m2-m2.T) < 1e-9*max(np.linalg.norm(m2),1.0),
                'Majorana flavour rotation destroyed symmetry')
    return a2,b2,m2,l2


def relative(a: Any, b: Any) -> float:
    x=np.asarray(a,complex);y=np.asarray(b,complex)
    return float(np.linalg.norm(x-y)/max(np.linalg.norm(x),np.linalg.norm(y),1e-30))


def deterministic_unitary(rng: np.random.Generator, n=3) -> np.ndarray:
    z = rng.normal(size=(n,n)) + 1j*rng.normal(size=(n,n))
    q, r = np.linalg.qr(z)
    diag = np.diag(r)
    return q @ np.diag(diag / np.abs(diag))


def test_flavour_covariance(y1: Any, y2: Any, mf: Any, lam: complex,
                            majorana: bool, seed=1234) -> dict[str, Any]:
    """Independent U(3) heavy and external rotations plus U(1) scalar phases."""
    rng = np.random.default_rng(seed)
    uc = deterministic_unitary(rng)
    ux = None if majorana else deterministic_unitary(rng)
    un = deterministic_unitary(rng)
    s1 = float(rng.uniform(-np.pi, np.pi))
    s2 = float(rng.uniform(-np.pi, np.pi))
    t = change_basis(y1,y2,mf,lam,uc,ux,un,s1,s2,majorana)
    before = source_flavour_numerator(y1,y2,mf,lam)
    after = source_flavour_numerator(*t)
    expected = un.T @ before @ un
    r = relative(after, expected)
    s1old = (before+before.T)/2
    s1new = (after+after.T)/2
    r_sym = relative(s1new, un.T @ s1old @ un)
    require(r < 1e-10 and r_sym < 1e-10,
            'Source spurion flavour rotation covariance failed')
    return {'seed': seed, 'majorana':majorana,
            'U3_flavour_covariance_residual':r,
            'symmetric_nu_pair_covariance_residual':r_sym,
            'source_spurion_numerator_frobenius':float(np.linalg.norm(before)),
            'conditional_only':True}


def sign_orbit_witness() -> dict[str, Any]:
    """Exact sign orbits from internal field rephasings by fourth roots of unity."""
    import itertools
    rows=[]
    for majorana in (False,True):
        q=exponent_matrix(majorana)[:,1:]
        allowed=set()
        for phases in itertools.product(range(4),repeat=q.cols):
            v=q*sp.Matrix(phases)
            mod4=tuple(int(z)%4 for z in v)
            if all(k in (0,2) for k in mod4):
                allowed.add(tuple(k//2 for k in mod4))
        # sign of conjugated coupling is the same -1 in Z2.
        require(len(allowed)==8, 'Internal rephasing orbit should contain 8 of 16 sign patterns')
        parity=lambda x: sum(x)%2
        require(all(parity(x)==0 for x in allowed), 'Odd spurion flip in field orbit')
        require(all(x in allowed for x in itertools.product((0,1),repeat=4)
                    if parity(x)==0), 'Not all even sign flips generated')
        rows.append({'fermion_case':'Majorana' if majorana else 'Dirac',
                     'allowed_internal_rephase_sign_patterns':len(allowed),
                     'even_sign_patterns_are_internal_rephasings':True,
                     'fourth_root_field_phases_sufficient':True,
                     'single_yukawa_flip_in_internal_orbit':(1,0,0,0) in allowed,
                     'both_yukawas_flip_in_internal_orbit':(1,1,0,0) in allowed,
                     'Yukawa_plus_MF_flip_in_internal_orbit':(1,0,1,0) in allowed,
                     'Yukawa_plus_lambda_flip_in_internal_orbit':(1,0,0,1) in allowed})
    return {'cases':rows, 'external_nu_convention_fixed':True}
