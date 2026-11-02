"""Simultaneous one-scale scalar-pair matching in the T3 real-scalar EFT1 basis.

EFT1: V4 = lambda[p,q,a,b] h_p h_q s_a s_b /4! as part of the
fully symmetric real quartic tensor, and C[fg,i,j,a,b] is the symmetric
LLss coefficient (no conjugation in the scalar-index contraction).

The sign/prefactor is anchored to the project's EFT1 RGE convention
  16*pi^2 d C_LLhh/d ln(mu) |_quartic = lambda_hhss C_LLss.
The heavy-leg bubble threshold is -lambda*C*B0bar/(2*16*pi^2).
It is a full finite one-loop contribution at O(C_LLss lambda) in the
unbroken-theory dimension-five truncation. It is NOT by itself a full
one-loop matching for arbitrary unprovided EFT1 operators or schemes.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Mapping

import numpy as np

LOOP = 16.0 * math.pi**2


def _finite(a: np.ndarray, label: str) -> np.ndarray:
    x = np.asarray(a)
    if not (np.isfinite(x.real).all() and np.isfinite(x.imag).all()):
        raise ValueError(f"Nonfinite {label}")
    return x


def b0bar(m2a: float, m2b: float, mu: float) -> float:
    """MSbar finite B0(0; m_a^2,m_b^2), mass squares in GeV^2; mu in GeV.

    Numerically stable for exactly and nearly degenerate positive masses.
    """
    a, b, mu = float(m2a), float(m2b), float(mu)
    if not (math.isfinite(a) and math.isfinite(b) and math.isfinite(mu)):
        raise ValueError("Masses and scale must be finite")
    if a <= 0 or b <= 0 or mu <= 0:
        raise ValueError("Integrated scalar masses and scale must be positive")
    c = 0.5 * (a + b)
    d = (a - b) / (a + b)
    if abs(d) < 1.0e-3:
        # (a ln a - b ln b)/(a-b) = ln(c)+1-d^2/6-d^4/20-...
        d2 = d * d
        return -math.log(c / (mu * mu)) + d2 / 6 + d2*d2/20 + d2**3/42
    # Log ratios prevent overflow when the common mass scale is large.
    return 1.0 - (a * math.log(a/(mu*mu)) - b * math.log(b/(mu*mu))) / (a - b)


def weinberg_tensor() -> np.ndarray:
    """Dimensionless Weinberg external tensor T[i,j,A,B].

    H=(h0+i*h1,h2+i*h3)/sqrt(2), epsilon_01=+1,
    O5=(L^i eps_ij H^j)(L^k eps_kl H^l), symmetric in real Higgs legs.
    The output normalization is C5*T; the source adapter must use the same
    Lorentz, flavour and Lagrangian prefactors. No hidden conjugation.
    """
    eps = np.array([[0.0, 1.0], [-1.0, 0.0]], complex)
    t = np.array([[1, 1j, 0, 0], [0, 0, 1, 1j]], complex)/math.sqrt(2)
    v = eps @ t
    out = np.einsum("ia,jb->ijab", v, v, optimize=True)
    return 0.5 * (out + out.swapaxes(2, 3))


@dataclass(frozen=True)
class EigenBasis:
    masses2: np.ndarray
    U: np.ndarray

    @classmethod
    def from_mass2(cls, mass2: np.ndarray) -> "EigenBasis":
        m = np.asarray(_finite(mass2, "real mass-squared tensor"))
        if m.ndim != 2 or m.shape[0] != m.shape[1] or m.shape[0] < 2:
            raise ValueError("mass2 must be square with at least two real scalar components")
        if np.any(m.imag != 0):
            raise ValueError("Real scalar components require a real symmetric mass2")
        m = m.real
        if not np.allclose(m, m.T, rtol=1e-11, atol=1e-10):
            raise ValueError("Non-symmetric scalar mass2")
        masses, U = np.linalg.eigh(m)
        if np.min(masses) <= 0:
            raise ValueError("Non-positive heavy scalar eigenvalue")
        if not np.allclose(U.T@m@U, np.diag(masses), rtol=1e-10, atol=1e-7):
            raise ValueError("Mass diagonalization failure")
        return cls(masses, U)


def check_tensors(C: np.ndarray, quartic: np.ndarray, n: int) -> tuple[np.ndarray, np.ndarray]:
    C = np.asarray(_finite(C, "C_LLSS"), complex)
    q = np.asarray(_finite(quartic, "lambda_HHSS"), complex)
    if C.shape != (3, 3, 2, 2, n, n):
        raise ValueError(f"C_LLSS shape should be (3,3,2,2,{n},{n}); got {C.shape}")
    if q.shape != (4, 4, n, n):
        raise ValueError(f"lambda_HHSS shape should be (4,4,{n},{n}); got {q.shape}")
    scale = max(float(np.max(abs(C))), 1e-30)
    if not np.allclose(C, C.transpose(1, 0, 3, 2, 4, 5), rtol=1e-9, atol=1e-10*scale):
        raise ValueError("C_LLSS violates combined lepton-generation/weak-index exchange")
    if not np.allclose(C, C.swapaxes(4, 5), rtol=1e-9, atol=1e-10*scale):
        raise ValueError("C_LLSS violates scalar-leg symmetry")
    scaleq = max(float(np.max(abs(q))), 1e-30)
    if not np.allclose(q, q.swapaxes(0, 1), rtol=1e-9, atol=1e-10*scaleq):
        raise ValueError("lambda violates external Higgs real-leg symmetry")
    if not np.allclose(q, q.swapaxes(2, 3), rtol=1e-9, atol=1e-10*scaleq):
        raise ValueError("lambda violates internal real-scalar leg symmetry")
    if np.any(q.imag != 0):
        raise ValueError("Quartic potential in a real scalar basis must be real")
    return C, q.real


def rotate_to_eigenbasis(C: np.ndarray, q: np.ndarray, U: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Real orthogonal rotation s_old = U s_eigen; no Hermitian conjugation."""
    c = np.einsum("fgijab,ak,bl->fgijkl", C, U, U, optimize=True)
    quartic = np.einsum("ABab,ak,bl->ABkl", q, U, U, optimize=True)
    return c, quartic


def project_weinberg(amplitude: np.ndarray, *, strict: bool = True,
                     tolerance: float = 1.0e-8) -> tuple[np.ndarray, float]:
    """Project full LLhh amplitude onto O5, rejecting non-Weinberg leftovers.

    The residual detects missing weak/Higgs tensor structures and, in a broken
    spectrum, possible higher-dimensional EWSB contamination. It cannot certify
    that the upstream CG normalization is correct.
    """
    amplitude = np.asarray(_finite(amplitude, "LLhh amplitude"), complex)
    if amplitude.shape != (3, 3, 2, 2, 4, 4):
        raise ValueError("Full amplitude must have shape (3,3,2,2,4,4)")
    T = weinberg_tensor()
    projection = np.einsum("fgijAB,ijAB->fg", amplitude, T.conj(), optimize=True)/np.vdot(T,T)
    leftover = amplitude - np.einsum("fg,ijAB->fgijAB", projection,T,optimize=True)
    residual = float(np.linalg.norm(leftover)/max(np.linalg.norm(amplitude),1e-30))
    flavor_residual = float(np.linalg.norm(projection - projection.T))
    if flavor_residual > tolerance * max(np.linalg.norm(projection), 1e-30):
        raise ValueError("Projected C5 is not complex symmetric")
    if strict and residual > tolerance:
        raise ValueError(f"LLHH tensor is not pure Weinberg: relative residual {residual:.4g}")
    return 0.5*(projection + projection.T), residual


def scalar_bubble_matching(
    C: np.ndarray,
    quartic: np.ndarray,
    mass_squared: np.ndarray,
    mu_s: float,
    *,
    reference_scale: float | None = None,
    group_names: list[str] | None = None,
) -> dict:
    """All heavy scalar eigenstates simultaneously at mu_s; no scalar EFT2.

    Each ordered eigenstate pair contributes exactly once to the tensor sum.
    The returned per-(unordered)-pair contribution includes its correct
    1 for a=a / 2 for a!=b multiplicity. Both scalar multiplets are kept.
    """
    basis = EigenBasis.from_mass2(mass_squared)
    n = len(basis.masses2)
    C, quartic = check_tensors(C, quartic, n)
    Cmass, qmass = rotate_to_eigenbasis(C, quartic, basis.U)
    mu_s = float(mu_s)
    mref = float(mu_s if reference_scale is None else reference_scale)
    if not (math.isfinite(mu_s) and math.isfinite(mref) and mu_s > 0 and mref > 0):
        raise ValueError("Invalid matching or reference scale")
    if group_names is not None and len(group_names) != n:
        raise ValueError("group_names must specify each real scalar component")

    amp = np.zeros((3,3,2,2,4,4),dtype=complex)
    amp_ref = np.zeros_like(amp)
    pair_amplitudes: dict[str, np.ndarray] = {}
    pair_data: list[dict] = []
    scalelog = 2*math.log(mu_s/mref)
    for a in range(n):
        for b in range(a,n):
            # Non-diagonal scalar vertices carry both ordered contractions.
            factor = 1 if a==b else 2
            contraction = np.einsum("fgij,AB->fgijAB", Cmass[...,a,b],qmass[...,a,b],optimize=True)
            b_ref = b0bar(basis.masses2[a],basis.masses2[b],mref)
            b_mu = b_ref + scalelog
            pair = -factor * contraction * b_mu/(2*LOOP)
            key = f"e{a:02d}_e{b:02d}"
            pair_amplitudes[key] = pair
            amp += pair
            amp_ref += -factor * contraction*b_ref/(2*LOOP)
            participation = None
            if group_names is not None:
                participation = {
                    label:[float(np.sum(basis.U[np.array(group_names)==label,k]**2)) for k in (a,b)]
                    for label in dict.fromkeys(group_names)
                }
            pair_data.append({"pair":key,"eigenstate_indices":[a,b],
                              "masses_gev":[float(math.sqrt(basis.masses2[k])) for k in (a,b)],
                              "multiplicity":factor,
                              "b0_msbar_at_mu_s":b_mu,
                              "b0_finite_at_reference":b_ref,
                              "b0_explicit_mu_log":scalelog,
                              "scalar_species_participation":participation,
                              "frobenius_norm_full_amplitude":float(np.linalg.norm(pair))})

    return {"amplitude":amp,"amplitude_finite_reference":amp_ref,
            "amplitude_scale_log":amp-amp_ref,
            "pair_amplitudes":pair_amplitudes,"pair_data":pair_data,
            "eigenvalues_sq_gev2":basis.masses2,
            "mass_rotation_U":basis.U,
            "mu_s_gev":mu_s,"reference_scale_gev":mref}


def match_c5(C: np.ndarray, quartic: np.ndarray, mass_squared: np.ndarray,
             mu_s: float, *, c5_eft1: np.ndarray | None = None,
             group_names: list[str] | None = None, reference_scale: float | None = None,
             strict_projection: bool = True, tolerance: float = 1e-8) -> dict:
    """Add transported EFT1 C5(mu_s) once and project the simultaneous scalar loop."""
    result = scalar_bubble_matching(C,quartic,mass_squared,mu_s,
                                    reference_scale=reference_scale,group_names=group_names)
    bubble, residual = project_weinberg(result["amplitude"],strict=strict_projection,tolerance=tolerance)
    ref, _ = project_weinberg(result["amplitude_finite_reference"],strict=False)
    logs, _ = project_weinberg(result["amplitude_scale_log"],strict=False)
    direct = np.zeros((3,3),complex) if c5_eft1 is None else np.asarray(_finite(c5_eft1, "C5_EFT1"),complex)
    if direct.shape != (3,3) or not np.allclose(direct,direct.T,rtol=1e-10,atol=1e-28):
        raise ValueError("C5_EFT1 must be a symmetric complex 3x3 matrix")
    pair_c5 = {}
    for key,amplitude in result["pair_amplitudes"].items():
        pair_c5[key],_ = project_weinberg(amplitude,strict=False)
    total = direct + bubble
    return {**result,"c5_eft1":direct,"c5_scalar_bubble":bubble,
            "c5_reference_finite":ref,"c5_explicit_mu_log":logs,
            "c5_total":total,"pair_c5":pair_c5,
            "projection_residual":residual,
            "projection_passed":residual<=tolerance,
            "rge_log_identity": "d(delta_C5)/d ln(mu_s) = -project(lambda*C)/(16*pi^2) at fixed masses and couplings"}
