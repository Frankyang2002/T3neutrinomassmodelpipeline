"""A(-4) fresh Matchete LLSS source -> RGBeta *coordinate* bridge.

Reconstructs the source T3MixCG Cartesian/spherical SU(2) triplet map from
its exact saved 2x2x3 SparseArray, then embeds the three-flavour tree LLSS
polynomial in the 45-Weyl / 12-real-scalar master-tensor coordinate basis.

The *coordinate map* is algebraic. The four-component GammaCC -> Weyl
signed normalisation, renormalisable Matchete/RGBeta quartic dictionary,
other dimension-five operators and finite MSbar matching are NOT certified.
The output is NOT a physical Wilson coefficient or C5.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Mapping

import numpy as np

from Numerical.parallel_t3.background.general_batch1_api import decode_native


class NativeWilsonBridgeError(ValueError):
    """The saved source or its algebraic representation changed."""


# Exact InputForm as published in the F-only A(-4) fresh Wilson seed.
# Deliberately fail on a different source CG rather than guessing its basis.
_SOURCE_CG = ("SparseArray[Automatic, {2, 2, 3}, 0, {1, {{0, 3, 6}, "
              "{{1, 1}, {1, 3}, {2, 2}, {1, 2}, {2, 1}, {2, 3}}}, "
              "{3^(-1/4), I/3^(1/4), -3^(-1/4), -3^(-1/4), "
              "-3^(-1/4), I/3^(1/4)}}]")


def _require(value: bool, message: str) -> None:
    if not value:
        raise NativeWilsonBridgeError(message)


def source_cg_from_seed(seed: Mapping[str, Any]) -> np.ndarray:
    """Read exactly the published Matchete T3MixCG CSR tensor, in its ordering.

    The first dimension is i=1,2; CSR row offsets are 0,3,6; each stored
    coordinate is (j,k). The three tensor columns are a Cartesian triplet.
    """
    _require(seed.get('IntegratedFields') == ['F'] and
             seed.get('ActiveHeavyFields') == ['S1', 'S2'] and
             seed.get('TreeWilsonTermCount') == 2,
             'Not a three-flavour F-only LLSS source with both scalars')
    terms = seed.get('TreeWilsonTerms')
    _require(isinstance(terms, list) and len(terms) == 2 and
             {t.get('Chirality') for t in terms} == {'PR', 'PL'},
             'Missing complete tree-level Matchete chiral Wilson pair')
    for t in terms:
        op = t.get('TermInputForm', '')
        _require('GammaCC' in op and 'T3MixCG' in op and
                 op.count('Coupling[y1,') == 1 and
                 op.count('Coupling[y2,') == 1 and
                 op.count('Coupling[MF,') == 1 and
                 op.count('Field[NewScalar1,') == 1 and
                 op.count('Field[NewScalar2,') == 1,
                 'Tree LLSS source operator inventory changed')
    reg = seed.get('CGRegistry')
    _require(isinstance(reg, list) and len(reg) == 1 and
             reg[0].get('Name') == 'T3MixCG' and
             reg[0].get('RepsInputForm') ==
             '{Bar[SU2L[fund]], Bar[SU2L[fund]], T3BSMd3}',
             'The Matchete CG representation/orientation is not the saved one')
    normalized = ''.join(str(reg[0].get('TensorInputForm', '')).split())
    expected = ''.join(_SOURCE_CG.split())
    _require(normalized == expected,
             'T3MixCG SparseArray changed; no unsafe source-basis inference')
    a = 3.0 ** (-0.25)
    cg = np.zeros((2, 2, 3), dtype=np.complex128)
    cg[0, 0, 0] = a
    cg[0, 0, 2] = 1j*a
    cg[0, 1, 1] = -a
    cg[1, 0, 1] = -a
    cg[1, 1, 0] = -a
    cg[1, 1, 2] = 1j*a
    _require(np.count_nonzero(cg) == 6 and
             np.max(np.abs(cg - cg.swapaxes(0, 1))) < 1e-14,
             'Source triplet tensor is no longer symmetric in doublet legs')
    return cg


def spherical_cg() -> np.ndarray:
    """Invariant LL S2 CG in the RGBeta highest-weight SU(2) basis.

    Aminus4EFT1Couplings._mix_tensor has the EQUIVARIANT HH->3 CG with
    entries (11,1), (12,2), (22,3), used with the conjugate S2 field.
    The PL Wilson source instead has L L S2 (unconjugated S2); it requires
    the SU(2)-invariant triplet metric J_(m,-m)=(-1)^(1-m), so the correct
    singlet CG is C_invariant[i,j,k] = C_HH[i,j,m] J[m,k].
    This choice is fixed by SU(2) Ward identities, not a free phase guess.
    """
    hh = np.zeros((2, 2, 3), complex)
    hh[0, 0, 0] = hh[1, 1, 2] = 1.
    hh[0, 1, 1] = hh[1, 0, 1] = 1/math.sqrt(2)
    metric = np.array([[0,0,1],[0,-1,0],[1,0,0]],complex)
    return np.einsum('ijm,mk->ijk', hh, metric)


def verify_su2_singlet_cg(cg: np.ndarray) -> float:
    """Check three infinitesimal SU(2) Ward identities for L L S2.

    The lepton SU(2) fundamental and RGBeta highest-weight triplet must
    transform together; this distinguishes LL*S2 from HH*S2dag.
    """
    from Numerical.running.Aminus4EFT1Couplings import su2_generators
    sigma = su2_generators(2)
    tau = su2_generators(3)
    errors=[]
    for d,t in zip(sigma,tau):
        ward = (np.einsum('pi,pjk->ijk',d,cg) +
                np.einsum('pj,ipk->ijk',d,cg) +
                np.einsum('pk,ijp->ijk',t,cg))
        errors.append(float(np.max(np.abs(ward))))
    error=max(errors)
    _require(error < 1e-13,
             'LLSS spherical CG fails fundamental/triplet SU2 Ward identities')
    return error


@dataclass(frozen=True)
class CGMap:
    cartesian_from_spherical: np.ndarray
    spherical_cg_scale: float
    real_s2_cartesian_from_spherical: np.ndarray
    max_reconstruction_residual: float
    max_unitarity_residual: float
    max_real_orthogonality_residual: float


def derive_triplet_map(seed: Mapping[str, Any]) -> CGMap:
    """Solve C_mat[ij,k] = b U[k,s] C_invariant[ij,s] from SOURCE components.

    Because the operator is C_mat[ij,k] * S2_mat[k] (not conjugate S2),
    choose S2_mat = conjugate(U) * S2_sph. This gives exactly
    C_mat * S2_mat = b C_invariant * S2_sph, with no CG projection.
    The SU(2) singlet condition removes the equivariant/conjugate ambiguity.
    """
    mat, sph = source_cg_from_seed(seed), spherical_cg()
    verify_su2_singlet_cg(sph)
    gram = np.einsum('ijs,ijt->st', sph.conj(), sph)
    _require(np.allclose(gram, np.eye(3), atol=2e-15),
             'Spherical triplet CG basis not orthonormal')
    b = math.sqrt(2) * 3**(-0.25)
    U = np.einsum('ijk,ijs->ks', mat, sph.conj()) / b
    unitary_res = float(np.max(np.abs(U.conj().T @ U - np.eye(3))))
    rebuilt = b*np.einsum('ks,ijs->ijk', U, sph)
    cg_res = float(np.max(np.abs(mat-rebuilt)))
    _require(unitary_res < 1e-13 and cg_res < 1e-13,
             'No unitary CG intertwiner in the saved source')
    Z = U.conj()
    R = np.zeros((6, 6))
    for i in range(3):
        for j in range(3):
            R[2*i, 2*j] = Z[i, j].real
            R[2*i, 2*j+1] = -Z[i, j].imag
            R[2*i+1, 2*j] = Z[i, j].imag
            R[2*i+1, 2*j+1] = Z[i, j].real
    ortho_res = float(np.max(np.abs(R @ R.T - np.eye(6))))
    _require(ortho_res < 1e-13, 'Complex-to-real triplet map not orthogonal')
    # Test the actual PL operator contraction under the S2 field map.
    contracted = np.einsum('ijk,ks->ijs', mat, Z)
    _require(np.allclose(contracted, b*sph, atol=1e-13),
             'PL CG * S2 field rotation failed')
    return CGMap(U, b, R, cg_res, unitary_res, ortho_res)


def three_flavour_wbar(f_endpoint: Mapping[str, Any]) -> tuple[np.ndarray, dict]:
    """Conjugated source W12[p,q] in the near-diagonal native MF frame.

    Full MF matrix basis rotation is not silently approximated when large.
    """
    _require(f_endpoint.get('model') == 'A_-4' and
             f_endpoint.get('all_heavy_fermions_to_integrate_out_together') is True and
             float(f_endpoint.get('scale_GeV', 0)) == 1e5,
             'Wrong F endpoint or integration order')
    v = f_endpoint['native_uv_rgbeta_couplings_at_muF']
    M, y1, y2 = (np.asarray(decode_native(v[k]), complex)
                 for k in ('MF','y1','y2'))
    _require(all(x.shape == (3,3) and np.isfinite(x).all()
                 for x in (M,y1,y2)), 'Not three finite heavy flavours')
    norm = np.linalg.norm(M)
    _require(norm > 0, 'Zero heavy mass matrix')
    rel = float(np.linalg.norm(M-np.diag(np.diag(M)))/norm)
    _require(rel <= 1e-8,
             'Heavy MF is non-diagonal; native chiral rotation required')
    diag = np.diag(M)
    _require(min(diag.real) > 0 and
             max(abs(diag.imag)) <= 1e-8*max(abs(diag.real)),
             'Complex/nonpositive heavy mass basis requires rotation')
    mr = diag.real
    Wbar = (y1.conj()/mr[None,:]) @ y2.conj().T
    return Wbar, {'heavy_fermions':3,'MF_offdiagonal_relative':rel,
                  'MF_diagonal_GeV':mr.tolist(),
                  'Wbar_ordered_frobenius_GeVinv':float(np.linalg.norm(Wbar)),
                  'MF_diagonal_approximation':True}


def formal_tree_master(wbar: np.ndarray, cgmap: CGMap, *,
                       signed_spinor_factor: float = 1.) -> np.ndarray:
    """Embed the exact polynomial orbit, no inferred physical spinor factor.

    Source monomial: Wbar[p,q] C_mat[i,j,k] (L[p,i] L[q,j]) S1bar S2_mat[k]
    = Wbar[p,q] b C_invariant[i,j,k] (L[p,i] L[q,j]) S1bar S2_sph[k].
    
    Each complex scalar uses (R+iI)/sqrt(2). The generic master polynomial
    C[ijab] psi_i psi_j phi_a phi_b is symmetric under i<->j, a<->b:
    divide an ordered monomial by four when distributing its two orbits.
    No extra 2, 1/2 or sign is attributed to GammaCC here.
    """
    wbar = np.asarray(wbar, complex)
    _require(wbar.shape == (3,3) and np.isfinite(wbar).all(),
             'W12 must be finite, complex 3x3 ordered flavour matrix')
    _require(np.isfinite(signed_spinor_factor) and signed_spinor_factor != 0,
             'Invalid FORMAL signed spinor factor')
    from Numerical.running.Aminus4EFT1Couplings import FI, S_LABELS, F_LABELS
    _require(len(F_LABELS) == 45 and len(S_LABELS) == 12 and
             S_LABELS[4:6] == ('S1[1,R]','S1[1,I]') and
             S_LABELS[6:] == tuple(f'S2[{m},{ri}]'
                 for m in (1,2,3) for ri in ('R','I')),
             'EFT1 real-scalar ordering changed')
    c = np.zeros((45,45,12,12),complex)
    sph = spherical_cg()
    # Expansion of S1bar * S2, including both 1/sqrt(2) factors.
    pair = ((4,0, .5), (4,1, .5j), (5,0, -.5j), (5,1, .5))
    for p in range(3):
        for q in range(3):
            for i in range(2):
                for j in range(2):
                    i0 = FI[f'L[{p+1},{i+1}]']
                    j0 = FI[f'L[{q+1},{j+1}]']
                    for m in range(3):
                        amp = (signed_spinor_factor * wbar[p,q] *
                               cgmap.spherical_cg_scale * sph[i,j,m])
                        if amp == 0:
                            continue
                        for a, bpar, phase in pair:
                            b = 6+2*m+bpar
                            val = amp*phase/4
                            c[i0,j0,a,b] += val
                            c[j0,i0,a,b] += val
                            c[i0,j0,b,a] += val
                            c[j0,i0,b,a] += val
    _require(np.all(np.isfinite(c)) and
             np.max(np.abs(c-c.swapaxes(0,1))) < 1e-22 and
             np.max(np.abs(c-c.swapaxes(2,3))) < 1e-22,
             'Symmetric real-scalar master orbit construction failed')
    return c


def source_validation(seed: Mapping[str, Any], tree_text: str,
                      status: Mapping[str, Any],
                      lag: Mapping[str, Any]) -> None:
    """Reject exploratory matching and an unverified complete F tree source."""
    source_cg_from_seed(seed)
    _require(status.get('model') == 'A_-4' and
             status.get('first_threshold_reused') is False and
             status.get('A_joint_scalar_match_performed') is True and
             status.get('A_physical_scalar_soft_Weinberg_terms') == 10 and
             status.get('exploratory_source_theory_used') is False,
             'Require successful fresh three-F, joint-two-S physical soft run')
    _require(lag.get('CheckLagrangian_original') == 'True' and
             lag.get('accepted_canonicalisation') == 'ExactSavedTree',
             'Exact saved EFT1 source failed native CheckLagrangian')
    _require(tree_text.count('GammaCC') == 2 and
             tree_text.count('Coupling[T3EFT1W12,') == 2 and
             'NewFermion' not in tree_text,
             'F tree must contain precisely the two LLSS chiral tree operators')
