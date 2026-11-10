#!/usr/bin/env python3
"""Read-only D(-2) common-canonical SU(2) Ward check, at retained MF.

Uses:
  rough_inputs/.../boundary.npz: C_LLSS_formal_muF
  rough_inputs/.../background.npz: lambda4, mu_gev
  rough_inputs/.../saved_native_trajectory.npz: eft1__lambdaT3, mu_eft1_gev

Wilson: source dual lepton indices -> fundamental via epsilon; source spin-1
triplet -> canonical spin-1 via unitary intertwiner; singlet slots left unchanged.
Quartic: replaces ONLY D-family lambdaT3 unit tensor's uncontracted triplet
CG by its spin-1-metric-contracted variant; same saved native coupling.
No production code, input, running or fitting is changed.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from Numerical.parallel_t3.wilson.basis import (
    build_basis, spin_generators, real_scalar_generator,
)
from Numerical.parallel_t3.wilson.coordinates import (
    embed_formal_llss, embed_c5,
)
from Numerical.parallel_t3.rough_inputs import builder

DEFAULT = Path(
    "output/parallel_t3/rough_inputs/fit16_cc2dcbec4a/a/D_m2"
)


def triplet_intertwiner():
    j = 1j
    m = [
        np.array([[0, 0, 0], [0, 0, -j], [0, j, 0]], complex),
        np.array([[0, -j, 0], [j, 0, 0], [0, 0, 0]], complex),
        np.array([[0, 0, j], [0, 0, 0], [-j, 0, 0]], complex),
    ]
    src = tuple(-x for x in m)
    can = spin_generators(3)
    a = np.vstack([
        np.kron(np.eye(3), t_can) - np.kron(t_src.T, np.eye(3))
        for t_can, t_src in zip(can, src)
    ])
    _, singular, vh = np.linalg.svd(a)
    u = vh[-1].conj().reshape((3, 3), order="F")
    u /= np.sqrt(np.trace(u.conj().T @ u).real / 3)
    u *= np.exp(-1j * np.angle(u[0, 0]))
    o = np.zeros((6, 6), float)
    for i in range(3):
        for k in range(3):
            z = u[i, k]
            o[2*i:2*i+2, 2*k:2*k+2] = [
                [z.real, -z.imag], [z.imag, z.real]
            ]
    return src, can, u, o, float(singular[-1])


def wilson_ward(c, tf, ts):
    d = (
        np.einsum("ik,kjab->ijab", tf, c, optimize=True)
        + np.einsum("jk,ikab->ijab", tf, c, optimize=True)
        + np.einsum("ac,ijcb->ijab", ts, c, optimize=True)
        + np.einsum("bc,ijac->ijab", ts, c, optimize=True)
    )
    return float(np.linalg.norm(d) / max(np.linalg.norm(c), 1e-300))


def quartic_ward(x, t):
    d = (
        np.einsum("ae,ebcd->abcd", t, x, optimize=True)
        + np.einsum("be,aecd->abcd", t, x, optimize=True)
        + np.einsum("ce,abed->abcd", t, x, optimize=True)
        + np.einsum("de,abce->abcd", t, x, optimize=True)
    )
    return float(np.linalg.norm(d) / max(np.linalg.norm(x), 1e-300))


def corrected_d_portal(ns):
    """Spherical triplet metric J contracts the symmetric HH CG index."""
    j = np.array([[0, 0, 1], [0, -1, 0], [1, 0, 0]], float)
    g = builder._cg()  # (triplet_a, H_i, H_j), verified current source API
    t = np.einsum("aij,ak->ijk", g, j, optimize=True)
    out = np.zeros((ns, ns, ns, ns), complex)
    # H H S1 S2*, where for D(-2) S2 is an SU(2) singlet.
    for i, k, a in np.argwhere(np.abs(t) > 1e-14):
        out += t[i, k, a] * builder._monomial(
            builder._u(ns, 0, int(i)),
            builder._u(ns, 0, int(k)),
            builder._u(ns, 4, int(a)),
            builder._u(ns, 10, 0, True),
        )
    return out


def load_arrays(folder):
    result = {}
    for name, keys in (
        ("boundary.npz", ("C_LLSS_formal_muF", "C5_hard")),
        ("background.npz", ("lambda4", "mu_gev")),
        ("saved_native_trajectory.npz",
         ("mu_eft1_gev", "eft1__lambdaT3")),
    ):
        path = folder / name
        if not path.is_file():
            raise FileNotFoundError(str(path))
        with np.load(path, allow_pickle=False) as f:
            missing = set(keys) - set(f.files)
            if missing:
                raise KeyError(f"{path}: missing {sorted(missing)}; keys={f.files}")
            for key in keys:
                result[key] = np.asarray(f[key]).copy()
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", type=Path, default=DEFAULT,
                   help="Saved D_m2 rough_inputs/a folder (read-only)")
    p.add_argument("--tolerance", type=float, default=1e-10)
    args = p.parse_args()
    if args.tolerance <= 0:
        p.error("tolerance must be positive")

    data = load_arrays(args.input)
    basis = build_basis("D_-2")
    raw = np.asarray(data["C_LLSS_formal_muF"], complex)
    if raw.shape != (3, 3, 2, 2, 8, 8):
        raise ValueError(f"Expected formal tensor (3,3,2,2,8,8), got {raw.shape}")

    src, can, u, o, sval = triplet_intertwiner()
    if np.linalg.norm(u.conj().T @ u - np.eye(3)) > args.tolerance:
        raise RuntimeError("Triplet map not unitary")
    o8 = np.eye(8)
    o8[:6, :6] = o
    # embed_formal_llss: phi_source = R phi_canonical, R=O.T
    r8 = o8.T
    eps = np.array([[0, 1], [-1, 0]], complex)
    canonical_weak = np.einsum(
        "iI,jJ,pqIJab->pqijab", eps, eps, raw, optimize=True
    )
    original = embed_formal_llss(raw, basis, np.eye(8), 1.0)
    fixed = embed_formal_llss(canonical_weak, basis, r8, 1.0)
    c5 = np.asarray(data["C5_hard"], complex)
    if c5.shape != (3, 3):
        raise ValueError(f"Expected C5_hard (3,3), got {c5.shape}")
    fixed += embed_c5(c5, basis)
    original += embed_c5(c5, basis)

    nf, ns = basis.nf, basis.ns
    inds = [i for pair in basis.leptons for i in pair]

    def ll_block(c):
        return c[np.ix_(inds, inds, range(ns), range(ns))]

    print("=== D(-2) canonical tensor check at MF ===")
    print("Source:", args.input)
    print("Triplet intertwiner singular value:", f"{sval:.6e}")
    print("Wilson original/fixed norms:",
          f"{np.linalg.norm(original):.12e}",
          f"{np.linalg.norm(fixed):.12e}")
    w0, w1 = [], []
    generators = basis.generators()["SU2"]
    for n, (tf, ts) in enumerate(generators, 1):
        old = wilson_ward(ll_block(original), tf[np.ix_(inds, inds)], ts)
        new = wilson_ward(ll_block(fixed), tf[np.ix_(inds, inds)], ts)
        w0.append(old)
        w1.append(new)
        print(f"Wilson T{n}: provisional={old:.12e}, transformed={new:.12e}")

    # The source D(-2) lambdaT3 portal is missing the spin-one invariant J.
    units = builder.quartic_unit_tensors("D_-2")
    portal_old = np.asarray(units["lambdaT3"], complex)
    portal_new = corrected_d_portal(ns)
    if portal_old.shape != portal_new.shape or portal_old.shape != (ns,) * 4:
        raise ValueError("Unexpected quartic unit tensor geometry")

    times = np.asarray(data["mu_eft1_gev"], float)
    bg_times = np.asarray(data["mu_gev"], float)
    lambdas = np.asarray(data["lambda4"], float)
    coeffs = np.asarray(data["eft1__lambdaT3"])
    if lambdas.shape != (len(bg_times), ns, ns, ns, ns):
        raise ValueError(f"Wrong lambda4 array shape: {lambdas.shape}")
    if len(times) != len(coeffs):
        raise ValueError("Native lambdaT3 trajectory length mismatch")
    if len(times) != len(bg_times) or not np.allclose(times, bg_times, rtol=1e-10, atol=1e-6):
        raise ValueError("Background and native mu grids disagree")
    mf_idx = int(np.argmax(times))
    if abs(times[mf_idx] - 1e5) > 1e-6 * 1e5:
        raise ValueError(f"Expected MF=1e5 GeV, got {times[mf_idx]}")
    z = complex(np.asarray(coeffs[mf_idx]).item())
    old_lambda = lambdas[mf_idx].copy()
    delta = 2 * np.real(z * (portal_new - portal_old))
    fixed_lambda = old_lambda + delta

    print("MF [GeV]:", f"{times[mf_idx]:.12g}")
    print("Native lambdaT3 at MF:", z)
    print("Quartic original/fixed norms:",
          f"{np.linalg.norm(old_lambda):.12e}",
          f"{np.linalg.norm(fixed_lambda):.12e}")
    q0, q1, qp0, qp1 = [], [], [], []
    for n, (_, ts) in enumerate(generators, 1):
        a = quartic_ward(old_lambda, ts)
        b = quartic_ward(fixed_lambda, ts)
        old_portal = 2 * np.real(z * portal_old)
        new_portal = 2 * np.real(z * portal_new)
        c = quartic_ward(old_portal, ts)
        d = quartic_ward(new_portal, ts)
        q0.append(a); q1.append(b); qp0.append(c); qp1.append(d)
        print(f"Quartic T{n}: original={a:.12e}, J-corrected={b:.12e}, "
              f"portal original={c:.12e}, portal J={d:.12e}")

    # Test the tensor structures independently from the coupling's magnitude.
    for name, portal in (("original unit portal", portal_old),
                         ("J-corrected unit portal", portal_new)):
        print(name, "SU2 residuals:",
              " ".join(f"{quartic_ward(portal, ts):.4e}"
                       for _, ts in generators))
    print("Wilson SU2:", "PASS" if max(w1) < args.tolerance else "FAIL")
    print("Quartic SU2:", "PASS" if max(q1) < args.tolerance else "FAIL")
    print("J portal unit SU2:",
          "PASS" if all(quartic_ward(portal_new, ts) < args.tolerance
                        for _, ts in generators) else "FAIL")
    print("No production or saved files modified.")
    print("NOT certified: source field meaning, singlet real/complex map, "
          "GammaCC/Weyl factor, threshold loop matching.")
    if max(w1) >= args.tolerance or max(q1) >= args.tolerance:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
