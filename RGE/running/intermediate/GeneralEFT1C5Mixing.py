"""Certified-boundary EFT1 -> Weinberg one-loop *source* mixing.

No existing matching/running module is modified.  The master RGE follows
Aebischer--Bresciani--Selimovic, arXiv:2512.16890v1, eq. (4.85), in the
repository's RGE/general/WilsonTensorRGE.py conventions:

 L_5 = C[i,j,a,b] psi_i psi_j phi_a phi_b + h.c.,
 V_4 = lambda[a,b,c,d] phi_a phi_b phi_c phi_d / 4!,
 H^0 = (h_R + i h_I)/sqrt(2),
 L_5^W = (1/2) C5[p,q] (L_p.H)(L_q.H) + h.c.

Consequently C[nu_p,nu_q,h_R,h_R] = C5[p,q]/4.

This is a fixed-overall-one-loop calculation.  The source C^(0) and
renormalizable couplings are frozen at mu_F.  Their one-loop evolution inside
an already-one-loop C5 source starts at two loops.  All seven RGE terms are
kept, including terms that vanish for a particular certified boundary.

The *only* permitted physical input is an explicitly normalized and fully
certified master-tensor boundary.  Legacy Matchete ordered tensors MUST NOT
be passed in directly: their spinor phase/factor is not yet certified.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import permutations
from math import log, pi
from typing import Any, Mapping

import numpy as np

MODELS = tuple(f"{family}_{a:+d}" for family, alphas in {
    "A": (-4, -2, 0), "B": (-3, -1, 1), "C": (-3, -1, 1),
    "D": (-2, 0, 2), "E": (-4, -2, 0, 2),
}.items() for a in alphas)
FORMAT = "t3-eft1-certified-master-weyl-v1"
RESULT_FORMAT = "t3-eft1-one-loop-c5-source-v1"
LOOP = 1.0 / (16.0 * pi**2)
TERMS = (
    "yukawa_wavefunction", "scalar_pair", "mixed_yukawa_gauge",
    "crossed_yukawa", "conjugate_coefficient_yukawa",
    "scalar_anomalous_dimension", "fermion_anomalous_dimension",
)
# Certification requires an external normalization/counterterm derivation, not
# simply a successful numerical Ward check.
CERTIFICATES = (
    "matchete_gamma_cc_to_weyl_factor_and_sign",
    "canonical_C_ijab_operator_normalization",
    "complete_tree_dim5_operator_inventory",
    "other_dim5_sources_audited_at_one_loop",
    "complete_EFT1_yukawa_and_gauge_dictionary",
    "complete_scalar_quartic_counterterms",
    "all_heavy_scalar_components_and_mixed_masses",
    "weinberg_C5_symmetric_convention",
    "MSbar_subtraction_and_matching_log_partition",
)


class BoundaryError(ValueError):
    """Absent/inconsistent physics input: no silent zeros or inferred factors."""


def _finite_float(x: Any, label: str) -> float:
    try:
        result = float(x)
    except (TypeError, ValueError, OverflowError) as exc:
        raise BoundaryError(f"{label}: expected finite float") from exc
    if not np.isfinite(result):
        raise BoundaryError(f"{label}: nonfinite float")
    return result


def _complex_value(value: Any, label: str) -> complex:
    if isinstance(value, Mapping):
        if set(value) != {"real", "imag"}:
            raise BoundaryError(f"{label}: complex object needs real/imag")
        z = complex(_finite_float(value["real"], label),
                    _finite_float(value["imag"], label))
    elif isinstance(value, (int, float)):
        z = complex(_finite_float(value, label))
    else:
        raise BoundaryError(f"{label}: complex entries must be real/imag objects")
    return z


def _complex_array(value: Any, shape: tuple[int, ...], label: str) -> np.ndarray:
    if isinstance(value, Mapping) and set(value) == {"real", "imag"}:
        try:
            arr = np.asarray(value["real"], dtype=float) + 1j * np.asarray(value["imag"], dtype=float)
        except (TypeError, ValueError) as exc:
            raise BoundaryError(f"{label}: invalid complex array") from exc
    else:
        try:
            arr = np.asarray(value, dtype=float).astype(complex)
        except (TypeError, ValueError) as exc:
            raise BoundaryError(f"{label}: invalid real array") from exc
    if arr.shape != shape or not np.isfinite(arr).all():
        raise BoundaryError(f"{label}: expected finite shape {shape}, got {arr.shape}")
    return arr


def encode_matrix(z: np.ndarray) -> dict[str, list]:
    z = np.asarray(z, dtype=complex)
    return {"real": z.real.tolist(), "imag": z.imag.tolist()}


def _index(v: Any, n: int, label: str) -> int:
    if type(v) is not int or not 1 <= v <= n:
        raise BoundaryError(f"{label} must be integer in 1..{n}")
    return v - 1


def _all_perms(key: tuple[int, ...]) -> set[tuple[int, ...]]:
    return set(permutations(key))


def _write_orbit(array: np.ndarray, members: set[tuple[int, ...]], value: complex,
                 occupied: dict[tuple[int, ...], complex], label: str) -> None:
    for key in members:
        if key in occupied and not np.isclose(occupied[key], value, rtol=1e-11, atol=1e-24):
            raise BoundaryError(f"{label}: conflicting symmetry assignments at {key}")
        array[key] = value
        occupied[key] = value


def _read_quartic(components: list, ns: int) -> np.ndarray:
    lam = np.zeros((ns,) * 4, dtype=float)
    used: dict[tuple[int, ...], complex] = {}
    if not isinstance(components, list):
        raise BoundaryError("quartic_components must be a list (possibly empty)")
    for k, row in enumerate(components):
        inds = row.get("indices")
        if not isinstance(inds, list) or len(inds) != 4:
            raise BoundaryError("quartic indices must have length four")
        key = tuple(_index(i, ns, f"quartic[{k}]") for i in inds)
        val = _complex_value(row["value"], f"quartic[{k}]")
        if abs(val.imag) > 1e-14 * max(1, abs(val.real)):
            raise BoundaryError("lambda_abcd is real in the real-scalar basis")
        _write_orbit(lam, _all_perms(key), val.real, used, "quartic")
    return lam


def _read_yukawa(components: list, nf: int, ns: int) -> np.ndarray:
    y = np.zeros((nf, nf, ns), dtype=complex)
    used: dict[tuple[int, ...], complex] = {}
    if not isinstance(components, list):
        raise BoundaryError("yukawa_components must be a list (possibly empty)")
    for k, row in enumerate(components):
        inds = row.get("indices")
        if not isinstance(inds, list) or len(inds) != 3:
            raise BoundaryError("Yukawa indices must have length three")
        i, j, a = (_index(inds[0], nf, "y.i"),
                   _index(inds[1], nf, "y.j"), _index(inds[2], ns, "y.a"))
        z = _complex_value(row["value"], f"yukawa[{k}]")
        _write_orbit(y, {(i, j, a), (j, i, a)}, z, used, "Yukawa")
    return y


def _read_wilson(components: list, nf: int, ns: int, label: str) -> np.ndarray:
    c = np.zeros((nf, nf, ns, ns), dtype=complex)
    used: dict[tuple[int, ...], complex] = {}
    if not isinstance(components, list) or not components:
        raise BoundaryError(f"{label}: nonempty canonical Wilson components required")
    for k, row in enumerate(components):
        inds = row.get("indices")
        if not isinstance(inds, list) or len(inds) != 4:
            raise BoundaryError("Wilson indices must have length four")
        i, j = (_index(inds[0], nf, "C.i"), _index(inds[1], nf, "C.j"))
        a, b = (_index(inds[2], ns, "C.a"), _index(inds[3], ns, "C.b"))
        z = _complex_value(row["value"], f"{label}.components[{k}]")
        orbit = {(fi, fj, sa, sb) for fi, fj in ((i, j), (j, i))
                 for sa, sb in ((a, b), (b, a))}
        _write_orbit(c, orbit, z, used, label)
    return c


@dataclass
class Context:
    model: str
    nf: int
    ns: int
    leptons: tuple[tuple[int, int], ...]
    higgs: tuple[int, int, int, int]
    mu_f: float
    mu_s: float
    yukawa: np.ndarray
    quartic: np.ndarray
    gauge: tuple[tuple[float, tuple[tuple[np.ndarray, np.ndarray], ...]], ...]
    operators: tuple[tuple[str, np.ndarray], ...]
    authority: str
    certificates: dict[str, bool]


def load_boundary(payload: Mapping[str, Any], *, allow_synthetic: bool = False) -> Context:
    if not isinstance(payload, Mapping) or payload.get("format") != FORMAT:
        raise BoundaryError("Expected certified master-Weyl boundary format; legacy ordered tensors forbidden")
    model = payload.get("model")
    if model not in MODELS:
        raise BoundaryError(f"Unknown model: {model}")
    auth = payload.get("authority")
    if auth not in ("certified", "synthetic"):
        raise BoundaryError("authority must be certified or synthetic")
    if auth == "synthetic" and not allow_synthetic:
        raise BoundaryError("Synthetic inputs require explicit allow_synthetic=True")
    certificates = payload.get("certificates")
    if not isinstance(certificates, Mapping):
        raise BoundaryError("Missing explicit certification dictionary")
    certified = {name: certificates.get(name) is True for name in CERTIFICATES}
    if auth == "certified" and not all(certified.values()):
        raise BoundaryError("Physical boundary missing certifications: " +
                            ", ".join(k for k, v in certified.items() if not v))
    if payload.get("scheme") != "MSbar":
        raise BoundaryError("Only explicit MSbar master-RGE convention supported")
    if payload.get("master_lagrangian_convention") != "C_ijab psi_i psi_j phi_a phi_b + hc":
        raise BoundaryError("Unrecognized master Wilson normalization")
    if payload.get("weinberg_lagrangian_convention") != "(1/2) C5_ij (L_i.H)(L_j.H) + hc":
        raise BoundaryError("Unrecognized symmetric C5 normalization")
    nf, ns = payload.get("n_fermions"), payload.get("n_real_scalars")
    if type(nf) is not int or nf < 6 or nf > 100:
        raise BoundaryError("n_fermions must be integer >= 6")
    if type(ns) is not int or ns < 4 or ns > 100:
        raise BoundaryError("n_real_scalars must be integer >= 4")
    raw_leptons = payload.get("lepton_components")
    if not isinstance(raw_leptons, list) or len(raw_leptons) != 3:
        raise BoundaryError("Three [nu,e] lepton flavour pairs required")
    leptons = tuple(tuple(_index(i, nf, "lepton index") for i in pair)
                    for pair in raw_leptons if isinstance(pair, list) and len(pair) == 2)
    if len(leptons) != 3 or len(set(sum((list(x) for x in leptons), []))) != 6:
        raise BoundaryError("Lepton Weyl indices must be six unique components")
    raw_h = payload.get("higgs_real_components")
    if not isinstance(raw_h, list) or len(raw_h) != 4:
        raise BoundaryError("Four [ReH+, ImH+, ReH0, ImH0] indices required")
    higgs = tuple(_index(i, ns, "Higgs index") for i in raw_h)
    if len(set(higgs)) != 4:
        raise BoundaryError("Higgs indices must be unique")
    blocks = payload.get("scalar_field_blocks")
    if not isinstance(blocks, list) or len(blocks) < 3:
        raise BoundaryError("Require explicit Higgs and BOTH heavy scalar field blocks")
    used_scalar_indices: set[int] = set()
    field_names: set[str] = set()
    heavy_count, higgs_count = 0, 0
    for block in blocks:
        if not isinstance(block, Mapping):
            raise BoundaryError("Malformed scalar field block")
        name, role, entries = block.get("name"), block.get("role"), block.get("indices")
        if not isinstance(name, str) or not name or name in field_names:
            raise BoundaryError("Scalar field names must be unique and nonempty")
        if role not in ("Higgs", "heavy"):
            raise BoundaryError("Scalar field role must be Higgs or heavy")
        if not isinstance(entries, list) or not entries:
            raise BoundaryError("Scalar field block needs nonempty component indices")
        inds = [_index(v, ns, name) for v in entries]
        if len(set(inds)) != len(inds) or set(inds) & used_scalar_indices:
            raise BoundaryError("Real scalar component repeated between field blocks")
        if role == "Higgs":
            higgs_count += 1
            if set(inds) != set(higgs):
                raise BoundaryError("Higgs scalar block does not match Higgs projector")
        else:
            heavy_count += 1
        used_scalar_indices.update(inds)
        field_names.add(name)
    if higgs_count != 1 or heavy_count < 2 or used_scalar_indices != set(range(ns)):
        raise BoundaryError("Higgs and BOTH heavy scalar fields must cover the real-scalar basis")
    mu_f = _finite_float(payload.get("mu_F_gev"), "mu_F")
    mu_s = _finite_float(payload.get("mu_S_gev"), "mu_S")
    if not (mu_f > mu_s > 0):
        raise BoundaryError("The fixed ordering requires 0 < mu_S < mu_F")
    y = _read_yukawa(payload.get("yukawa_components"), nf, ns)
    lam = _read_quartic(payload.get("quartic_components"), ns)
    gauge = []
    for k, sector in enumerate(payload.get("gauge_sectors", [])):
        g = _finite_float(sector.get("coupling"), f"gauge[{k}].coupling")
        generators = sector.get("generators")
        if not isinstance(generators, list) or not generators:
            raise BoundaryError("Every declared gauge sector requires complete generator set")
        gens = []
        for n, pair in enumerate(generators):
            f = _complex_array(pair["fermion"], (nf, nf), "fermion generator")
            s = _complex_array(pair["scalar"], (ns, ns), "scalar generator")
            if not np.allclose(f, f.conj().T, atol=1e-10, rtol=1e-10):
                raise BoundaryError(f"Non-Hermitian fermion generator {k}/{n}")
            if not np.allclose(s, s.conj().T, atol=1e-10, rtol=1e-10):
                raise BoundaryError(f"Non-Hermitian scalar generator {k}/{n}")
            gens.append((f, s))
        gauge.append((g, tuple(gens)))
    if not isinstance(payload.get("gauge_sectors"), list):
        raise BoundaryError("gauge_sectors must be explicitly enumerated (empty only if synthetic)")
    if auth == "certified":
        if not gauge:
            raise BoundaryError("Physical SM EFT1 needs explicit SU(2) and U(1) generators")
        raw_sectors = payload["gauge_sectors"]
        named = {row.get("name"): k for k, row in enumerate(raw_sectors)}
        if len(raw_sectors) != 2 or set(named) != {"SU2L", "U1Y"}:
            raise BoundaryError("Certified EFT1 needs exactly named SU2L and U1Y gauge sectors")
        su2 = gauge[named["SU2L"]]
        u1 = gauge[named["U1Y"]]
        if len(su2[1]) != 3 or len(u1[1]) != 1 or su2[0] <= 0 or u1[0] <= 0:
            raise BoundaryError("SU2L needs 3 generators, U1Y needs 1, both with positive g")
        for which in (0, 1):
            T = [x[which] for x in su2[1]]
            U = u1[1][0][which]
            for a, b, z in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
                if not np.allclose(T[a]@T[b]-T[b]@T[a], 1j*T[z], rtol=1e-10, atol=1e-10):
                    raise BoundaryError("SU2L generators fail [T1,T2]=iT3 algebra")
            for g in T:
                if not np.allclose(g@U-U@g, 0, rtol=1e-10, atol=1e-10):
                    raise BoundaryError("U1Y must commute with SU2L")
        # Validate the SM Higgs and lepton blocks against the checked repository
        # representation conventions, not just abstract gauge invariance.
        sigma = [np.array([[0, 1], [1, 0]], complex)/2,
                 np.array([[0, -1j], [1j, 0]], complex)/2,
                 np.diag([0.5, -0.5]).astype(complex)]
        for a, f in enumerate(sigma):
            for nu, elec in leptons:
                sub = su2[1][a][0][np.ix_([nu, elec], [nu, elec])]
                if not np.allclose(sub, f, rtol=1e-10, atol=1e-10):
                    raise BoundaryError("Lepton SU2L generators violate source doublet convention")
        uf = u1[1][0][0]
        for nu, elec in leptons:
            if not np.allclose(uf[np.ix_([nu, elec], [nu, elec])],
                               -0.5*np.eye(2), rtol=1e-10, atol=1e-10):
                raise BoundaryError("Lepton U1Y must be -1/2 in the source convention")
        h = list(higgs)
        uh = u1[1][0][1]
        expected = np.zeros((4,4), complex)
        expected[0,1] = expected[2,3] = 0.5j
        expected[1,0] = expected[3,2] = -0.5j
        if not np.allclose(uh[np.ix_(h,h)], expected, rtol=1e-10, atol=1e-10):
            raise BoundaryError("Higgs U1Y=+1/2 real-component convention violated")
        for a,t in enumerate(sigma):
            theta = np.zeros((4,4), complex)
            for p in range(2):
                for q in range(2):
                    x = t[p,q]
                    theta[2*p,2*q] += 1j*x.imag
                    theta[2*p,2*q+1] += 1j*x.real
                    theta[2*p+1,2*q] -= 1j*x.real
                    theta[2*p+1,2*q+1] += 1j*x.imag
            if not np.allclose(su2[1][a][1][np.ix_(h,h)], theta, rtol=1e-10, atol=1e-10):
                raise BoundaryError("Higgs SU2L real-component convention violated")
    raw_ops = payload.get("operators")
    if not isinstance(raw_ops, list) or not raw_ops:
        raise BoundaryError("At least one tree-level operator source required")
    names: set[str] = set()
    ops = []
    for op in raw_ops:
        name = op.get("source_id")
        if not isinstance(name, str) or not name or name in names:
            raise BoundaryError("Unique nonblank operator source_ids required")
        names.add(name)
        ops.append((name, _read_wilson(op.get("components"), nf, ns, name)))
    return Context(str(model), nf, ns, leptons, higgs, mu_f, mu_s,
                   y, lam, tuple(gauge), tuple(ops), str(auth), certified)


def _collinear(context: Context) -> tuple[np.ndarray, np.ndarray]:
    """Exact repository gamma_c conventions, with real scalars and Weyl fermions."""
    y = context.yukawa
    gs = np.zeros((context.ns, context.ns), dtype=complex)
    gf = np.zeros((context.nf, context.nf), dtype=complex)
    for coupling, generators in context.gauge:
        for f, s in generators:
            gs -= 4 * coupling**2 * (s @ s)
            gf -= 3 * coupling**2 * (f @ f)
    # gamma_s[a,b] = 1/2 sum_ij (y[j,i,a] y*[j,i,b] + y[j,i,b] y*[j,i,a]).
    gram = np.einsum("ija,ijb->ab", y, y.conj(), optimize=True)
    gs += 0.5 * (gram + gram.T)
    gf += 0.5 * np.einsum("ika,jka->ij", y, y.conj(), optimize=True)
    return gs, gf


def calculate_master_beta(context: Context, c: np.ndarray) -> dict[str, np.ndarray]:
    """All seven terms of the current general WilsonTensorRGE, eq. (4.85).

    Returns 16*pi**2 dC/dln(mu) for *every* tensor component.  This also
    exposes non-Weinberg components for the one-loop closure validation.
    """
    y, lam = context.yukawa, context.quartic
    cc = c.conj()
    dst = {}
    ycontract = np.einsum("kld,klab->dab", y.conj(), c, optimize=True)
    dst["yukawa_wavefunction"] = np.einsum("ijd,dab->ijab", y, ycontract, optimize=True)

    s = np.einsum("abcd,ijcd->ijab", lam, c, optimize=True)
    for coupling, generators in context.gauge:
        for _, theta in generators:
            s -= 2 * coupling**2 * np.einsum("ac,bd,ijcd->ijab", theta, theta, c,
                                             optimize=True)
    dst["scalar_pair"] = 0.5 * (s + s.swapaxes(2, 3))

    base = (2 * np.einsum("jkc,klb,ilac->ijab", y, y.conj(), c, optimize=True)
            + np.einsum("jkb,klc,ilac->ijab", y, y.conj(), c, optimize=True))
    for coupling, generators in context.gauge:
        for tf, ts in generators:
            base += 4 * coupling**2 * np.einsum("lj,bc,ilac->ijab", tf, ts, c,
                                                  optimize=True)
    dst["mixed_yukawa_gauge"] = (base + base.swapaxes(0, 1)
                                  + base.swapaxes(2, 3)
                                  + base.swapaxes(0, 1).swapaxes(2, 3))

    crossed = np.einsum("jla,ikc,lkbc->ijab", y, y, c, optimize=True)
    dst["crossed_yukawa"] = crossed + crossed.transpose(1, 0, 3, 2)

    dst["conjugate_coefficient_yukawa"] = (
        2 * np.einsum("ild,jkd,klab->ijab", y, y, cc, optimize=True)
        + np.einsum("ijd,kld,klab->ijab", y, y, cc, optimize=True))

    gamma_s, gamma_f = _collinear(context)
    gam_s = np.einsum("bd,ijad->ijab", gamma_s, c, optimize=True)
    dst["scalar_anomalous_dimension"] = gam_s + gam_s.swapaxes(2, 3)
    gam_f = np.einsum("jk,ikab->ijab", gamma_f, c, optimize=True)
    dst["fermion_anomalous_dimension"] = gam_f + gam_f.swapaxes(0, 1)
    return dst


def _weinberg_template(c5: np.ndarray, context: Context) -> np.ndarray:
    """Master C_ijab equivalent to (1/2) C5_ij (L.H)(L.H) + h.c."""
    n = context.nf
    ns = context.ns
    t = np.zeros((n, n, ns, ns), dtype=complex)
    hp_r, hp_i, h0_r, h0_i = context.higgs
    mapping = {}
    for p, (nu, elec) in enumerate(context.leptons):
        vnu = np.zeros(ns, dtype=complex)
        ve = np.zeros(ns, dtype=complex)
        vnu[h0_r], vnu[h0_i] = 1 / np.sqrt(2), 1j / np.sqrt(2)
        ve[hp_r], ve[hp_i] = -1 / np.sqrt(2), -1j / np.sqrt(2)
        mapping[nu] = (p, vnu)
        mapping[elec] = (p, ve)
    for i, (p, u) in mapping.items():
        for j, (q, v) in mapping.items():
            t[i, j] = 0.25 * c5[p, q] * (np.outer(u, v) + np.outer(v, u))
    return t


def project_c5(beta: np.ndarray, context: Context) -> tuple[np.ndarray, float, float]:
    """Extract complex symmetric C5, validate its entire LLHH SU(2) orbit.

    The projection is only physical when the total beta is Weinberg-shaped.
    Independent individual diagram terms need *not* pass this test.
    """
    hr = context.higgs[2]
    matrix = 4 * np.array([[beta[p[0], q[0], hr, hr]
                             for q in context.leptons] for p in context.leptons])
    symmetric_residual = float(np.max(np.abs(matrix - matrix.T)))
    matrix = 0.5 * (matrix + matrix.T)
    expected = _weinberg_template(matrix, context)
    ext = [i for pair in context.leptons for i in pair]
    hi = list(context.higgs)
    actual = beta[np.ix_(ext, ext, hi, hi)]
    target = expected[np.ix_(ext, ext, hi, hi)]
    residual = float(np.max(np.abs(actual - target)))
    return matrix, residual, symmetric_residual


def solve_boundary(payload: Mapping[str, Any], *, allow_synthetic: bool = False,
                   ward_rtol: float = 2e-9, ward_atol: float = 1e-26,
                   trajectory_points: int = 21) -> dict:
    """Fixed-order one-loop log, retaining every source and every master term.

    Does not include finite simultaneous-scalar matching or claim a neutrino
    prediction without that matching.  No self-running LLSS iterates are made.
    """
    context = load_boundary(payload, allow_synthetic=allow_synthetic)
    if type(trajectory_points) is not int or not 2 <= trajectory_points <= 1001:
        raise BoundaryError("trajectory_points must be integer in 2..1001")
    log_ratio = log(context.mu_s / context.mu_f)
    factor = LOOP * log_ratio
    operator_results = []
    beta_sum = np.zeros((3, 3), dtype=complex)
    beta_tensor_sum = np.zeros((context.nf, context.nf, context.ns, context.ns), dtype=complex)
    per_term_sum = {name: np.zeros((3, 3), dtype=complex) for name in TERMS}
    for source_id, c in context.operators:
        beta_terms = calculate_master_beta(context, c)
        term_matrices = {}
        total_tensor = np.zeros_like(c)
        for name in TERMS:
            b = beta_terms[name]
            if not np.all(np.isfinite(b)):
                raise BoundaryError(f"{source_id}/{name}: nonfinite beta")
            projected, _, _ = project_c5(b, context)
            per_term_sum[name] += projected
            term_matrices[name] = {
                "beta_C5_16pi2": encode_matrix(projected),
                "delta_C5_log_gev_inverse": encode_matrix(factor * projected),
                "order": "tree_D5_times_one_loop_renormalizable_couplings",
            }
            total_tensor += b
        projected, ward, symmetry = project_c5(total_tensor, context)
        term_error = float(np.max(np.abs(projected - sum((
            _complex_array(term_matrices[name]["beta_C5_16pi2"], (3, 3), name)
            for name in TERMS), np.zeros((3, 3), complex)))))
        if term_error > 1e-9 * max(1e-20, float(np.max(abs(projected)))):
            raise RuntimeError("Internal diagram decomposition fails recomposition")
        beta_sum += projected
        beta_tensor_sum += total_tensor
        operator_results.append({
            "source_id": source_id,
            "beta_C5_16pi2": encode_matrix(projected),
            "delta_C5_log_gev_inverse": encode_matrix(factor * projected),
            "weinberg_subspace_absolute_residual": ward,
            "flavour_symmetric_absolute_residual": symmetry,
            "per_rge_term": term_matrices,
        })
    ward_tolerance = ward_atol + ward_rtol * max(float(np.max(abs(beta_sum))), 1e-30)
    total_matrix_again, total_ward_error, total_symmetry_error = project_c5(beta_tensor_sum, context)
    if not np.allclose(beta_sum, total_matrix_again, rtol=1e-12, atol=1e-28):
        raise RuntimeError("Total tensor and per-source projections disagree")
    physical = (context.authority == "certified" and all(context.certificates.values())
                and total_ward_error <= ward_tolerance)
    if context.authority == "certified" and not physical:
        raise BoundaryError("Certified source does not close into gauge-invariant Weinberg tensor; "
                            f"max Ward residual {total_ward_error:g} > {ward_tolerance:g}")
    expected_sum = sum((per_term_sum[name] for name in TERMS), np.zeros((3, 3), complex))
    if not np.allclose(beta_sum, expected_sum, rtol=1e-11, atol=1e-28):
        raise RuntimeError("Summed source and term channels inconsistent")
    return {
        "format": RESULT_FORMAT, "model": context.model,
        "status": "CertifiedLeadingOneLoopLog_RequiresFiniteJointScalarMatching" if physical
                  else "SyntheticDiagnostic_NotPhysicalC5",
        "authority": context.authority,
        "physics_completeness": {
            "psi2phi2_master_anomalous_dimension_terms_implemented": list(TERMS),
            "external_boundary_certificates": context.certificates,
            "physical_eft1_source_verified": physical,
            "joint_two_scalar_finite_matching_included": False,
            "complete_one_loop_neutrino_mass": False,
            "LLSS_self_running_iterated": False,
        },
        "scheme": "MSbar", "scale_convention": "mu*d/dmu = d/dln(mu)",
        "scales_gev": {"mu_F": context.mu_f, "mu_S": context.mu_s},
        "ln_muS_over_muF": log_ratio,
        "one_over_16pi2": LOOP,
        "operator_results": operator_results,
        "summed_by_rge_term_beta_C5_16pi2": {
            name: encode_matrix(per_term_sum[name]) for name in TERMS},
        "total_beta_C5_16pi2": encode_matrix(beta_sum),
        "total_delta_C5_log_gev_inverse": encode_matrix(factor * beta_sum),
        "trajectory": [{
            "mu_gev": float(context.mu_f * np.exp(log_ratio * k / (trajectory_points - 1))),
            "delta_C5_log_gev_inverse": encode_matrix(
                LOOP * (log_ratio * k / (trajectory_points - 1)) * beta_sum),
        } for k in range(trajectory_points)],
        "matching_interface": {
            "EFT1_C5_muS_equals_C5_muF_plus_log_piece": True,
            "d_delta_C5_log_dln_muS": encode_matrix(LOOP * beta_sum),
            "matching_must_subtract_EFT1_IR_log_once": True,
            "finite_joint_scalar_matching_provided_by_this_module": False,
            "finite_matching_scheme": "must use same MSbar master operator conventions",
            "note": "No finite term is reconstructed from the UV amplitude. "
                    "To avoid double-counting, the finite scalar matching provider "
                    "must specify its EFT1 infrared-log subtraction explicitly.",
        },
        "ward_residual_max": total_ward_error,
        "flavour_symmetry_residual_max": total_symmetry_error,
        "ward_validation_tolerance": ward_tolerance,
        "C5_to_neutrino_mass": "m_nu = -(v^2/2) C5",
    }
