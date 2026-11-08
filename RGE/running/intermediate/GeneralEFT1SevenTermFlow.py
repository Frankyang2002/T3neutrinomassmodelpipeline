"""Matrix-free, real-linear seven-term psi^2 phi^2 EFT1 Wilson evolution.

Source baseline: WilsonTensorRGE.py / GeneralEFT1C5Mixing.py at
12128de602ac6bf109dbad50dbac8446ed629434.  The four independent crossed exchanges transcribe the unstarred-C formula
in arXiv:2512.16890v1 Eq.(4.85), but this formula FAILS gauge covariance
for general complex Weyl representations. The separately selectable
``ward_covariant_candidate`` variant uses C* and transposed fermion Casimir.
``fonseca_eq22_six_term`` omits the crossed term and retains the latter
Casimir fix. All three modes are explicitly UNVERIFIED as complete physical
A(-4) anomalous dimensions pending the version-of-record reconciliation.

L5 = C_ijab psi_i psi_j phi_a phi_b + h.c. (no factorial),
V4 = lambda_abcd phi_a phi_b phi_c phi_d / 4! ;
16*pi^2 dC/dln(mu) = sum_{seven terms} beta^(n).

All arrays are in the canonical left-Weyl/real-scalar basis, 0-based in memory;
JSON index lists are 1-based. Fermion/scalar pairs are *separately* symmetric.
No scalar matching and no physical certification are performed by this engine.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import log, pi
from typing import Any, Mapping, Sequence

import numpy as np
from scipy.integrate import solve_ivp

LOOP = 1.0 / (16 * pi**2)
TERMS = ("yukawa_wavefunction", "scalar_pair", "mixed_yukawa_gauge",
         "crossed_yukawa", "conjugate_coefficient_yukawa",
         "scalar_anomalous_dimension", "fermion_anomalous_dimension")
COUPLING_FORMAT = "t3-eft1-seven-term-coupling-samples-v1"
RESULT_FORMAT = "t3-a-m4-full-seven-term-wilson-trajectory-v1"


def _finite(z: np.ndarray, description: str) -> np.ndarray:
    if not np.isfinite(z).all():
        raise ValueError(description + " contains non-finite entries")
    return z


@dataclass(frozen=True)
class Couplings:
    """Complete master-RGE tensors at one logarithmic scale.

    Each gauge sector is (g, ((t_A,theta_A),...)); t and theta are
    Hermitian generators in the Weyl and *real* scalar representations.
    SU(3) must be included for a physical quark-containing EFT1 basis.
    """
    yukawa: np.ndarray             # (nf,nf,ns), y_ija=y_jia
    quartic: np.ndarray            # (ns,ns,ns,ns), fully symmetric real
    gauge: tuple                   # ((g, ((t_A,theta_A), ...)), ...)

    def __post_init__(self):
        y, lam = np.asarray(self.yukawa), np.asarray(self.quartic)
        if y.ndim != 3 or y.shape[0] != y.shape[1] or lam.shape != (y.shape[2],)*4:
            raise ValueError("Yukawa/quartic tensor shapes incompatible")
        if not np.allclose(y, y.swapaxes(0, 1), atol=1e-13, rtol=1e-11):
            raise ValueError("Yukawa fermion exchange symmetry violated")
        if not np.allclose(lam, lam.real, atol=1e-13, rtol=1e-11):
            raise ValueError("Real-scalar quartic must be real")
        for axis in (1, 2, 3):
            if not np.allclose(lam, np.swapaxes(lam, 0, axis), atol=1e-13, rtol=1e-11):
                raise ValueError("Quartic S4 symmetry violated")
        _finite(y, "yukawa"); _finite(lam, "quartic")
        for g, pairs in self.gauge:
            if not np.isfinite(g) or not pairs:
                raise ValueError("Gauge sector requires finite coupling and generators")
            for t, theta in pairs:
                if np.shape(t) != (y.shape[0],)*2 or np.shape(theta) != (y.shape[2],)*2:
                    raise ValueError("Gauge generator dimensions disagree with basis")
                if not np.allclose(t, np.asarray(t).conj().T, atol=1e-11):
                    raise ValueError("Non-Hermitian fermion generator")
                if not np.allclose(theta, np.asarray(theta).conj().T, atol=1e-11):
                    raise ValueError("Non-Hermitian real-scalar generator")


@dataclass(frozen=True)
class TensorBasis:
    nf: int
    ns: int
    leptons: tuple[tuple[int, int], ...]  # 0-based (nu,charged), 3 generations
    higgs: tuple[int, int, int, int]       # 0-based (Re+, Im+, Re0, Im0)

    def __post_init__(self):
        if self.nf < 6 or self.ns < 4 or len(self.leptons) != 3 or len(self.higgs) != 4:
            raise ValueError("Require 3-generation leptons and four Higgs real components")
        fi = [x for pair in self.leptons for x in pair]
        if len(set(fi)) != 6 or not all(0 <= x < self.nf for x in fi):
            raise ValueError("Invalid lepton component map")
        if len(set(self.higgs)) != 4 or not all(0 <= x < self.ns for x in self.higgs):
            raise ValueError("Invalid Higgs component map")

    @property
    def independent_size(self):
        return (self.nf * (self.nf + 1)//2) * (self.ns * (self.ns + 1)//2)


class SymmetricCoordinates:
    """Complete, not lepton-truncated, separate-pair-symmetric coordinates."""
    def __init__(self, basis: TensorBasis):
        self.basis = basis
        self.fi, self.fj = np.triu_indices(basis.nf)
        self.sa, self.sb = np.triu_indices(basis.ns)
        self.shape = (len(self.fi), len(self.sa))
        self.size = int(np.prod(self.shape))

    def pack(self, c: np.ndarray) -> np.ndarray:
        return np.asarray(c)[self.fi[:, None], self.fj[:, None],
                             self.sa[None, :], self.sb[None, :]].reshape(-1)

    def unpack(self, z: np.ndarray) -> np.ndarray:
        mat = np.asarray(z, dtype=np.complex128).reshape(self.shape)
        n, s = self.basis.nf, self.basis.ns
        result = np.zeros((n, n, s, s), dtype=np.complex128)
        for f, j in ((self.fi, self.fj), (self.fj, self.fi)):
            for a, b in ((self.sa, self.sb), (self.sb, self.sa)):
                result[f[:, None], j[:, None], a[None, :], b[None, :]] = mat
        return result

    def labels(self) -> np.ndarray:
        f, j = np.repeat(self.fi, len(self.sa)), np.repeat(self.fj, len(self.sa))
        a, b = np.tile(self.sa, len(self.fi)), np.tile(self.sb, len(self.fi))
        return np.stack((f, j, a, b), axis=1).astype(np.int32) + 1


def symmetry_residual(c: np.ndarray) -> float:
    denom = max(float(np.max(np.abs(c))), 1e-300)
    return float(max(np.max(np.abs(c-c.swapaxes(0, 1))),
                     np.max(np.abs(c-c.swapaxes(2, 3)))) / denom)


def gauge_ward_residual(c: np.ndarray, t: np.ndarray, theta: np.ndarray) -> float:
    """Infinitesimal gauge-invariance residual for C_ijab in master basis.

    Gauge action on the four fields is t^T on each Weyl index and theta^T
    on each real-scalar index. Returns relative max-norm residual.
    """
    variation = (np.einsum("ki,kjab->ijab",t,c,optimize=True)
                 + np.einsum("kj,ikab->ijab",t,c,optimize=True)
                 + np.einsum("ca,ijcb->ijab",theta,c,optimize=True)
                 + np.einsum("cb,ijac->ijab",theta,c,optimize=True))
    return float(np.max(np.abs(variation))/max(np.max(np.abs(c)),1e-300))


def collinear(couplings: Couplings):
    """Exactly the AnomalousDimensions.py (A.2,A.3) contractions."""
    y = couplings.yukawa
    n, _, ns = y.shape
    gs = np.zeros((ns, ns), complex)
    gf = np.zeros((n, n), complex)
    for g, pairs in couplings.gauge:
        for t, theta in pairs:
            gs -= 4*g*g*(theta @ theta)
            gf -= 3*g*g*(t @ t)
    gram = np.einsum("ija,ijb->ab", y, y.conj(), optimize=True)
    gs += 0.5*(gram + gram.T)
    gf += 0.5*np.einsum("ika,jka->ij", y, y.conj(), optimize=True)
    return gs, gf


def master_beta(c: np.ndarray, couplings: Couplings, *, corrected_crossed: bool = True) -> dict[str, np.ndarray]:
    """Matrix-free full seven-term master tensor; no projected-only shortcut.

    corrected_crossed=False reproduces the repository's original two-exchange
    symbolic/numerical source *for audit only*, not for physical evolution.
    """
    y, lam = couplings.yukawa, couplings.quartic
    nf, _, ns = y.shape
    if np.shape(c) != (nf, nf, ns, ns):
        raise ValueError("Full Wilson tensor shape does not match couplings")
    _finite(c, "Wilson")
    conj = c.conj()
    contract = np.einsum
    dst = {}
    ycon = contract("kld,klab->dab", y.conj(), c, optimize=True)
    dst[TERMS[0]] = contract("ijd,dab->ijab", y, ycon, optimize=True)
    pair = contract("abcd,ijcd->ijab", lam, c, optimize=True)
    for g, pairs in couplings.gauge:
        for _, theta in pairs:
            pair -= 2*g*g*contract("ac,bd,ijcd->ijab", theta, theta, c, optimize=True)
    dst[TERMS[1]] = 0.5*(pair+pair.swapaxes(2, 3))
    base = (2*contract("jkc,klb,ilac->ijab", y, y.conj(), c, optimize=True)
            + contract("jkb,klc,ilac->ijab", y, y.conj(), c, optimize=True))
    for g, pairs in couplings.gauge:
        for tf, ts in pairs:
            base += 4*g*g*contract("lj,bc,ilac->ijab", tf, ts, c, optimize=True)
    dst[TERMS[2]] = (base + base.swapaxes(0, 1) + base.swapaxes(2, 3)
                     + base.swapaxes(0, 1).swapaxes(2, 3))
    crossed = contract("jla,ikc,lkbc->ijab", y, y, c, optimize=True)
    if corrected_crossed:
        dst[TERMS[3]] = (crossed + crossed.swapaxes(0, 1) + crossed.swapaxes(2, 3)
                         + crossed.swapaxes(0, 1).swapaxes(2, 3))
    else:
        dst[TERMS[3]] = crossed + crossed.transpose(1, 0, 3, 2)
    # NumPy greedy cannot find a useful path for the unreduced seven-index
    # contraction: explicit two-stage contraction cuts ~2e10 scalar products
    # at nf=45 to a sequence of BLAS-friendly tensor operations.
    anti_tmp = contract("ild,klab->ikabd", y, conj, optimize=True)
    dst[TERMS[4]] = (2*contract("jkd,ikabd->ijab", y, anti_tmp, optimize=True)
                     + contract("ijd,kld,klab->ijab", y, y, conj, optimize=True))
    gs, gf = collinear(couplings)
    base_s = contract("bd,ijad->ijab", gs, c, optimize=True)
    dst[TERMS[5]] = base_s + base_s.swapaxes(2, 3)
    base_f = contract("jk,ikab->ijab", gf, c, optimize=True)
    dst[TERMS[6]] = base_f + base_f.swapaxes(0, 1)
    return dst


def total_beta(terms: Mapping[str, np.ndarray]) -> np.ndarray:
    return sum((terms[n] for n in TERMS), np.zeros_like(terms[TERMS[0]]))


def c5_projection(c: np.ndarray, basis: TensorBasis) -> np.ndarray:
    """4 C[nu_i,nu_j,ReH0,ReH0], NOT final scalar-matched SMEFT C5."""
    nu = [p[0] for p in basis.leptons]
    h = basis.higgs[2]
    return 4*c[np.ix_(nu, nu, [h], [h])][:, :, 0, 0]


def _encode(z: np.ndarray) -> dict:
    return {"real": np.asarray(z).real.tolist(), "imag": np.asarray(z).imag.tolist()}


def interpolate(samples: Sequence[Couplings], log_scales: np.ndarray, t: float) -> Couplings:
    hi = max(1, min(len(samples)-1, int(np.searchsorted(log_scales, t, side="right"))))
    lo = hi-1
    w = float(np.clip((t-log_scales[lo])/(log_scales[hi]-log_scales[lo]), 0, 1))
    a, b = samples[lo], samples[hi]
    gauge = []
    for (ga, pa), (gb, pb) in zip(a.gauge, b.gauge):
        if len(pa) != len(pb) or any(not np.array_equal(x, y)
                                      for xpair, ypair in zip(pa, pb)
                                      for x, y in zip(xpair, ypair)):
            raise ValueError("Gauge basis must remain fixed across coupling samples")
        gauge.append(((1-w)*ga+w*gb, pa))
    return Couplings((1-w)*a.yukawa+w*b.yukawa,
                     (1-w)*a.quartic+w*b.quartic, tuple(gauge))


def evolve(basis: TensorBasis, initial: np.ndarray, mu_f: float, mu_s: float,
           sample_scales: Sequence[float], samples: Sequence[Couplings], *,
           points: int = 11, rtol: float = 1e-8, atol: float = 1e-18,
           method: str = "DOP853",
           beta_mode: str = "source_v1_unstarred") -> dict:
    """Run the complete ODE with log(mu)-interpolated running Y,lambda,g.

    Trajectory is complete in unique symmetric coordinates, not sparse support.
    The returned coordinate vector plus coordinates() reconstructs all C_ijab.
    """
    supported = ("source_v1_unstarred", "ward_covariant_candidate",
                 "fonseca_eq22_six_term")
    if beta_mode not in supported:
        raise ValueError("Unknown beta_mode: " + ", ".join(supported))
    if beta_mode == "ward_covariant_candidate":
        from RGE.running.intermediate.GeneralEFT1CrossedWardCandidate import candidate_beta
        evaluate_beta = candidate_beta
    elif beta_mode == "fonseca_eq22_six_term":
        from RGE.running.intermediate.GeneralEFT1SixTermReferenceCandidate import six_term_beta
        evaluate_beta = six_term_beta
    else:
        evaluate_beta = master_beta
    scales = np.asarray(sample_scales, dtype=float)
    if not np.isfinite(scales).all() or not np.isfinite([mu_f, mu_s]).all():
        raise ValueError("Nonfinite threshold or sample")
    if not 0 < mu_s < mu_f or not 2 <= points <= 201:
        raise ValueError("Require 0 < mu_S < mu_F and 2..201 points")
    if len(samples) != len(scales) or len(scales) < 2 or np.any(np.diff(scales) <= 0):
        raise ValueError("At least two strictly ascending running-coupling samples required")
    if scales[0] > mu_s*(1+1e-12) or scales[-1] < mu_f*(1-1e-12):
        raise ValueError("Coupling samples must bracket both thresholds")
    for s in samples:
        if np.shape(s.yukawa) != (basis.nf,basis.nf,basis.ns):
            raise ValueError("Couplings do not cover the *full* EFT1 fermion basis")
    coords = SymmetricCoordinates(basis)
    initial = np.asarray(initial, dtype=complex)
    if initial.shape != (basis.nf,basis.nf,basis.ns,basis.ns):
        raise ValueError("Canonical boundary must contain entire tensor")
    if symmetry_residual(initial) > 1e-11:
        raise ValueError("Canonical boundary lacks separate pair symmetries")
    z0 = coords.pack(_finite(initial, "initial"))
    t0, t1, x = log(mu_f), log(mu_s), np.log(scales)

    def rhs(t, z):
        c = coords.unpack(z)
        cb = interpolate(samples, x, t)
        return LOOP*coords.pack(total_beta(evaluate_beta(c, cb)))

    # Solve exactly up to every interpolation knot. An adaptive solver that
    # crosses a kink in dY/dln(mu), d(lambda)/dln(mu), etc. can otherwise
    # produce step-size-dependent source-recomposition errors.
    ts = np.linspace(t0,t1,points)
    knots = [t0]+[float(q) for q in x if t1 < q < t0][::-1]+[t1]
    if np.any(np.diff(knots) >= 0):
        raise ValueError("Coupling interpolation knots must be descending")
    trajectory = np.empty((points,coords.size),dtype=complex)
    filled = np.zeros(points,dtype=bool)
    zz = z0.copy()
    evaluations=0
    messages=[]
    for upper,lower in zip(knots[:-1],knots[1:]):
        mask = (ts <= upper+1e-12) & (ts >= lower-1e-12) & (~filled)
        selected = np.flatnonzero(mask)
        requested = ts[selected].tolist()
        # A final lower-knot evaluation is essential even when no output
        # time lies in the current interpolation interval.
        if not requested or abs(requested[-1]-lower)>1e-12:
            requested.append(lower)
        requested = [min(upper,max(lower,float(q))) for q in requested]
        requested = sorted(set(requested),reverse=True)
        sub = solve_ivp(rhs,(upper,lower),zz,t_eval=requested,
                        rtol=rtol,atol=atol,method=method)
        if not sub.success or not np.isfinite(sub.y).all():
            raise RuntimeError("Seven-term Wilson ODE failed: " + sub.message)
        zz = sub.y[:,-1].copy()
        evaluations += sub.nfev
        messages.append(sub.message)
        for idx in selected:
            hits = np.flatnonzero(np.abs(sub.t-ts[idx])<1e-11)
            if len(hits)!=1:
                raise RuntimeError("Requested output sample was not integrated")
            trajectory[idx] = sub.y[:,hits[0]]
            filled[idx] = True
    if not np.all(filled):
        raise RuntimeError("Internal trajectory grid incomplete")
    _finite(trajectory, "trajectory")
    end = coords.unpack(trajectory[-1])
    terms0 = evaluate_beta(initial, interpolate(samples, x, t0))
    terms1 = evaluate_beta(end, interpolate(samples, x, t1))
    linear = initial + LOOP*(t1-t0)*total_beta(terms0)
    projections, norms, symmetry = [], [], []
    for t, z in zip(ts, trajectory):
        c = coords.unpack(z)
        terms = evaluate_beta(c, interpolate(samples, x, t))
        projections.append({name: _encode(c5_projection(terms[name], basis)) for name in TERMS} |
                           {"total": _encode(c5_projection(total_beta(terms), basis)),
                            "induced_LLHH": _encode(c5_projection(c, basis))})
        norms.append({name: float(np.linalg.norm(terms[name])) for name in TERMS} |
                     {"total": float(np.linalg.norm(total_beta(terms)))})
        symmetry.append(symmetry_residual(total_beta(terms)))
    zero_start = np.flatnonzero(z0 == 0)
    threshold = max(float(np.max(np.abs(z0))), 1e-300)*1e-10
    generated = zero_start[np.abs(trajectory[-1,zero_start]) > threshold]
    order = generated[np.argsort(-np.abs(trajectory[-1, generated]))]
    labels = coords.labels()
    return {
        "beta_mode": beta_mode,
        "beta_physics_certified": False,
        "anti_linear_term_names": ([TERMS[3], TERMS[4]] if beta_mode == "ward_covariant_candidate"
                                   else [TERMS[4]]),
        "basis": basis, "coordinates": labels, "scales_gev": np.exp(ts),
        "trajectory": trajectory, "initial": initial, "final": end,
        "beta_initial": terms0, "beta_final": terms1,
        "projections": projections, "norms": norms,
        "symmetry_residuals": symmetry, "solver": {
            "method": method, "rtol": rtol, "atol": atol, "nfev": int(evaluations),
            "interpolation_segments":len(knots)-1,
            "success": True, "message": "; ".join(messages)},
        "mu_F_gev": float(mu_f), "mu_S_gev": float(mu_s),
        "fixed_order_full": linear,
        "fixed_order_difference_relative": float(np.linalg.norm(end-linear) /
                   max(np.linalg.norm(initial), 1e-300)),
        "generated_indices": generated, "generated_sorted": order,
        "generated_threshold": threshold,
        "initial_support_count": int(np.count_nonzero(z0)),
    }


def _complex_matrix(raw, shape, label):
    if isinstance(raw, Mapping):
        if set(raw) != {"real", "imag"}:
            raise ValueError(f"{label}: complex matrix requires real/imag")
        arr = np.asarray(raw["real"], dtype=float)+1j*np.asarray(raw["imag"], dtype=float)
    else:
        arr = np.asarray(raw, dtype=complex)
    if arr.shape != shape:
        raise ValueError(f"{label} expected {shape}, got {arr.shape}")
    return _finite(arr, label)


def load_coupling_samples(payload: Mapping[str,Any], boundary: Any, *,
                          allow_synthetic: bool = False) -> tuple[np.ndarray, list[Couplings]]:
    """Strict Chat B handoff: independent (running) coupling tensors.

    The static generators live in payload['gauge_generators'] and are never
    interpolated. Samples contain independent y,lambda components and g's.
    Input is deliberately separate from Chat A's canonical C boundary.
    """
    from RGE.running.intermediate.GeneralEFT1C5Mixing import _read_yukawa, _read_quartic

    if payload.get("format") != COUPLING_FORMAT or payload.get("model") != boundary.model:
        raise ValueError("Missing Batch34 coupling schema or model mismatch")
    if payload.get("scheme") != "MSbar" or payload.get("basis") != "canonical_left_Weyl_real_scalar":
        raise ValueError("Different renormalization/basis convention")
    if payload.get("n_fermions") != boundary.nf or payload.get("n_real_scalars") != boundary.ns:
        raise ValueError("Couplings do not match Chat A full basis")
    if payload.get("authority") not in ("synthetic", "certified"):
        raise ValueError("Couplings require explicit authority")
    if payload["authority"] == "synthetic" and not allow_synthetic:
        raise ValueError("Synthetic running must be explicitly allowed")
    if payload["authority"] == "certified":
        if boundary.authority != "certified" or payload.get("complete_general_EFT1_RGE_verified") is not True:
            raise ValueError("Certified coupling input missing independent verification")
    sectors = payload.get("gauge_generators")
    if not isinstance(sectors, list):
        raise ValueError("Static gauge_generators list is mandatory")
    static = []
    for sector in sectors:
        pairs = []
        for raw in sector.get("generators", []):
            f = _complex_matrix(raw["fermion"], (boundary.nf,boundary.nf), "fermion generator")
            s = _complex_matrix(raw["scalar"], (boundary.ns,boundary.ns), "scalar generator")
            pairs.append((f, s))
        if not pairs or not isinstance(sector.get("name"), str):
            raise ValueError("Each gauge group requires name and complete generator set")
        static.append((sector["name"], tuple(pairs)))
    if len({name for name,_ in static}) != len(static):
        raise ValueError("Duplicate gauge sectors")
    if payload["authority"] == "certified":
        if boundary.nf != 45 or boundary.ns != 12:
            raise ValueError("Physical A(-4) must use complete 45-Weyl, 12-real-scalar basis")
        if {name for name,_ in static} != {"SU3C", "SU2L", "U1Y"}:
            raise ValueError("Complete quark-containing EFT1 requires SU3C, SU2L, U1Y")
        for name,pairs in static:
            if len(pairs) != {"SU3C":8, "SU2L":3, "U1Y":1}[name]:
                raise ValueError("Wrong generator inventory for " + name)
            for tf,ts in pairs:
                if not np.allclose(ts, -ts.T, rtol=1e-10, atol=1e-10) or not np.allclose(ts.real, 0, atol=1e-10):
                    raise ValueError("Real-scalar gauge generators must be imaginary antisymmetric")
                if name == "SU3C" and np.max(np.abs(ts)) > 1e-10:
                    raise ValueError("A(-4) scalars must be SU(3)C singlets")
        # Check that distinct gauge factors commute. This is a necessary
        # representation test, not a substitute for Chat B's Ward proof.
        for k,(name,pp) in enumerate(static):
            for other,qq in static[k+1:]:
                for t0,s0 in pp:
                    for t1,s1 in qq:
                        if not np.allclose(t0@t1,t1@t0,rtol=1e-10,atol=1e-10) or not np.allclose(s0@s1,s1@s0,rtol=1e-10,atol=1e-10):
                            raise ValueError("Different gauge factors must commute")
    rows = payload.get("samples")
    if not isinstance(rows, list) or len(rows) < 2:
        raise ValueError("Coupling samples must be at least two scales; fixed high-scale freezing forbidden")
    scales, samples = [], []
    for row in rows:
        mu = float(row["mu_gev"])
        if not np.isfinite(mu) or mu <= 0:
            raise ValueError("Invalid coupling sample scale")
        g = row.get("gauge_couplings")
        if not isinstance(g, Mapping) or set(g) != {name for name,_ in static}:
            raise ValueError("Gauge couplings must list every static group")
        sectors_at_mu = []
        for name, pairs in static:
            value = float(g[name])
            if not np.isfinite(value) or value < 0:
                raise ValueError("Gauge coupling must be finite/nonnegative")
            sectors_at_mu.append((value, pairs))
        y = _read_yukawa(row["yukawa_components"], boundary.nf, boundary.ns)
        lam = _read_quartic(row["quartic_components"], boundary.ns)
        scales.append(mu)
        samples.append(Couplings(y, lam, tuple(sectors_at_mu)))
    scales = np.asarray(scales)
    if np.any(np.diff(scales) <= 0):
        raise ValueError("Coupling samples must be strictly ascending")
    if scales[0] > boundary.mu_s or scales[-1] < boundary.mu_f:
        raise ValueError("Running-coupling interval does not cover both thresholds")
    if payload.get("couplings_scale_dependent") is not True:
        raise ValueError("Running samples must explicitly declare scale dependence")
    return scales, samples


def run_contracts(boundary_payload: Mapping[str, Any], coupling_payload: Mapping[str,Any], *,
                  allow_synthetic: bool = False, **kwargs) -> dict:
    """Chat A boundary + Chat B trajectory contract, fail-closed on provenance.

    Existing certified-master-Weyl boundary parser supplies the independent
    Matchete-to-Weyl normalization preflight; this module does not redo it.
    """
    from RGE.running.intermediate.GeneralEFT1C5Mixing import load_boundary
    context = load_boundary(boundary_payload, allow_synthetic=allow_synthetic)
    decompose = bool(kwargs.pop("source_decomposition", False))
    beta_mode = kwargs.get("beta_mode", "source_v1_unstarred")
    if context.authority == "certified" or coupling_payload.get("authority") == "certified":
        raise ValueError("Physical A(-4) input is blocked: crossed term loop convention "
                         "is not independently certified; use synthetic diagnostics only")
    if context.model != "A_-4":
        raise ValueError("This Batch34 runner is A(-4)-scoped")
    basis = TensorBasis(context.nf, context.ns, context.leptons, context.higgs)
    # The canonical boundary (Chat A) and the full coupling basis (Chat B)
    # must refer to EXACTLY the same fields. A certified 45-Weyl handoff
    # must never silently replace 39 quark/charge-conjugate coordinates with 0.
    fermion_blocks = coupling_payload.get("fermion_field_blocks")
    scalar_blocks = coupling_payload.get("scalar_field_blocks")
    if fermion_blocks is not None:
        seen, names = [], set()
        for entry in fermion_blocks:
            if entry.get("name") in names or not entry.get("name"):
                raise ValueError("Duplicated or empty fermion block name")
            names.add(entry["name"])
            seen.extend(entry.get("indices",[]))
        if sorted(seen) != list(range(1,context.nf+1)):
            raise ValueError("Fermion block map must cover ALL Weyl components once")
    if scalar_blocks is not None and scalar_blocks != boundary_payload.get("scalar_field_blocks"):
        raise ValueError("Chat B scalar index map must match Chat A exactly")
    if coupling_payload.get("authority") == "certified":
        if fermion_blocks is None or scalar_blocks is None:
            raise ValueError("Certified Chat B handoff requires complete fermion/scalar index maps")
    c0 = sum((a for _,a in context.operators), np.zeros((context.nf,context.nf,context.ns,context.ns), complex))
    mu, samples = load_coupling_samples(coupling_payload, context, allow_synthetic=allow_synthetic)
    result = evolve(basis, c0, context.mu_f, context.mu_s, mu, samples, **kwargs)
    result["authority"] = "certified_inputs_unvalidated_physics" if (
        context.authority == "certified" and coupling_payload["authority"] == "certified") else "synthetic"
    result["source_ids"] = [name for name,_ in context.operators]
    result["published_formula_reference"] = (
        "Fonseca-Olgoso-Santiago JHEP07(2026)094 Eq.(22), six non-crossed terms"
        if result["beta_mode"] == "fonseca_eq22_six_term"
        else "arXiv:2512.16890v1 Eq.(4.85), crossed-term convention unresolved")
    result["unverified_loop_convention_blocker"] = True
    result["fermion_field_blocks"] = fermion_blocks or []
    result["scalar_field_blocks"] = boundary_payload.get("scalar_field_blocks",[])
    result["input_certificates"] = dict(context.certificates)
    result["coupling_source_verified"] = bool(coupling_payload.get("complete_general_EFT1_RGE_verified"))
    # This is a direct executable-source comparison with the unmodified
    # GeneralEFT1C5Mixing.calculate_master_beta, not a claim that its crossed
    # permutation convention is physically correct.
    from copy import copy
    from RGE.running.intermediate.GeneralEFT1C5Mixing import calculate_master_beta
    old_context = copy(context)
    at_mf = interpolate(samples,np.log(mu), log(context.mu_f))
    old_context.yukawa, old_context.quartic, old_context.gauge = (
        at_mf.yukawa, at_mf.quartic, at_mf.gauge)
    upstream = calculate_master_beta(old_context,c0)
    calculated = result["beta_initial"]
    if context.authority == "certified":
        checks = {"at_MF":{},"at_MS":{}}
        for target,beta_c,point in (("at_MF",c0,at_mf),("at_MS",result["final"],samples[0])):
            for sector_idx,(g,pairs) in enumerate(point.gauge):
                for index,(tf,theta) in enumerate(pairs):
                    value=gauge_ward_residual(beta_c,tf,theta)
                    checks[target][f"sector{sector_idx}_generator{index}"]=value
                    if value > 2e-7:
                        raise ValueError("Certified Wilson tensor fails gauge Ward closure")
        result["gauge_ward_residuals"] = checks
    result["upstream_reference"] = {
        name: {"max_abs_full_tensor_difference":float(np.max(np.abs(calculated[name]-upstream[name]))),
               "direct_LLHH_projection_upstream":_encode(c5_projection(upstream[name],basis)),
               "direct_LLHH_projection_selected_variant":_encode(c5_projection(calculated[name],basis))}
        for name in TERMS}
    if decompose:
        source_runs = []
        for name,source in context.operators:
            channel = evolve(basis,source,context.mu_f,context.mu_s,mu,samples,
                             **{**kwargs,"points":2})
            source_runs.append((name,channel))
        assembled = sum((run["final"] for _,run in source_runs),np.zeros_like(c0))
        residual = float(np.linalg.norm(assembled-result["final"])/
                         max(np.linalg.norm(result["final"]),1e-300))
        if residual > 10*max(kwargs.get("rtol",1e-8),1e-11):
            raise RuntimeError("Independent source-channel ODEs fail real-linearity check")
        result["source_runs"] = source_runs
        result["source_recomposition_residual"] = residual
    return result


def save_results(result: dict, output_npz: str, summary_json: str) -> None:
    """Persist complete independent coordinates, all seven endpoint beta tensors."""
    import json
    from pathlib import Path
    p, s = Path(output_npz), Path(summary_json)
    if p.exists() or s.exists():
        raise FileExistsError("Refusing to overwrite Batch34 result")
    p.parent.mkdir(parents=True, exist_ok=True)
    s.parent.mkdir(parents=True, exist_ok=True)
    source_arrays = {f"source_{k:03d}_full_C_at_MS_gev_inv":run["final"]
                     for k,(_,run) in enumerate(result.get("source_runs",[]))}
    np.savez_compressed(p, coordinates_1based=result["coordinates"],
        mu_gev=result["scales_gev"],
        complete_symmetric_C_trajectory_gev_inv=result["trajectory"],
        full_C_at_MF_gev_inv=result["initial"], full_C_at_MS_gev_inv=result["final"],
        generated_zero_boundary_indices=result["generated_indices"],
        generated_zero_boundary_coordinates_1based=result["coordinates"][result["generated_indices"]],
        fixed_order_full_C_at_MS=result["fixed_order_full"],
        **{f"beta_{key}_at_MF_16pi2":result["beta_initial"][key] for key in TERMS},
        **{f"beta_{key}_at_MS_16pi2":result["beta_final"][key] for key in TERMS},
        **source_arrays)
    basis = result["basis"]
    order = result["generated_sorted"]
    zz = result["trajectory"][-1]
    report = {
        "format": RESULT_FORMAT,
        "beta_mode": result["beta_mode"],
        "beta_physics_certified": False,
        "anti_linear_term_names": result["anti_linear_term_names"],
        "source_equation_conflict_unresolved": True,
        "status": "SyntheticNumericalClosure_NotPhysicalC5" if result.get("authority") == "synthetic"
                  else "ConditionalVerifiedInputNumericalTransport_RequiresScalarMatching",
        "model": "A_-4", "authority": result.get("authority", "synthetic"),
        "source_ids": result.get("source_ids", []),
        "scheme": "MSbar", "RGE": "16pi2*dC/dlnmu=sum(seven beta terms)",
        "normalization": "L5=C_ijab psi_i psi_j phi_a phi_b + hc; V4=lambda_abcd phi_a phi_b phi_c phi_d/4!",
        "basis": {"n_fermions":basis.nf, "n_real_scalars":basis.ns,
            "lepton_components_1based":(np.asarray(basis.leptons)+1).tolist(),
            "higgs_real_components_1based":(np.asarray(basis.higgs)+1).tolist(),
            "independent_symmetric_coordinates":basis.independent_size,
            "fermion_symmetry":"i<->j", "scalar_symmetry":"a<->b",
            "coordinate_order":"np.triu_indices fermions, then np.triu_indices scalars",
            "fermion_field_blocks":result.get("fermion_field_blocks",[]),
            "scalar_field_blocks":result.get("scalar_field_blocks",[])},
        "mu_F_gev":result["mu_F_gev"], "mu_S_gev":result["mu_S_gev"],
        "input_certificates":result.get("input_certificates",{}),
        "coupling_source_verified":result.get("coupling_source_verified",False),
        "upstream_fixed_order_term_comparison":result.get("upstream_reference",{}),
        "upstream_fixed_order_reference_module":"GeneralEFT1C5Mixing.calculate_master_beta",
        "upstream_crossed_term_is_two_exchange": True,
        "source_v1_four_exchange_is_not_gauge_covariant": True,
        "candidate_crossed_Cstar_is_not_loop_certified": True,
        "source_decomposition": [{"source_id":name,
            "initial_norm":float(np.linalg.norm(run["initial"])),
            "final_norm":float(np.linalg.norm(run["final"])),
            "induced_LLHH_at_MS":_encode(c5_projection(run["final"],basis)),
            "seven_beta_projections_at_MS":{
                term:_encode(c5_projection(run["beta_final"][term],basis)) for term in TERMS}}
            for name,run in result.get("source_runs",[])],
        "source_recomposition_residual":result.get("source_recomposition_residual"),
        "input_couplings_interpolated_in":"ln(mu)",
        "sample_scales_gev":result["scales_gev"].tolist(),
        "all_seven_terms":list(TERMS),
        "per_sample_beta_projected_16pi2_and_induced_LLHH":result["projections"],
        "per_sample_full_beta_frobenius_norms":result["norms"],
        "per_sample_total_beta_symmetry_residual":result["symmetry_residuals"],
        "gauge_ward_residuals_when_certified":result.get("gauge_ward_residuals",{}),
        "at_MS_induced_LLHH_not_matched_SMEFT_C5":_encode(c5_projection(result["final"],basis)),
        "at_MF_induced_LLHH":_encode(c5_projection(result["initial"],basis)),
        "initial_nonzero_independent_coordinates":result["initial_support_count"],
        "generated_initially_zero_coordinate_count":len(result["generated_indices"]),
        "generated_threshold_gev_inv":result["generated_threshold"],
        "largest_generated": [{"indices_1based":result["coordinates"][i].tolist(),
                                "coefficient_gev_inv": {"real":float(zz[i].real), "imag":float(zz[i].imag)}}
                               for i in order[:30]],
        "anti_linear_per_term_beta_norm_at_MF":{
            key:float(np.linalg.norm(result["beta_initial"][key]))
            for key in result["anti_linear_term_names"]},
        "anti_linear_per_term_beta_norm_at_MS":{
            key:float(np.linalg.norm(result["beta_final"][key]))
            for key in result["anti_linear_term_names"]},
        "conjugate_term_beta_norm_at_MF":float(np.linalg.norm(result["beta_initial"][TERMS[4]])),
        "conjugate_term_beta_norm_at_MS":float(np.linalg.norm(result["beta_final"][TERMS[4]])),
        "fixed_order_full_relative_difference":result["fixed_order_difference_relative"],
        "solver":result["solver"],
        "raw_arrays_npz":p.name,
        "physical_C5_certified":False,
        "hard_scalar_threshold_matching_performed":False,
        "note":"The LLHH projection in EFT1 is not the complete matched SMEFT Wilson coefficient."
    }
    s.write_text(json.dumps(report, indent=2, allow_nan=False)+"\n", encoding="utf-8")
