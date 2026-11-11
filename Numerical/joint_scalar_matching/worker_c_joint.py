"""Worker C: strict joint EFT1 -> SMEFT scalar threshold interface.

New additive module; delegates the verified LLSS*HHSS zero-momentum finite
MSbar bubble to Numerical.joint_scalar_matching.core. No production solver edits.

IMPORTANT: this is the complete ONE-INSERTION LLSS*HHSS bubble, not the full
scalar one-loop matching for a general EFT1 operator inventory. Other finite
operator, counterterm and wavefunction contributions are NOT silently set to
zero as a claim of physical completeness.

Input convention: T3_REAL_LLSS_HHSS_WEINBERG_v1; same as core.py and adapter.py.
EFT1 Wilsons and scalar background MUST BOTH be evaluated at mu_s.
"""
from __future__ import annotations

from collections.abc import Mapping
import math
import numpy as np

from Numerical.joint_scalar_matching.core import match_c5, check_tensors, EigenBasis

CONVENTION = "T3_REAL_LLSS_HHSS_WEINBERG_v1"
SCHEME = "MSbar"


def _matrix(value: object, key: str) -> np.ndarray:
    a = np.asarray(value, dtype=complex)
    if a.shape != (3, 3) or not (np.isfinite(a.real).all() and np.isfinite(a.imag).all()):
        raise ValueError(f"{key} must be finite complex (3,3)")
    if not np.allclose(a, a.T, rtol=1e-10, atol=1e-28):
        raise ValueError(f"{key} must be flavour-symmetric")
    return a


def _require_metadata(data: Mapping, kind: str, mu: float) -> None:
    if data.get("tensor_convention") != CONVENTION:
        raise ValueError(f"{kind} needs explicit tensor_convention={CONVENTION}")
    if data.get("renormalization_scheme") != SCHEME:
        raise ValueError(f"{kind} needs renormalization_scheme=MSbar")
    scale = float(data.get("mu_s_gev", float("nan")))
    if not math.isfinite(scale) or not math.isclose(scale, mu, rel_tol=1e-10):
        raise ValueError(f"{kind} must be evaluated at mu_s_gev={mu}")
    if not isinstance(data.get("provenance"), Mapping) or not data["provenance"].get("source"):
        raise ValueError(f"{kind} requires provenance.source")


def match_joint_scalars(eft1_wilsons: Mapping, scalar_background: Mapping,
                        mu_s: float) -> dict:
    """Match both S1 and S2 simultaneously, returning distinct C5 sectors.

    eft1_wilsons:
      C5_EFT1: (3,3) already evolved, INCLUDING all EFT1-generated C5;
      C_LLSS: (3,3,2,2,n,n) at mu_s;
      mu_s_gev, tensor_convention, renormalization_scheme, provenance.
    scalar_background:
      mass_squared_gev2: (n,n) full unbroken real scalar quadratic form;
      lambda_HHSS: (4,4,n,n) complete HHSS quartic tensor;
      scalar_groups: list length n, at least one S1 and one S2;
      scalar_mass_origin: 'unbroken_electroweak' or 'EWSB_resummed';
      optional lambda_HHSS_channels and quartic_channel_names;
      mu_s_gev, tensor_convention, renormalization_scheme, provenance.

    Results are numerical numpy arrays, not JSON. The output has *no*
    independently computed finite classes beyond the scalar bubble.
    """
    if not isinstance(eft1_wilsons, Mapping) or not isinstance(scalar_background, Mapping):
        raise TypeError("Both interfaces must be mappings")
    mu = float(mu_s)
    if not math.isfinite(mu) or mu <= 0:
        raise ValueError("mu_s must be finite and positive GeV")
    _require_metadata(eft1_wilsons, "eft1_wilsons", mu)
    _require_metadata(scalar_background, "scalar_background", mu)
    c5 = _matrix(eft1_wilsons["C5_EFT1"], "C5_EFT1")
    mass2 = np.asarray(scalar_background["mass_squared_gev2"])
    basis = EigenBasis.from_mass2(mass2)
    n = len(basis.masses2)
    groups = list(scalar_background["scalar_groups"])
    if len(groups) != n or set(groups) != {"S1", "S2"}:
        raise ValueError("scalar_groups must label all real components of both S1 and S2")
    # Preserve canonical input-field ordering as in the project's adapter.
    if groups != ["S1"] * groups.count("S1") + ["S2"] * groups.count("S2"):
        raise ValueError("scalar basis must order all S1 real components before S2")
    if groups.count("S1") % 2 or groups.count("S2") % 2:
        raise ValueError("each complex scalar multiplet has 2*d real components")
    mass_origin = scalar_background.get("scalar_mass_origin")
    if mass_origin not in {"unbroken_electroweak", "EWSB_resummed"}:
        raise ValueError("scalar_mass_origin must be declared")
    C, q = check_tensors(eft1_wilsons["C_LLSS"],
                         scalar_background["lambda_HHSS"], n)
    channels = scalar_background.get("lambda_HHSS_channels")
    channel_breakdown = {}
    if channels is not None:
        channels = np.asarray(channels)
        names = list(scalar_background.get("quartic_channel_names", []))
        if channels.ndim != 5 or channels.shape[1:] != q.shape:
            raise ValueError("lambda_HHSS_channels must be (k,4,4,n,n)")
        if len(names) != len(channels) or len(set(names)) != len(names) or not all(isinstance(x, str) and x for x in names):
            raise ValueError("quartic channels need unique names")
        for ch in channels:
            check_tensors(C, ch, n)
        if not np.allclose(channels.sum(axis=0), q, rtol=1e-11,
                           atol=1e-13 * max(float(np.max(np.abs(q))), 1e-300)):
            raise ValueError("Quartic channel sum differs from lambda_HHSS")
    elif scalar_background.get("quartic_channel_names"):
        raise ValueError("Channel names supplied without channel tensors")

    result = match_c5(C, q, mass2, mu, c5_eft1=c5,
                      group_names=groups, strict_projection=True)
    if channels is not None:
        for name, ch in zip(names, channels):
            part = match_c5(C, ch, mass2, mu, c5_eft1=np.zeros((3, 3), complex),
                            group_names=groups, strict_projection=False)
            channel_breakdown[name] = {"c5_finite_scalar": part["c5_scalar_bubble"],
                                       "projection_residual": part["projection_residual"]}
        summed = sum((d["c5_finite_scalar"] for d in channel_breakdown.values()),
                     np.zeros((3,3), complex))
        if not np.allclose(summed, result["c5_scalar_bubble"], rtol=1e-9, atol=1e-28):
            raise ValueError("Quartic-channel decomposition does not close")

    return {
        "c5_inherited": c5.copy(),
        "c5_finite_scalar": result["c5_scalar_bubble"],
        "c5_total": result["c5_total"],
        "contribution_breakdown": {
            "inherited_eft1_already_run": c5.copy(),
            "finite_LLSS_HHSS_joint_bubble": result["c5_scalar_bubble"],
            "eigenstate_pairs": result["pair_c5"],
            "quartic_channels": channel_breakdown,
            "additional_scalar_one_loop_terms": "NOT_CALCULATED",
        },
        "scalar_eigenvalues_sq_gev2": result["eigenvalues_sq_gev2"],
        "scalar_rotation_U": result["mass_rotation_U"],
        "weinberg_projection_residual": result["projection_residual"],
        "mu_s_gev": mu,
        "assumptions": [
            "Dimension-five truncation, zero external momentum, one C_LLSS*lambda_HHSS insertion",
            "Finite MSbar B0; all real components of S1 and S2 jointly integrated out",
            "C5_EFT1 is already evolved to mu_s and is added once; no new EFT1 RG log inserted",
            "Input C_LLSS, C5 and HHSS tensors have consistent Weyl/CG/operator normalization",
            "Input scalar mass matrix and quartic include all desired running scalar directions",
            "EWSB-resummed masses can contain dimension-seven-and-higher effects" if mass_origin == "EWSB_resummed" else "Unbroken-electroweak scalar masses",
        ],
        "provenance": {
            "wilson": dict(eft1_wilsons["provenance"]),
            "scalar": dict(scalar_background["provenance"]),
            "kernel": "Numerical.joint_scalar_matching.core.match_c5 (public repository main inspected 2026-11-12)",
            "tensor_convention": CONVENTION,
            "scheme": SCHEME,
        },
        "computed_terms": ["LLSS_HHSS_finite_scalar_bubble"],
        "omitted_terms": ["Other EFT1 dimension-five operator insertions if present",
                          "Finite scalar threshold wavefunction and counterterm corrections if present",
                          "Other scalar one-loop topologies not reducible to the C_LLSS*lambda_HHSS bubble"],
        "full_one_loop_scalar_matching_completed": False,
    }
