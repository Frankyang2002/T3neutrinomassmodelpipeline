"""T3 Batch 7: source-preserving pole and operator-completeness diagnostics.

This is an audit of saved Batches 1/2/4/5/6, not an independently
calculated Matchete UV amplitude. The source-level C5 projection need NOT
be symmetric before quartic-channel reassembly; only the *physical total*
C5 is asserted to be symmetric. No antisymmetric parts are discarded from
saved source arrays. No counterterm is fitted or injected into C5.

Matching order: UV -> integrate out ALL F -> EFT1 (BOTH S) -> integrate
out BOTH scalars -> SMEFT. No one-scalar intermediate threshold.
"""
from __future__ import annotations

from collections.abc import Mapping, Sequence
from math import isfinite
from typing import Any

import numpy as np

SCHEMA = "t3-native-source-pole-closure-batch7-v1"
STATUS = "NUMERICAL_SOURCE_RESOLVED_COUNTERTERM_AND_POLE_AUDIT_CONDITIONAL"
MODELS = (
    'A_-4', 'A_-2', 'A_+0', 'B_-3', 'B_-1', 'B_+1',
    'C_-3', 'C_-1', 'C_+1', 'D_-2', 'D_+0', 'D_+2',
    'E_-4', 'E_-2', 'E_+0', 'E_+2',
)
COMPLETE_MATCHING_GATES = (
    "native_expression_proven_bare_UV_1PI",
    "native_UV_counterterms_explicit",
    "complete_hard_topology_inventory",
    "signed_GammaCC_Weyl_normalisation_verified",
    "UV_minus_EFT1_soft_matching_verified",
    "all_scalar_threshold_matching_complete",
)


def _matrix(x: Any, label: str) -> np.ndarray:
    a = np.asarray(x, dtype=np.complex128)
    if a.shape != (3, 3) or not np.isfinite(a).all():
        raise ValueError(f"{label}: expected finite complex 3x3")
    return a


def relative(a: Any, b: Any) -> float:
    x, y = np.asarray(a), np.asarray(b)
    if x.shape != y.shape or not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError("Cannot compute relative norm for different/nonfinite tensors")
    return float(np.linalg.norm(x - y) / max(float(np.linalg.norm(y)), 1e-100))


def symmetry_split(x: Any, label: str = "matrix") -> dict[str, Any]:
    """Keep complex transpose, NEVER Hermitian conjugation for Majorana C5."""
    a = _matrix(x, label)
    s = (a + a.T) / 2
    antisym = (a - a.T) / 2
    return {
        "raw": a, "symmetric": s, "antisymmetric": antisym,
        "raw_norm": float(np.linalg.norm(a)),
        "symmetric_norm": float(np.linalg.norm(s)),
        "antisymmetric_norm": float(np.linalg.norm(antisym)),
        "antisymmetric_fraction": float(np.linalg.norm(antisym) /
                                        max(float(np.linalg.norm(a)), 1e-100)),
        "transpose_residual": relative(a, a.T),
    }


def physical_total(x: Any, label: str, tolerance: float = 1e-9) -> np.ndarray:
    """Only total Wilson coefficients or certified physical tensors use this gate."""
    if not isfinite(tolerance) or tolerance <= 0:
        raise ValueError("Invalid physical symmetry tolerance")
    part = symmetry_split(x, label)
    if part["transpose_residual"] > tolerance:
        raise ValueError(f"{label}: physical C5 violates flavour symmetry, "
                         f"relative={part['transpose_residual']:.6g}")
    return part["symmetric"]


def matrix_record(x: Any, label: str, *, physical: bool = False) -> dict[str, Any]:
    """JSON-ready exact raw matrix, and a NON-DESTRUCTIVE symmetry diagnostic.

    For a source summand with nonzero antisymmetry there is deliberately no
    mass interpretation and no symmetrisation. The raw matrix is unchanged.
    """
    p = symmetry_split(x, label)
    if physical:
        physical_total(x, label)
    a = p["raw"]
    return {
        "real": a.real.tolist(), "imag": a.imag.tolist(),
        "frobenius_GeV_inverse": p["raw_norm"],
        "symmetric_norm": p["symmetric_norm"],
        "antisymmetric_norm": p["antisymmetric_norm"],
        "transpose_residual": p["transpose_residual"],
        "antisymmetric_fraction": p["antisymmetric_fraction"],
        "physical_total_validation_applied": bool(physical),
    }


def _source_stack(arr: Any, names: Sequence[str], label: str) -> np.ndarray:
    x = np.asarray(arr, dtype=complex)
    if x.ndim != 3 or x.shape[1:] != (3, 3) or not len(names) == len(x):
        raise ValueError(f"{label}: require one finite complex 3x3 per source")
    if not np.isfinite(x).all() or not names or len(set(names)) != len(names):
        raise ValueError(f"{label}: invalid finite source tensor/name inventory")
    if any(not isinstance(n, str) or not n for n in names):
        raise ValueError(f"{label}: source names must be nonempty strings")
    return x


def source_cancellation_ledger(source_matrices: Any, names: Sequence[str],
                               expected_total: Any,
                               *, label: str, tolerance: float = 1e-9) -> dict[str, Any]:
    """Check complete source summation, preserving antisymmetric cancellation.

    A per-source antisymmetric contribution is NOT an independently physical
    Wilson coefficient. If their sum is symmetric they can cancel; if their
    sum is not symmetric the aggregate remains blocked.
    """
    arrays = _source_stack(source_matrices, names, label)
    expected = physical_total(expected_total, f"{label} expected total", tolerance)
    summed_raw = arrays.sum(axis=0)
    rel = relative(summed_raw, expected)
    if rel > tolerance:
        raise ValueError(f"{label}: source sum fails total reassembly {rel:.6g}")
    combined = physical_total(summed_raw, f"{label} source sum", tolerance)
    anti_sources = np.asarray([symmetry_split(m)["antisymmetric"] for m in arrays])
    total_source_antisym_norm = float(np.sum(np.linalg.norm(anti_sources, axis=(1, 2))))
    residual_anti_norm = float(np.linalg.norm(anti_sources.sum(axis=0)))
    return {
        "source_count": len(names),
        "source_names": list(names),
        "per_source": {name: matrix_record(value, name) for name, value in zip(names, arrays)},
        "total": matrix_record(combined, label + " total", physical=True),
        "reassembly_residual": rel,
        "sum_abs_antisymmetric_norm": total_source_antisym_norm,
        "net_antisymmetric_norm": residual_anti_norm,
        "antisymmetric_cancellation_fraction": (
            float(1 - residual_anti_norm / total_source_antisym_norm)
            if total_source_antisym_norm > 0 else 1.0),
        "source_matrices_preserved_unsymmetrised": True,
    }


def source_by_native_label(kernels: Any, labels: Sequence[str]) -> dict[str, np.ndarray]:
    """Retain *raw* native quartic grouping before any flavour symmetry gate.

    Native input geometry is [Matchete kernel, 3 heavy F, 3, 3]. A native
    source may not be separately symmetric, while the sum must be.
    """
    z = np.asarray(kernels, complex)
    if z.ndim != 4 or z.shape[1:] != (3, 3, 3) or len(labels) != len(z):
        raise ValueError("Native Matchete source tensor/label mismatch")
    if not np.isfinite(z).all() or not labels:
        raise ValueError("Nonfinite/empty native source tensor")
    grouped: dict[str, np.ndarray] = {}
    for label, piece in zip(labels, z):
        if not isinstance(label, str) or not label.startswith("lambda"):
            raise ValueError("Native quartic label absent")
        if label not in grouped:
            grouped[label] = np.zeros((3, 3), complex)
        grouped[label] += piece.sum(axis=0)
    total = z.sum(axis=(0, 1))
    if relative(sum(grouped.values(), np.zeros((3, 3), complex)), total) > 1e-10:
        raise ValueError("Native labelled source aggregation failed")
    return grouped


def linear_coefficient_diagnostic(reference: Any, candidate: Any) -> dict[str, Any]:
    """Complex least-squares scalar mapping, DIAGNOSTIC ONLY.

    Even a tiny residual is not a derivation of a spinor/Matchete convention;
    fitted coefficients MUST NOT feed the matching or neutrino mass.
    """
    x = _matrix(reference, "reference")
    y = _matrix(candidate, "candidate")
    norm_sq = float(np.vdot(x, x).real)
    if norm_sq <= 1e-200:
        return {"fit_defined": False, "coefficient_real": None,
                "coefficient_imag": None, "postfit_residual_relative": None,
                "fit_used_in_physics": False}
    z = np.vdot(x, y) / norm_sq
    fitted = z * x
    return {"fit_defined": True, "coefficient_real": float(z.real),
            "coefficient_imag": float(z.imag),
            "postfit_residual_relative": relative(y, fitted),
            "candidate_minus_reference_relative": relative(y, x),
            "fit_used_in_physics": False}


def native_pole_evidence(*, native: Any, eft1_ct: Any,
                         candidate_adm: Any,
                         candidate_scalar_adm: Any) -> dict[str, Any]:
    """Quantify, never artificially cancel, all known pole mismatches.

    Matchete MatchedEFT PL expressions have NOT been shown to represent the
    bare UV 1PI amplitude. Thus differences here are comparison diagnostics,
    NOT a sum of renormalized UV + EFT1 counterterms.
    """
    n = physical_total(native, "total native Laurent pole")
    c = physical_total(eft1_ct, "EFT1 scalar bubble CT")
    adm = physical_total(candidate_adm, "seven-channel EFT1 ADM CT")
    scalar = physical_total(candidate_scalar_adm, "EFT1 scalar-pair ADM CT")
    return {
        "native_pole": matrix_record(n, "native pole", physical=True),
        "EFT1_scalar_bubble_CT": matrix_record(c, "EFT1 scalar CT", physical=True),
        "native_minus_EFT1_scalar_CT": matrix_record(n - c, "native - EFT1 scalar CT", physical=True),
        "candidate_EFT1_seven_channel_CT": matrix_record(adm, "EFT1 candidate ADM", physical=True),
        "native_minus_candidate_ADM": matrix_record(n - adm, "native - candidate ADM", physical=True),
        "native_vs_scalar_CT_relative": relative(n, c),
        "native_vs_full_EFT1_ADM_relative": relative(n, adm),
        "explicit_scalar_CT_vs_candidate_scalar_ADM_relative": relative(c, scalar),
        "formal_normalisation_fit_native_vs_scalar_CT": linear_coefficient_diagnostic(c, n),
        "formal_normalisation_fit_native_vs_full_ADM": linear_coefficient_diagnostic(adm, n),
        "native_is_proven_bare_UV_amplitude": False,
        "native_UV_counterterms_calculated": False,
        "UV_pole_cancellation_verified": False,
        "counterterm_fitted_to_cancel_native": False,
    }


def scalar_integral_checks(*, explicit_ct: Any,
                           batch5_soft_pole: Any,
                           b0_scale_ct: Any,
                           atol_relative: float = 1e-8) -> dict[str, Any]:
    c = physical_total(explicit_ct, "explicit scalar CT")
    neg_soft = physical_total(-_matrix(batch5_soft_pole, "soft bare pole"), "negative soft pole")
    slope = physical_total(b0_scale_ct, "finite B0 scale slope")
    r1 = relative(c, neg_soft)
    r2 = relative(c, slope)
    if max(r1, r2) > atol_relative:
        raise ValueError(f"Explicit scalar CT disagrees with B0 or saved soft pole: {r1}, {r2}")
    return {"CT_vs_negative_soft_pole_relative": r1,
            "CT_vs_B0_scale_slope_relative": r2,
            "EFT1_scalar_bubble_self_pole_cancelled": True,
            "complete_UV_pole_cancelled": False}


def certify_matching(evidence: Mapping[str, bool] | None = None) -> dict[str, Any]:
    """Fail closed: no batch is allowed to imply external UV validation."""
    evidence = dict(evidence or {})
    extra = set(evidence) - set(COMPLETE_MATCHING_GATES)
    if extra:
        raise ValueError(f"Unknown certification gate(s): {sorted(extra)}")
    gates = {name: bool(evidence.get(name, False)) for name in COMPLETE_MATCHING_GATES}
    return {"gates": gates, "complete_matching_claim_allowed": all(gates.values()),
            "missing_gates": [name for name, ok in gates.items() if not ok],
            "no_synthetic_UV_counterterm": True}
