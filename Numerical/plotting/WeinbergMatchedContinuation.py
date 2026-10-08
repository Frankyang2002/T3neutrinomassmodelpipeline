"""Checked display-only continuation of T3 Weinberg matching diagnostics.

Between MF and MS the displayed 'hard + direct' curve is a threshold-anchored
reconstruction, NOT an intermediate-EFT Wilson coefficient.  The hard term is
held fixed at its MS value while the direct LLSS -> Weinberg term varies.

For newer diagnostics, prefer full complex delta_c5_real/delta_c5_imag arrays.
Older four-scenario comparison diagnostics store only |delta_c5(mu)|.  A
complex continuation cannot in general be recovered from magnitudes alone.
For the project's *real* benchmark couplings only, recover the unique signed
log-affine component curve and verify every stored magnitude.  Reject complex
or ambiguous cases rather than making up phases or sums of magnitudes.

No RGEs, loop matching, or physical data are recomputed here.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np


@dataclass(frozen=True)
class MatchedContinuation:
    mu_gev: np.ndarray
    hard: np.ndarray                 # complex (3, 3) at MS
    direct: np.ndarray               # complex (n, 3, 3)
    combined: np.ndarray             # complex (n, 3, 3)
    method: str
    maximum_relative_magnitude_residual: float


def _matrix(block: dict[str, Any], name: str) -> np.ndarray:
    real = np.asarray(block["real"], dtype=float)
    imag = np.asarray(block["imag"], dtype=float)
    if real.shape != (3, 3) or imag.shape != real.shape:
        raise ValueError(f"{name} must contain real and imag 3x3 matrices")
    result = real + 1j * imag
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{name} contains nonfinite entries")
    return result


def _max_relative_error(actual: np.ndarray, expected: np.ndarray) -> float:
    scale = float(max(np.max(np.abs(expected)), np.max(np.abs(actual)), 1e-40))
    return float(np.max(np.abs(actual - expected)) / scale)


def reconstruct_matched_continuation(payload: dict[str, Any]) -> MatchedContinuation:
    """Build a threshold-matched continuation and validate the full curves.

    Raises ValueError when existing saved diagnostics do not justify an exact
    continuation (missing complex data, non-real legacy data, failed matching,
    or non-log-affine legacy magnitudes).
    """
    block = payload["intermediate_direct_weinberg"]
    threshold = payload["c5_threshold_contributions"]
    scales = payload["scales_gev"]
    mu_f = float(scales["mu_fermion_threshold"])
    mu_s = float(scales["mu_scalar_threshold"])
    if not (mu_f > mu_s > 0):
        raise ValueError("Expected a descending fermion-to-scalar interval")
    mu = np.asarray(block["mu_gev"], dtype=float)
    observed = np.asarray(block["delta_c5_abs"], dtype=float)
    if mu.ndim != 1 or mu.size < 3 or observed.shape != (mu.size, 3, 3):
        raise ValueError("Intermediate direct C5 magnitude trajectory is incomplete")
    if (not np.all(np.isfinite(mu)) or not np.all(np.isfinite(observed))
            or np.any(mu <= 0) or np.any(observed < 0)):
        raise ValueError("Intermediate C5 scales or magnitudes are invalid")
    if np.any(np.diff(mu) >= 0) or abs(mu[0] / mu_f - 1) > 1e-7 or abs(mu[-1] / mu_s - 1) > 1e-7:
        raise ValueError("Intermediate trajectory must cover MF down to MS")

    hard = _matrix(threshold["hard"], "hard")
    direct_s = _matrix(threshold["direct_running"], "direct_running")
    combined_s = _matrix(threshold["combined"], "combined")
    if _max_relative_error(hard + direct_s, combined_s) > 1e-8:
        raise ValueError("Hard + direct does not reproduce the stored complex threshold C5")
    if _max_relative_error(observed[-1], np.abs(direct_s)) > 1e-7:
        raise ValueError("Intermediate direct endpoint disagrees with threshold direct C5")

    final = payload["final_running"]
    final_mu = np.asarray(final["mu_gev"], dtype=float)
    final_abs = np.asarray(final["c5_abs"], dtype=float)
    if final_mu.size < 1 or final_abs.shape != (final_mu.size, 3, 3):
        raise ValueError("Invalid final C5 diagnostics")
    if abs(final_mu[0] / mu_s - 1) > 1e-7:
        raise ValueError("Final C5 trajectory does not start at the scalar threshold")
    if _max_relative_error(final_abs[0], np.abs(combined_s)) > 1e-7:
        raise ValueError("Final C5 magnitude is inconsistent with matching at MS")

    if "delta_c5_real" in block and "delta_c5_imag" in block:
        real = np.asarray(block["delta_c5_real"], dtype=float)
        imag = np.asarray(block["delta_c5_imag"], dtype=float)
        if real.shape != observed.shape or imag.shape != observed.shape:
            raise ValueError("Stored complex intermediate C5 has incorrect dimensions")
        direct = real + 1j * imag
        method = "saved_complex_intermediate_trajectory"
    else:
        # A real one-loop amplitude is log-affine in each component.  The
        # magnitudes determine its real sign at MF if the complete sampled
        # trajectory selects one unambiguous signed branch.
        physical_scale = max(float(np.max(np.abs(hard))), float(np.max(np.abs(direct_s))), 1e-40)
        if (float(np.max(np.abs(np.imag(hard)))) > 1e-10 * physical_scale
                or float(np.max(np.abs(np.imag(direct_s)))) > 1e-10 * physical_scale):
            raise ValueError("Complex legacy direct C5 cannot be reconstructed from magnitudes")
        log_fraction = np.log(mu / mu_s) / np.log(mu_f / mu_s)
        direct = np.zeros_like(observed, dtype=complex)
        for i in range(3):
            for j in range(i, 3):
                ds = float(np.real(direct_s[i, j]))
                at_f = float(observed[0, i, j])
                sampled = observed[:, i, j]
                candidates = [ds + (sign * at_f - ds) * log_fraction for sign in (+1, -1)]
                errors = [
                    _max_relative_error(np.abs(values), sampled)
                    for values in candidates
                ]
                best = int(np.argmin(errors))
                if errors[best] > 2e-6:
                    raise ValueError(
                        f"Real legacy C5[{i+1},{j+1}] fails fixed-order log-affine validation: {errors[best]:.3g}"
                    )
                # Reject an unresolved sign unless the two candidate curves
                # are numerically equivalent over this component's full scale.
                other = 1 - best
                if errors[other] <= 2e-6 and _max_relative_error(candidates[0], candidates[1]) > 2e-6:
                    raise ValueError(f"Ambiguous sign for real legacy C5[{i+1},{j+1}]")
                direct[:, i, j] = candidates[best]
                direct[:, j, i] = candidates[best]
        method = "validated_real_log_affine_reconstruction"

    if not np.all(np.isfinite(direct)):
        raise ValueError("Intermediate complex reconstruction contains nonfinite values")
    if _max_relative_error(np.abs(direct), observed) > 2e-6:
        raise ValueError("Reconstructed intermediate curve does not match all saved magnitudes")
    if _max_relative_error(direct[-1], direct_s) > 1e-7:
        raise ValueError("Reconstructed intermediate complex C5 does not match MS endpoint")
    if _max_relative_error(direct, np.swapaxes(direct, 1, 2)) > 1e-7:
        raise ValueError("Reconstructed intermediate C5 is not symmetric")

    combined = hard[None, :, :] + direct
    if _max_relative_error(combined[-1], combined_s) > 1e-7:
        raise ValueError("Matched total does not join the final complex threshold C5")

    return MatchedContinuation(
        mu_gev=mu,
        hard=hard,
        direct=direct,
        combined=combined,
        method=method,
        maximum_relative_magnitude_residual=_max_relative_error(np.abs(direct), observed),
    )
