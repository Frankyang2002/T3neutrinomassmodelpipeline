"""Fixed-order diagnostic trajectory for scalar-only LLSS Wilson coefficients.

This module samples the *existing* symbolic component transport implemented by
``RGE.running.intermediate.ScalarOnlyWilsonFlow``.  It does not introduce a new
RGE, does not RG-improve the Wilson coefficient, and does not modify the
authoritative final Weinberg coefficient.

For every requested scale mu between MF and MS it evaluates the symbolic
fixed-order expression

    C_LLSS(mu) = C_LLSS^(0)
               + hbar * log(mu/MF) * beta^(1)[C_LLSS^(0)]
               + O(hbar^2),

with hbar = 1/(16*pi^2).

The returned trajectory is symbolic at component level.  This is intentional:
the component Wilson pipeline is a one-generation/group-theory object, whereas
the production numerical benchmark uses full 3x3 lepton-heavy Yukawa matrices.
No arbitrary scalar reduction of those matrices is made here.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Sequence

import numpy as np

from RGE.running.intermediate.ScalarOnlyWilsonFlow import (
    run_component_wilson_transport,
)


def sample_fixed_order_llss_transport(
    *,
    wilson_seed_path: Path,
    wilson_rge_path: Path,
    rgbeta_path: Path,
    mu_high_gev: float,
    scales_gev: Sequence[float],
) -> dict[str, Any]:
    """Return symbolic fixed-order LLSS component transport at saved scales."""

    scales = np.asarray(scales_gev, dtype=float)
    mu_high = float(mu_high_gev)

    if scales.ndim != 1 or scales.size < 2:
        raise ValueError(
            "scales_gev must be one-dimensional with at least two entries."
        )
    if not np.isfinite(mu_high) or mu_high <= 0.0:
        raise ValueError("mu_high_gev must be finite and positive.")
    if not np.all(np.isfinite(scales)) or np.any(scales <= 0.0):
        raise ValueError("All scales_gev entries must be finite and positive.")

    # The production hierarchy runs downward from the fermion threshold.
    tolerance = 1.0e-12 * max(1.0, abs(mu_high))
    if np.any(scales > mu_high + tolerance):
        raise ValueError(
            "LLSS diagnostic scales cannot lie above the upper threshold."
        )

    points: list[dict[str, Any]] = []
    reference_counts: tuple[int, int] | None = None

    for mu in scales:
        transported = run_component_wilson_transport(
            wilson_seed_path=Path(wilson_seed_path),
            wilson_rge_path=Path(wilson_rge_path),
            rgbeta_path=Path(rgbeta_path),
            mu_high=mu_high,
            mu_low=float(mu),
            output_path=None,
        )

        counts = (
            int(transported["tree_boundary_component_count"]),
            int(transported["beta_component_count"]),
        )
        if reference_counts is None:
            reference_counts = counts
        elif counts != reference_counts:
            raise RuntimeError(
                "LLSS component basis changed while sampling one trajectory."
            )

        points.append(
            {
                "mu_gev": float(mu),
                "log_mu_over_mu_high": float(np.log(float(mu) / mu_high)),
                "tree_boundary_components": transported[
                    "tree_boundary_components"
                ],
                "one_loop_running_components": transported[
                    "one_loop_running_components"
                ],
                "full_fixed_order_components": transported[
                    "full_fixed_order_components"
                ],
                "running_correction_component_count": int(
                    transported["running_correction_component_count"]
                ),
                "generated_by_running_component_count": int(
                    transported["generated_by_running_component_count"]
                ),
            }
        )

    first = points[0]
    equal_scale_sampled = bool(
        np.isclose(first["mu_gev"], mu_high, rtol=1.0e-12, atol=0.0)
    )
    equal_scale_running_vanishes = None
    if equal_scale_sampled:
        equal_scale_running_vanishes = (
            len(first["one_loop_running_components"]) == 0
            or all(
                str(value).strip() in {"0", "0.0"}
                for value in first["one_loop_running_components"].values()
            )
        )

    return {
        "status": "Success",
        "diagnostic_only": True,
        "quantity": "scalar-only LLSS component Wilson coefficient",
        "representation": "symbolic one-generation/group-component",
        "order": "fixed overall one loop",
        "rg_improved": False,
        "used_in_authoritative_final_c5": False,
        "convention": {
            "hbar": "1/(16*pi^2)",
            "rge": "16*pi^2*dC/dln(mu)=beta^(1)",
            "transport": (
                "C_LLSS(mu)=C_LLSS^(0)"
                "+hbar*log(mu/MF)*beta^(1)[C_LLSS^(0)]+O(hbar^2)"
            ),
        },
        "power_counting": {
            "llss_self_running": "O(hbar)",
            "llss_self_running_inserted_into_scalar_loop": "O(hbar^2)",
            "authoritative_one_loop_final_c5": (
                "does not include the O(hbar^2) insertion"
            ),
        },
        "flavor_note": (
            "No numerical norm is reported: this component transport is a "
            "one-generation/group-theory object, while production Yukawas are "
            "full 3x3 flavor matrices. No arbitrary matrix-to-scalar reduction "
            "is introduced."
        ),
        "mu_high_gev": mu_high,
        "mu_gev": scales.tolist(),
        "tree_boundary_component_count": (
            reference_counts[0] if reference_counts is not None else 0
        ),
        "beta_component_count": (
            reference_counts[1] if reference_counts is not None else 0
        ),
        "equal_scale_sampled": equal_scale_sampled,
        "equal_scale_running_vanishes": equal_scale_running_vanishes,
        "points": points,
    }
