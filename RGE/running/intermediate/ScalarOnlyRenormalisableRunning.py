"""Renormalisable running after the T3 fermion has been removed."""

from __future__ import annotations

from RGE.running.rgbeta.RGBetaT3Running import run_rgbeta_t3_scalar_only


def run_scalar_only_renormalisable_rge(
    d_s1: int,
    d_s2: int,
    d_f: int,
    alpha: int,
    *,
    shared_scalar: bool = False,
):
    """Run RGBeta for the active scalar sector between thresholds."""
    return run_rgbeta_t3_scalar_only(
        d_s1,
        d_s2,
        d_f,
        alpha,
        shared_scalar=shared_scalar,
    )
