"""Renormalisable running for the scalar-only EFT after the T3 fermion is removed."""

from __future__ import annotations

from RGE.running.intermediate.LegacyEFT1Compatibility import (
    legacy_scalar_only_renormalisable_rge,
)


def run_scalar_only_renormalisable_rge(
    d_s1: int,
    d_s2: int,
    d_f: int,
    alpha: int,
    *,
    shared_scalar: bool = False,
):
    """Run the verified RGBeta calculation for the active T3 scalar sector.

    This descriptive interface is independent of ordinal EFT numbering.  The
    historical implementation is reached only through the dedicated legacy
    compatibility boundary so existing JSON metadata and reports remain
    unchanged during the refactor.
    """
    return legacy_scalar_only_renormalisable_rge(
        d_s1,
        d_s2,
        d_f,
        alpha,
        shared_scalar=shared_scalar,
    )
