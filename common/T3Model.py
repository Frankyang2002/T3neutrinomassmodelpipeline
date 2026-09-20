# This file gives the python definition for T3 representations 
# and is used across different files, we standardise the convention here


from __future__ import annotations

FormalT3Dimensions = tuple[int, int, int] # For separate scalar fields
SharedScalarDimensions = tuple[int, int] # For same scalar fields

T3_CLASSES: dict[str, FormalT3Dimensions] = {
    "A": (1, 3, 2),
    "B": (2, 2, 1),
    "C": (2, 2, 3),
    "D": (3, 1, 2),
    "E": (3, 3, 2),
}

SMOKE = [("B", -1), ("C", -1), ("A", 0), ("D", -2), ("E", 0)]


def encode_alpha(alpha: int) -> str:
    """Convert negative integers to ones that work in wolfram
    Eg: -5 = m5, 5 = p5."""
    return f"m{abs(alpha)}" if alpha < 0 else f"p{alpha}"


def _supported_su2_dimension(dimension: int) -> bool:
    """Return whether one SU(2) irrep dimension is in current pipeline scope."""
    return 1 <= dimension <= 3


def _scalar_can_couple_to_lepton_and_fermion(
    scalar_dimension: int,
    fermion_dimension: int,
) -> bool:
    """Check the SU(2) dimension condition for an L-F-S Yukawa interaction.
    dS = dF pm 1 condition"""
    return abs(scalar_dimension - fermion_dimension) == 1


def _scalar_pair_contains_triplet(
    d_s1: int,
    d_s2: int,
) -> bool:
    """Check whether S1 x S2 contains the integer-spin J=1 channel to make triplets"""
    j1 = (d_s1 - 1) / 2
    j2 = (d_s2 - 1) / 2

    return (
        abs(j1 - j2) <= 1 <= j1 + j2
        and float(j1 + j2).is_integer()
    )


def valid_t3_dimensions(d_s1: int, d_s2: int, d_f: int) -> bool:
    """Validate formal T3 topology dimensions (S1, S2, F) 
    and requires conditions to form weinberg operator
    """
    dimensions = (d_s1, d_s2, d_f)

    if not all(_supported_su2_dimension(d) for d in dimensions):
        return False

    if not _scalar_can_couple_to_lepton_and_fermion(d_s1, d_f):
        return False

    if not _scalar_can_couple_to_lepton_and_fermion(d_s2, d_f):
        return False

    return _scalar_pair_contains_triplet(d_s1, d_s2)


def valid_shared_scalar_dimensions(d_s: int, d_f: int) -> bool:
    """Validate physical dimensions for the supported shared-scalar branch.
    """
    return d_s == 2 and d_f in (1, 3)


def shared_scalar_formal_dimensions(
    d_s: int,
    d_f: int,
) -> FormalT3Dimensions:
    """Map physical shared-scalar dimensions (S, F) to formal (S1, S2, F).
    We do not treat the two of them separately but we represent them as separate for calculations
    """
    if not valid_shared_scalar_dimensions(d_s, d_f):
        raise ValueError(
            "Shared-scalar mode currently supports dS=2 with dF=1 or 3."
        )

    return d_s, d_s, d_f


def physical_dimensions_from_formal(
    d_s1: int,
    d_s2: int,
    d_f: int,
    *,
    shared_scalar: bool = False,
) -> FormalT3Dimensions | SharedScalarDimensions:
    """Return the physical dimension tuple represented by formal T3 data.

    Ordinary T3 has two physical scalars and returns (dS1, dS2, dF).
    Shared-scalar mode has one physical scalar and returns (dS, dF)
    
    We distinguish between the two
    """
    if not shared_scalar:
        return d_s1, d_s2, d_f

    if d_s1 != d_s2:
        raise ValueError(
            "Shared-scalar formal dimensions require dS1 == dS2."
        )

    if not valid_shared_scalar_dimensions(d_s1, d_f):
        raise ValueError(
            "Shared-scalar mode currently supports dS=2 with dF=1 or 3."
        )

    return d_s1, d_f


def identify_t3_class(d_s1: int, d_s2: int, d_f: int) -> str | None:
    """Return the A-E label for matching formal T3 dimensions, if any."""
    dimensions = (d_s1, d_s2, d_f)

    for model_class, known_dimensions in T3_CLASSES.items():
        if dimensions == known_dimensions:
            return model_class

    return None
