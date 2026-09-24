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
    """Return whether one SU(2) irrep dimension is in current production scope."""
    return 1 <= dimension <= 3


def _positive_su2_dimension(dimension: int) -> bool:
    """Return whether a dimension can label a finite SU(2) irrep."""
    return dimension >= 1


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


def valid_t3_topology_dimensions(d_s1: int, d_s2: int, d_f: int) -> bool:
    """Validate the representation-theory conditions required by the T3 topology.

    Unlike :func:`valid_t3_dimensions`, this does not impose the current
    singlet/doublet/triplet production-support restriction. It is therefore the
    validator used by the explicit ``--force`` path.
    """
    dimensions = (d_s1, d_s2, d_f)

    if not all(_positive_su2_dimension(d) for d in dimensions):
        return False

    if not _scalar_can_couple_to_lepton_and_fermion(d_s1, d_f):
        return False

    if not _scalar_can_couple_to_lepton_and_fermion(d_s2, d_f):
        return False

    return _scalar_pair_contains_triplet(d_s1, d_s2)


def valid_t3_dimensions(d_s1: int, d_s2: int, d_f: int) -> bool:
    """Validate a T3 representation inside the current production support."""
    dimensions = (d_s1, d_s2, d_f)
    return (
        all(_supported_su2_dimension(d) for d in dimensions)
        and valid_t3_topology_dimensions(d_s1, d_s2, d_f)
    )


def valid_shared_scalar_topology_dimensions(d_s: int, d_f: int) -> bool:
    """Validate the topology for one physical scalar used in both T3 scalar roles."""
    return valid_t3_topology_dimensions(d_s, d_s, d_f)


def valid_shared_scalar_dimensions(d_s: int, d_f: int) -> bool:
    """Validate physical dimensions for the supported shared-scalar branch."""
    return d_s == 2 and d_f in (1, 3)


def _su2_irrep_contains_neutral_component(
    dimension: int,
    twice_hypercharge: int,
) -> bool:
    """Return whether an SU(2) irrep contains a state with electric charge Q=0.

    The project uses the convention Q = T3 + Y. For an irrep of dimension
    ``d`` the allowed values of ``2*T3`` are ``-(d-1), -(d-3), ..., d-1``.
    A neutral component therefore exists exactly when ``-2Y`` is one of those
    weights.
    """
    if not _positive_su2_dimension(dimension):
        return False

    highest_twice_t3 = dimension - 1
    return (
        abs(twice_hypercharge) <= highest_twice_t3
        and (highest_twice_t3 - twice_hypercharge) % 2 == 0
    )


def t3_neutral_component_fields(
    d_s1: int,
    d_s2: int,
    d_f: int,
    alpha: int,
) -> tuple[str, ...]:
    """Return the BSM T3 multiplets that contain an electrically neutral state.

    The T3 hypercharges are
      2 Y(S1) = alpha,
      2 Y(S2) = alpha + 2,
      2 Y(F)  = alpha + 1.
    """
    candidates = (
        ("S1", d_s1, alpha),
        ("S2", d_s2, alpha + 2),
        ("F", d_f, alpha + 1),
    )
    return tuple(
        name
        for name, dimension, twice_hypercharge in candidates
        if _su2_irrep_contains_neutral_component(dimension, twice_hypercharge)
    )


def t3_has_neutral_bsm_component(
    d_s1: int,
    d_s2: int,
    d_f: int,
    alpha: int,
) -> bool:
    """Return whether at least one BSM T3 multiplet has a neutral component."""
    return bool(t3_neutral_component_fields(d_s1, d_s2, d_f, alpha))


def neutral_t3_class_points(
    alphas,
    model_classes: tuple[str, ...] = ("A", "B", "C", "D", "E"),
) -> tuple[tuple[str, int], ...]:
    """Return class/alpha points with at least one neutral BSM component."""
    points: list[tuple[str, int]] = []
    for model_class in model_classes:
        if model_class not in T3_CLASSES:
            raise ValueError(f"Unknown T3 model class: {model_class}")
        d_s1, d_s2, d_f = T3_CLASSES[model_class]
        for alpha in alphas:
            if t3_has_neutral_bsm_component(d_s1, d_s2, d_f, alpha):
                points.append((model_class, alpha))
    return tuple(points)


def shared_scalar_formal_dimensions(
    d_s: int,
    d_f: int,
    *,
    force: bool = False,
) -> FormalT3Dimensions:
    """Map physical shared-scalar dimensions (S, F) to formal (S1, S2, F).
    We do not treat the two of them separately but we represent them as separate for calculations
    """
    validator = (
        valid_shared_scalar_topology_dimensions
        if force
        else valid_shared_scalar_dimensions
    )
    if not validator(d_s, d_f):
        if force:
            raise ValueError(
                "Forced shared-scalar dimensions must still form a T3 topology: "
                "positive dimensions with dS=dF±1 and S⊗S containing the triplet."
            )
        raise ValueError(
            "Shared-scalar mode currently supports dS=2 with dF=1 or 3. "
            "Use --force to bypass only this production-support restriction."
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
