from __future__ import annotations

T3_CLASSES = {
    "A": (1, 3, 2),
    "B": (2, 2, 1),
    "C": (2, 2, 3),
    "D": (3, 1, 2),
    "E": (3, 3, 2),
}

INTERESTING = [("A", 0), ("B", -1), ("C", -1), ("D", -2), ("E", 0)]
SMOKE = [("B", -1), ("C", -1), ("A", 0), ("D", -2), ("E", 0)]
EXTENDED = [
    ("A", 0), ("A", -2), ("B", -1), ("C", -1),
    ("D", -2), ("E", 0), ("E", -2),
]


def encode_alpha(alpha: int) -> str:
    return f"m{abs(alpha)}" if alpha < 0 else f"p{alpha}"


def valid_t3_dimensions(d_s1: int, d_s2: int, d_f: int) -> bool:
    if any(d > 3 for d in (d_s1, d_s2, d_f)):
        return False
    if min(d_s1, d_s2, d_f) < 1:
        return False
    if abs(d_s1 - d_f) != 1 or abs(d_s2 - d_f) != 1:
        return False
    j1 = (d_s1 - 1) / 2
    j2 = (d_s2 - 1) / 2
    return abs(j1 - j2) <= 1 <= j1 + j2 and float(j1 + j2).is_integer()


def valid_shared_scalar_dimensions(d_s: int, d_f: int) -> bool:
    """Supported one-physical-scalar branch.

    Current production scope is the inert scalar doublet with a neutral real
    singlet/triplet fermion: dS=2, dF=1 or 3.  This covers the usual type-I
    and type-III scotogenic branches and avoids the known neutral pseudoreal
    dF=2 ambiguity in the current RGBeta setup.
    """
    return d_s == 2 and d_f in (1, 3)


def shared_scalar_formal_dimensions(d_s: int, d_f: int) -> tuple[int, int, int]:
    if not valid_shared_scalar_dimensions(d_s, d_f):
        raise ValueError(
            "Shared-scalar mode currently supports dS=2 with dF=1 or 3."
        )
    return d_s, d_s, d_f


def identify_t3_class(d_s1: int, d_s2: int, d_f: int) -> str | None:
    dimensions = (d_s1, d_s2, d_f)
    for model_class, known_dimensions in T3_CLASSES.items():
        if dimensions == known_dimensions:
            return model_class
    return None
