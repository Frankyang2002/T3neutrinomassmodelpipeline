from __future__ import annotations


# Known T3 models from the original classification with specific dimensions
T3_CLASSES = {
    "A": (1, 3, 2),
    "B": (2, 2, 1),
    "C": (2, 2, 3),
    "D": (3, 1, 2),
    "E": (3, 3, 2),
}


# Some interesting models to look into.
# The first entry is the known T3 class and the second is alpha
# Hypercharges are fixed by gauge invariance of the Yukawa interactions
# and the four-scalar interaction.  The historical classification uses
# doubled hypercharge labels; in our physical convention Q = T3 + Y:
#   Y(S1) = alpha / 2
#   Y(F)  = (alpha + 1) / 2
#   Y(S2) = (alpha + 2) / 2
# The A-E labels only specify the SU(2) representations.
INTERESTING = [("A", 0), ("B", -1), ("C", -1), ("D", -2), ("E", 0)]

SMOKE = [("B", -1), ("C", -1), ("A", 0), ("D", -2), ("E", 0)]

EXTENDED = [
    ("A", 0),
    ("A", -2),
    ("B", -1),
    ("C", -1),
    ("D", -2),
    ("E", 0),
    ("E", -2),
]


def encode_alpha(alpha: int) -> str:
    """For Wolfram, recode +5 -> p5 and -5 -> m5."""
    return f"m{abs(alpha)}" if alpha < 0 else f"p{alpha}"

def valid_t3_dimensions(d_s1: int, d_s2: int, d_f: int) -> bool:
    """Check whether the SU(2) dimensions are supported by the current T3 pipeline."""

    # Current production scope: singlet, doublet and triplet only.
    # Higher representations are intentionally deferred for possible future work.
    if any(d > 3 for d in (d_s1, d_s2, d_f)):
        return False

    if min(d_s1, d_s2, d_f) < 1:
        return False

    # Each Yukawa contains a lepton doublet.  For SU(2), a scalar that can
    # couple to L and F must therefore have dS = dF +/- 1.
    if abs(d_s1 - d_f) != 1 or abs(d_s2 - d_f) != 1:
        return False

    # The identical Higgs pair is symmetric, so H x H contributes through
    # the triplet (J=1) channel.  S1 x S2 must therefore contain J=1.
    j1 = (d_s1 - 1) / 2
    j2 = (d_s2 - 1) / 2

    return abs(j1 - j2) <= 1 <= j1 + j2 and float(j1 + j2).is_integer()

def identify_t3_class(d_s1: int, d_s2: int, d_f: int) -> str | None:
    """Return the known A-E label if these dimensions match one."""

    dimensions = (d_s1, d_s2, d_f)

    for model_class, known_dimensions in T3_CLASSES.items():
        if dimensions == known_dimensions:
            return model_class

    return None

