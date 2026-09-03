from __future__ import annotations

from collections import Counter
from itertools import combinations_with_replacement
from typing import Iterable


# We store an SU(2) irrep by twice its spin:
# singlet -> 0, doublet -> 1, triplet -> 2.
#
# For the current project scope d <= 3:
#     twice_spin = dimension - 1.
T3_CLASSES = {
    "A": ((1, 3, 2), 0),
    "B": ((2, 2, 1), -1),
    "C": ((2, 2, 3), -1),
    "D": ((3, 1, 2), -2),
    "E": ((3, 3, 2), 0),
}

SPECIES = ("H", "H†", "S1", "S1†", "S2", "S2†")


def compositions(total: int, parts: int, prefix: tuple[int, ...] = ()) -> Iterable[tuple[int, ...]]:
    if parts == 1:
        yield prefix + (total,)
        return

    for value in range(total + 1):
        yield from compositions(total - value, parts - 1, prefix + (value,))


def symmetric_power_decomposition(twice_spin: int, power: int) -> dict[int, int]:
    """Decompose Sym^power(j) into SU(2) irreps.

    Repeated identical scalar fields commute, so repeated copies of the same
    species live in the symmetric tensor power rather than the unrestricted
    tensor product.
    """

    if power == 0:
        return {0: 1}

    weights = list(range(-twice_spin, twice_spin + 1, 2))
    weight_multiplicity: Counter[int] = Counter()

    for chosen_weights in combinations_with_replacement(weights, power):
        weight_multiplicity[sum(chosen_weights)] += 1

    decomposition: Counter[int] = Counter()

    while weight_multiplicity:
        # Remove zero entries created by previous subtractions.
        weight_multiplicity += Counter()
        if not weight_multiplicity:
            break

        highest_weight = max(weight_multiplicity)
        multiplicity = weight_multiplicity[highest_weight]
        decomposition[highest_weight] += multiplicity

        for weight in range(-highest_weight, highest_weight + 1, 2):
            weight_multiplicity[weight] -= multiplicity
            if weight_multiplicity[weight] == 0:
                del weight_multiplicity[weight]

    return dict(decomposition)


def tensor_product(
    left: dict[int, int],
    right: dict[int, int],
) -> dict[int, int]:
    result: Counter[int] = Counter()

    for j1, multiplicity1 in left.items():
        for j2, multiplicity2 in right.items():
            for j in range(abs(j1 - j2), j1 + j2 + 1, 2):
                result[j] += multiplicity1 * multiplicity2

    return dict(result)


def su2_singlet_multiplicity(
    counts: tuple[int, ...],
    d_s1: int,
    d_s2: int,
) -> int:
    reps = (
        1,          # H
        1,          # H†
        d_s1 - 1,   # S1
        d_s1 - 1,   # S1†
        d_s2 - 1,   # S2
        d_s2 - 1,   # S2†
    )

    decomposition = {0: 1}

    for rep, count in zip(reps, counts):
        decomposition = tensor_product(
            decomposition,
            symmetric_power_decomposition(rep, count),
        )

    return decomposition.get(0, 0)


def conjugate_counts(counts: tuple[int, ...]) -> tuple[int, ...]:
    return (
        counts[1],
        counts[0],
        counts[3],
        counts[2],
        counts[5],
        counts[4],
    )


def canonical_under_hc(counts: tuple[int, ...]) -> tuple[int, ...]:
    conjugate = conjugate_counts(counts)
    return min(counts, conjugate)


def monomial_name(counts: tuple[int, ...]) -> str:
    factors: list[str] = []

    for species, count in zip(SPECIES, counts):
        if count == 0:
            continue
        factors.append(species if count == 1 else f"{species}^{count}")

    return " ".join(factors)


def current_potential_structures(alpha: int) -> dict[tuple[int, ...], int]:
    """Field multisets currently represented in T3RGBetaModel.wl.

    The value is how many independent SU(2) contractions are explicitly
    represented by the current model for that field multiset.
    """

    structures = {
        canonical_under_hc((2, 2, 0, 0, 0, 0)): 1,  # lambdaH
        canonical_under_hc((0, 0, 2, 2, 0, 0)): 1,  # lambdaS1
        canonical_under_hc((0, 0, 0, 0, 2, 2)): 1,  # lambdaS2
        canonical_under_hc((1, 1, 1, 1, 0, 0)): 1,  # lambdaH1
        canonical_under_hc((1, 1, 0, 0, 1, 1)): 1,  # lambdaH2
        canonical_under_hc((0, 0, 1, 1, 1, 1)): 1,  # lambda12
        canonical_under_hc((2, 0, 1, 0, 0, 1)): 1,  # lambdaT3 + h.c.
    }

    return structures


def allowed_quartics(
    d_s1: int,
    d_s2: int,
    alpha: int,
) -> list[tuple[tuple[int, ...], int, bool]]:
    # Twice the physical hypercharge, so every charge is integral:
    #
    #   2Y(H)  = 1
    #   2Y(S1) = alpha
    #   2Y(S2) = alpha + 2
    hypercharges = (
        1,
        -1,
        alpha,
        -alpha,
        alpha + 2,
        -(alpha + 2),
    )

    found: list[tuple[tuple[int, ...], int, bool]] = []
    seen: set[tuple[int, ...]] = set()

    for counts in compositions(4, 6):
        # H is Z2-even; S1 and S2 are Z2-odd.
        bsm_count = sum(counts[2:])
        if bsm_count % 2 != 0:
            continue

        total_hypercharge = sum(
            count * hypercharge
            for count, hypercharge in zip(counts, hypercharges)
        )
        if total_hypercharge != 0:
            continue

        multiplicity = su2_singlet_multiplicity(
            counts,
            d_s1,
            d_s2,
        )
        if multiplicity == 0:
            continue

        canonical = canonical_under_hc(counts)
        if canonical in seen:
            continue
        seen.add(canonical)

        has_distinct_hc = conjugate_counts(canonical) != canonical
        found.append((canonical, multiplicity, has_distinct_hc))

    return found


def main() -> None:
    for label, ((d_s1, d_s2, _d_f), alpha) in T3_CLASSES.items():
        allowed = allowed_quartics(d_s1, d_s2, alpha)
        current = current_potential_structures(alpha)

        print()
        print("=" * 78)
        print(f"T3-{label}: dS1={d_s1}, dS2={d_s2}, alpha={alpha}")
        print("=" * 78)

        missing_total = 0

        for counts, allowed_multiplicity, has_hc in allowed:
            present_multiplicity = current.get(counts, 0)
            missing = max(allowed_multiplicity - present_multiplicity, 0)

            if missing == 0:
                status = "covered"
            elif present_multiplicity:
                status = f"MISSING {missing} extra SU(2) contraction(s)"
            else:
                status = f"MISSING field structure ({allowed_multiplicity} invariant(s))"

            hc = " + h.c." if has_hc else ""

            print(
                f"{monomial_name(counts):34s}"
                f" multiplicity={allowed_multiplicity}"
                f"  current={present_multiplicity}"
                f"  -> {status}{hc}"
            )

            missing_total += missing

        total_allowed = sum(item[1] for item in allowed)
        total_current = sum(
            min(current.get(counts, 0), multiplicity)
            for counts, multiplicity, _ in allowed
        )

        print()
        print(
            f"Independent gauge+Z2 quartic invariants: {total_allowed}; "
            f"currently represented: {total_current}; "
            f"missing: {missing_total}"
        )


if __name__ == "__main__":
    main()
