"""Exhaustive scalar-operator counting for the 16 supported ordinary T3 points.

This is a *static audit* of the source at Git commit
41194fcd78fc7fe3c9378c9f8ee41228ad6473e7 (24 October 2026).
It does not run Matchete, RGBeta, matching, or numerical evolution.

Count SU(2) singlets in tensor products of symmetric powers of identical
commuting scalar species. Apply U(1)_Y and Z2, then identify a monomial and
its Hermitian conjugate. Report independent real coefficient dimensions
before quotienting by possible field-phase redefinitions.

Project path: tests/python/check_t3_scalar_quartics.py
Run: python tests/python/check_t3_scalar_quartics.py --output output/uv_operator_audit.json
"""
from __future__ import annotations

import argparse
from collections import Counter
from itertools import combinations_with_replacement
import json
from pathlib import Path

HEAD_SHA = "41194fcd78fc7fe3c9378c9f8ee41228ad6473e7"
MODEL_POINTS = {
    "A": ((1, 3, 2), (-4, -2, 0)),
    "B": ((2, 2, 1), (-3, -1, 1)),
    "C": ((2, 2, 3), (-3, -1, 1)),
    "D": ((3, 1, 2), (-2, 0, 2)),
    "E": ((3, 3, 2), (-4, -2, 0, 2)),
}
SPECIES = ("H", "Hdag", "S1", "S1dag", "S2", "S2dag")


def compositions(total: int, parts: int):
    if parts == 1:
        yield (total,)
    else:
        for n in range(total + 1):
            for rest in compositions(total - n, parts - 1):
                yield (n, *rest)


def symmetric_power(twice_j: int, count: int) -> dict[int, int]:
    """SU(2) Sym^count(j), with integer weights represented as 2j."""
    if count == 0:
        return {0: 1}
    weights = range(-twice_j, twice_j + 1, 2)
    multiplicities = Counter(
        sum(chosen)
        for chosen in combinations_with_replacement(weights, count)
    )
    irreps: Counter[int] = Counter()
    while multiplicities:
        maximum = max(multiplicities)
        n = multiplicities[maximum]
        irreps[maximum] += n
        for weight in range(-maximum, maximum + 1, 2):
            multiplicities[weight] -= n
            if multiplicities[weight] == 0:
                del multiplicities[weight]
        if any(v < 0 for v in multiplicities.values()):
            raise AssertionError("Invalid SU(2) symmetric-power decomposition")
    return dict(irreps)


def singlet_multiplicity(counts: tuple[int, ...], d1: int, d2: int) -> int:
    twice_spins = (1, 1, d1 - 1, d1 - 1, d2 - 1, d2 - 1)
    current: dict[int, int] = {0: 1}
    for twice_j, power in zip(twice_spins, counts):
        nxt: Counter[int] = Counter()
        for left, left_n in current.items():
            for right, right_n in symmetric_power(twice_j, power).items():
                for joined in range(abs(left - right), left + right + 1, 2):
                    nxt[joined] += left_n * right_n
        current = dict(nxt)
    return current.get(0, 0)


def conjugate(counts: tuple[int, ...]) -> tuple[int, ...]:
    return (counts[1], counts[0], counts[3], counts[2], counts[5], counts[4])


def allowed(d1: int, d2: int, alpha: int, degree: int) -> dict[tuple[int, ...], dict]:
    twice_y = (1, -1, alpha, -alpha, alpha + 2, -alpha - 2)
    found: dict[tuple[int, ...], dict] = {}
    for n in compositions(degree, 6):
        if (sum(n[2:]) % 2) != 0:  # exact Z2
            continue
        if sum(a * y for a, y in zip(n, twice_y)) != 0:
            continue
        number = singlet_multiplicity(n, d1, d2)
        if not number:
            continue
        key = min(n, conjugate(n))
        paired = key != conjugate(key)
        found[key] = {
            "counts": list(key),
            "field_structure": " ".join(
                name + (f"^{power}" if power > 1 else "")
                for name, power in zip(SPECIES, key) if power
            ),
            "invariants": number,
            "complex_coefficient": paired,
            "real_parameters": number * (2 if paired else 1),
        }
    return found


def expected_builder_structures() -> tuple[tuple[int, ...], ...]:
    """Structures inspected in Lagrangian/interactions/{ScalarPotential,T3Topology}.wl.

    H's own quartic comes from LoadModel[\"SM\"]. Each has all its independent
    SU(2) contractions (not necessarily the same coupling basis as RGBeta).
    """
    return (
        (2, 2, 0, 0, 0, 0),  # H^2 Hdag^2
        (0, 0, 2, 2, 0, 0),  # S1^2 S1dag^2
        (0, 0, 0, 0, 2, 2),  # S2^2 S2dag^2
        (1, 1, 1, 1, 0, 0),  # portal 1
        (1, 1, 0, 0, 1, 1),  # portal 2
        (0, 0, 1, 1, 1, 1),  # cross scalar
        (2, 0, 1, 0, 0, 1),  # T3 topology
    )


def extra_rgbeta_b_c_alpha_minus_one() -> tuple[tuple[int, ...], ...]:
    """Six additional, named complex quartics in T3RGBetaAddRGClosedQuartics."""
    return (
        (0, 2, 0, 0, 2, 0),  # lambdaHHdagS2S2
        (0, 2, 0, 2, 0, 0),  # lambdaHHdagS1barS1bar
        (0, 0, 0, 2, 0, 2),  # lambdaS1bar2S2bar2
        (0, 0, 0, 1, 1, 2),  # lambdaS1barS2S2bar2
        (0, 0, 1, 2, 0, 1),  # lambdaS1S1bar2S2bar
        (1, 1, 0, 1, 0, 1),  # lambdaHHdagS1barS2barCross
    )


def canonical(mon: tuple[int, ...]) -> tuple[int, ...]:
    return min(mon, conjugate(mon))


def inspection(d1: int, d2: int, alpha: int) -> dict:
    all_q4 = allowed(d1, d2, alpha, 4)
    all_q2 = allowed(d1, d2, alpha, 2)
    representation = {canonical(m) for m in expected_builder_structures()}
    rgbeta_special = (
        {canonical(m) for m in extra_rgbeta_b_c_alpha_minus_one()}
        if (d1, d2, alpha) == (2, 2, -1)
        else set()
    )

    # RGBeta lists all generic two-pair SU(2) basis contractions; the
    # topology source guarantees that the topology contraction is one.
    represented_matchete = {
        k: all_q4[k]["invariants"] for k in representation
        if k in all_q4
    }
    represented_rgbeta = dict(represented_matchete)
    for k in rgbeta_special:
        represented_rgbeta[k] = min(
            all_q4[k]["invariants"], represented_rgbeta.get(k, 0) + 1
        )

    def count_real(mapping: dict) -> int:
        return sum(
            n * (2 if all_q4[k]["complex_coefficient"] else 1)
            for k, n in mapping.items()
        )

    total_q4 = sum(entry["real_parameters"] for entry in all_q4.values())
    missing_mat = []
    missing_rg = []
    for key, entry in all_q4.items():
        for out, represented in (
            (missing_mat, represented_matchete),
            (missing_rg, represented_rgbeta),
        ):
            delta = entry["invariants"] - represented.get(key, 0)
            if delta:
                out.append({
                    **entry,
                    "missing_invariants": delta,
                    "missing_real_parameters": delta * (
                        2 if entry["complex_coefficient"] else 1
                    ),
                })

    # Scalar masses in both ordinary implementations currently contain only
    # S1dag S1, S2dag S2: the remaining bilinears are genuine omissions from
    # the unconstrained complex-field potential.
    standard_q2 = {
        canonical((0, 0, 1, 1, 0, 0)),
        canonical((0, 0, 0, 0, 1, 1)),
    }
    missing_q2 = [v for k, v in all_q2.items() if k not in standard_q2
                  and k != canonical((1, 1, 0, 0, 0, 0))]

    return {
        "dimension_two_missing": missing_q2,
        "quartic_total_real_parameters_including_SM_Higgs": total_q4,
        "quartic_matchete_real_parameters": count_real(represented_matchete),
        "quartic_rgbeta_real_parameters": count_real(represented_rgbeta),
        "quartic_missing_matchete": missing_mat,
        "quartic_missing_rgbeta": missing_rg,
        "quartic_missing_matchete_real_parameters": total_q4 - count_real(represented_matchete),
        "quartic_missing_rgbeta_real_parameters": total_q4 - count_real(represented_rgbeta),
    }


def audit():
    output = {}
    for name, (dims, alphas) in MODEL_POINTS.items():
        for a in alphas:
            key = f"{name}_{a:+d}"
            output[key] = {
                "class": name,
                "alpha": a,
                "dimensions": list(dims),
                **inspection(dims[0], dims[1], a),
            }
    return output


def self_test() -> None:
    # Explicit representation-theoretic crosschecks, and regression against
    # the inspected 2026-10-24 source: these are counts, not physics runs.
    expected = {
        "A_-4": (10, 10, 10, 0),
        "A_-2": (20, 10, 10, 1),
        "A_+0": (20, 10, 10, 1),
        "B_-3": (13, 11, 11, 0),
        "B_-1": (25, 11, 23, 1),
        "B_+1": (13, 11, 11, 0),
        "C_-3": (13, 11, 11, 0),
        "C_-1": (25, 11, 23, 1),
        "C_+1": (13, 11, 11, 0),
        "D_-2": (20, 10, 10, 1),
        "D_+0": (20, 10, 10, 1),
        "D_+2": (10, 10, 10, 0),
        "E_-4": (14, 14, 14, 0),
        "E_-2": (26, 14, 14, 1),
        "E_+0": (26, 14, 14, 1),
        "E_+2": (14, 14, 14, 0),
    }
    data = audit()
    assert len(data) == len(expected) == 16
    for key, (all_c, mat_c, rg_c, n_missing_mass) in expected.items():
        value = data[key]
        actual = (
            value["quartic_total_real_parameters_including_SM_Higgs"],
            value["quartic_matchete_real_parameters"],
            value["quartic_rgbeta_real_parameters"],
            len(value["dimension_two_missing"]),
        )
        assert actual == (all_c, mat_c, rg_c, n_missing_mass), (key, actual)
    assert singlet_multiplicity((0, 0, 2, 2, 0, 0), 3, 1) == 2
    assert singlet_multiplicity((0, 0, 1, 1, 1, 1), 3, 3) == 3


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--require-complete", action="store_true",
                        help="Return exit code 1 when any allowed operator is absent")
    args = parser.parse_args()
    self_test()
    data = audit()
    print("T3 scalar completeness: independent real quartic parameters, INCLUDING SM Higgs")
    print("Case     full  Matchete  RGBeta  missing Mat  missing RG  extra dim-2")
    incomplete = False
    for value in data.values():
        missing_mat = value["quartic_missing_matchete_real_parameters"]
        missing_rg = value["quartic_missing_rgbeta_real_parameters"]
        n2 = len(value["dimension_two_missing"])
        incomplete |= bool(missing_mat or missing_rg or n2)
        print(f"{value['class']} {value['alpha']:+3d}   "
              f"{value['quartic_total_real_parameters_including_SM_Higgs']:>4d}  "
              f"{value['quartic_matchete_real_parameters']:>8d}  "
              f"{value['quartic_rgbeta_real_parameters']:>6d}  "
              f"{missing_mat:>11d}  {missing_rg:>10d}  {n2:>11d}")
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps({"source_commit": HEAD_SHA,
                           "scope": "static complex-scalar operator count",
                           "models": data}, indent=2) + "\n", encoding="utf-8")
        print("Wrote", args.output)
    print("PASS: self-test; no Matchete, RGBeta or numerical runs were performed")
    return int(args.require_complete and incomplete)


if __name__ == "__main__":
    raise SystemExit(main())
