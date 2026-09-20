from __future__ import annotations

"""Weinberg-normalization, flavor-symmetry, and EFT1 Wilson-prefactor regressions.

Default mode
------------
Run the normalization-cleanup regression plus the mixed-C12 flavor-symmetry
checks consolidated from ``RGE/running/C12FlavorSymmetryValidation.py``:

* ``build_flavor_c5_matrix`` must produce the physical symmetric coefficient
  C5 = A_pq + A_qp, giving 2*A_ordered in the one-generation symbolic test.
* The mixed C12 flavor kernel must be symmetric for non-diagonal complex
  Yukawas and equal 1/2(K + K^T), where K is the ordered y1* y2* kernel.
* ``build_neutrino_mass_matrix`` must preserve the previous physical mass
  bookkeeping.

Optional Wilson-adapter mode
----------------------------
If a Matchete EFT1 Wilson-seed JSON is supplied, also check that
``EFT1TensorAdapters.build_eft1_wilson_tensor`` preserves an arbitrary outer
numerical prefactor in every PL tree Wilson term.  This is the former
``tests/check_wilson_prefactor_preservation.py`` regression, consolidated here
without changing the adapter or any physics convention.
"""

import argparse
import copy
import json
from pathlib import Path
import sys

import sympy as sp

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from RGE.matching.FlavorC5Matching import build_flavor_c5_matrix
from RGE.running.eft1.EFT1TensorAdapters import build_eft1_wilson_tensor
from RGE.running.weinberg.FlavorMatchedWeinbergStage import build_neutrino_mass_matrix


def canonical_items(tensor) -> dict[tuple[int, ...], sp.Expr]:
    items: dict[tuple[int, ...], sp.Expr] = {}

    for key, raw in tensor.nonzero_items():
        value = sp.simplify(raw)
        if value != 0:
            items[tuple(key)] = value

    return items


def check_scaled_term(
    seed_payload: dict,
    term: dict,
    d_s1: int,
    d_s2: int,
    *,
    scale: sp.Integer,
) -> dict:
    one = copy.deepcopy(seed_payload)
    one["TreeWilsonTerms"] = [copy.deepcopy(term)]
    one["TreeWilsonTermCount"] = 1

    baseline = build_eft1_wilson_tensor(
        one,
        d_s1=d_s1,
        d_s2=d_s2,
        chirality="PL",
        project_operator_symmetry=True,
    )
    baseline_items = canonical_items(baseline)

    scaled = copy.deepcopy(one)
    scaled_term = scaled["TreeWilsonTerms"][0]
    scaled_term["TermInputForm"] = (
        f"({sp.sstr(scale)})*({scaled_term['TermInputForm']})"
    )

    scaled_tensor = build_eft1_wilson_tensor(
        scaled,
        d_s1=d_s1,
        d_s2=d_s2,
        chirality="PL",
        project_operator_symmetry=True,
    )
    scaled_items = canonical_items(scaled_tensor)

    keys = sorted(set(baseline_items) | set(scaled_items))
    failures = []
    observed_ratios = set()

    for key in keys:
        baseline_value = sp.simplify(
            baseline_items.get(key, 0)
        )
        scaled_value = sp.simplify(
            scaled_items.get(key, 0)
        )

        if baseline_value == 0:
            if scaled_value != 0:
                failures.append(
                    (
                        key,
                        baseline_value,
                        scaled_value,
                        "new component appeared",
                    )
                )
            continue

        ratio = sp.simplify(
            scaled_value / baseline_value
        )
        observed_ratios.add(
            sp.sstr(ratio)
        )

        if sp.simplify(
            scaled_value - scale * baseline_value
        ) != 0:
            failures.append(
                (
                    key,
                    baseline_value,
                    scaled_value,
                    f"ratio={ratio}",
                )
            )

    return {
        "term_index": int(term["Index"]),
        "cg_names": term.get("CGNames", []),
        "baseline_nonzero": len(baseline_items),
        "scaled_nonzero": len(scaled_items),
        "expected_scale": sp.sstr(scale),
        "observed_ratios": sorted(observed_ratios),
        "passes_exact_scaling": not failures,
        "failure_count": len(failures),
        "first_failures": [
            {
                "component": list(key),
                "baseline": sp.sstr(baseline_value),
                "scaled": sp.sstr(scaled_value),
                "detail": detail,
            }
            for (
                key,
                baseline_value,
                scaled_value,
                detail,
            ) in failures[:10]
        ],
    }


def _ordered_c12_kernel(
    y1: sp.MatrixBase,
    y2: sp.MatrixBase,
    masses: list[sp.Expr] | tuple[sp.Expr, ...],
) -> sp.Matrix:
    """Return K_pq = sum_r y1^*_{pr} y2^*_{qr}/M_r."""
    if y1.shape != y2.shape:
        raise ValueError("y1 and y2 must have the same shape.")

    n_lepton, n_heavy = y1.shape
    if len(masses) != n_heavy:
        raise ValueError("Need one heavy-fermion mass per heavy generation.")

    return sp.Matrix(
        n_lepton,
        n_lepton,
        lambda p, q: sp.simplify(
            sum(
                sp.conjugate(y1[p, r])
                * sp.conjugate(y2[q, r])
                / masses[r]
                for r in range(n_heavy)
            )
        ),
    )


def check_c12_flavor_symmetry() -> bool:
    """Preserve the former standalone non-diagonal C12 flavor validation."""
    y1 = sp.Matrix(
        [
            [1 + sp.I, 2],
            [3, 1 - 2 * sp.I],
            [2 - sp.I, -1],
        ]
    )
    y2 = sp.Matrix(
        [
            [2, -sp.I],
            [1 + sp.I, 4],
            [-2, 3 + sp.I],
        ]
    )
    masses = [sp.Integer(5), sp.Integer(7)]

    ordered = _ordered_c12_kernel(y1, y2, masses)

    # build_flavor_c5_matrix substitutes MF -> M_r inside the supplied
    # kernel for each heavy generation.  Therefore use the explicit 1/MF
    # kernel here so the result is the physical symmetric sum K + K^T.
    # The mixed C12 convention itself is 1/2(K + K^T), so compare after the
    # explicit factor of 1/2.
    MF = sp.Symbol("MF")
    physical_c5 = build_flavor_c5_matrix(
        1 / MF,
        y1,
        y2,
        fermion_mass_symbol=MF,
        heavy_masses=masses,
    )
    mixed_c12 = sp.Rational(1, 2) * physical_c5

    symmetry_residual = (mixed_c12 - mixed_c12.T).applyfunc(sp.simplify)
    reconstruction = (
        mixed_c12
        - sp.Rational(1, 2) * (ordered + ordered.T)
    ).applyfunc(sp.simplify)

    ordered_is_nonsymmetric = (
        (ordered - ordered.T).applyfunc(sp.simplify)
        != sp.zeros(3)
    )

    ok = (
        symmetry_residual == sp.zeros(3)
        and reconstruction == sp.zeros(3)
        and ordered_is_nonsymmetric
    )

    print()
    print("=" * 72)
    print("MIXED C12 FLAVOR-SYMMETRY REGRESSION")
    print("=" * 72)
    print(
        "PASS" if ok else "FAIL",
        ": C12 = 1/2 (K + K^T) for non-diagonal complex Yukawas",
    )

    return ok


def check_weinberg_normalization_cleanup() -> bool:
    A, y1, y2, v = sp.symbols(
        "A y1 y2 v"
    )

    kernel = A
    Y1 = sp.Matrix([[y1]])
    Y2 = sp.Matrix([[y2]])

    C5 = build_flavor_c5_matrix(
        kernel,
        Y1,
        Y2,
        heavy_masses=None,
    )

    ordered = (
        A
        * sp.conjugate(y1)
        * sp.conjugate(y2)
    )
    c5_difference = sp.simplify(
        C5[0, 0] - 2 * ordered
    )

    m_new = build_neutrino_mass_matrix(
        C5,
        vev=v,
    )
    m_old = -v**2 * ordered
    mass_difference = sp.simplify(
        m_new[0, 0] - m_old
    )

    print("=" * 72)
    print(
        "WEINBERG NORMALIZATION CLEANUP REGRESSION"
    )
    print("=" * 72)
    print()
    print(
        "Difference C5_physical - "
        "2*A_ordered:"
    )
    print(c5_difference)
    print()
    print(
        "Physical neutrino-mass difference "
        "old vs new bookkeeping:"
    )
    print(mass_difference)
    print()

    if c5_difference != 0:
        print(
            "FAIL: physical C5 normalization"
        )
        return False

    if mass_difference != 0:
        print(
            "FAIL: physical neutrino mass changed"
        )
        return False

    print("PASS")
    return True


def check_wilson_prefactor_preservation(
    seed_path: Path,
    *,
    d_s1: int,
    d_s2: int,
    scale: int,
) -> dict:
    seed = json.loads(
        seed_path.read_text(
            encoding="utf-8"
        )
    )

    pl_terms = [
        term
        for term in seed.get(
            "TreeWilsonTerms",
            [],
        )
        if term.get("Chirality") == "PL"
    ]

    if not pl_terms:
        raise RuntimeError(
            "No PL terms found in the Wilson seed."
        )

    results = [
        check_scaled_term(
            seed,
            term,
            d_s1=d_s1,
            d_s2=d_s2,
            scale=sp.Integer(scale),
        )
        for term in pl_terms
    ]

    return {
        "status": "Success",
        "seed": str(seed_path),
        "all_terms_preserve_outer_prefactor": all(
            item["passes_exact_scaling"]
            for item in results
        ),
        "terms": results,
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Run the Weinberg-normalization and C12 flavor-symmetry "
            "regressions and, optionally, the EFT1 Wilson outer-prefactor "
            "preservation regression."
        )
    )
    parser.add_argument(
        "--wilson-seed",
        type=Path,
        default=None,
        help=(
            "Optional Matchete EFT1 Wilson seed JSON. "
            "If supplied, run the exact outer-prefactor "
            "preservation check as well."
        ),
    )
    parser.add_argument(
        "--d-s1",
        type=int,
        default=None,
        help=(
            "S1 SU(2) dimension for --wilson-seed."
        ),
    )
    parser.add_argument(
        "--d-s2",
        type=int,
        default=None,
        help=(
            "S2 SU(2) dimension for --wilson-seed."
        ),
    )
    parser.add_argument(
        "--scale",
        type=int,
        default=7,
        help=(
            "Arbitrary outer numerical scaling used by "
            "the Wilson-prefactor regression."
        ),
    )
    parser.add_argument(
        "--prefactor-output",
        type=Path,
        default=None,
        help=(
            "Optional JSON output for the Wilson-prefactor "
            "regression."
        ),
    )
    args = parser.parse_args()

    normalization_ok = (
        check_weinberg_normalization_cleanup()
    )
    flavor_symmetry_ok = check_c12_flavor_symmetry()

    prefactor_ok = True

    if args.wilson_seed is not None:
        if args.d_s1 is None or args.d_s2 is None:
            parser.error(
                "--d-s1 and --d-s2 are required "
                "with --wilson-seed."
            )

        print()
        print("=" * 72)
        print(
            "EFT1 WILSON OUTER-PREFACTOR REGRESSION"
        )
        print("=" * 72)

        payload = (
            check_wilson_prefactor_preservation(
                args.wilson_seed,
                d_s1=args.d_s1,
                d_s2=args.d_s2,
                scale=args.scale,
            )
        )

        print(
            json.dumps(
                payload,
                indent=2,
            )
        )

        if args.prefactor_output is not None:
            args.prefactor_output.parent.mkdir(
                parents=True,
                exist_ok=True,
            )
            args.prefactor_output.write_text(
                json.dumps(
                    payload,
                    indent=2,
                ),
                encoding="utf-8",
            )

        prefactor_ok = bool(
            payload[
                "all_terms_preserve_outer_prefactor"
            ]
        )

    return (
        0
        if normalization_ok and flavor_symmetry_ok and prefactor_ok
        else 1
    )


if __name__ == "__main__":
    raise SystemExit(main())
