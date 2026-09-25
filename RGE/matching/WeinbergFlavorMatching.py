"""Flavor lifting of a one-generation matched Weinberg coefficient.

The Matchete matching result uses one lepton Yukawa at each T3 vertex. Gauge
and SU(2) contractions are generation independent, so the one-generation
coefficient can be factorized as

    kappa_1g = F_loop * conjugate(y1) * conjugate(y2).

This module owns the symbolic lift of that result to the physical symmetric
Majorana coefficient C5_pq. It does not parse the hierarchical final-C5 JSON.
"""

from __future__ import annotations

from pathlib import Path

import sympy as sp

from RGE.matching.MatcheteC5Parsing import parse_matchete_c5


def extract_one_generation_loop_kernel(
    kappa_1g: sp.Expr,
    *,
    y1_name: str = "y1",
    y2_name: str = "y2",
) -> sp.Expr:
    """Remove the one-generation Yukawa product from the matched C5."""

    y1 = sp.Symbol(y1_name)
    y2 = sp.Symbol(y2_name)

    yukawa_product = sp.conjugate(y1) * sp.conjugate(y2)
    kernel = sp.cancel(kappa_1g / yukawa_product)

    if sp.simplify(kernel * yukawa_product - kappa_1g) != 0:
        raise ValueError(
            "Could not factor the matched C5 into "
            "F_loop * conjugate(y1) * conjugate(y2)."
        )

    return sp.factor(kernel)


def symbolic_t3_yukawas(
    n_lepton: int = 3,
    n_heavy: int = 3,
) -> tuple[sp.Matrix, sp.Matrix]:
    """Return symbolic y1_pr and y2_pr matrices."""

    y1 = sp.Matrix(
        n_lepton,
        n_heavy,
        lambda p, r: sp.Symbol(f"y1_{p + 1}{r + 1}"),
    )
    y2 = sp.Matrix(
        n_lepton,
        n_heavy,
        lambda p, r: sp.Symbol(f"y2_{p + 1}{r + 1}"),
    )

    return y1, y2


def build_majorana_c5_flavor_matrix(
    kernel: sp.Expr,
    y1: sp.MatrixBase,
    y2: sp.MatrixBase,
    *,
    fermion_mass_symbol: sp.Symbol | None = None,
    heavy_masses: list[sp.Expr] | tuple[sp.Expr, ...] | None = None,
) -> sp.Matrix:
    """Build the symmetric lepton-flavor Majorana C5 matrix."""

    if y1.shape != y2.shape:
        raise ValueError("y1 and y2 must have the same shape.")

    n_lepton, n_heavy = y1.shape

    if heavy_masses is not None and len(heavy_masses) != n_heavy:
        raise ValueError(
            "heavy_masses must contain one mass per heavy generation."
        )

    if heavy_masses is not None and fermion_mass_symbol is None:
        fermion_mass_symbol = sp.Symbol("MF")

    def entry(p: int, q: int) -> sp.Expr:
        value = sp.S.Zero

        for r in range(n_heavy):
            kernel_r = kernel

            if heavy_masses is not None:
                kernel_r = kernel_r.subs(
                    fermion_mass_symbol,
                    heavy_masses[r],
                )

            value += kernel_r * (
                sp.conjugate(y1[p, r]) * sp.conjugate(y2[q, r])
                + sp.conjugate(y2[p, r]) * sp.conjugate(y1[q, r])
            )

        return value

    c5 = sp.MutableDenseMatrix.zeros(n_lepton, n_lepton)

    for p in range(n_lepton):
        for q in range(p, n_lepton):
            value = entry(p, q)
            c5[p, q] = value
            c5[q, p] = value

    return sp.Matrix(c5)


def build_flavor_c5_from_matchete(
    c5_path: Path,
    *,
    n_lepton: int = 3,
    n_heavy: int = 3,
    split_heavy_masses: bool = True,
) -> dict:
    """Read a one-generation Matchete C5 and lift it to full flavor."""

    c5_path = Path(c5_path)

    kappa_1g = parse_matchete_c5(
        c5_path.read_text(encoding="utf-8")
    )
    kernel = extract_one_generation_loop_kernel(kappa_1g)

    y1, y2 = symbolic_t3_yukawas(
        n_lepton=n_lepton,
        n_heavy=n_heavy,
    )

    if split_heavy_masses:
        heavy_masses = [
            sp.Symbol(f"MF{r + 1}")
            for r in range(n_heavy)
        ]
    else:
        heavy_masses = None

    c5 = build_majorana_c5_flavor_matrix(
        kernel,
        y1,
        y2,
        heavy_masses=heavy_masses,
    )

    return {
        "kappa_1g": kappa_1g,
        "kernel": kernel,
        "y1": y1,
        "y2": y2,
        "heavy_masses": heavy_masses,
        "K": c5,
    }


def write_flavor_matching_outputs(
    output_dir: Path,
    result: dict,
    *,
    debug_outputs: bool = False,
) -> dict:
    """Write the existing full-flavor debug artifacts."""

    output_dir = Path(output_dir)
    debug_dir = output_dir / "debug"

    matrix_path = debug_dir / "c5_flavor_matrix.txt"
    kernel_path = debug_dir / "c5_loop_kernel.txt"

    if debug_outputs:
        debug_dir.mkdir(parents=True, exist_ok=True)
        matrix_path.write_text(
            sp.sstr(result["K"]) + "\n",
            encoding="utf-8",
        )
        kernel_path.write_text(
            sp.sstr(result["kernel"]) + "\n",
            encoding="utf-8",
        )

    return {
        "C5FlavorMatrixFile": (
            matrix_path.relative_to(output_dir).as_posix()
            if debug_outputs
            else ""
        ),
        "C5LoopKernelFile": (
            kernel_path.relative_to(output_dir).as_posix()
            if debug_outputs
            else ""
        ),
        "HeavyFlavorGenerations": result["y1"].cols,
        "LeptonFlavorGenerations": result["y1"].rows,
    }


# Historical function aliases. New code should use the descriptive names above.
extract_t3_loop_kernel = extract_one_generation_loop_kernel
build_flavor_c5_matrix = build_majorana_c5_flavor_matrix
match_c5_flavor_from_file = build_flavor_c5_from_matchete
