"""Flavor lifting of matched T3 Weinberg coefficients.

The original two-vertex restricted matching remains unchanged. Additional
Yukawa channels are supported by ``build_general_majorana_c5_flavor_matrix``;
each channel requires a genuinely matched one-generation loop kernel.
This file does not derive scalar mixing, SU(2) factors or loop integrals.
"""

from __future__ import annotations

from pathlib import Path
from typing import Mapping

import numpy as np
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
    y1 = sp.Matrix(n_lepton, n_heavy,
                   lambda p, r: sp.Symbol(f"y1_{p + 1}{r + 1}"))
    y2 = sp.Matrix(n_lepton, n_heavy,
                   lambda p, r: sp.Symbol(f"y2_{p + 1}{r + 1}"))
    return y1, y2


def build_majorana_c5_flavor_matrix(
    kernel: sp.Expr,
    y1: sp.MatrixBase,
    y2: sp.MatrixBase,
    *,
    fermion_mass_symbol: sp.Symbol | None = None,
    heavy_masses: list[sp.Expr] | tuple[sp.Expr, ...] | None = None,
) -> sp.Matrix:
    """Build the original symmetric lepton-flavor Majorana C5 matrix."""
    if y1.shape != y2.shape:
        raise ValueError("y1 and y2 must have the same shape.")
    n_lepton, n_heavy = y1.shape
    if heavy_masses is not None and len(heavy_masses) != n_heavy:
        raise ValueError("heavy_masses must contain one mass per heavy generation.")
    if heavy_masses is not None and fermion_mass_symbol is None:
        fermion_mass_symbol = sp.Symbol("MF")

    def entry(p: int, q: int) -> sp.Expr:
        value = sp.S.Zero
        for r in range(n_heavy):
            kernel_r = kernel
            if heavy_masses is not None:
                kernel_r = kernel_r.subs(fermion_mass_symbol, heavy_masses[r])
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


def build_general_majorana_c5_flavor_matrix(
    yukawas: Mapping[str, sp.MatrixBase],
    channel_kernels: Mapping[tuple[str, str], sp.Expr],
    *,
    fermion_mass_symbol: sp.Symbol | None = None,
    heavy_masses: list[sp.Expr] | tuple[sp.Expr, ...] | None = None,
) -> sp.Matrix:
    """Sum *supplied, matched* one-loop channels into a symmetric C5.

    Keys like ("y1", "y2Conj") identify external Yukawa vertices, NOT a
    prescription for which channels are physically allowed. Only pass a
    channel after deriving its full loop kernel including its SU(2) factor,
    masses, scalar interactions and relative normalization.

    Each kernel multiplies conjugate(Y_a[p,r])*conjugate(Y_b[q,r]) plus the
    external-lepton interchange, exactly as in the restricted implementation.
    The formula assumes Yukawa matrices are in the diagonal heavy-fermion mass
    basis. Off-diagonal fermion masses require a separate basis rotation.
    """
    if not channel_kernels:
        raise ValueError("channel_kernels is empty: no matched C5 calculation")
    if fermion_mass_symbol is None:
        fermion_mass_symbol = sp.Symbol("MF")
    shapes = {tuple(matrix.shape) for matrix in yukawas.values()}
    if len(shapes) != 1:
        raise ValueError("All Yukawa matrices must have matching dimensions")
    if not shapes:
        raise ValueError("No Yukawa matrices supplied")
    n_lepton, n_heavy = next(iter(shapes))
    if heavy_masses is not None and len(heavy_masses) != n_heavy:
        raise ValueError("Heavy mass list length does not match Yukawa columns")
    result = sp.zeros(n_lepton)
    for channel, kernel in channel_kernels.items():
        if not isinstance(channel, tuple) or len(channel) != 2:
            raise ValueError(f"Invalid matching channel {channel!r}")
        first, second = channel
        if first not in yukawas or second not in yukawas:
            raise ValueError(f"Missing Yukawa matrix for channel {channel!r}")
        y_first, y_second = yukawas[first], yukawas[second]
        for r in range(n_heavy):
            coefficient = sp.sympify(kernel)
            if heavy_masses is not None:
                coefficient = coefficient.subs(fermion_mass_symbol, heavy_masses[r])
            for p in range(n_lepton):
                for q in range(p, n_lepton):
                    term = coefficient * (
                        sp.conjugate(y_first[p, r]) * sp.conjugate(y_second[q, r])
                        + sp.conjugate(y_second[p, r]) * sp.conjugate(y_first[q, r])
                    )
                    result[p, q] += term
                    if q != p:
                        result[q, p] += term
    return sp.Matrix(result)


def build_flavor_c5_from_matchete(
    c5_path: Path,
    *,
    n_lepton: int = 3,
    n_heavy: int = 3,
    split_heavy_masses: bool = True,
) -> dict:
    """Read a restricted one-generation Matchete C5 and lift it to flavor."""
    c5_path = Path(c5_path)
    kappa_1g = parse_matchete_c5(c5_path.read_text(encoding="utf-8"))
    kernel = extract_one_generation_loop_kernel(kappa_1g)
    y1, y2 = symbolic_t3_yukawas(n_lepton=n_lepton, n_heavy=n_heavy)
    heavy_masses = ([sp.Symbol(f"MF{r + 1}") for r in range(n_heavy)]
                    if split_heavy_masses else None)
    c5 = build_majorana_c5_flavor_matrix(kernel, y1, y2,
                                         heavy_masses=heavy_masses)
    return {"kappa_1g": kappa_1g, "kernel": kernel, "y1": y1,
            "y2": y2, "heavy_masses": heavy_masses, "K": c5}


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
        matrix_path.write_text(sp.sstr(result["K"]) + "\n", encoding="utf-8")
        kernel_path.write_text(sp.sstr(result["kernel"]) + "\n", encoding="utf-8")
    return {
        "C5FlavorMatrixFile": matrix_path.relative_to(output_dir).as_posix()
        if debug_outputs else "",
        "C5LoopKernelFile": kernel_path.relative_to(output_dir).as_posix()
        if debug_outputs else "",
        "HeavyFlavorGenerations": result["y1"].cols,
        "LeptonFlavorGenerations": result["y1"].rows,
    }


extract_t3_loop_kernel = extract_one_generation_loop_kernel
build_flavor_c5_matrix = build_majorana_c5_flavor_matrix
match_c5_flavor_from_file = build_flavor_c5_from_matchete
