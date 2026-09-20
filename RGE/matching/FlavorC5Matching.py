from __future__ import annotations

"""
Lift a flavor-suppressed T3 Matchete coefficient into a full lepton-flavor
Weinberg matrix.

The Matchete matching currently uses one symbolic lepton Yukawa at each T3
vertex, y1 and y2.  Gauge and SU(2) contractions are independent of lepton
generation, so the one-generation result can be factorized as

    kappa_1g = F_loop * conjugate(y1) * conjugate(y2).

For three lepton flavors and n heavy-fermion generations, the symmetric
Weinberg matrix is then

    K_pq =
        sum_r F_loop(MF -> MF_r)
        [ conjugate(y1_pr) conjugate(y2_qr)
        + conjugate(y2_pr) conjugate(y1_qr) ].

Matchete supplies one ordered flavor coefficient A_pq. In the project convention
L_EFT contains (1/2) C5_pq O_pq + h.c., so the physical symmetric coefficient is
C5_pq = A_pq + A_qp. The one-generation physical C5 is therefore twice the
ordered Matchete coefficient.

Assumptions
-----------
- The heavy-fermion mass basis is diagonal.
- S1 and S2 masses are common across heavy-fermion generations.
- Flavor enters only through the two T3 Yukawa vertices.
"""

import json
from pathlib import Path

import sympy as sp

from RGE.matching.MatchedWeinbergRGE import parse_matchete_c5


def extract_t3_loop_kernel(
    kappa_1g: sp.Expr,
    *,
    y1_name: str = "y1",
    y2_name: str = "y2",
) -> sp.Expr:
    """Remove the one-generation T3 Yukawa product from the matched C5."""

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


def build_flavor_c5_matrix(
    kernel: sp.Expr,
    y1: sp.MatrixBase,
    y2: sp.MatrixBase,
    *,
    fermion_mass_symbol: sp.Symbol | None = None,
    heavy_masses: list[sp.Expr] | tuple[sp.Expr, ...] | None = None,
) -> sp.Matrix:
    """
    Build the symmetric lepton-flavor C5 matrix from the scalar loop kernel.

    If heavy_masses is supplied, MF in the kernel is replaced by MF_r for each
    heavy generation before summing.
    """

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

    K = sp.MutableDenseMatrix.zeros(n_lepton, n_lepton)

    for p in range(n_lepton):
        for q in range(p, n_lepton):
            value = entry(p, q)
            K[p, q] = value
            K[q, p] = value

    return sp.Matrix(K)


def match_c5_flavor_from_file(
    c5_path: Path,
    *,
    n_lepton: int = 3,
    n_heavy: int = 3,
    split_heavy_masses: bool = True,
) -> dict:
    """Read Matchete C5 and construct its full symbolic flavor matrix."""

    c5_path = Path(c5_path)

    kappa_1g = parse_matchete_c5(
        c5_path.read_text(encoding="utf-8")
    )

    kernel = extract_t3_loop_kernel(kappa_1g)

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

    K = build_flavor_c5_matrix(
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
        "K": K,
    }


def write_flavor_matching_outputs(
    output_dir: Path,
    result: dict,
    *,
    debug_outputs: bool = False,
) -> dict:
    """Write the full matched C5 matrix and loop kernel."""

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


# ---------------------------------------------------------------------------
# Final hierarchical Weinberg-coefficient adapter
# Consolidated from the former standalone adapter module.
# ---------------------------------------------------------------------------

def _load_payload(path: Path) -> dict:
    path = Path(path)
    payload = json.loads(path.read_text(encoding="utf-8"))

    if payload.get("status") != "Success":
        raise ValueError(
            "final_weinberg_coefficient.json is not successful."
        )

    combined = payload.get("combined", {})
    if not combined.get("ready_for_physical_majorana_numerics", False):
        raise ValueError(
            "Final Weinberg coefficient is not marked ready for physical "
            "Majorana numerics."
        )

    return payload


def _parse_direct_c12_prefactor(raw_expression: str) -> sp.Expr:
    """
    Parse A in

        A * C12[p,q]

    from the direct-running expression.

    Since
        C12[p,q] = 1/(2 MF) sum_r (y1* y2* + y2* y1*),
    the kernel passed to build_flavor_c5_matrix is A/MF.
    """
    marker = "C12[p,q]"

    if marker not in raw_expression:
        raise ValueError(
            "Direct running expression does not contain C12[p,q]."
        )

    text = raw_expression.replace(marker, "1")
    text = text.replace("^", "**")

    prefactor = sp.sympify(
        text,
        locals={
            "log": sp.log,
            "sqrt": sp.sqrt,
        },
    )

    return sp.simplify(prefactor)


def load_final_weinberg_flavor_matrix(
    final_weinberg_path: Path,
    *,
    n_lepton: int = 3,
    n_heavy: int = 3,
    split_heavy_masses: bool = True,
) -> dict:
    """
    Build the symbolic physical Majorana C5 matrix from the final JSON.

    The hard part is read from the renormalized ordered Matchete expression.
    Its scalar one-generation kernel is extracted, then the already validated
    symmetric flavor lift is rebuilt explicitly.

    The direct-running part is read from its C12 coefficient.  The stored
    running object is order-one, so the physical correction receives one
    explicit factor of hbar here.
    """
    final_weinberg_path = Path(final_weinberg_path)
    payload = _load_payload(final_weinberg_path)

    hard = payload.get("renormalized_hard_threshold", {})
    hard_expression = hard.get("expression")

    if hard.get("status") != "Success" or not hard_expression:
        raise ValueError(
            "Renormalized hard threshold is unavailable."
        )

    # Matchete's ordered hard expression reduces to the same scalar kernel
    # used by the legacy flavor constructor.
    hard_1g = parse_matchete_c5(hard_expression)
    hard_kernel = extract_t3_loop_kernel(hard_1g)

    direct = payload.get("direct_running_contribution", {})
    raw_direct = direct.get("raw_expression", "")

    if not direct.get("manifestly_symmetric_under_pq", False):
        raise ValueError(
            "Direct running contribution is not validated as p<->q symmetric."
        )

    direct_prefactor = _parse_direct_c12_prefactor(raw_direct)
    MF = sp.Symbol("MF")
    hbar = sp.Symbol("hbar")

    # build_flavor_c5_matrix contributes 1/2 * kernel * symmetric Yukawas.
    # C12 itself contributes 1/(2 MF), hence kernel = prefactor/MF.
    running_kernel = sp.simplify(hbar * direct_prefactor / MF)

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

    K_hard = build_flavor_c5_matrix(
        hard_kernel,
        y1,
        y2,
        fermion_mass_symbol=MF,
        heavy_masses=heavy_masses,
    )

    K_running = build_flavor_c5_matrix(
        running_kernel,
        y1,
        y2,
        fermion_mass_symbol=MF,
        heavy_masses=heavy_masses,
    )

    K = sp.Matrix(
        n_lepton,
        n_lepton,
        lambda p, q: sp.simplify(K_hard[p, q] + K_running[p, q]),
    )

    symmetry_difference = K - K.T
    if any(sp.simplify(entry) != 0 for entry in symmetry_difference):
        raise RuntimeError(
            "Final Weinberg adapter produced a non-symmetric C5 matrix."
        )

    return {
        "source_kind": "final_weinberg_json",
        "payload": payload,
        "hard_kernel": hard_kernel,
        "running_kernel": running_kernel,
        "y1": y1,
        "y2": y2,
        "heavy_masses": heavy_masses,
        "K_hard": K_hard,
        "K_running": K_running,
        "K": K,
    }


def is_final_weinberg_json(path: Path) -> bool:
    path = Path(path)

    if path.suffix.lower() != ".json":
        return False

    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return False

    return (
        payload.get("scheme") == "fixed_order_one_loop_MSbar_hierarchical"
        and "physical_majorana_hard_threshold" in payload
        and "direct_running_contribution" in payload
    )
