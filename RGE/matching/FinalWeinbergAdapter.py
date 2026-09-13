from __future__ import annotations

"""
Adapter from final_weinberg_coefficient.json to the downstream symbolic
full-flavor C5 matrix.

The hierarchical matching pipeline now produces a physical Majorana-symmetric
coefficient directly.  This module consumes that object without reconstructing
flavor from the old one-generation c5_coefficient.txt file.

Legacy/common-threshold paths remain in FlavorMatchedC5.py.
"""

import json
import re
from pathlib import Path

import sympy as sp

from RGE.matching.FlavorMatchedC5 import (
    build_flavor_matched_c5,
    extract_t3_loop_kernel,
    symbolic_t3_yukawas,
)
from RGE.matching.MatchedEFTRGE import parse_matchete_c5


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
    the kernel passed to build_flavor_matched_c5 is A/MF.
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

    # build_flavor_matched_c5 contributes 1/2 * kernel * symmetric Yukawas.
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

    K_hard = build_flavor_matched_c5(
        hard_kernel,
        y1,
        y2,
        fermion_mass_symbol=MF,
        heavy_masses=heavy_masses,
    )

    K_running = build_flavor_matched_c5(
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
