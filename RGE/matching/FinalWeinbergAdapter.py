"""Adapter from the hierarchical final-C5 JSON to symbolic flavor C5.

The final JSON already combines the renormalized hard matching contribution and
the direct intermediate-EFT running contribution. This module converts that
serialized matching result into the symbolic physical Majorana C5 matrix used
by downstream RGE and neutrino calculations.
"""

from __future__ import annotations

import json
from pathlib import Path

import sympy as sp

from RGE.matching.MatcheteC5Parsing import parse_matchete_c5
from RGE.matching.WeinbergFlavorMatching import (
    build_majorana_c5_flavor_matrix,
    extract_one_generation_loop_kernel,
    symbolic_t3_yukawas,
)


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
    """Extract A from a direct-running expression A*C12[p,q]."""

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


def load_hierarchical_majorana_c5(
    final_weinberg_path: Path,
    *,
    n_lepton: int = 3,
    n_heavy: int = 3,
    split_heavy_masses: bool = True,
) -> dict:
    """Build the symbolic physical Majorana C5 from the final matching JSON."""

    final_weinberg_path = Path(final_weinberg_path)
    payload = _load_payload(final_weinberg_path)

    hard = payload.get("renormalized_hard_threshold", {})
    hard_expression = hard.get("expression")

    if hard.get("status") != "Success" or not hard_expression:
        raise ValueError(
            "Renormalized hard threshold is unavailable."
        )

    hard_1g = parse_matchete_c5(hard_expression)
    hard_kernel = extract_one_generation_loop_kernel(hard_1g)

    direct = payload.get("direct_running_contribution", {})
    raw_direct = direct.get("raw_expression", "")

    if not direct.get("manifestly_symmetric_under_pq", False):
        raise ValueError(
            "Direct running contribution is not validated as p<->q symmetric."
        )

    direct_prefactor = _parse_direct_c12_prefactor(raw_direct)
    MF = sp.Symbol("MF")
    hbar = sp.Symbol("hbar")

    # C12[p,q] carries 1/(2 MF) times the symmetric Yukawa structure.
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

    c5_hard = build_majorana_c5_flavor_matrix(
        hard_kernel,
        y1,
        y2,
        fermion_mass_symbol=MF,
        heavy_masses=heavy_masses,
    )

    c5_running = build_majorana_c5_flavor_matrix(
        running_kernel,
        y1,
        y2,
        fermion_mass_symbol=MF,
        heavy_masses=heavy_masses,
    )

    c5 = sp.Matrix(
        n_lepton,
        n_lepton,
        lambda p, q: sp.simplify(c5_hard[p, q] + c5_running[p, q]),
    )

    symmetry_difference = c5 - c5.T
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
        "K_hard": c5_hard,
        "K_running": c5_running,
        "K": c5,
    }


def is_hierarchical_final_c5(path: Path) -> bool:
    """Return whether path is the authoritative hierarchical final-C5 JSON."""

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


# Historical names retained for callers that have not migrated yet.
load_final_weinberg_flavor_matrix = load_hierarchical_majorana_c5
is_final_weinberg_json = is_hierarchical_final_c5
