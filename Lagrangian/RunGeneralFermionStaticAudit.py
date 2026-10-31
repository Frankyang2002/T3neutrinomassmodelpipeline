#!/usr/bin/env python3
"""Independent, no-Matchete catalogue/Yukawa hypercharge audit for ordinary T3.

Reads no RGE/matching/numerical files. Produces the 16-model candidate table;
this cannot prove Matchete tensor normalization or physical mass completeness.
"""
from __future__ import annotations

import argparse
import csv
import json
from fractions import Fraction
from pathlib import Path

# Class dimensions are copied from T3ModelCatalog.wl at GitHub SHA in report.
CLASS_DIMS = {
    "A": (1, 3, 2), "B": (2, 2, 1), "C": (2, 2, 3),
    "D": (3, 1, 2), "E": (3, 3, 2),
}
# The 16 explicitly registered alpha points, from the scalar trial manifest.
ALPHAS = {
    "A": (-4, -2, 0), "B": (-3, -1, 1),
    "C": (-3, -1, 1), "D": (-2, 0, 2),
    "E": (-4, -2, 0, 2),
}
EXTRA_SCALAR_QUARTICS = {
    "A_-4": 0, "A_-2": 5, "A_0": 5,
    "B_-3": 1, "B_-1": 7, "B_1": 1,
    "C_-3": 1, "C_-1": 7, "C_1": 1,
    "D_-2": 5, "D_0": 5, "D_2": 0,
    "E_-4": 0, "E_-2": 6, "E_0": 6, "E_2": 0,
}
EXTRA_SCALAR_BILINEARS = {
    "A_-2", "A_0", "B_-1", "C_-1",
    "D_-2", "D_0", "E_-2", "E_0",
}


def su2_yukawa_allowed(d_f: int, d_s: int) -> bool:
    # SU(2): 2 x d_f = (d_f - 1) + (d_f + 1), each once.
    return abs(d_f - d_s) == 1


def get_yukawas(d1: int, d2: int, df: int, alpha: int):
    y1 = Fraction(alpha, 2)
    y2 = Fraction(alpha + 2, 2)
    yf = Fraction(alpha + 1, 2)
    # 4-component shorthand: Lbar F^(c optional) S^(dagger optional).
    raw = []
    for species, d_s, ys in (("S1", d1, y1), ("S2", d2, y2)):
        if not su2_yukawa_allowed(df, d_s):
            continue
        for fermion, fy in (("F", yf), ("Fc", -yf)):
            for scalar_conj, sy in ((False, ys), (True, -ys)):
                # Hypercharge of Lbar is +1/2.
                if Fraction(1, 2) + fy + sy == 0:
                    raw.append((fermion, species, scalar_conj))
    # If F has Majorana reality (Y=0 and odd-dimensional irrep), F=F^c.
    majorana = yf == 0 and df % 2 == 1
    normalized = sorted({
        ("F" if majorana else f, s, dagger) for (f, s, dagger) in raw
    })
    return normalized, majorana


def baseline_real_scalar_quartics(d1: int, d2: int) -> int:
    # Identical complex scalars: dim Inv[Sym^2(R*) x Sym^2(R)]
    # = 1 for singlet/doublet; 2 for triplet.
    selfs = {1: 1, 2: 1, 3: 2}
    # H*H x S*S channels: singlet and, for non-singlet S, triplet.
    portals = {1: 1, 2: 2, 3: 2}
    # S1*S1 x S2*S2: common intermediate isospins.
    cross = len(set(range(d1)).intersection(range(d2)))
    return selfs[d1] + selfs[d2] + portals[d1] + portals[d2] + cross


def model_rows():
    rows = []
    for cls, alphas in ALPHAS.items():
        d1, d2, df = CLASS_DIMS[cls]
        for alpha in alphas:
            name = f"{cls}_{alpha}"
            y1, y2, yf = (Fraction(alpha, 2), Fraction(alpha + 2, 2),
                          Fraction(alpha + 1, 2))
            invariant_yukawas, majorana = get_yukawas(d1, d2, df, alpha)
            standard = {
                ("F" if majorana else "Fc", "S1", False),
                ("F", "S2", True),
            }
            assert standard.issubset(set(invariant_yukawas)), name
            extras = sorted(set(invariant_yukawas) - standard)
            # The six special models gain one scalar-conjugate Yukawa.
            assert len(extras) == int(alpha in (0, -2)), (name, extras)
            assert not (y1 == 0 and y2 == 0)
            quartic_extra = EXTRA_SCALAR_QUARTICS[name]
            rows.append({
                "model": name,
                "class": cls,
                "alpha": alpha,
                "su2_s1": d1,
                "su2_s2": d2,
                "su2_f": df,
                "y_s1": str(y1),
                "y_s2": str(y2),
                "y_f": str(yf),
                "z2_sm": "+1",
                "z2_s1_s2_f": "-1",
                "fermion_type": "Majorana" if majorana else "Dirac",
                "yukawa_matrices_independent": len(invariant_yukawas),
                "yukawa_couplings_complex_entries": 9 * len(invariant_yukawas),
                "extra_yukawa": ";".join(
                    f"Lbar {f} {s}{'dag' if bar else ''}" for f, s, bar in extras
                ),
                "extra_yukawa_matrices": len(extras),
                "baseline_scalar_quartics_real": baseline_real_scalar_quartics(d1, d2),
                "baseline_scalar_quartics_complex": 1,
                "additional_scalar_quartics_complex": quartic_extra,
                "holomorphic_scalar_bilinears_complex": int(name in EXTRA_SCALAR_BILINEARS),
                "diagonal_scalar_mass_squared_real": 2,
                "heavy_fermion_generations": 3,
                "heavy_fermion_mass_parameters_expected_after_diagonalization": 3,
                "heavy_fermion_mass_parameter_symbols_in_current_fields": 1,
                "live_matchete_fermion_extension": "NOT RUN",
                "full_uv_symmetry_certified": False,
            })
    assert len(rows) == 16
    assert sum(r["extra_yukawa_matrices"] for r in rows) == 6
    assert sum(r["additional_scalar_quartics_complex"] for r in rows) == 50
    assert sum(r["holomorphic_scalar_bilinears_complex"] for r in rows) == 8
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", type=Path, required=True)
    args = ap.parse_args()
    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    rows = model_rows()
    with (outdir / "t3_uv_completeness_matrix.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    with (outdir / "t3_uv_completeness_matrix.json").open("w", encoding="utf-8") as f:
        json.dump({
            "kind": "STATIC_SOURCE_AND_SYMMETRY_PREDICTION_NOT_PHYSICAL_CERTIFICATE",
            "github_main_sha_observed": "305035e9df8b90e34fc74c7a3f57453ea6c6aae3",
            "new_fermion_z2": "odd",
            "sm_z2": "even",
            "hypercharge_normalization": "Q=T3+Y",
            "convention": "Lbar F(c) S(dag)",
            "counts": {"models": 16, "extra_yukawa_matrices": 6,
                       "extra_scalar_quartics_complex": 50,
                       "extra_scalar_bilinears_complex": 8},
            "notes": [
                "A successful CheckLagrangian is not operator-basis proof.",
                "SU(2) tensor normalization and all Matchete checks are unrun here.",
                "Z2 parity is enforced in enumeration, not registered as a Matchete charge.",
                "Mass MF is a single symbol for 3 flavours; general flavour masses unverified.",
                "Real/complex counts describe couplings, not real dof; SM-only terms excluded.",
            ],
            "models": rows,
        }, f, indent=2)
    print("16-model static audit passed; 6 extra Yukawa matrices, 50 extra scalar quartics, 8 bilinears")
    print(outdir / "t3_uv_completeness_matrix.json")
    print(outdir / "t3_uv_completeness_matrix.csv")


if __name__ == "__main__":
    main()
