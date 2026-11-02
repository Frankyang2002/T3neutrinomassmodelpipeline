"""Construct exactly 16 UNPHYSICAL schema/linear-algebra test boundaries.

No file generated here is a phenomenological T3 prediction.  Only the
representation *sizes* mimic the five T3 families; holomorphic quartic and
CG coefficients are synthetic, and full gauge fields are deliberately absent.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from RGE.running.intermediate.GeneralEFT1C5Mixing import MODELS, FORMAT, CERTIFICATES

DIMENSIONS = {"A": (1, 3), "B": (2, 2), "C": (2, 2),
              "D": (3, 1), "E": (3, 3)}


def fixture(model: str) -> dict:
    if model not in MODELS:
        raise ValueError("Unknown model")
    d1, d2 = DIMENSIONS[model[0]]
    ns = 4 + 2*d1 + 2*d2
    s1_r, s1_i = 5, 6
    s2_r, s2_i = 5 + 2*d1, 6 + 2*d1
    # Complex symmetric 3x3 matrix with all three families and CP phases.
    mat = [[complex(2.3e-11/(1+p+q), 2.1e-12*(p+q+1))
            for q in range(3)] for p in range(3)]
    comps_1 = [
        {"indices": [2*p+1, 2*q+1, s1_r, s2_r],
         "value": {"real": mat[p][q].real, "imag": mat[p][q].imag}}
        for p in range(3) for q in range(p, 3)
    ]
    comps_2 = [
        {"indices": [2*p+1, 2*q+1, s1_i, s2_i],
         "value": {"real": -0.31*mat[p][q].real, "imag": 0.2*mat[p][q].imag}}
        for p in range(3) for q in range(p, 3)
    ]
    quartic = [
        {"indices": [3, 3, s1_r, s2_r], "value": 0.09},
        {"indices": [4, 4, s1_r, s2_r], "value": -0.09},
        {"indices": [1, 1, s1_r, s2_r], "value": 0.09},
        {"indices": [2, 2, s1_r, s2_r], "value": -0.09},
        {"indices": [3, 3, s1_i, s2_i], "value": -0.04},
        {"indices": [4, 4, s1_i, s2_i], "value": 0.04},
    ]
    return {
        "format": FORMAT, "model": model, "authority": "synthetic",
        "synthetic_notice": "NOT PHYSICAL: not gauge invariant or Matchete normalized",
        "scheme": "MSbar", "mu_F_gev": 100000.0, "mu_S_gev": 2000.0,
        "n_fermions": 6, "n_real_scalars": ns,
        "lepton_components": [[1, 2], [3, 4], [5, 6]],
        "higgs_real_components": [1, 2, 3, 4],
        "scalar_field_blocks": [
            {"name": "H", "role": "Higgs", "indices": [1, 2, 3, 4]},
            {"name": "toy_S1", "role": "heavy", "indices": list(range(5, 5+2*d1))},
            {"name": "toy_S2", "role": "heavy", "indices": list(range(5+2*d1, ns+1))},
        ],
        "master_lagrangian_convention": "C_ijab psi_i psi_j phi_a phi_b + hc",
        "weinberg_lagrangian_convention": "(1/2) C5_ij (L_i.H)(L_j.H) + hc",
        "certificates": {name: False for name in CERTIFICATES},
        "yukawa_components": [], "gauge_sectors": [],
        "quartic_components": quartic,
        "operators": [
            {"source_id": "synthetic_LLSS_1", "components": comps_1},
            {"source_id": "synthetic_LLSS_2", "components": comps_2},
        ],
    }


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output-root", type=Path, required=True)
    args = p.parse_args(argv)
    for model in MODELS:
        tag = model.replace("+", "p").replace("-", "m")
        dest = args.output_root / tag / "boundary.json"
        dest.parent.mkdir(parents=True, exist_ok=True)
        with dest.open("x", encoding="utf-8") as f:
            json.dump(fixture(model), f, indent=2, allow_nan=False)
            f.write("\n")
    print(f"Generated {len(MODELS)} synthetic-only fixtures in {args.output_root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
