"""T3-A-m2 exact symbolic C12/C5 readiness audit, without modifying production.

Install at Numerical/analysis/CheckT3AC12Readiness.py
Run: python -m Numerical.analysis.CheckT3AC12Readiness
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from RGE.group_factors.validation.C12RepresentationGroupFactors import derive_from_saved_data
from RGE.group_factors.validation.C12ReverseMixingGroupFactors import derive_reverse_mixing

ROOT = Path(__file__).resolve().parents[2]
MODEL = "T3_dS1_1_dS2_3_dF_2_alpha_m2"
ALIAS = "T3_A_alpha_m2"
FILES = ("eft1_after_F_wilson_seed.json",
         "eft1_after_F_scalar_quartic_seed.json",
         "eft1_rgbeta_rge.json", "eft1_wilson_rge.json")

def locate_data():
    root = ROOT / "output/full/comparison/smallY_smallL" / MODEL
    candidates = sorted(root.glob(f"**/{ALIAS}/data"))
    hits = [d for d in candidates if all((d/n).is_file() for n in FILES)]
    if len(hits) != 1:
        raise FileNotFoundError(
            f"Expected one complete T3-A-m2 symbolic data folder under {root}; "
            f"found {len(hits)}: {hits}")
    return hits[0]

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=None)
    parser.add_argument("--out", type=Path, default=(
        ROOT / "output/studies/t3a_coupled_readiness"))
    args = parser.parse_args()
    data = args.data_dir or locate_data()
    forward = derive_from_saved_data(data)
    reverse = derive_reverse_mixing(data)
    output = {
        "model": MODEL,
        "kind": "exact_symbolic_C12_C5_group_factor_audit",
        "interpretation": ("Does not constitute an end-to-end coupled numerical "
                           "prediction; complete operator-basis closure remains to be validated."),
        "source_data_dir": str(data),
        "c12_diagonal_and_forward_seed_projection": forward,
        "c5_to_c12_reverse_projection": reverse,
    }
    args.out.mkdir(parents=True, exist_ok=True)
    path = args.out / "t3a_m2_symbolic_readiness.json"
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite previous audit: {path}")
    path.write_text(json.dumps(output, indent=2, default=str)+"\n",
                    encoding="utf-8")
    print("Saved:", path)
    print("Reverse status:", reverse.get("status"))
    print("Reverse factor:", reverse.get("factor"))
    print("Reverse residual count:", reverse.get("residual_count"))
    print("C12 sector closures:")
    for name, sector in forward.get("sectors", {}).items():
        print(" ", name, "coefficient=", sector.get("coefficient"),
              "residuals=", sector.get("residual_count", "inspect JSON"))

if __name__ == "__main__":
    main()
