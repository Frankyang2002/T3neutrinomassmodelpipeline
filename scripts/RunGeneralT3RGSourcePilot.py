"""Evaluate RGBeta one-loop source terms as an illustrative leading-log step."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from Numerical.experimental.GeneralT3RGSourcePilot import one_step_source


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--proof", type=Path, default=Path("output/general_rgbeta_source_proof.json"))
    p.add_argument("--model", choices=("B_m1", "C_m1"), required=True)
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    ns = p.parse_args()
    try:
        config = json.loads(ns.input.read_text(encoding="utf-8-sig"))
        result = one_step_source(ns.proof, ns.model, config)
        ns.output.parent.mkdir(parents=True, exist_ok=True)
        ns.output.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
        for name, entry in result["results"].items():
            print(ns.model, name, "delta_lambda:", entry["delta_lambda_leading_log"])
        print("Report:", ns.output)
        return 0
    except (ValueError, KeyError, OSError, TypeError) as exc:
        print("BLOCKED:", exc, file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
