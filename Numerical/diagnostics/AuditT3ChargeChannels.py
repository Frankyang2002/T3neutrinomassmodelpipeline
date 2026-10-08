"""Charge-conservation-only audit of T3 neutrino-mass internal channels.

The UV interactions are bar(L) F^c S1, bar(L) F S2^dagger,
and H H S1 S2^dagger. For neutral external legs all three
internal particle charges must agree: QF=QS1=QS2.

These charge intersections are NECESSARY conditions, not a verification
of nonzero CG coefficients, equal loop weights, or multiplicity factors.
"""
from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path
from typing import Any


def _parse_charge(value: str) -> Fraction:
    return Fraction(str(value))


def charge_channels(case: dict[str, Any]) -> dict[str, Any]:
    fields = case["fields"]
    charges = {
        name: {_parse_charge(x["Q"]) for x in fields[name]["weights_and_charges"]}
        for name in ("F", "S1", "S2")
    }
    common = sorted(charges["F"] & charges["S1"] & charges["S2"])
    return {
        "model": case["model"],
        "internal_charge_channels": [str(c) for c in common],
        "charge_allowed_channel_count": len(common),
        "cg_contraction_verified": False,
        "equal_channel_weights_verified": False,
        "verified_loop_multiplicity": None,
        "compensation_factor": None,
    }


def audit(payload: dict[str, Any]) -> dict[str, Any]:
    cases = payload.get("models", [])
    if not cases:
        raise ValueError("Charge-spectrum input contains no models")
    rows = [charge_channels(case) for case in cases]
    return {
        "status": "charge_allowed_channels_only",
        "count": len(rows),
        "selection_rule": "Q(F)=Q(S1)=Q(S2) for external neutrinos and neutral Higgs",
        "models": rows,
        "warning": "Charge matching is necessary, not sufficient. Do not set lambdaT3 compensation factors until transformed CG weights and propagator contractions are verified.",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path,
                        default=Path("Reports/output/full/multiplicity/charge_spectrum_audit.json"))
    parser.add_argument("--output", type=Path,
                        default=Path("Reports/output/full/multiplicity/charge_channels.json"))
    args = parser.parse_args()
    result = audit(json.loads(args.input.read_text(encoding="utf-8-sig")))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"Models: {result['count']}")
    for row in result["models"]:
        m = row["model"]
        print(f"  ({m['dS1']},{m['dS2']},{m['dF']}), alpha={m['alpha']}: "
              f"charges={row['internal_charge_channels']}; "
              f"charge-allowed channels={row['charge_allowed_channel_count']}")
    print(f"JSON: {args.output}")
    print("Physical loop multiplicity: not yet verified")


if __name__ == "__main__":
    main()
