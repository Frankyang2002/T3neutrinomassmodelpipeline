from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
from typing import Any

import sympy as sp
from sympy.printing.mathematica import mathematica_code

SAFE_SCALAR_COUPLINGS = {
    "MF", "gY", "g2", "g3",
    "lambdaH", "lambdaS1", "lambdaS2",
    "lambdaH1", "lambdaH2", "lambda12", "lambdaT3",
    "lambdaS1Inv1", "lambdaS1Inv2",
    "lambdaS2Inv1", "lambdaS2Inv2",
    "lambdaH1Inv1", "lambdaH1Inv2",
    "lambdaH2Inv1", "lambdaH2Inv2",
    "lambda12Inv1", "lambda12Inv2",
}
ALIASES = {"g2": "gL", "g3": "gs"}


def _expr(raw: Any) -> sp.Expr:
    return sp.sympify(str(raw), locals={
        "conjugate": sp.conjugate,
        "sqrt": sp.sqrt,
        "log": sp.log,
        "I": sp.I,
        "pi": sp.pi,
    })


def _to_wolfram_scalar(raw: Any) -> str:
    text = mathematica_code(sp.simplify(_expr(raw)))
    for name in sorted(SAFE_SCALAR_COUPLINGS, key=len, reverse=True):
        target = ALIASES.get(name, name)
        text = re.sub(
            rf"\b{re.escape(name)}\b",
            f"Coupling[{target}, {{}}, 0]",
            text,
        )
    text = text.replace(
        "Log[MS/Coupling[MF, {}, 0]]",
        "Coupling[EFT1RunLog, {}, 0]",
    )
    return text


def _c12_definition(info: dict[str, Any]) -> str:
    one_gen = _expr(info["one_generation_reduction"])
    y1 = sp.Symbol("y1")
    y2 = sp.Symbol("y2")
    mass_factor = sp.simplify(sp.conjugate(y1) * sp.conjugate(y2) / one_gen)
    denom = _to_wolfram_scalar(2 * mass_factor)
    return (
        'EFT1C12[p_, q_] := Module[{r = Unique["nfl$"]},\n'
        '  (\n'
        '    Bar[Coupling[y1, {p, Index[r, NFlavor]}, 0]] *\n'
        '    Bar[Coupling[y2, {q, Index[r, NFlavor]}, 0]] +\n'
        '    Bar[Coupling[y2, {p, Index[r, NFlavor]}, 0]] *\n'
        '    Bar[Coupling[y1, {q, Index[r, NFlavor]}, 0]]\n'
        f'  ) / ({denom})\n'
        '];'
    )


def export_direct_weinberg_matchete(
    transport_path: Path,
    output_path: Path,
) -> dict[str, Any]:
    payload = json.loads(Path(transport_path).read_text(encoding="utf-8"))
    if payload.get("status") != "Success":
        raise ValueError("Direct Weinberg transport is not successful.")
    if not payload.get("one_generation_reduction_matches"):
        raise ValueError("One-generation regression failed.")
    if not payload.get("equal_scale_running_vanishes"):
        raise ValueError("Equal-scale running check failed.")

    boundary = payload["boundary_tensors"]["C12"]
    correction = payload["running_corrections"]["Weinberg"]
    terms = correction.get("terms") or []
    if len(terms) != 1 or terms[0].get("source") != "C12":
        raise ValueError("Expected one direct C12 -> Weinberg running term.")

    running_coeff = _to_wolfram_scalar(terms[0]["running_coefficient"])
    c12_def = _c12_definition(boundary)

    wolfram = f'''(* Auto-generated direct one-loop Weinberg transport.

   The heavy LLSS self-running is O(hbar).  Feeding that correction through
   the scalar loop would be O(hbar^2), so the authoritative fixed-one-loop
   path sets the heavy running insertion to zero and carries only the direct
   LLSS -> Weinberg mixing.
*)

ClearAll[
  EFT1C12,
  EFT1DeltaWeinberg,
  EFT1WilsonRunningHeavyInsertion,
  EFT1WilsonRunningWeinbergCoefficient,
  EFT1WilsonRunningLogReplacement,
  EFT1WilsonFlavorRunningMetadata
];

If[!KeyExistsQ[GetCouplings[], EFT1RunLog],
  DefineCoupling[EFT1RunLog, SelfConjugate -> True];
];

EFT1WilsonRunningLogReplacement =
  Coupling[EFT1RunLog, {{}}, 0] ->
    Log[MS/Coupling[MF, {{}}, 0]];

{c12_def}

EFT1DeltaWeinberg[p_, q_] :=
  ({running_coeff}) * EFT1C12[p, q];

EFT1WilsonRunningHeavyInsertion = 0;

EFT1WilsonRunningWeinbergCoefficient[p_, q_] :=
  EFT1DeltaWeinberg[p, q];

EFT1WilsonFlavorRunningMetadata = <|
  "MuHigh" -> "{payload['scales']['mu_high']}",
  "MuLow" -> "{payload['scales']['mu_low']}",
  "LogRatio" -> "{payload['scales']['log_ratio']}",
  "HeavyOperatorCount" -> 0,
  "WeinbergOperatorCount" -> 1,
  "DirectWeinbergOnlyAtOneLoop" -> True,
  "MixedC12FlavorSymmetrized" -> True,
  "OneGenerationRegression" -> True,
  "EqualScaleRunningVanishes" -> True
|>;
'''

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(wolfram, encoding="utf-8")

    return {
        "status": "Success",
        "output": str(output_path),
        "heavy_operator_count": 0,
        "weinberg_operator_count": 1,
        "full_flavor": True,
        "direct_weinberg_only_at_one_loop": True,
        "one_generation_reduction_matches": True,
        "equal_scale_running_vanishes": True,
        "log_ratio": payload["scales"]["log_ratio"],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("transport_json", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    result = export_direct_weinberg_matchete(args.transport_json, args.output)
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
