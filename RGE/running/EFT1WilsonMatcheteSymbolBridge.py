from __future__ import annotations

"""Build the symbol bridge needed before EFT1 running can enter Matchete.

The Python EFT1 Wilson RGE is presently a one-generation component calculation.
Its projected running multipliers therefore contain plain symbols such as
``y1``, ``y2`` and ``ye``, whereas the Matchete model defines those Yukawas as
indexed ``Coupling`` objects.

Unindexed parameters (gauge couplings, scalar couplings and masses) can be
mapped back mechanically. Indexed Yukawas cannot: their flavor contractions
must first be lifted from the one-generation beta coefficient to the full
Matchete flavor tensor.
"""

import argparse
import json
from pathlib import Path
from typing import Any


SAFE_MATCHETE_COUPLINGS = {
    "MF", "gY", "g2", "g3",
    "lambdaH", "lambdaS1", "lambdaS2",
    "lambdaH1", "lambdaH2", "lambda12", "lambdaT3",
    "lambdaS1Inv1", "lambdaS1Inv2",
    "lambdaS2Inv1", "lambdaS2Inv2",
    "lambdaH1Inv1", "lambdaH1Inv2",
    "lambdaH2Inv1", "lambdaH2Inv2",
    "lambda12Inv1", "lambda12Inv2",
}

FLAVOR_INDEXED_COUPLINGS = {"yu", "yd", "ye", "y1", "y2"}
SCALE_SYMBOLS = {"MS", "MS1", "MS2", "muHigh", "muLow"}


def _load_symbols(path: Path) -> list[str]:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))

    if "referenced_plain_symbols" in payload:
        return sorted(set(payload["referenced_plain_symbols"]))

    if "operator_coefficients" not in payload:
        raise ValueError(
            "Expected exporter metadata JSON or operator-projection JSON."
        )

    import sympy as sp

    symbols: set[str] = set()
    for entry in payload["operator_coefficients"]:
        expr = sp.sympify(
            str(entry.get("running_multiplier", "0")),
            locals={
                "conjugate": sp.conjugate,
                "sqrt": sp.sqrt,
                "log": sp.log,
                "I": sp.I,
                "pi": sp.pi,
            },
        )
        symbols.update(str(symbol) for symbol in expr.free_symbols)

    return sorted(symbols)


def classify_symbols(symbols: list[str]) -> dict[str, list[str]]:
    safe: list[str] = []
    flavor: list[str] = []
    scales: list[str] = []
    unknown: list[str] = []

    for symbol in sorted(set(symbols)):
        if symbol in SAFE_MATCHETE_COUPLINGS:
            safe.append(symbol)
        elif symbol in FLAVOR_INDEXED_COUPLINGS:
            flavor.append(symbol)
        elif symbol in SCALE_SYMBOLS:
            scales.append(symbol)
        else:
            unknown.append(symbol)

    return {
        "safe_scalar_couplings": safe,
        "flavor_indexed_couplings": flavor,
        "scale_symbols": scales,
        "unknown_symbols": unknown,
    }


def write_wolfram_bridge(
    classification: dict[str, list[str]],
    output_path: Path,
) -> None:
    safe = classification["safe_scalar_couplings"]
    flavor = classification["flavor_indexed_couplings"]
    scales = classification["scale_symbols"]
    unknown = classification["unknown_symbols"]

    rules = ",\n  ".join(
        f"{name} -> Coupling[{name}, {{}}, 0]"
        for name in safe
    ) or "Nothing"

    flavor_list = ", ".join(f'"{name}"' for name in flavor)
    scale_list = ", ".join(f'"{name}"' for name in scales)
    unknown_list = ", ".join(f'"{name}"' for name in unknown)

    wolfram = f'''(* Auto-generated safe symbol bridge for EFT1 Wilson running.

   Only unindexed Matchete couplings are converted here.
   Flavor-indexed Yukawas are intentionally NOT replaced.
*)

ClearAll[
  EFT1WilsonRunningSafeParameterRules,
  EFT1WilsonRunningFlavorSymbols,
  EFT1WilsonRunningScaleSymbols,
  EFT1WilsonRunningUnknownSymbols,
  EFT1WilsonRunningSymbolBridgeReadyQ
];

EFT1WilsonRunningSafeParameterRules = {{
  {rules}
}};

EFT1WilsonRunningFlavorSymbols = {{{flavor_list}}};
EFT1WilsonRunningScaleSymbols = {{{scale_list}}};
EFT1WilsonRunningUnknownSymbols = {{{unknown_list}}};

EFT1WilsonRunningSymbolBridgeReadyQ =
  EFT1WilsonRunningFlavorSymbols === {{}} &&
  EFT1WilsonRunningUnknownSymbols === {{}};
'''

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(wolfram, encoding="utf-8")


def build_symbol_bridge(
    metadata_or_projection_path: Path,
    output_path: Path | None = None,
) -> dict[str, Any]:
    symbols = _load_symbols(metadata_or_projection_path)
    classification = classify_symbols(symbols)

    if output_path is not None:
        write_wolfram_bridge(classification, Path(output_path))

    flavor = classification["flavor_indexed_couplings"]
    unknown = classification["unknown_symbols"]

    if unknown:
        status = "UnknownSymbols"
    elif flavor:
        status = "NeedsFlavorLift"
    else:
        status = "Ready"

    return {
        "status": status,
        "ready_for_matchete_stage2": status == "Ready",
        **classification,
        "safe_rule_count": len(classification["safe_scalar_couplings"]),
        "flavor_lift_required": bool(flavor),
        "output": str(output_path) if output_path is not None else None,
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Classify Python-side EFT1 running symbols and emit safe "
            "Matchete replacement rules."
        )
    )
    parser.add_argument("metadata_or_projection_json", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    result = build_symbol_bridge(
        args.metadata_or_projection_json,
        output_path=args.output,
    )
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
