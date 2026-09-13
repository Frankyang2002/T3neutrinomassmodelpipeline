from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
from typing import Any

import sympy as sp
from sympy.printing.mathematica import mathematica_code


MATCHETE_COUPLING_ALIASES = {
    # Python/RGE names -> actual Matchete SM coupling symbols.
    "g2": "gL",
    "g3": "gs",
    "ye": "Ye",
    "yu": "Yu",
    "yd": "Yd",
}

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

BOUNDARY_RE = re.compile(
    r"conjugate\((?P<y1>[A-Za-z0-9_]+)\[p,r\]\).*?"
    r"conjugate\((?P<y2>[A-Za-z0-9_]+)\[q,r\]\)"
)

MIXED_BOUNDARY_RE = re.compile(
    r"conjugate\((?P<y1>[A-Za-z0-9_]+)\[p,r\]\).*?"
    r"conjugate\((?P<y2>[A-Za-z0-9_]+)\[q,r\]\).*?\+.*?"
    r"conjugate\((?P=y2)\[p,r\]\).*?"
    r"conjugate\((?P=y1)\[q,r\]\)"
)


def _load(path: Path) -> dict[str, Any]:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    if payload.get("status") != "Success":
        raise ValueError(
            f"Flavor transport is not successful: {payload.get('status')!r}"
        )
    if not payload.get("one_generation_reduction_matches"):
        raise ValueError("One-generation regression is not successful.")
    if not payload.get("equal_scale_running_vanishes"):
        raise ValueError("Equal-scale running check is not successful.")
    return payload


def _sympy_expr(raw: Any) -> sp.Expr:
    return sp.sympify(
        str(raw),
        locals={
            "conjugate": sp.conjugate,
            "sqrt": sp.sqrt,
            "log": sp.log,
            "I": sp.I,
            "pi": sp.pi,
        },
    )


def _to_wolfram_scalar(raw: Any) -> str:
    expr = sp.simplify(_sympy_expr(raw))
    text = mathematica_code(expr)

    for name in sorted(SAFE_SCALAR_COUPLINGS, key=len, reverse=True):
        matchete_name = MATCHETE_COUPLING_ALIASES.get(name, name)
        text = re.sub(
            rf"\b{re.escape(name)}\b",
            f"Coupling[{matchete_name}, {{}}, 0]",
            text,
        )

    # Matchete rejects ordinary, undefined symbols inside a Lagrangian.
    # The running logarithm is therefore represented during matching by a
    # registered real, dimensionless dummy coupling.  RunThresholdStage.wl
    # restores the physical Log[MS/MF] only after Match has finished.
    text = text.replace(
        "Log[MS/Coupling[MF, {}, 0]]",
        "Coupling[EFT1RunLog, {}, 0]",
    )

    return text


def _extract_boundary_definition(name: str, info: dict[str, Any]) -> str:
    text = str(info["definition_text"])
    mixed_match = MIXED_BOUNDARY_RE.search(text)
    match = mixed_match or BOUNDARY_RE.search(text)
    if not match:
        raise ValueError(f"Could not parse flavor kernel for {name}: {text}")
    y_a = match.group("y1")
    y_b = match.group("y2")
    one_gen = _sympy_expr(info["one_generation_reduction"])
    yuk = sp.conjugate(sp.Symbol(y_a)) * sp.conjugate(sp.Symbol(y_b))
    denom = _to_wolfram_scalar(sp.simplify(yuk / one_gen))
    if mixed_match:
        denom = f"2*({denom})"
        return (
            f'{name}[p_, q_] := Module[{{r = Unique["nfl$"]}},\n'
            f'  (\n'
            f'    Bar[Coupling[{y_a}, {{p, Index[r, NFlavor]}}, 0]] *\n'
            f'    Bar[Coupling[{y_b}, {{q, Index[r, NFlavor]}}, 0]] +\n'
            f'    Bar[Coupling[{y_b}, {{p, Index[r, NFlavor]}}, 0]] *\n'
            f'    Bar[Coupling[{y_a}, {{q, Index[r, NFlavor]}}, 0]]\n'
            f'  ) / ({denom})\n'
            f'];'
        )
    return (
        f'{name}[p_, q_] := Module[{{r = Unique["nfl$"]}},\n'
        f'  Bar[Coupling[{y_a}, {{p, Index[r, NFlavor]}}, 0]] *\n'
        f'  Bar[Coupling[{y_b}, {{q, Index[r, NFlavor]}}, 0]] /\n'
        f'  ({denom})\n'
        f'];'
    )


BARRED_TERM_COUPLING_RE = re.compile(
    r"Bar\[Coupling\[(?P<name>[A-Za-z0-9_]+),\s*"
    r"\{Index\[(?P<flavor>[^,\]]+),\s*Flavor\],\s*"
    r"Index\[(?P<nflavor>[^,\]]+),\s*NFlavor\]\},\s*0\]\]"
)


def _strip_tree_flavor_coefficient(
    term: str,
    *,
    raw_mass_denominator_factor: str,
    flavor_label: str,
) -> tuple[str, str, str]:
    """Divide a raw Matchete term by its canonical flavor coefficient.

    The raw term may contain arbitrary representation factors before the
    Yukawas and/or multiplying MF, e.g. -Sqrt[3]/2 or 1/Sqrt[3].  Those are
    operator-geometry factors and must remain in the skeleton.
    """
    outer = term.strip()
    matches = list(BARRED_TERM_COUPLING_RE.finditer(outer))
    if len(matches) != 2:
        raise ValueError(
            "Expected exactly two barred T3 Yukawas in PL operator term."
        )

    p_dummy = matches[0].group("flavor").strip()
    q_dummy = matches[1].group("flavor").strip()

    # Remove the MF denominator completely; its non-flavor prefactor is
    # restored below relative to the canonical C11/C12/C22 normalization.
    numerator, count = re.subn(
        r"/\((?:.*?)\*Coupling\[MF,\s*\{\},\s*0\]\)\s*$",
        "",
        outer,
    )
    if count == 0:
        numerator, count = re.subn(
            r"/Coupling\[MF,\s*\{\},\s*0\]\s*$",
            "",
            outer,
        )
    if count != 1:
        raise ValueError(
            "Could not strip the final MF denominator from PL operator term."
        )

    # Replace only the two flavor Yukawa factors, preserving all CG/numerical
    # prefactors in the raw operator geometry.
    for match in reversed(list(BARRED_TERM_COUPLING_RE.finditer(numerator))[:2]):
        numerator = numerator[:match.start()] + "1" + numerator[match.end():]

    canonical_mass_factor = 2 if flavor_label in {"C11", "C22"} else 1
    raw_factor = str(raw_mass_denominator_factor).strip()

    if raw_factor == str(canonical_mass_factor):
        skeleton = numerator
    else:
        skeleton = (
            f"({canonical_mass_factor})/({raw_factor})*({numerator})"
        )

    return skeleton, p_dummy, q_dummy


def _he_definition() -> str:
    # Matchete's SM model names the charged-lepton Yukawa coupling Ye
    # (capital Y), not the lowercase Python/RGE symbol ye.
    #
    # H_e[p,q] = (Y_e Y_e^dagger)[p,q]
    return (
        'EFT1He[p_, q_] := Module[{r = Unique["fl$"]},\n'
        '  Coupling[Ye, {p, Index[r, Flavor]}, 0] *\n'
        '  Bar[Coupling[Ye, {q, Index[r, Flavor]}, 0]]\n'
        '];'
    )


def _tensor_term_wl(record: dict[str, Any], p: str, q: str) -> str:
    source = record["source"]
    coeff = _to_wolfram_scalar(record["running_coefficient"])

    if record["kind"] == "flavor_blind":
        tensor = f"EFT1{source}[{p}, {q}]"
    elif record["kind"] == "He_two_leg_action":
        tensor = (
            'Module[{s = Unique["fl$"]}, '
            f'EFT1He[{p}, Index[s, Flavor]]*'
            f'EFT1{source}[Index[s, Flavor], {q}] + '
            f'EFT1{source}[{p}, Index[s, Flavor]]*'
            f'EFT1He[{q}, Index[s, Flavor]]]'
        )
    else:
        raise ValueError(f"Unknown running tensor kind: {record['kind']!r}")

    return f"({coeff})*({tensor})"


def _running_function_definition(target: str, correction: dict[str, Any]) -> str:
    pieces = [
        _tensor_term_wl(record, "p", "q")
        for record in correction.get("terms", [])
    ]
    body = "0" if not pieces else " +\n    ".join(pieces)
    return f"EFT1Delta{target}[p_, q_] :=\n  {body};"



def _validate_wolfram_delimiters(text: str) -> None:
    """Catch exporter-generated bracket/parenthesis corruption early.

    This is not a full Wolfram parser.  It is sufficient for the generated
    file because strings/comments are simple and the previous failure was an
    unmatched delimiter introduced while stripping the original coefficient.
    """
    pairs = {")": "(", "]": "[", "}": "{"}
    stack: list[tuple[str, int]] = []
    in_string = False
    escaped = False
    i = 0

    while i < len(text):
        ch = text[i]

        if in_string:
            if escaped:
                escaped = False
            elif ch == "\\\\":
                escaped = True
            elif ch == '"':
                in_string = False
            i += 1
            continue

        if ch == '"':
            in_string = True
            i += 1
            continue

        # Skip Mathematica comments (* ... *), including the generated header.
        if ch == "(" and i + 1 < len(text) and text[i + 1] == "*":
            end = text.find("*)", i + 2)
            if end == -1:
                raise ValueError("Unterminated Wolfram comment.")
            i = end + 2
            continue

        if ch in "([{":
            stack.append((ch, i))
        elif ch in ")]}":
            if not stack or stack[-1][0] != pairs[ch]:
                raise ValueError(
                    f"Unbalanced Wolfram delimiter {ch!r} at character {i}."
                )
            stack.pop()

        i += 1

    if in_string:
        raise ValueError("Unterminated Wolfram string.")
    if stack:
        opener, position = stack[-1]
        raise ValueError(
            f"Unclosed Wolfram delimiter {opener!r} at character {position}."
        )


def export_flavor_matchete(
    transport_path: Path,
    output_path: Path,
) -> dict[str, Any]:
    payload = _load(transport_path)

    boundaries = payload["boundary_tensors"]
    operator_directions = payload.get("operator_directions") or {}
    corrections = payload["running_corrections"]

    required_heavy = tuple(operator_directions)
    if not required_heavy:
        raise ValueError("No heavy operator directions were exported.")
    for name in required_heavy:
        if name not in corrections:
            raise ValueError(f"Missing running correction {name}.")
    if "Weinberg" not in corrections:
        raise ValueError("Missing Weinberg running correction.")

    boundary_names = tuple(boundaries)
    boundary_defs = [
        _extract_boundary_definition(f"EFT1{name}", boundaries[name])
        for name in boundary_names
    ]
    delta_defs = [
        _running_function_definition(name, corrections[name])
        for name in (*required_heavy, "Weinberg")
    ]

    heavy_terms = []
    for name in required_heavy:
        info = operator_directions[name]
        skeleton, p_dummy, q_dummy = _strip_tree_flavor_coefficient(
            info["original_term_input_form"],
            raw_mass_denominator_factor=info[
                "raw_mass_denominator_factor"
            ],
            flavor_label=info["flavor_label"],
        )
        p_index = f"Index[{p_dummy}, Flavor]"
        q_index = f"Index[{q_dummy}, Flavor]"
        heavy_terms.append(
            f"EFT1Delta{name}[{p_index}, {q_index}] * ({skeleton})"
        )

    heavy_body = " +\n  ".join(heavy_terms)

    boundary_symbol_list = ", ".join(f"EFT1{name}" for name in boundary_names)
    delta_symbol_list = ", ".join(
        f"EFT1Delta{name}" for name in (*required_heavy, "Weinberg")
    )

    wolfram = f'''(* Auto-generated by EFT1WilsonFlavorMatcheteExporter.py

   Full-flavor fixed-order EFT1 running insertion at the S1,S2 threshold.

   Heavy-field running terms are tree-matched at threshold 2.
   The Weinberg running coefficient is carried directly to the final EFT.

   IMPORTANT:
   - hbar multiplies the whole running correction.
   - The stage-1 one-loop matching boundary is separate and is NOT evolved.
*)

ClearAll[
  {boundary_symbol_list},
  EFT1He,
  {delta_symbol_list},
  EFT1WilsonRunningHeavyInsertion,
  EFT1WilsonRunningWeinbergCoefficient,
  EFT1WilsonRunningLogReplacement,
  EFT1WilsonFlavorRunningMetadata
];

(* Matchete-valid stand-in for Log[MS/MF] during threshold matching. *)
If[!KeyExistsQ[GetCouplings[], EFT1RunLog],
  DefineCoupling[EFT1RunLog, SelfConjugate -> True];
];

EFT1WilsonRunningLogReplacement =
  Coupling[EFT1RunLog, {{}}, 0] ->
    Log[MS/Coupling[MF, {{}}, 0]];

{chr(10).join(boundary_defs)}

{_he_definition()}

{chr(10).join(delta_defs)}

(* The Python flavor seed contains the holomorphic PL sector.  A physical
   Lagrangian must also contain its Hermitian conjugate.  PlusHc is Matchete's
   native wrapper and is accepted directly by Match. *)
EFT1WilsonRunningHeavyInsertion =
  hbar * PlusHc[
  {heavy_body}
  ];

EFT1WilsonRunningWeinbergCoefficient[p_, q_] :=
  EFT1DeltaWeinberg[p, q];

EFT1WilsonFlavorRunningMetadata = <|
  "MuHigh" -> "{payload["scales"]["mu_high"]}",
  "MuLow" -> "{payload["scales"]["mu_low"]}",
  "LogRatio" -> "{payload["scales"]["log_ratio"]}",
  "HeavyOperatorCount" -> {len(required_heavy)},
  "WeinbergOperatorCount" -> 1,
  "MixedC12FlavorSymmetrized" -> True,
  "OneGenerationRegression" -> True,
  "EqualScaleRunningVanishes" -> True
|>;
'''

    # Fail on malformed generated WL here rather than much later inside the
    # fresh Matchete threshold kernel.
    _validate_wolfram_delimiters(wolfram)

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(wolfram, encoding="utf-8")

    return {
        "status": "Success",
        "output": str(output_path),
        "heavy_operator_count": len(required_heavy),
        "weinberg_operator_count": 1,
        "full_flavor": True,
        "one_generation_reduction_matches": True,
        "equal_scale_running_vanishes": True,
        "log_ratio": payload["scales"]["log_ratio"],
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Export full-flavor EFT1 leading-log running as a Matchete-ready "
            "Wolfram insertion for threshold 2."
        )
    )
    parser.add_argument("flavor_transport_json", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    result = export_flavor_matchete(
        args.flavor_transport_json,
        args.output,
    )
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
