'''Get EFT1 RGE components into full flavour running 
We do not work with real components now, we use flavour tensors
We can get a flavour blind factor for beta to wilson coefficient
In here we have
C_{12,pq} for 12 with our yukawa flavour structure and
p,q as lepton SM flavours electron, mu, tau 
C{12,pq}∝∑_r[y*1,pry*2,qr+y*2,pry*1,qr]
'''
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any
import sympy as sp
import re
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

# ---------------------------------------------------------------------------
# Direct Weinberg flavor transport in EFT1
# ---------------------------------------------------------------------------

def _expr(raw: Any) -> sp.Expr:
    return sp.sympify(
        str(raw),
        locals={
            "conjugate": sp.conjugate,
            "sqrt": sp.sqrt,
            "log": sp.log,
            "I": sp.I,
            "pi": sp.pi,
            "MF": sp.Symbol("MF", nonzero=True),
            "y1": sp.Symbol("y1"),
            "y2": sp.Symbol("y2"),
        },
    )


def _text(expr: sp.Expr) -> str:
    return sp.sstr(sp.factor(sp.simplify(expr)))


def _load_success(path: Path) -> dict[str, Any]:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    if payload.get("status") != "Success":
        raise ValueError(f"{path} is not successful.")
    return payload


def _mixed_c12_boundary(flavor_seed: dict[str, Any]) -> dict[str, Any]:
    '''Search for flavor in mixed y1-y2 coefficient
    We obtain the flavor tensor form it C_{12,pq}'''
    matches = [
        term
        for term in flavor_seed.get("terms", [])
        if str(term.get("flavor_label", "")) == "C12"
        or {item.get("name") for item in term.get("yukawas", [])} == {"y1", "y2"}
    ]
    if not matches:
        raise ValueError("No mixed y1-y2 flavor boundary was found.")

    reference = matches[0]
    for term in matches[1:]:
        for key in (
            "flavor_kernel_text",
            "flavor_kernel_wolfram",
            "one_generation_reduction",
        ):
            if str(term[key]) != str(reference[key]):
                raise ValueError(
                    "Mixed contractions do not share one canonical C12 flavor boundary."
                )
    return reference


def build_direct_weinberg_flavor_transport(
    flavor_seed_path: Path,
    component_rge_path: Path,
    mu_high: str | sp.Expr,
    mu_low: str | sp.Expr,
    output_path: Path | None = None,
) -> dict[str, Any]:
    '''
    We obtain weinberg tensor
    We check that it is proportional to Weinberg tensor
    We obtain delta kappa as the change in Weinberg coefficient running'''
    flavor_seed = _load_success(flavor_seed_path)
    component_rge = _load_success(component_rge_path)

    if not flavor_seed.get("one_generation_reduction_matches"):
        raise ValueError("Flavor-seed one-generation regression failed.")

    validation = component_rge.get("weinberg_subspace_validation") or {}
    if not validation.get("matches_weinberg_subspace"):
        raise ValueError("Component RGE did not validate the Weinberg subspace.")

    beta_kappa = _expr(validation["beta_kappa_16pi2"])
    c12 = _mixed_c12_boundary(flavor_seed)
    c12_1g = _expr(c12["one_generation_reduction"])

    prefactor = sp.simplify(beta_kappa / c12_1g)

    # These y1,y2 and MF should not be in the prefactor, we have them in our C_12
    y1 = sp.Symbol("y1")
    y2 = sp.Symbol("y2")
    MF = sp.Symbol("MF", nonzero=True)
    if prefactor.has(y1, y2, sp.conjugate(y1), sp.conjugate(y2), MF):
        raise ValueError(
            "Could not factor beta_kappa into a representation coefficient times C12."
        )

    # L = ln(mulow/muhigh)
    hi = _expr(mu_high)
    lo = _expr(mu_low)
    log_ratio = sp.simplify(sp.log(lo / hi))
    running_prefactor = sp.simplify(log_ratio * prefactor)

    one_gen_residual = sp.simplify(prefactor * c12_1g - beta_kappa)
    one_gen_ok = one_gen_residual == 0
    equal_scale_ok = sp.simplify(running_prefactor.subs(lo, hi)) == 0

    payload = {
        "status": "Success" if one_gen_ok and equal_scale_ok else "RegressionFailed",
        "scheme": "fixed_order_one_loop_direct_weinberg_only",
        "purpose": "full-flavor direct EFT1 mixing of the tree LLSS operator into the Weinberg operator",
        "power_counting": {
            "tree_LLSS_boundary_order": "hbar^0",
            "direct_LLSS_to_Weinberg_running_order": "hbar^1",
            "heavy_LLSS_self_running_order": "hbar^1",
            "heavy_running_then_scalar_loop_order": "hbar^2",
            "authoritative_one_loop_includes_heavy_self_running": False,
        },
        "scales": {
            "mu_high": str(mu_high),
            "mu_low": str(mu_low),
            "log_ratio": _text(log_ratio),
        },
        "boundary_tensors": {
            "C12": {
                "flavor_label": "C12",
                "definition_text": c12["flavor_kernel_text"],
                "definition_wolfram": c12["flavor_kernel_wolfram"],
                "one_generation_reduction": c12["one_generation_reduction"],
                "manifestly_symmetric_under_pq": True,
            }
        },
        "direct_weinberg_beta": {
            "beta_kappa_16pi2_one_generation": _text(beta_kappa),
            "canonical_c12_one_generation": _text(c12_1g),
            "flavor_blind_prefactor": _text(prefactor),
            "full_flavor_beta_tensor_text": f"({_text(prefactor)})*C12[p,q]",
        },
        "running_corrections": {
            "Weinberg": {
                "terms": [
                    {
                        "source": "C12",
                        "kind": "flavor_blind",
                        "beta_coefficient": _text(prefactor),
                        "running_coefficient": _text(running_prefactor),
                        "tensor_text": "C12[p,q]",
                        "tensor_wolfram": "C12[p, q]",
                    }
                ],
                "beta_tensor_text": f"({_text(prefactor)})*C12[p,q]",
                "running_tensor_text": f"({_text(running_prefactor)})*C12[p,q]",
                "running_tensor_wolfram": f"({_text(running_prefactor)})*C12[p, q]",
            }
        },
        "one_generation_regression": {
            "matches": one_gen_ok,
            "residual": _text(one_gen_residual),
        },
        "one_generation_reduction_matches": one_gen_ok,
        "equal_scale_running_vanishes": equal_scale_ok,
        "manifestly_symmetric_under_pq": True,
    }

    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    return payload


def _to_wolfram_scalar(raw: Any) -> str:
    '''SymPy scalar converted to mathematica'''
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
    '''Generate full flavour wolfram definition of C_12'''
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
    """Export the validated direct EFT1 C12 -> Weinberg transport to Wolfram."""
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
    parser.add_argument("flavor_seed_json", type=Path)
    parser.add_argument("component_rge_json", type=Path)
    parser.add_argument("--mu-high", default="MF")
    parser.add_argument("--mu-low", default="MS")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    result = build_direct_weinberg_flavor_transport(
        args.flavor_seed_json,
        args.component_rge_json,
        args.mu_high,
        args.mu_low,
        output_path=args.output,
    )
    print(json.dumps({
        "status": result["status"],
        "prefactor": result["direct_weinberg_beta"]["flavor_blind_prefactor"],
        "one_generation_reduction_matches": result["one_generation_reduction_matches"],
        "equal_scale_running_vanishes": result["equal_scale_running_vanishes"],
    }, indent=2))
    return 0 if result["status"] == "Success" else 1
