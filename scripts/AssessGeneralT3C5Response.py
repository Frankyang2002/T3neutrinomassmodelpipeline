"""Validate/export the C5 response to the independent extra quartic atoms.

Uses the *formal* Mathematica results from CompareGeneralT3C5.wl. No physical
numeric C5 or RGBeta coupling conversion is inferred. Fails closed on missing
models, unproven limits, nonlinearity or altered unrelated model coefficients.

    python scripts/AssessGeneralT3C5Response.py \\
       --root output/general_t3_trial \\
       --output output/general_c5_response.json --strict
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any

EXPECTED_MODELS = frozenset({
    "A_m4", "A_m2", "A_p0", "B_m3", "B_m1", "B_p1",
    "C_m3", "C_m1", "C_p1", "D_m2", "D_p0", "D_p2",
    "E_m4", "E_m2", "E_p0", "E_p2",
})
ACTIVE_MODELS = frozenset({"B_m1", "C_m1"})
ACTIVE = frozenset({"lambdaT3ExtraQ020020I1", "lambdaT3ExtraQ020200I1"})
NAME_RE = re.compile(r"lambdaT3ExtraQ[0-9]+I[0-9]+")
YUKAWA_RE = re.compile(r"Coupling\[(y[12]),")


def _read_json(path: Path) -> dict[str, Any]:
    obj = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(obj, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return obj


def classify_atom(atom: dict[str, Any]) -> dict[str, Any]:
    name = str(atom.get("coupling_atom", ""))
    coeff = str(atom.get("linear_coefficient_inputform", ""))
    symbols = sorted(set(NAME_RE.findall(name)))
    yukawas = sorted(set(YUKAWA_RE.findall(coeff)))
    return {
        "coupling_atom": name,
        "matchete_couplings": symbols,
        "holomorphic_or_conjugate": (
            "conjugate" if name.startswith(("Bar[", "Conjugate["))
            else "direct"
        ),
        "yukawa_species_in_coefficient": yukawas,
        "contains_regulator": any(s in coeff for s in ("ϵ", "\\[Epsilon]", "Epsilon")),
        "contains_matching_log": any(s in coeff for s in ("Log[", "log[")),
        "depends_on_other_extra_couplings": bool(NAME_RE.search(coeff)),
        "coefficient_inputform": coeff,
    }


def assess(root: Path) -> dict[str, Any]:
    report_path = root / "c5_algebraic_report.json"
    formal = _read_json(report_path)
    rows = formal.get("results")
    if not isinstance(rows, list):
        raise ValueError("Formal report has no results list")
    by_model = {item.get("model"): item for item in rows if isinstance(item, dict)}
    errors = []
    if set(by_model) != EXPECTED_MODELS:
        errors.append("Formal C5 report is missing or has unexpected model IDs: "
                      f"missing={sorted(EXPECTED_MODELS-set(by_model))}, "
                      f"unexpected={sorted(set(by_model)-EXPECTED_MODELS)}")
    results = []
    for model in sorted(EXPECTED_MODELS):
        item = by_model.get(model, {})
        status = item.get("status", "Missing")
        if item.get("zero_extra_limit_proven") is not True:
            errors.append(f"{model}: zero-extra limit not proven")
        response = []
        if model in ACTIVE_MODELS:
            if status != "ZeroExtraLimitProven":
                errors.append(f"{model}: expected nontrivial proven delta; got {status}")
            if item.get("linear_in_extra_coupling_atoms") is not True or (
                item.get("coefficient_independent_of_other_extra_couplings") is not True
            ):
                errors.append(f"{model}: delta nonlinear or not certified linear in independent extra atoms")
            if item.get("response_reconstruction_residual_proven_zero") is not True:
                errors.append(f"{model}: linear reconstruction residual was not proven zero")
            if item.get("response_at_zero_proven_zero") is not True:
                errors.append(f"{model}: extracted coupling response is nonzero at zero couplings")
            path = root / model / "general_c5_algebraic" / "delta_c5_atoms.json"
            if not path.is_file():
                errors.append(f"{model}: missing derivative export {path}")
            else:
                data = _read_json(path)
                if data.get("model") != model or not isinstance(data.get("atoms"), list):
                    errors.append(f"{model}: malformed derivative export")
                elif (
                    data.get("response_reconstruction_residual_proven_zero") is not True
                    or data.get("response_at_zero_proven_zero") is not True
                    or data.get("polynomial_in_extra_atoms") is not True
                ):
                    errors.append(f"{model}: derivative export lacks a complete algebraic proof")
                else:
                    response = [classify_atom(row) for row in data["atoms"]]
                    present = {nm for entry in response for nm in entry["matchete_couplings"]}
                    if present != ACTIVE:
                        errors.append(f"{model}: expected exactly two C5-active quartics, "
                                      f"got {sorted(present)}")
                    if not response or any(
                        x["depends_on_other_extra_couplings"] for x in response
                    ):
                        errors.append(f"{model}: nonlinear or empty extra C5 derivative export")
                    if any(str(row.get("linear_coefficient_inputform", "")).strip() == "0"
                           or row.get("coefficient_formally_zero") is True
                           for row in data["atoms"]):
                        errors.append(f"{model}: zero-response atoms should be filtered")
        elif status != "TextIdentical":
            errors.append(f"{model}: baseline and extended C5 unexpectedly differ: {status}")
        results.append({
            "model": model,
            "formal_status": status,
            "zero_extra_limit_proven": item.get("zero_extra_limit_proven") is True,
            "extra_c5_derivatives": response,
        })
    return {
        "status": "VerifiedFormalResponse" if not errors else "ReviewRequired",
        "scope": "16 ordinary models; one-loop C5 from extra scalar quartics only",
        "c5_changed_models": sorted(ACTIVE_MODELS),
        "numeric_matching_implemented": False,
        "rgbeta_coupling_normalisations_verified": False,
        "warnings": [
            "Formal coefficient derivatives are not numerical neutrino masses.",
            "The exported C5 may contain poles or matching logarithms; renormalisation and threshold conventions must be applied before numerical use.",
            "Quadratic scalar mixing remains outside this calculation.",
        ],
        "errors": errors,
        "results": results,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("output/general_t3_trial"))
    parser.add_argument("--output", type=Path, default=Path("output/general_c5_response.json"))
    parser.add_argument("--strict", action="store_true")
    ns = parser.parse_args()
    root = ns.root.resolve()
    result = assess(root)
    ns.output.parent.mkdir(parents=True, exist_ok=True)
    ns.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print("C5 formal-response status:", result["status"])
    for item in result["results"]:
        if item["extra_c5_derivatives"]:
            for atom in item["extra_c5_derivatives"]:
                print(item["model"], atom["matchete_couplings"],
                      "Yukawas", atom["yukawa_species_in_coefficient"],
                      "pole", atom["contains_regulator"],
                      "log", atom["contains_matching_log"])
    for error in result["errors"]:
        print("ERROR:", error)
    print("Report:", ns.output)
    return 2 if ns.strict and result["errors"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
