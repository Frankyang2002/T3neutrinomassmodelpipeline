"""Validate a fresh A(-4) EFT1 installed-RGBeta operator-level witness.

This checks RGBeta's live $quartics registry and formal Lam tensor, with
source hashing, against the *declared* A(-4) invariants. It never invents a
numerical tGen or delS2 basis, changes beta functions or issues a physical
Matchete-to-RGBeta certificate. Native numerical axes are a separate gate.

Usage:
 python scripts/CheckAminus4RGBetaOperatorWitness.py \
   --native output/a_minus4_eft1_batch2/A_m4_native_operator.json \
   --report output/a_minus4_eft1_batch2/A_m4_operator_check.json
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import re
import sys
from typing import Mapping, Any

ROOT = Path(__file__).resolve().parents[1]
FORMAT = "a-minus4-rgbeta-operator-witness-v1"
STATUS = "NATIVE_SYMBOLIC_OPERATOR_TENSOR_CAPTURED_COMPONENT_MAP_PENDING"
MODEL_SOURCE = Path("RGE/running/rgbeta/T3RGBetaModel.wl")
QUARTIC_FIELDS = {
    "lambdaH": "{Bar[H],H,Bar[H],H}",
    "lambdaS1": "{Bar[S1],S1,Bar[S1],S1}",
    "lambdaS2": "{Bar[S2],S2,Bar[S2],S2}",
    "lambdaH1": "{Bar[H],H,Bar[S1],S1}",
    "lambdaH2": "{Bar[H],H,Bar[S2],S2}",
    "lambda12": "{Bar[S1],S1,Bar[S2],S2}",
    "lambdaT3": "{H,H,S1,Bar[S2]}",
    "lambdaH2Adj": "{Bar[H],H,Bar[S2],S2}",
    "lambdaS2Adj": "{Bar[S2],S2,Bar[S2],S2}",
}
EXPECTED_GROUP_HEADS = {
    "lambdaH": {"del[SU2L[fund]": 2},
    "lambdaS1": {},
    "lambdaS2": {"del[SU2L[S2]": 2},
    "lambdaH1": {"del[SU2L[fund]": 1},
    "lambdaH2": {"del[SU2L[fund]": 1, "del[SU2L[S2]": 1},
    "lambda12": {"del[SU2L[S2]": 1},
    "lambdaT3": {"delS2[SU2L": 1},
    "lambdaH2Adj": {"tGen[SU2L[fund]": 1, "tGen[SU2L[S2]": 1},
    "lambdaS2Adj": {"tGen[SU2L[S2]": 2},
}
HALVED_INVARIANTS = frozenset({"lambdaH", "lambdaS1", "lambdaS2", "lambdaS2Adj"})


class OperatorWitnessError(ValueError):
    """Unsafe, stale, inconsistent or incomplete native operator record."""


def require(test: Any, message: str) -> None:
    if not test:
        raise OperatorWitnessError(message)


def clean(expression: str) -> str:
    require(isinstance(expression, str), "Native InputForm must be a string")
    return re.sub(r"\s+", "", expression)


def check_invariants(rows: list[Mapping[str, Any]]) -> list[dict[str, Any]]:
    require(isinstance(rows, list) and len(rows) == 9,
            "Native quartic registry must contain exactly nine entries")
    found = set()
    checked = []
    for item in rows:
        require(isinstance(item, Mapping), "Native quartic operator is not a mapping")
        name = item.get("name")
        require(name in QUARTIC_FIELDS and name not in found,
                f"Unexpected/duplicate native quartic {name!r}")
        found.add(name)
        require(item.get("registration_complete") is True,
                f"Incomplete native registry: {name}")
        require(clean(item.get("fields")) == QUARTIC_FIELDS[name],
                f"Incorrect field order for {name}")
        expected_real = name != "lambdaT3"
        require(item.get("self_conjugate") is expected_real,
                f"Incorrect self-conjugation for {name}")
        require(item.get("native_unique_arrangements", 0) > 0
                and isinstance(item["native_unique_arrangements"], int),
                f"Missing unique field arrangements: {name}")
        require(clean(item.get("coupling_inputform")) == f"6*{name}",
                f"Unexpected native quartic internal tensor ratio: {name}")
        bar_expected = f"6*{name}" if expected_real else "6*Bar[lambdaT3]"
        require(clean(item.get("coupling_bar_inputform")) == bar_expected,
                f"Unexpected conjugate quartic tensor ratio: {name}")
        inv = clean(item.get("invariant_inputform"))
        require(inv and not any(bad in inv for bad in ("$Failed", "Missing[")),
                f"Native invariant unavailable for {name}")
        for head, count in EXPECTED_GROUP_HEADS[name].items():
            require(inv.count(head) == count,
                    f"{name}: expected {count} occurrence(s) of {head}")
        if name in HALVED_INVARIANTS:
            # Both Mathematica 1/2 and the trivially equivalent 0.5 are
            # accepted as diagnostics; no numerical index conventions are
            # deduced from this literal source-factor check.
            require("/2" in inv or "0.5" in inv or inv == "1/2",
                    f"Native invariant factor 1/2 missing: {name}")
        checked.append({
            "name": name,
            "field_order_verified": True,
            "symbolic_invariant_group_structure_verified": True,
            "native_internal_coupling_prefactor": 6,
            "physical_action_coefficient_verified": False,
        })
    require(found == set(QUARTIC_FIELDS), "Native registry incomplete")
    return checked


def analyze(native: Mapping[str, Any], root: Path) -> dict[str, Any]:
    require(isinstance(native, Mapping), "Native witness must be a JSON object")
    require(native.get("format") == FORMAT and native.get("status") == STATUS,
            "Native witness missing/blocked or wrong format")
    require(native.get("model") == "A_-4" and native.get("stage") == "EFT1"
            and native.get("dimensions_alpha") == [1, 3, 2, -4],
            "Wrong representation or EFT stage")
    require(native.get("source_model_path") == MODEL_SOURCE.as_posix(),
            "Unexpected model source")
    model_path = Path(root) / MODEL_SOURCE
    require(model_path.is_file(), f"Missing source file: {model_path}")
    local_sha = sha256(model_path.read_bytes()).hexdigest()
    require(native.get("source_model_sha256") == local_sha,
            "Native witness is stale: T3RGBetaModel.wl hash differs")
    package = native.get("rgbeta_package_file")
    require(isinstance(package, str) and "RGBeta" in package,
            "Installed RGBeta package provenance missing")
    symbols = native.get("native_quartic_symbols")
    require(isinstance(symbols, list) and len(symbols) == 9
            and sorted(symbols) == sorted(QUARTIC_FIELDS),
            "Live native quartic symbol inventory mismatch")
    checked = check_invariants(native.get("quartic_operators"))
    species = native.get("scalar_species_indices")
    require(isinstance(species, Mapping) and set(species) ==
            {"H", "Bar[H]", "S1", "Bar[S1]", "S2", "Bar[S2]"},
            "Native real/complex scalar species inventory mismatch")
    indices = list(species.values())
    require(all(isinstance(x, int) and not isinstance(x, bool) and x > 0
                for x in indices) and len(set(indices)) == 6,
            "Native field index map missing/ambiguous")
    require(native.get("native_Lam_symbol_names") == ["RGBeta`PackageScope`Lam"],
            "Internal RGBeta Lam tensor symbol not uniquely resolved")
    require(native.get("native_Lam_formal_tensor_available") is True,
            "The native package did not evaluate Lam[a,b,c,d]")
    formal = native.get("native_Lam_formal_tensor_inputform")
    require(sorted(native.get("native_Lam_visible_quartic_symbols", [])) ==
            sorted(QUARTIC_FIELDS),
            "Native Lam tensor did not visibly contain every quartic coupling")
    require(isinstance(formal, str) and "TStructure[" in formal
            and "SparseArray[" in formal and len(formal) > 50,
            "Native Lam tensor lacks TStructure/SparseArray evidence")
    for name in QUARTIC_FIELDS:
        require(re.search(r"(?<![A-Za-z0-9_$])" + re.escape(name) +
                          r"(?![A-Za-z0-9_$])", formal) is not None,
                f"Native Lam tensor does not contain independent symbol {name}")
    require(native.get("native_Lam_string_length", -1) >= len(formal),
            "Malformed Lam tensor length metadata")
    require(native.get("native_Lam_truncated") is False,
            "Native Lam tensor was truncated; cannot validate full expression")
    # Don't allow a diagnostic JSON to silently state an unsupported physical
    # mapping, irrespective of how convincing the native symbolic material is.
    for field in ("matchete_physical_coupling_map_certified",
                  "native_to_reference_component_axes_certified",
                  "rgbeta_quartic_betas_certified", "production_modified"):
        require(native.get(field) is False,
                f"Unsafe or missing certification/production flag: {field}")
    require(native.get("native_operator_registry_complete") is True,
            "Incomplete source-registered native operator set")
    numeric_triplet = native.get("native_triplet_generator_components_numeric") is True
    numeric_cg = native.get("native_triplet_CG_components_numeric") is True
    numeric_fund = native.get("native_fund_generator_components_numeric") is True
    for key, available in (("native_triplet_generator_components_inputform", numeric_triplet),
                           ("native_triplet_CG_components_inputform", numeric_cg),
                           ("native_fund_generator_components_inputform", numeric_fund)):
        val = native.get(key)
        require(isinstance(val, str) and
                ((val != "unavailable") if available else val == "unavailable"),
                f"Missing native numeric-evidence field: {key}")
    # A numeric matrix in the installed package still does not identify its
    # SU(2) axes with Matchete without an independent CG/generator map.
    return {
        "format": "a-minus4-rgbeta-operator-check-v1",
        "status": "SYMBOLIC_OPERATOR_REGISTRY_VERIFIED_NUMERIC_AXES_PENDING",
        "model": "A_-4", "stage": "EFT1",
        "native_source_sha256_verified": local_sha,
        "installed_rgbeta_source": package,
        "registered_operators_verified": len(checked),
        "operator_results": checked,
        "native_formal_Lam_tensor_recorded": True,
        "native_numeric_fund_generators": numeric_fund,
        "native_numeric_triplet_generators": numeric_triplet,
        "native_numeric_triplet_CG": numeric_cg,
        "native_to_reference_component_axes_certified": False,
        "physical_matchete_quartic_map_certified": False,
        "scalar_beta_polynomials_certified": False,
        "next_gate": (
            "Compare native symbolic Lam quartic kernels and projector "
            "against independent invariant tensors; obtain a native-to-Matchete "
            "component-axis witness before any physical beta certification"
        ),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args(argv)
    if args.native.resolve() == args.report.resolve():
        parser.error("Native and validation JSON files must differ")
    try:
        raw = json.loads(args.native.read_text(encoding="utf-8-sig"))
        result = analyze(raw, ROOT)
    except (ValueError, OSError, KeyError, TypeError) as exc:
        print(f"BLOCKED: {exc}", file=sys.stderr)
        return 2
    if args.report.exists():
        print(f"BLOCKED: refusing to overwrite {args.report}", file=sys.stderr)
        return 2
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print("Native A(-4) operator registry: 9/9 VERIFIED (symbolic)")
    print("Installed RGBeta Lam tensor: RECORDED (native formal tensor)")
    print("Native numeric triplet generator / CG:",
          "AVAILABLE" if result["native_numeric_triplet_generators"] and
          result["native_numeric_triplet_CG"] else "NOT EXPOSED")
    print("Physical Matchete component/coupling map certified: NO")
    print("Report:", args.report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
