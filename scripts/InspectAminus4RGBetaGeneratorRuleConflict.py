"""Validate a fresh RGBeta v1.2.0 generator-rule conflict witness.

Read-only. No modification of RGBeta, source, scalar betas, or solver.
The uploaded installed GroupsAndIndices.m SHA-256 is pinned to prevent
confusing this particular source audit with a different package version.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
FORMAT = "a-minus4-rgbeta-generator-rule-conflict-v1"
EXPECTED_GROUPS_SHA256 = "213a312fae28ad72f738584bacd7985dd1179883523be2f1559723f8c1f8c92a"
STATUSES = frozenset({
    "NATIVE_SAME_ORIENTATION_FIERZ_PATH_DISAGREEMENT",
    "NATIVE_SAME_ORIENTATION_PATHS_AGREE",
    "NATIVE_SAME_ORIENTATION_PATH_COMPARISON_UNRESOLVED",
})


class AuditBlocked(ValueError):
    """Native witness lacks required evidence or has changed source."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AuditBlocked(message)


def sha256_file(path: Path) -> str:
    require(path.is_file(), f"missing required source file: {path}")
    return hashlib.sha256(path.read_bytes()).hexdigest()


def installed_groups_default() -> Path:
    roaming = os.environ.get("APPDATA")
    if roaming:
        return (Path(roaming) / "Wolfram" / "Applications" / "RGBeta"
                / "Package" / "GroupsAndIndices.m")
    return (Path.home() / "AppData" / "Roaming" / "Wolfram" /
            "Applications" / "RGBeta" / "Package" / "GroupsAndIndices.m")


def classify_native(data: dict[str, Any]) -> dict[str, Any]:
    require(isinstance(data, dict), "native witness is not a JSON object")
    require(data.get("format") == FORMAT, "unknown native witness format")
    require(data.get("model") == "A_-4" and data.get("stage") == "EFT1",
            "native model/stage mismatch")
    require(data.get("status") in STATUSES, "invalid native status")
    require(data.get("test_kind") ==
            "same-orientation generator contraction; no beta functions",
            "unexpected native test")
    for key in (
        "fund_direct_reduction_inputform",
        "fund_four_distinct_reduction_inputform",
        "fund_fierz_then_contract_reduction_inputform",
        "fund_direct_minus_fierz_path_inputform",
    ):
        require(isinstance(data.get(key), str) and data[key].strip(),
                f"missing native algebra evidence: {key}")
    for key in (
        "native_direct_matches_generic",
        "native_direct_matches_fierz",
        "native_fierz_path_matches_reference",
        "native_two_paths_disagree",
        "physical_scalar_beta_error_proven",
        "physical_matchete_coupling_map_certified",
        "source_modified",
    ):
        require(type(data.get(key)) is bool, f"invalid native flag: {key}")
    require(data["physical_scalar_beta_error_proven"] is False and
            data["physical_matchete_coupling_map_certified"] is False and
            data["source_modified"] is False,
            "native source or certification safety flag violated")
    status = data["status"]
    conflict = data["native_two_paths_disagree"]
    if status == "NATIVE_SAME_ORIENTATION_FIERZ_PATH_DISAGREEMENT":
        require(conflict and data["native_fierz_path_matches_reference"],
                "contradiction claim not supported by two native paths")
        category = "NATIVE_GROUP_REDUCTION_CONFLICT_OBSERVED"
    elif status == "NATIVE_SAME_ORIENTATION_PATHS_AGREE":
        require(not conflict, "native status disagrees with conflict flag")
        category = "NATIVE_PATHS_AGREE_FOR_TESTED_CONTRACTION"
    else:
        require(not conflict, "unresolved status disagrees with conflict flag")
        category = "INCONCLUSIVE_NATIVE_REDUCTION"
    return {
        "native_algebra_status": category,
        "native_direct_matches_generic_rule": data["native_direct_matches_generic"],
        "native_direct_matches_su2_fierz": data["native_direct_matches_fierz"],
        "native_fierz_path_matches_su2_reference":
            data["native_fierz_path_matches_reference"],
        "native_two_paths_disagree": conflict,
        "beta_function_causality_established": False,
        "matchete_physical_basis_certified": False,
    }


def inspect(native_path: Path, model_path: Path,
            installed_source: Path) -> dict[str, Any]:
    require(native_path.is_file(), f"missing native witness: {native_path}")
    data = json.loads(native_path.read_text(encoding="utf-8-sig"))
    require(isinstance(data.get("model_source_sha256"), str),
            "missing model-source digest")
    actual_model = sha256_file(model_path)
    require(actual_model == data["model_source_sha256"],
            "local T3RGBetaModel.wl differs from native witness")
    actual_groups = sha256_file(installed_source)
    require(actual_groups == EXPECTED_GROUPS_SHA256,
            "installed GroupsAndIndices.m differs from uploaded audited source")
    require(actual_groups == data.get("installed_groups_source_sha256"),
            "native kernel and Python observed different GroupsAndIndices.m")
    result = classify_native(data)
    return {
        "format": "a-minus4-rgbeta-generator-rule-conflict-audit-v1",
        "status": result["native_algebra_status"],
        "model": "A_-4",
        "stage": "EFT1",
        "native_file": str(native_path),
        "installed_groups_source": str(installed_source),
        "installed_groups_sha256": actual_groups,
        "model_source_sha256": actual_model,
        "native_algebra": result,
        "known_source_rules": {
            "generic_same_orientation":
                "GroupsAndIndices.m:190-191; -C2(R) delta",
            "fundamental_SU2_fierz":
                "GroupsAndIndices.m:429-431; +1/4 delta for summed same-orientation fundamental pair",
            "one_loop_scalar_contraction":
                "TensorCalculations.m:979 and RGCoefficients.m:54",
        },
        "physical_quartic_beta_certified": False,
        "matchete_coupling_map_certified": False,
        "production_modified": False,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--installed-source", type=Path,
                        default=installed_groups_default())
    a = parser.parse_args(argv)
    try:
        result = inspect(a.native,
                         ROOT / "RGE" / "running" / "rgbeta" / "T3RGBetaModel.wl",
                         a.installed_source)
        require(not a.report.exists(), f"refusing to overwrite report {a.report}")
        a.report.parent.mkdir(parents=True, exist_ok=True)
        a.report.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    except (OSError, UnicodeError, ValueError, json.JSONDecodeError) as exc:
        print(f"BLOCKED: {exc}", file=sys.stderr)
        return 2
    print("Native algebra status:", result["status"])
    print("Native direct/Fierz path difference:",
          result["native_algebra"]["native_two_paths_disagree"])
    print("Quartic-beta causality certified: NO")
    print("Report:", a.report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
