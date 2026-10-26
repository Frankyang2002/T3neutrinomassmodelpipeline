"""Analyse experimental Matchete scalar-quartic audit outputs.

This reports *observed* Matchete validation separately from static operator
counting. It does not imply RGBeta or numerical consistency.

Example:
 python scripts/ReportGeneralScalarTrial.py --root output/general_t3_trial --all
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[1]
MANIFEST = PROJECT / "docs/general_scalar_trial_operator_manifest.json"


def evaluate(root: Path, *, all_models: bool = False) -> dict:
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))["models"]
    results: dict[str, dict] = {}
    for key, entry in manifest.items():
        model, alpha = key.split("_")
        token = ("m" + alpha[1:]) if alpha.startswith("-") else ("p" + alpha[1:])
        label = model + "_" + token
        path = root / label / "audit.json"
        if not all_models and not path.exists():
            continue
        expected = entry["extra_quartic_complex_couplings"]
        result = {"expected_extra_quartics": expected, "file": str(path)}
        if path.exists():
            try:
                data = json.loads(path.read_text(encoding="utf-8-sig"))
                observed = data.get("ExtraQuarticCandidateCount")
                accepted = data.get("AcceptedExtraQuarticCount")
                status = data.get("Status", "Unspecified")
                embedded_expected = data.get('ExpectedExtraQuarticCount')
                result.update(
                    status=status,
                    embedded_expected_extra_quartics=embedded_expected,
                    observed_extra_quartics=observed,
                    accepted_extra_quartics=accepted,
                    baseline_valid=data.get("BaselineLagrangianValid"),
                    combined_valid=data.get("CombinedLagrangianValid"),
                    success=(status == "Success" and observed == expected
                             and accepted == expected
                             and (embedded_expected is None or embedded_expected == expected)
                             and data.get("BaselineLagrangianValid") is True
                             and data.get("CombinedLagrangianValid") is True),
                )
                if status == "FrontendUnavailable":
                    result["diagnostic"] = "Front end unavailable; operator validity has not been tested"
                elif observed != expected and observed is not None:
                    result["diagnostic"] = "Generated invariant count disagrees with static enumeration"
                elif embedded_expected is not None and embedded_expected != expected:
                    result["diagnostic"] = "Matchete-embedded and Python operator manifests disagree"
            except (OSError, ValueError, TypeError) as err:
                result.update(status="UnreadableJSON", success=False, diagnostic=str(err))
        else:
            result.update(status="NotRun", success=False)
        results[key] = result
    return {
        "model_count_checked": len(results),
        "passed": sum(result.get("success", False) for result in results.values()),
        "failed_or_unrun": sum(not result.get("success", False) for result in results.values()),
        "models": results,
        "note": "Symbolic scalar-quartic verification only; no validated general RGBeta, scalar mass mixing, or numerical prediction",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("output/general_t3_trial"))
    parser.add_argument("--all", action="store_true", help="report all 16 expected models, including ones not run")
    parser.add_argument("--output", type=Path, help="write machine-readable report JSON")
    parser.add_argument("--require-complete", action="store_true", help="fail if any inspected model failed or is missing")
    args = parser.parse_args()
    root = args.root if args.root.is_absolute() else PROJECT / args.root
    report = evaluate(root, all_models=args.all)
    print(f"Audits: {report['model_count_checked']} | passed: {report['passed']} | failed/unrun: {report['failed_or_unrun']}")
    if report["model_count_checked"] == 0:
        print("No audit JSON files found. Run scripts/RunGeneralScalarTrial.ps1 first; "
              "use --all to list all 16 expected model points.")
    for key, item in report["models"].items():
        print(f"{key:7} {item['status']:34} "
              f"expected={item['expected_extra_quartics']:2} "
              f"observed={str(item.get('observed_extra_quartics', '-')):>4} "
              f"accepted={str(item.get('accepted_extra_quartics', '-')):>4}")
        if item.get("diagnostic"):
            print("  ", item["diagnostic"])
    if args.output:
        output = args.output if args.output.is_absolute() else PROJECT / args.output
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        print(f"JSON report: {output}")
    if report["model_count_checked"] == 0:
        return 1
    return 1 if (args.require_complete and report["failed_or_unrun"] > 0) else 0


if __name__ == "__main__":
    raise SystemExit(main())
