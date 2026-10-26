"""Consolidate 16-model Matchete audit, textual C5 and formal algebraic C5.

A missing or unverified algebraic result is NEVER labelled a passed proof.
Usage (project root):
    python scripts/SummarizeGeneralT3Physics.py --root output/general_t3_trial \
        --output output/general_t3_physics_summary.json --require-proof
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[1]
EXPECTED_CLASSES = {
    "A_m4", "A_m2", "A_p0", "B_m3", "B_m1", "B_p1",
    "C_m3", "C_m1", "C_p1", "D_m2", "D_p0", "D_p2",
    "E_m4", "E_m2", "E_p0", "E_p2",
}
C5_CLASSES = {"B_m1", "C_m1"}


def _json(path: Path) -> dict | None:
    if not path.is_file():
        return None
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    return value if isinstance(value, dict) else None


def compile_summary(root: Path) -> dict:
    symbolic = _json(root.parent / "general_matching_symbolic_report.json")
    if symbolic is None:
        symbolic = _json(root / "general_matching_symbolic_report.json")
    audits = _json(root.parent / "general_scalar_audit_report.json")
    if audits is None:
        audits = _json(root / "general_scalar_audit_report.json")
    proof = _json(root / "c5_algebraic_report.json")
    if proof is None:
        proof = _json(root.parent / "c5_algebraic_report.json")

    symbolic_rows = {
        item.get("model"): item for item in (symbolic or {}).get("results", [])
        if isinstance(item, dict)
    }
    proof_rows = {
        item.get("model"): item for item in (proof or {}).get("results", [])
        if isinstance(item, dict)
    }
    # Keep field names minimal and grounded in these reports. The complete
    # operator audit is separately checked by --require-complete there.
    rows = []
    for name in sorted(EXPECTED_CLASSES):
        text = symbolic_rows.get(name)
        algebra = proof_rows.get(name)
        observed = (
            (text.get("variants", {}).get("extended", {})
             .get("extra_couplings_observed", [])) if text else None
        )
        if name in C5_CLASSES:
            symbolic_status = ("NewCouplingsObserved" if observed else
                               "MissingOrNoNewCouplings")
        else:
            symbolic_status = ("NoExtraCouplingsObserved" if observed == [] else
                               "UnexpectedOrMissing")
        proof_status = algebra.get("status", "NotRun") if algebra else "NotRun"
        reduced = bool(algebra and algebra.get("zero_extra_limit_proven") is True)
        rows.append({
            "model": name,
            "expected_additional_one_loop_c5": name in C5_CLASSES,
            "symbolic_status": symbolic_status,
            "extra_c5_couplings": observed,
            "baseline_extended_text_equal": (
                text.get("baseline_extended_c5_text_equal_ignoring_whitespace")
                if text else None
            ),
            "zero_extra_proof_status": proof_status,
            "zero_extra_limit_proven": reduced,
            "formal_delta_file": algebra.get("delta_inputform_file") if algebra else None,
            "atom_coefficients_file": algebra.get("atom_coefficients_file") if algebra else None,
        })
    symbolically_consistent = all(
        row["symbolic_status"] ==
        ("NewCouplingsObserved" if row["expected_additional_one_loop_c5"]
         else "NoExtraCouplingsObserved") for row in rows
    )
    all_proven = all(row["zero_extra_limit_proven"] for row in rows)
    return {
        "scope": "16 separate supported ordinary T3 models; added scalar quartics only",
        "expected_model_count": len(EXPECTED_CLASSES),
        "symbolically_consistent": symbolically_consistent,
        "all_zero_extra_limits_formally_proven": all_proven,
        "c5_changed_models": [row["model"] for row in rows
                              if row["symbolic_status"] == "NewCouplingsObserved"],
        "assumptions": [
            "Quadratic scalar mixing is not implemented.",
            "Matchete-to-RGBeta coupling-normalisation mapping is not yet verified.",
            "No extended numerical RG evolution or neutrino observables claimed.",
            "Algebraic reduction, when proven, uses zero new quartic coefficients at the matching boundary, not zero values at all RG scales.",
        ],
        "audits_report_available": audits is not None,
        "proof_report_available": proof is not None,
        "models": rows,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("output/general_t3_trial"))
    parser.add_argument("--output", type=Path,
                        default=Path("output/general_t3_physics_summary.json"))
    parser.add_argument("--require-proof", action="store_true")
    ns = parser.parse_args()
    root = ns.root if ns.root.is_absolute() else PROJECT / ns.root
    destination = ns.output if ns.output.is_absolute() else PROJECT / ns.output
    data = compile_summary(root)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    csv_path = destination.with_suffix(".csv")
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=[
            "model", "expected_additional_one_loop_c5", "symbolic_status",
            "extra_c5_couplings", "baseline_extended_text_equal",
            "zero_extra_proof_status", "zero_extra_limit_proven",
        ])
        writer.writeheader()
        for row in data["models"]:
            writer.writerow({k: (";".join(v) if isinstance(v, list) else v)
                             for k, v in row.items() if k in writer.fieldnames})
    print(f"Symbolically consistent: {data['symbolically_consistent']}")
    print(f"Zero-extra limits proven: {data['all_zero_extra_limits_formally_proven']}")
    print(f"Models with new C5 terms: {', '.join(data['c5_changed_models'])}")
    print(f"JSON: {destination}")
    print(f"CSV: {csv_path}")
    if ns.require_proof and not (data["symbolically_consistent"] and
                                data["all_zero_extra_limits_formally_proven"]):
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
