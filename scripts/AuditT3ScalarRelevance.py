"""Batch 1: auditable scalar-operator relevance / fixed-order power-counting inventory.

Read-only: consumes the existing 16-model scalar inventories and saved Matchete
comparison reports.  It does NOT calculate Feynman diagrams, beta functions,
C5 coefficients, or change production inputs.

Run from the project root:
    python -m scripts.AuditT3ScalarRelevance
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
MODELS = (
    "A_-4", "A_-2", "A_+0", "B_-3", "B_-1", "B_+1", "C_-3", "C_-1",
    "C_+1", "D_-2", "D_+0", "D_+2", "E_-4", "E_-2", "E_+0", "E_+2",
)
SPECIES = ("H", "Hdag", "S1", "S1dag", "S2", "S2dag")
INPUTS = {
    "manifest": "docs/general_scalar_trial_operator_manifest.json",
    "matching": "output/general_matching_symbolic_report.json",
    "physics": "output/general_t3_physics_summary.json",
    "audit": "output/general_scalar_audit_report.json",
    "quadratic": "output/general_t3_quadratic_eft1/quadratic_basis_audit.json",
}


class AuditError(ValueError):
    """Input missing, inconsistent, or not sufficiently verified for this audit."""


def require(test: bool, reason: str) -> None:
    if not test:
        raise AuditError(reason)


def tag(model: str) -> str:
    family, alpha = model.split("_")
    return family + "_" + ("m" + alpha[1:] if alpha.startswith("-") else "p" + alpha[1:])


def counts_key(counts: list[int]) -> str:
    require(len(counts) == 6 and all(type(n) is int and n >= 0 for n in counts),
            f"Invalid field count vector: {counts!r}")
    return "".join(map(str, counts))


def load_sources(root: Path) -> tuple[dict[str, Any], dict[str, str]]:
    sources: dict[str, Any] = {}
    hashes: dict[str, str] = {}
    for name, rel in INPUTS.items():
        path = root / rel
        require(path.is_file(), f"Missing required source: {rel}")
        blob = path.read_bytes()
        hashes[rel] = hashlib.sha256(blob).hexdigest()
        try:
            sources[name] = json.loads(blob)
        except (ValueError, UnicodeError) as exc:
            raise AuditError(f"Invalid JSON: {rel}: {exc}") from exc
        require(isinstance(sources[name], dict), f"Expected JSON object: {rel}")
    return sources, hashes


def validate_sources(src: dict[str, Any], strict: bool = True) -> tuple[dict[str, dict], dict[str, dict]]:
    manifest = src["manifest"]
    models: dict[str, dict] = manifest.get("models", {})
    expected = set(MODELS) if strict else set(models)
    require(bool(expected), "No model records in operator manifest")
    require(set(models) == expected, "Model set mismatch in operator manifest")
    if strict:
        require(manifest.get("model_count") == 16, "16-model manifest count incorrect")

    audits: dict[str, dict] = src["audit"].get("models", {})
    physics_list = src["physics"].get("models", [])
    matching_list = src["matching"].get("results", [])
    quadratic_models: dict[str, dict] = src["quadratic"].get("models", {})
    require(isinstance(physics_list, list) and isinstance(matching_list, list),
            "Unexpected physics or matching summary format")
    physics = {r["model"]: r for r in physics_list}
    matching = {r["model"]: r for r in matching_list}
    require(len(physics) == len(physics_list) and len(matching) == len(matching_list),
            "Duplicated model records in C5 reports")
    tags = {tag(model) for model in expected}
    require(set(audits) == expected, "Quartic audit model set mismatch")
    require(set(quadratic_models) == expected, "Quadratic audit model set mismatch")
    require(set(physics) == tags, "Physics summary model set mismatch")
    require(set(matching) == tags, "Matching comparison model set mismatch")
    require(src["physics"].get("all_zero_extra_limits_formally_proven") is True,
            "All-model zero-extra algebraic limit is not certified")
    require(src["physics"].get("symbolically_consistent") is True,
            "Physics source not marked symbolically consistent")
    if strict:
        require(src["audit"].get("passed") == 16 and src["audit"].get("failed_or_unrun") == 0,
                "16/16 saved quartic audits are not successful")
        require(src["quadratic"].get("models_checked") == 16 and
                src["quadratic"].get("quadratic_slots") == 8,
                "Quadratic inventory not certified 8/16")
    for model in expected:
        a = audits[model]
        require(a.get("success") is True and a.get("baseline_valid") is True
                and a.get("combined_valid") is True, f"Quartic Matchete audit failed: {model}")
        p = physics[tag(model)]
        require(p.get("zero_extra_limit_proven") is True,
                f"Zero-extra limit unproven: {model}")
        m = matching[tag(model)]
        require(m.get("extended_log_operator_count_matches_manifest") is True,
                f"Extended operator count not supported: {model}")
        for variant in ("baseline", "extended"):
            require(m.get("variants", {}).get(variant, {}).get("matched_c5_artifact_available") is True,
                    f"Missing {variant} C5 artifact: {model}")
        n = models[model].get("extra_quartic_complex_couplings")
        require(a.get("accepted_extra_quartics") == n and
                m.get("accepted_extra_quartics_in_extended_log") == n,
                f"Quartic invariant counts differ: {model}")
        q_num = models[model].get("quadratic_complex_couplings_not_enabled")
        qm = quadratic_models[model]
        require(q_num == qm.get("independent_complex_mass_parameters_in_declared_basis"),
                f"Quadratic counts differ: {model}")
        require(q_num == len(models[model].get("quadratic_structures_not_enabled", [])),
                f"Quadratic field-structure count mismatch: {model}")
    return physics, matching


def _entry(model: str, kind: str, name: str, structure: str, counts: list[int] | None,
           index: int | None, relevance: str, evidence: str, rationale: str,
           omission: str, **extra: Any) -> dict[str, Any]:
    return {
        "model": model, "kind": kind, "coupling": name, "field_structure": structure,
        "counts_H_Hdag_S1_S1dag_S2_S2dag": counts, "invariant_index": index,
        "higgs_legs": counts[0] + counts[1] if counts is not None else None,
        "classification": relevance, "evidence_level": evidence,
        "physics_reason": rationale, "omission_rule": omission, **extra,
    }


def classify_sources(src: dict[str, Any], *, strict: bool = True) -> dict[str, Any]:
    physics, matching = validate_sources(src, strict=strict)
    models = src["manifest"]["models"]
    rows: list[dict[str, Any]] = []
    for model in sorted(models):
        data = models[model]
        code = tag(model)
        observed = set(matching[code]["variants"]["extended"]["extra_couplings_observed"])
        proven_record = set(physics[code].get("extra_c5_couplings", []))
        require(observed == proven_record,
                f"Saved C5 coupling reports disagree: {model}")
        all_couplings: set[str] = set()

        # Topology-relevant interaction is present in T3Topology.wl.
        rows.append(_entry(
            model, "baseline_topology", "lambdaT3", "H H S1 S2dag + h.c.",
            [2, 0, 1, 0, 0, 1], None,
            "leading_one_loop_topology", "source_code_and_topology",
            "The two Yukawa vertices and HH-S1-S2dag vertex form the selected T3 loop.",
            "RETAIN for the stated T3 contribution; not evidence that other leading couplings vanish.",
        ))
        for mass, scalar in (("MS1", "S1"), ("MS2", "S2")):
            rows.append(_entry(
                model, "baseline_quadratic", mass, f"{scalar}dag {scalar}", None, None,
                "tree_level_scalar_spectrum", "source_code",
                "Enters scalar propagators already at tree level.",
                "RETAIN; cannot be dismissed by loop-order counting.",
            ))
        for name, structure, higgs in (
            ("lambdaS1Inv*", "(S1dag S1)(S1dag S1)", False),
            ("lambdaS2Inv*", "(S2dag S2)(S2dag S2)", False),
            ("lambdaH1Inv*", "(Hdag H)(S1dag S1)", True),
            ("lambdaH2Inv*", "(Hdag H)(S2dag S2)", True),
            ("lambda12Inv*", "(S1dag S1)(S2dag S2)", False),
        ):
            rows.append(_entry(
                model, "baseline_quartic_family", name, structure, None, None,
                "one_loop_higgs_vertex_candidate" if higgs else "not_direct_single_vertex_LLSS_to_LLHH",
                "source_code_plus_graph_power_counting",
                ("Contains two Higgs fields and may affect the selected LLHH amplitude; topology/flavour "
                 "contractions require a dedicated check." if higgs else
                 "No Higgs external fields: cannot alone saturate LLSS -> LLHH with one quartic vertex. "
                 "Can affect UV running, mass renormalisation and other amplitudes."),
                ("DO NOT DROP solely from Higgs-leg counting." if higgs else
                 "Can omit its DIRECT single-vertex LLSS->LLHH insertion at one loop, "
                 "not its independent value or effects on other quantities."),
            ))

        quartic_total = 0
        for struct in data["extra_quartic_structures"]:
            counts = struct["counts"]
            require(sum(counts) == 4, f"Nonquartic count in {model}: {counts}")
            hlegs = counts[0] + counts[1]
            require(hlegs in (0, 2), f"Unexpected quartic Higgs-leg pattern {model}: {counts}")
            n = struct["invariants"]
            require(type(n) is int and n >= 1, f"Invalid SU2 multiplicity: {model}")
            for index in range(1, n + 1):
                coupling = "lambdaT3ExtraQ" + counts_key(counts) + f"I{index}"
                require(coupling not in all_couplings, f"Repeated extra coupling: {coupling}")
                all_couplings.add(coupling)
                quartic_total += 1
                if coupling in observed:
                    require(hlegs == 2, f"Unexpected zero-Higgs source directly present in C5: {model}")
                    cl = "observed_in_saved_one_loop_C5"
                    ev = "symbolic_matching_artifact"
                    rationale = ("The named coupling occurs in the saved extended one-loop C5; "
                                 "the algebraic zero-extra limit is proven. This is symbolic evidence, "
                                 "not independent full physical normalization.")
                    omission = "RETAIN if nonzero; zero is an explicit boundary restriction, not a higher-order argument."
                elif hlegs == 2:
                    cl = "higgs_vertex_not_seen_in_saved_C5"
                    ev = "saved_matching_absence_only"
                    rationale = ("Two Higgs legs make a direct one-vertex scalar-threshold diagram possible "
                                 "on power counting, but the named coupling is absent from the available C5 artifact. "
                                 "No claim for other mass bases, thresholds or orders.")
                    omission = "Cannot categorically omit; Batch 2 must prove absence/cancellation for the claimed setup."
                else:
                    cl = "no_direct_single_vertex_LLSS_to_LLHH"
                    ev = "field_leg_counting"
                    rationale = ("Four BSM scalar legs and zero Higgs legs: with one LLSS insertion, "
                                 "this quartic alone cannot supply the two external H fields in "
                                 "the selected one-loop LLSS -> LLHH amplitude. Indirect RG/mass effects remain.")
                    omission = ("May omit direct contribution to this specific one-loop amplitude; "
                                "not a declaration that the physical coupling vanishes.")
                rows.append(_entry(model, "additional_quartic", coupling, struct["field_structure"],
                                   counts, index, cl, ev, rationale, omission,
                                   observed_in_saved_C5=(coupling in observed)))
        require(quartic_total == data["extra_quartic_complex_couplings"],
                f"Extra quartic inventory mismatch: {model}")
        require(all(c in all_couplings for c in observed),
                f"A coupling in saved C5 is absent from the inventory: {model}: {observed-all_couplings}")

        quadratic_total = 0
        for struct in data["quadratic_structures_not_enabled"]:
            counts = struct["counts"]
            require(sum(counts) == 2 and counts[0] == counts[1] == 0,
                    f"Unexpected quadratic content: {model}: {counts}")
            n = struct["invariants"]
            for index in range(1, n + 1):
                quadratic_total += 1
                rows.append(_entry(
                    model, "additional_quadratic", "mu2_extra_" + counts_key(counts) + f"_I{index}",
                    struct["field_structure"], counts, index,
                    "unsuppressed_tree_scalar_mass_effect", "symmetry_inventory_and_mass_matrix",
                    "An allowed quadratic operator changes physical masses or scalar mixing at tree level; "
                    "it is not suppressed by an additional loop factor. Name here is an audit identifier, "
                    "not a declared Wolfram coupling symbol.",
                    "RETAIN for arbitrary nonzero masses/mixing or explicitly set coefficient to zero "
                    "as a benchmark assumption; verify any RG regeneration.",
                ))
        require(quadratic_total == data["quadratic_complex_couplings_not_enabled"],
                f"Extra quadratic inventory mismatch: {model}")

    new_q = [r for r in rows if r["kind"] == "additional_quartic"]
    new_m = [r for r in rows if r["kind"] == "additional_quadratic"]
    status_count = dict(sorted(Counter(row["classification"] for row in rows).items()))
    summary = {
        "scope": "fixed-order LLSS->LLHH operator relevance, 16 supported ordinary T3 models",
        "strict_complete": strict,
        "models_checked": len(models),
        "extra_quartic_slots": len(new_q),
        "extra_quadratic_slots": len(new_m),
        "observed_additional_quartic_C5_slots": sum(r.get("observed_in_saved_C5", False) for r in new_q),
        "observed_models": sorted({r["model"] for r in new_q if r.get("observed_in_saved_C5")}),
        "classifications": status_count,
        "warnings": [
            "This is a graph-leg and saved-artifact classification, not a newly calculated C5 or a proof of generic zero couplings.",
            "Omitting a quartic from the specific direct one-loop LLSS->LLHH graph does NOT omit it from the full UV scalar potential or scalar RGEs.",
            "Higgs-bearing interactions not seen in saved matching are NOT proven irrelevant in another scheme, scalar basis or mass spectrum.",
            "Allowed quadratic coefficients affect tree-level scalar masses/mixing and cannot be dropped by one-loop power counting.",
            "One-loop RG generation followed by one-loop mass matching is formally higher order for zero-boundary couplings, but logarithmic enhancement and boundary-scale conventions must be checked.",
        ],
    }
    if strict:
        require(len(new_q) == 50 and len(new_m) == 8, "Expected 50 new quartic and 8 quadratic model slots")
        require(summary["observed_models"] == ["B_-1", "C_-1"] and
                summary["observed_additional_quartic_C5_slots"] == 4,
                "Changed saved C5 dependence: expected B/C(-1), two named couplings each")
    return {"summary": summary, "operators": rows}


def write_reports(result: dict[str, Any], hashes: dict[str, str], dest: Path) -> None:
    dest.mkdir(parents=True, exist_ok=True)
    result = {"input_sha256": hashes, **result}
    (dest / "operator_classification.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    with (dest / "operator_classification.csv").open("w", newline="", encoding="utf-8") as fh:
        columns = ["model", "kind", "coupling", "field_structure", "invariant_index",
                   "higgs_legs", "classification", "evidence_level", "physics_reason",
                   "omission_rule", "observed_in_saved_C5"]
        w = csv.DictWriter(fh, columns, extrasaction="ignore")
        w.writeheader()
        w.writerows(result["operators"])
    summary = result["summary"]
    lines = [
        "# T3 scalar relevance — Batch 1 source audit", "",
        f"- Models: {summary['models_checked']}",
        f"- Extra quartic slots: {summary['extra_quartic_slots']}",
        f"- Extra quadratic slots: {summary['extra_quadratic_slots']}",
        f"- Extra quartic slots observed in saved C5: {summary['observed_additional_quartic_C5_slots']}",
        "- Affected models: " + ", ".join(summary["observed_models"]), "",
        "## Classification counts", "",
        "| Classification | Rows |", "|---|---:|",
        *[f"| {k} | {v} |" for k, v in summary["classifications"].items()], "",
        "## Caveats", "",
        *[f"- {w}" for w in summary["warnings"]], "",
        "## Input hashes", "",
        *[f"- `{path}` SHA256 `{sha}`" for path, sha in sorted(hashes.items())], "",
        "See `operator_classification.csv` for a per-operator reason and omission criterion.",
    ]
    (dest / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT, help="Project repository root")
    parser.add_argument("--out-dir", type=Path, default=Path("output/scalar_relevance"))
    args = parser.parse_args(argv)
    root = args.root.resolve()
    dest = args.out_dir if args.out_dir.is_absolute() else root / args.out_dir
    try:
        sources, hashes = load_sources(root)
        result = classify_sources(sources, strict=True)
        write_reports(result, hashes, dest)
    except (AuditError, OSError) as exc:
        parser.exit(1, f"Scalar relevance audit FAILED (no verified report): {exc}\n")
    brief = result["summary"]
    print(f"Batch 1: {brief['models_checked']} models, "
          f"{brief['extra_quartic_slots']} extra quartic slots, "
          f"{brief['extra_quadratic_slots']} extra quadratic slots, "
          f"{brief['observed_additional_quartic_C5_slots']} directly observed extra C5 couplings")
    print(f"Reports: {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
