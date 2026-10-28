"""Batch 2: necessary single-quartic-vertex Wick-contraction audit of T3 C5.

This proves a *necessary* field-flow condition in the unbroken SU(2)xU(1),
Z2-odd, diagonal complex-scalar-propagator basis. It does NOT perform Matchete
matching, claim complete one-loop matching, or include nonzero holomorphic
quadratic masses/mixed scalar propagators.

Source Yukawa conventions (Lagrangian/interactions/T3Topology.wl, a02e04c6):
    y1 * Bar[l] * CConj[F] * S1
    y2 * Bar[l] * F        * S2dag
Two lepton-type vertices therefore source S1 and/or S2dag. With diagonal
complex scalar propagators, the quartic end must be the conjugate:
    S1dag and/or S2.
The holomorphic Weinberg projection is calibrated to the h.c. of the project's
HH S1 S2dag topology vertex, i.e. Hdag Hdag S1dag S2.

Only one insertion of a quartic vertex and the two original Yukawa vertices
is tested. Exclusion from this graph is NOT exclusion from a complete scalar
potential or arbitrary scalar-mixing basis.
"""
from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
EXPECTED_MODELS = 16
EXPECTED_EXTRA_QUARTICS = 50
EXPECTED_EXTRA_QUADRATICS = 8
EXPECTED_HIGGS_EXTRAS = 24
EXPECTED_OBSERVED = 4
BASIS_SOURCE = "Lagrangian/interactions/T3Topology.wl @ a02e04c6def5cf8522b2f612ccd6bff7d9963b91"
COUNTS_KEY = "counts_H_Hdag_S1_S1dag_S2_S2dag"

# Fields required at the quartic end to contract two original Yukawa scalar
# ends, using only <S1 S1dag> and <S2 S2dag> propagators.
SCALAR_ENDS = {
    (0, 2, 0, 0): "y1_y1",
    (0, 1, 1, 0): "y1_y2",
    (0, 0, 2, 0): "y2_y2",
}


class AuditError(ValueError):
    """Fail rather than issuing a misleading certified relevance statement."""


def conjugate_counts(counts: tuple[int, ...]) -> tuple[int, ...]:
    if len(counts) != 6:
        raise AuditError(f"expected six bosonic counts, received {counts!r}")
    return tuple(counts[i] for i in (1, 0, 3, 2, 5, 4))


def canonical_vertex_test(counts: tuple[int, ...]) -> dict[str, Any]:
    """Test both Hermitian orientations against the same lepton-side source."""
    if len(counts) != 6 or any(not isinstance(x, int) or x < 0 for x in counts):
        raise AuditError(f"invalid six-species counts: {counts}")
    if sum(counts) != 4 or sum(counts[:2]) != 2:
        raise AuditError(f"not a two-Higgs/two-BSM quartic: {counts}")

    orientations = (("as_named", counts), ("hermitian_conjugate", conjugate_counts(counts)))
    checks: list[dict[str, Any]] = []
    for label, oriented in orientations:
        h, hdag, s1, s1dag, s2, s2dag = oriented
        higgs_ok = h == 0 and hdag == 2
        ends = (s1, s1dag, s2, s2dag)
        channel = SCALAR_ENDS.get(ends) if higgs_ok else None
        checks.append({
            "orientation": label,
            "higgs_pair_correct": higgs_ok,
            "scalar_pair_compatible": ends in SCALAR_ENDS,
            "yukawa_pair_channel": channel,
        })
    matched = [c for c in checks if c["yukawa_pair_channel"]]
    if len(matched) > 1:
        raise AuditError(f"unexpected duplicated holomorphic orientation: {counts}")
    if matched:
        status = "single_vertex_wick_eligible"
        reason = "Exactly two Hdag and the conjugate of two original Yukawa scalar ends; this is a necessary, not sufficient, nonzero-amplitude condition."
    elif not any(c["higgs_pair_correct"] for c in checks):
        status = "higgs_chirality_mismatch"
        reason = "Both orientations contain H Hdag, not the Hdag Hdag vertex used by the lepton-side holomorphic Weinberg projection."
    else:
        status = "scalar_wick_orientation_mismatch"
        reason = "The Hdag Hdag orientation does not supply two scalar legs conjugate to S1/S2dag Yukawa sources with diagonal scalar propagators."
    return {"status": status, "reason": reason,
            "channel": matched[0]["yukawa_pair_channel"] if matched else None,
            "orientation": matched[0]["orientation"] if matched else None,
            "orientations_tested": checks}


def validate_batch1(report: dict[str, Any]) -> None:
    summ = report.get("summary", {})
    ops = report.get("operators")
    if not isinstance(ops, list) or not summ.get("strict_complete"):
        raise AuditError("input must be a strict-complete Batch 1 operator report")
    n_models = {o.get("model") for o in ops}
    if len(n_models) != EXPECTED_MODELS or summ.get("models_checked") != EXPECTED_MODELS:
        raise AuditError("Batch 1 report does not have exactly 16 models")
    kinds = collections.Counter(o.get("kind") for o in ops)
    expected = {"baseline_topology": 16, "baseline_quadratic": 32,
                "baseline_quartic_family": 80,
                "additional_quartic": 50, "additional_quadratic": 8}
    if kinds != expected:
        raise AuditError(f"unexpected Batch 1 kind inventory: {kinds}")
    if (summ.get("extra_quartic_slots") != EXPECTED_EXTRA_QUARTICS or
            summ.get("extra_quadratic_slots") != EXPECTED_EXTRA_QUADRATICS or
            summ.get("observed_additional_quartic_C5_slots") != EXPECTED_OBSERVED):
        raise AuditError("Batch 1 headline totals have changed")
    uniq = set()
    for o in ops:
        label = (o["model"], o["kind"], o["coupling"])
        if label in uniq:
            raise AuditError(f"duplicate model/kind/coupling: {label}")
        uniq.add(label)
        if o["kind"] == "additional_quartic":
            counts = o.get(COUNTS_KEY)
            if not isinstance(counts, list) or len(counts) != 6:
                raise AuditError(f"missing counts for {label}")
            # The symbolic coupling suffix is a six-digit field-count signature.
            expected_prefix = "lambdaT3ExtraQ" + "".join(str(x) for x in counts) + "I"
            if not o["coupling"].startswith(expected_prefix):
                raise AuditError(f"coupling/field-count mismatch: {label}")
        if o["kind"] == "additional_quadratic" and sum(o.get(COUNTS_KEY) or []) != 2:
            raise AuditError(f"non-bilinear additional quadratic: {label}")
    for model in n_models:
        own = collections.Counter(o["kind"] for o in ops if o["model"] == model)
        for k, v in (("baseline_topology", 1), ("baseline_quadratic", 2),
                     ("baseline_quartic_family", 5)):
            if own[k] != v:
                raise AuditError(f"{model}: unexpected {k} count {own[k]}")


def verify_source_hashes(report: dict[str, Any], root: Path) -> None:
    hashes = report.get("input_sha256", {})
    if len(hashes) < 5:
        raise AuditError("missing Batch 1 input provenance hashes")
    for path_string, expected_hash in hashes.items():
        path = root / path_string
        if not path.is_file():
            raise AuditError(f"missing Batch 1 evidence file: {path}")
        observed_hash = hashlib.sha256(path.read_bytes()).hexdigest()
        if observed_hash != expected_hash:
            raise AuditError(f"stale Batch 1 evidence: {path_string}; rerun Batch 1 first")
    topology_path = root / "Lagrangian/interactions/T3Topology.wl"
    if not topology_path.is_file():
        raise AuditError(f"missing source Yukawa convention: {topology_path}")
    topology = topology_path.read_text(encoding="utf-8")
    signatures = (
        "useCConj = which === 1;",
        "scalarBar = which === 2;",
        "BuildT3MixingCandidates",
        "s1 = If[d1 === 1, NewScalar1[]",
        "s2 = If[d2 === 1, Bar[NewScalar2[]]",
    )
    for snippet in signatures:
        if snippet not in topology:
            raise AuditError(f"Yukawa/topology source contract changed: {snippet}; review assumptions")


def analyze_batch1(report: dict[str, Any], *, source_hashes_verified: bool) -> dict[str, Any]:
    validate_batch1(report)
    out: list[dict[str, Any]] = []
    hits = collections.Counter()
    detail = collections.Counter()
    per_model: dict[str, dict[str, int]] = {}
    for o in report["operators"]:
        model = o["model"]
        m = per_model.setdefault(model, {"new_two_H_vertices": 0,
                                         "wick_eligible": 0,
                                         "chirality_mismatch": 0,
                                         "scalar_orientation_mismatch": 0,
                                         "observed_in_saved_C5": 0})
        row = {
            "model": model, "kind": o["kind"], "coupling": o["coupling"],
            "field_structure": o["field_structure"],
            "source_Batch1_classification": o["classification"],
            "observed_in_saved_C5": o.get("observed_in_saved_C5"),
            "counts_H_Hdag_S1_S1dag_S2_S2dag": o.get(COUNTS_KEY),
        }
        k = o["kind"]
        if k == "baseline_topology":
            test = canonical_vertex_test(tuple(o[COUNTS_KEY]))
            if test["status"] != "single_vertex_wick_eligible" or test["channel"] != "y1_y2":
                raise AuditError(f"baseline T3 calibrator unexpectedly fails in {model}")
            row.update({"status": "baseline_topology_reproduced", "channel": test["channel"],
                        "basis_assumption": "diagonal_complex_scalar_propagators",
                        "reason": "Hermitian conjugate of HH S1 S2dag saturates the Yukawa scalar ends."})
        elif k == "additional_quartic":
            hlegs = sum(o[COUNTS_KEY][:2])
            if hlegs == 2:
                m["new_two_H_vertices"] += 1
                test = canonical_vertex_test(tuple(o[COUNTS_KEY]))
                row.update(test)
                row["basis_assumption"] = "diagonal_complex_scalar_propagators"
                if test["status"] == "single_vertex_wick_eligible":
                    m["wick_eligible"] += 1
                    if o.get("observed_in_saved_C5"):
                        m["observed_in_saved_C5"] += 1
                        row["status"] = "single_vertex_wick_eligible_and_saved_C5_observed"
                    else:
                        row["status"] = "single_vertex_wick_eligible_not_verified"
                elif test["status"] == "higgs_chirality_mismatch":
                    m["chirality_mismatch"] += 1
                elif test["status"] == "scalar_wick_orientation_mismatch":
                    m["scalar_orientation_mismatch"] += 1
                else:
                    raise AuditError(f"unsupported result: {test}")
                if o.get("observed_in_saved_C5") and test["status"] != "single_vertex_wick_eligible":
                    raise AuditError(f"saved C5 contradicts Wick prerequisite: {model}/{o['coupling']}")
                hits[row["status"]] += 1
            elif hlegs == 0:
                row.update({"status": "no_H_fields_at_one_vertex",
                            "reason": "Cannot furnish the two external Higgs fields with a single quartic vertex; indirect effects remain."})
            else:
                raise AuditError(f"unexpected quartic Higgs leg count: {model}/{o['coupling']}")
        elif k == "baseline_quartic_family":
            if o["classification"] == "one_loop_higgs_vertex_candidate":
                row.update({"status": "baseline_portal_higgs_chirality_mismatch",
                            "reason": "The Hdag H portal cannot furnish the Hdag Hdag vertex of the unbroken holomorphic Weinberg projection alone; v/M mass shifts and RG effects remain."})
            else:
                row.update({"status": "baseline_no_direct_two_H_vertex",
                            "reason": "No Higgs pair for the selected single-quartic vertex; scalar dynamics/RG effects remain."})
        elif k == "additional_quadratic":
            row.update({"status": "unsuppressed_tree_level_mass_effect",
                        "reason": "Can generate anomalous/off-diagonal propagators or physical mass splittings at tree level. Nonzero coupling INVALIDATES this Wick-screening assumption."})
        elif k == "baseline_quadratic":
            row.update({"status": "required_diagonal_mass_input",
                        "reason": "Determines the propagators used by the test."})
        else:
            raise AuditError(f"unexpected kind: {k}")
        out.append(row)
        detail[row["status"]] += 1
    if (sum(m["new_two_H_vertices"] for m in per_model.values()) != EXPECTED_HIGGS_EXTRAS or
            sum(m["wick_eligible"] for m in per_model.values()) != EXPECTED_OBSERVED or
            sum(m["observed_in_saved_C5"] for m in per_model.values()) != EXPECTED_OBSERVED or
            sum(m["chirality_mismatch"] for m in per_model.values()) != 10 or
            sum(m["scalar_orientation_mismatch"] for m in per_model.values()) != 10):
        raise AuditError(f"new-vertex totals do not reproduce the known report: {per_model}")
    if sorted(m for m, v in per_model.items() if v["wick_eligible"]) != ["B_-1", "C_-1"]:
        raise AuditError("unexpected one-loop extra-quartic models")
    return {
        "format": "t3-one-loop-scalar-single-quartic-wick-audit-v1",
        "scope": "unbroken electroweak, LLSS -> LLHH / holomorphic Weinberg projection at fixed one-loop order",
        "status": "WickNecessaryConditionsReproduceSavedC5Inventory",
        "new_matching_calculation_performed": False,
        "source_provenance_verified_on_local_checkout": source_hashes_verified,
        "source_basis_convention": BASIS_SOURCE,
        "batch1_report_sha256": None,
        "assumptions": [
            "Use both lepton-side Yukawa vertices as in T3Topology.wl: y1 couples S1, y2 couples S2dag.",
            "Free scalar propagators connect S1 to S1dag and S2 to S2dag only; all eight additional holomorphic quadratic coefficients are zero for this test.",
            "A single quartic with two Higgs legs is inserted, and both external Higgs fields follow the conjugate orientation of the benchmark HH S1 S2dag topology coupling.",
            "Work in the unbroken gauge phase and at leading dimension-five, fixed-order one-loop matching; no additional Higgs VEV insertions or quartic vertices.",
            "Fermion chirality/Majorana contractions, nonzero CG coefficients and coefficient normalization are NOT proven by Wick-leg counting; saved C5 observations are a separate check.",
        ],
        "counts": {
            "models_checked": EXPECTED_MODELS,
            "additional_quartics_total": EXPECTED_EXTRA_QUARTICS,
            "new_two_H_vertices": EXPECTED_HIGGS_EXTRAS,
            "new_direct_wick_eligible": EXPECTED_OBSERVED,
            "saved_C5_observed_new": EXPECTED_OBSERVED,
            "excluded_H_Hdag_by_higgs_orientation": 10,
            "excluded_Hdag2_by_scalar_orientation": 10,
            "new_no_H_vertex": 26,
            "extra_quadratics_not_screened": EXPECTED_EXTRA_QUADRATICS,
            "baseline_portal_families_H_Hdag": 32,
        },
        "classification_counts": dict(sorted(detail.items())),
        "per_model": dict(sorted(per_model.items())),
        "operators": out,
        "limitations": [
            "This is a necessary contraction test in the DIAGONAL scalar basis, NOT a proof of irrelevance with nonzero extra quadratic mixing.",
            "Fermion contractions and group tensors can make a Wick-eligible vertex vanish; the four observed terms are supported separately by saved symbolic Matchete matching.",
            "Non-observed 20 are excluded only from the defined direct one-quartic-vertex amplitude, not from the general UV scalar potential.",
            "Higher-loop, RG, v^2/M^2, vacuum stability, and general two-threshold Wilson matching are not addressed here.",
        ],
    }


def markdown_summary(result: dict[str, Any]) -> str:
    c = result["counts"]
    lines = ["# T3 scalar relevance — Batch 2: Wick-connection audit", "",
        f"- Models checked: **{c['models_checked']}**",
        f"- New Higgs-bearing quartic slots: **{c['new_two_H_vertices']}**",
        f"- Wick-eligible and present in the saved C5: **{c['new_direct_wick_eligible']}** (two each in B_-1/C_-1)",
        f"- Excluded from the specific single-vertex graph: **20** (10 H Hdag, 10 scalar-orientation mismatch)",
        f"- Four-BSM/no-H quartics: **{c['new_no_H_vertex']}**; no direct single-quartic insertion",
        f"- Additional quadratic mass parameters: **{c['extra_quadratics_not_screened']}**, not excluded", "",
        "## Essential physical derivation", "",
        "Original lepton Yukawas source `S1` and `S2dag`. With diagonal complex scalar propagators, the Higgs-bearing quartic must supply `S1dag` or `S2`. A single quartic feeding the matched holomorphic Weinberg sector must supply `Hdag^2`, using the conjugate orientation of the baseline `H^2 S1 S2dag` coupling. Therefore the only scalar field pairs admitted at the quartic end are `S1dag^2` (`y1y1`), `S1dag S2` (`y1y2`) and `S2^2` (`y2y2`). Both hermitian orientations of every quartic were tested.", "",
        "## Model breakdown", "",
        "| Model | New 2-H quartics | Wick-eligible | H Hdag mismatch | Scalar mismatch |",
        "|---|---:|---:|---:|---:|",]
    for model, m in result["per_model"].items():
        lines.append(f"| {model} | {m['new_two_H_vertices']} | {m['wick_eligible']} | {m['chirality_mismatch']} | {m['scalar_orientation_mismatch']} |")
    lines += ["", "## Direct new vertices in the saved one-loop C5", ""]
    for o in result["operators"]:
        if o["status"] == "single_vertex_wick_eligible_and_saved_C5_observed":
            lines.append(f"- `{o['model']}` — `{o['coupling']}`: `{o['field_structure']}` ({o['channel']})")
    lines += ["", "## Limitations", ""]
    lines.extend("- " + x for x in result["limitations"])
    lines += ["", f"- Local Batch 1 source hash verification: **{'yes' if result['source_provenance_verified_on_local_checkout'] else 'not performed (uploaded report only)'}**", ""]
    return "\n".join(lines)


def write_result(result: dict[str, Any], directory: Path) -> list[Path]:
    directory.mkdir(parents=True, exist_ok=True)
    js = directory / "leading_c5_vertex_audit.json"
    csv_file = directory / "leading_c5_vertex_audit.csv"
    md = directory / "leading_c5_vertex_summary.md"
    js.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    keys = ["model", "kind", "coupling", "field_structure", "status", "channel",
            "reason", "observed_in_saved_C5"]
    with csv_file.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=keys)
        writer.writeheader()
        for o in result["operators"]:
            writer.writerow({k: o.get(k) for k in keys})
    md.write_text(markdown_summary(result), encoding="utf-8")
    return [js, csv_file, md]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path("output/scalar_relevance/operator_classification.json"))
    parser.add_argument("--output-dir", type=Path, default=Path("output/scalar_relevance"))
    parser.add_argument("--skip-source-hash-check", action="store_true",
                        help="for offline uploaded-report reanalysis only; report is marked as not source-verified")
    args = parser.parse_args()
    try:
        input_path = args.input if args.input.is_absolute() else ROOT / args.input
        output_dir = args.output_dir if args.output_dir.is_absolute() else ROOT / args.output_dir
        report = json.loads(input_path.read_text(encoding="utf-8"))
        validate_batch1(report)
        if not args.skip_source_hash_check:
            verify_source_hashes(report, ROOT)
        result = analyze_batch1(report, source_hashes_verified=not args.skip_source_hash_check)
        result["batch1_report_sha256"] = hashlib.sha256(input_path.read_bytes()).hexdigest()
        paths = write_result(result, output_dir)
        print("Batch 2 audit succeeded:", result["counts"])
        print("Source evidence: " + ("hash-verified" if not args.skip_source_hash_check else "NOT hash-verified"))
        for path in paths:
            print(path)
        return 0
    except (AuditError, OSError, ValueError, KeyError, TypeError) as error:
        parser.exit(status=1, message=f"Batch 2 audit FAILED: {error}\n")


if __name__ == "__main__":
    raise SystemExit(main())
