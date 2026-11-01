"""Large read-only physics/report batch for the complete 16-model T3 symbolic run.

Runs independent, exact one-loop *gauge* factor projections on native RGBeta UV
and EFT1 expressions for all allowed scalar quartics and bilinears, plus UV
Yukawas and gY/g2/g3. Audits the inventory using the existing provenanced
structural checker. Emits a coupling-basis ledger and consolidated reports.

Crucial scope limit: exact gauge-factor equality inside RGBeta is NOT proof
of the physical Matchete<->RGBeta coupling dictionary, including SU(2) tensor
intertwiners, signs, phases, conjugations, and Lagrangian-vs-potential factors.
Nor is it a check of all non-gauge scalar/Yukawa beta coefficients.

Run from repository root:
  python -X utf8 -m scripts.RunT3GeneralBatchVerification --self-test
  python -X utf8 -m scripts.RunT3GeneralBatchVerification --model A_-2
  python -X utf8 -m scripts.RunT3GeneralBatchVerification
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from fractions import Fraction
import json
from pathlib import Path
import re
from typing import Any

from scripts import AuditFull16GeneralSymbolicPhysics as audit
from scripts import RunT3GeneralSymbolicStudy as general
from scripts import RunT3CompleteGeneralScalarSymbolic as scalar
from scripts.T3GeneralGroupFactorProjection import (
    C2, T2, ProjectionError, expected_gauge_betas,
    expected_scalar_gauge_factors, expected_yukawa_gauge_factors,
    project_native_beta, _top_level_terms,
)
from Numerical.experimental.GeneralT3ScalarQuadraticMasses import REP, load_manifest, quadratic_kind

ROOT = Path(__file__).resolve().parents[1]
REPORT_NAME = "GeneralBatchVerification"
FORMAT = "t3-full16-exact-gauge-projection-and-physical-dictionary-ledger-v1"
STAGE_NAMES = ("UV", "EFT1")
GAUGES = ("gY", "g2", "g3")
BASE_SCALAR: dict[str, tuple[int, ...]] = {
    "lambdaH": (2, 2, 0, 0, 0, 0),
    "lambdaS1": (0, 0, 2, 2, 0, 0),
    "lambdaS2": (0, 0, 0, 0, 2, 2),
    "lambdaH1": (1, 1, 1, 1, 0, 0),
    "lambdaH2": (1, 1, 0, 0, 1, 1),
    "lambda12": (0, 0, 1, 1, 1, 1),
    "lambdaT3": (2, 0, 1, 0, 0, 1),
}
BASE_MASS: dict[str, tuple[int, ...]] = {
    "mS1Sq": (0, 0, 1, 1, 0, 0),
    "mS2Sq": (0, 0, 0, 0, 1, 1),
}
QUADRATIC_SYMBOLS = {
    "self_s1": "muT3GeneralS1Sq",
    "self_s2": "muT3GeneralS2Sq",
    "mixed": "muT3GeneralS1S2",
}
SCALAR_FIELDS = ("H", "Hdag", "S1", "S1dag", "S2", "S2dag")
STATUS_PASS = "VerifiedRGBetaGaugeCoefficient"
STATUS_MISMATCH = "GaugeCoefficientMismatch"
STATUS_UNPROJECTABLE = "NotProjected"
STATUS_UNCERTIFIED = "NotCertified"


def require(test: bool, msg: str) -> None:
    if not test:
        raise ValueError(msg)


def save_csv(path: Path, rows: list[dict[str, Any]], columns: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def model_numbers(model: str) -> tuple[tuple[int, int, int], int]:
    family, alpha = model.split("_")
    return REP[family], int(alpha)


def count_name(counts: tuple[int, ...]) -> str:
    return " ".join(f"{SCALAR_FIELDS[i]}^{n}" for i, n in enumerate(counts) if n) or "1"


def scalar_interactions(model: str, entry: dict, kind: str) -> list[dict]:
    """Provide every renormalisable scalar beta key and its exact field counts.

    Non-singlet baseline invariants share the same field multiplicities and
    hence the same universal gauge-linear coefficient. Neither a single CG
    normalisation nor the cross-invariant tensor basis is assumed here.
    """
    data: list[dict] = []
    for symbol, counts in BASE_SCALAR.items():
        data.append({"name": symbol, "counts": counts, "kind": "baseline_quartic"})
    for symbol, counts in BASE_MASS.items():
        data.append({"name": symbol, "counts": counts, "kind": "baseline_mass"})
    for item in entry["new_rgbeta_symbols"]:
        data.append({"name": item["name"], "counts": tuple(item["counts"]),
                     "kind": "additional_quartic_new_rgbeta"})
    for item in entry["preexisting_rgbeta_symbols"]:
        data.append({"name": item["existing_name"], "counts": tuple(item["counts"]),
                     "kind": "additional_quartic_existing_rgbeta"})
    # Standard representation-dependent independent SU(2) contractions.
    d1, d2, _ = REP[model[0]]
    dependent = [
        ("lambdaH1Adj", BASE_SCALAR["lambdaH1"], d1 > 1),
        ("lambdaH2Adj", BASE_SCALAR["lambdaH2"], d2 > 1),
        ("lambdaS1Adj", BASE_SCALAR["lambdaS1"], d1 == 3),
        ("lambdaS2Adj", BASE_SCALAR["lambdaS2"], d2 == 3),
        ("lambda12Adj", BASE_SCALAR["lambda12"], d1 > 1 and d2 > 1),
        ("lambda12Cross", BASE_SCALAR["lambda12"], d1 == d2 == 3),
    ]
    for symbol, counts, allowed in dependent:
        if allowed:
            data.append({"name": symbol, "counts": counts, "kind": "baseline_nonsinglet_quartic"})
    if kind != "none":
        mass_counts = {"self_s1": (0, 0, 0, 2, 0, 0),
                       "self_s2": (0, 0, 0, 0, 0, 2),
                       "mixed": (0, 0, 0, 1, 0, 1)}[kind]
        data.append({"name": QUADRATIC_SYMBOLS[kind], "counts": mass_counts,
                     "kind": "additional_bilinear"})
    # Other baseline adjoint couplings are not mapped here: their field
    # counts can be inferred from source, but this batch does not silently
    # claim a complete physical adjoint tensor dictionary.
    keys = [x["name"] for x in data]
    require(len(keys) == len(set(keys)), f"{model}: duplicate scalar beta key")
    return data


def check_one_coeff(model: str, stage: str, beta_key: str, raw: str,
                    gauge: str, expected: Fraction, section: str,
                    *, matrix: bool = False, field_counts: str = "") -> dict:
    base = {"model": model, "stage": stage, "sector": section,
            "beta": beta_key, "gauge": gauge,
            "field_counts": field_counts, "expected": str(expected),
            "extracted": "", "residual": "", "status": STATUS_UNPROJECTABLE,
            "detail": ""}
    try:
        actual = project_native_beta(raw, gauge=gauge,
                                     matrix_coupling=beta_key if matrix else None,
                                     coupling=None if section == "gauge" or matrix else beta_key,
                                     zero_is_valid=(expected == 0))
        base["extracted"] = str(actual)
        base["residual"] = str(actual - expected)
        base["status"] = STATUS_PASS if actual == expected else STATUS_MISMATCH
    except (ProjectionError, ValueError, TypeError, ArithmeticError) as exc:
        base["detail"] = str(exc)[:250]
    return base


def project_model(model: str, entry: dict, kind: str, folder: Path) -> tuple[list[dict], list[dict], list[dict]]:
    dims, alpha = model_numbers(model)
    scalars = scalar_interactions(model, entry, kind)
    scalar_rows: list[dict] = []
    coefficient_rows: list[dict] = []
    native_term_rows: list[dict] = []
    for stage in STAGE_NAMES:
        data = audit.read(folder / f"{stage.lower()}_rgbeta.json")
        beta = data["report_betas"]
        require(isinstance(beta, dict), f"{model}/{stage}: no native report_betas")
        # Full symbolic term index. This preserves native coefficients verbatim;
        # unlike the gauge projections these entries are NOT independently checked.
        for coupling, raw in beta.items():
            try:
                terms = _top_level_terms(str(raw))
                split_status = "TopLevelTermsIndexed"
            except ProjectionError:
                terms = [str(raw)]
                split_status = "RawExpressionOnly"
            for number, term in enumerate(terms, 1):
                native_term_rows.append({"model": model, "stage": stage,
                                         "beta": coupling, "term": number,
                                         "native_expression": term,
                                         "source_format": split_status,
                                         "coefficient_validation": STATUS_UNCERTIFIED})
        for gauge, expected in expected_gauge_betas(dims, alpha, stage).items():
            require(gauge in beta, f"{model}/{stage}: gauge beta {gauge} absent")
            coefficient_rows.append(check_one_coeff(model, stage, gauge, beta[gauge],
                                                     gauge, expected, "gauge"))
        for item in scalars:
            key = item["name"]
            require(key in beta, f"{model}/{stage}: expected scalar beta absent: {key}")
            pred = expected_scalar_gauge_factors(item["counts"], dims, alpha)
            scalar_rows.append({"model": model, "stage": stage, "beta": key,
                                "sector": item["kind"], "fields": count_name(item["counts"]),
                                "g2_factor_expected": str(pred["g2"]),
                                "gY_factor_expected": str(pred["gY"]),
                                "physical_coupling_dictionary": STATUS_UNCERTIFIED})
            for gauge, expected in pred.items():
                coefficient_rows.append(check_one_coeff(model, stage, key, beta[key],
                                                        gauge, expected, item["kind"],
                                                        field_counts=count_name(item["counts"])))
        if stage == "UV":
            y = expected_yukawa_gauge_factors(dims, alpha)
            ykeys = ["y1", "y2"]
            additional = general.expected_yukawa(model)
            if additional:
                ykeys.append(additional)
            for key in ykeys:
                require(key in beta, f"{model}: missing Yukawa beta {key}")
                for gauge, expected in y.items():
                    coefficient_rows.append(check_one_coeff(model, stage, key,
                                                            beta[key], gauge, expected,
                                                            "yukawa", matrix=True))
    return coefficient_rows, scalar_rows, native_term_rows


def coupling_ledger(model: str, entry: dict, manifest: dict, kind: str) -> list[dict]:
    """Inventory the cross-package operator mapping, without inventing factors."""
    dims, alpha = model_numbers(model)
    rows: list[dict] = []
    counts_sources: dict[tuple[int, ...], list[tuple[str, str]]] = {}
    for item in entry["new_rgbeta_symbols"]:
        key = tuple(item["counts"])
        counts_sources.setdefault(key, []).append((item["name"], "new"))
    for item in entry["preexisting_rgbeta_symbols"]:
        key = tuple(item["counts"])
        counts_sources.setdefault(key, []).append((item["existing_name"], "existing"))
    for family in manifest["models"][model]["extra_quartic_structures"]:
        counts = tuple(family["counts"])
        n = family["invariants"]
        symbols = counts_sources.get(counts, [])
        # Reused RGBeta symbols are *not* one-to-one with Matchete's
        # independent contraction labels, especially in B/C(-1).
        require(len(symbols) == n,
                f"{model}: missing or extra RGBeta symbols in quartic block {counts}")
        matchete_base = "lambdaT3ExtraQ" + "".join(str(x) for x in counts)
        for i in range(1, n + 1):
            rows.append({"model": model, "sector": "extra_quartic",
                         "field_counts": count_name(counts),
                         "matchete": f"{matchete_base}I{i}",
                         "rgbeta_candidates": ";".join(s[0] for s in symbols),
                         "basis_rank": n,
                         "correspondence": "BasisBlockUndetermined" if n > 1 else "CandidateOnly",
                         "relative_sign": STATUS_UNCERTIFIED,
                         "normalization": STATUS_UNCERTIFIED,
                         "complex_conjugation": STATUS_UNCERTIFIED,
                         "group_tensor_map": STATUS_UNCERTIFIED})
    for scalar_mass in manifest["models"][model]["quadratic_structures_not_enabled"]:
        counts = tuple(scalar_mass["counts"])
        rows.append({"model": model, "sector": "extra_bilinear",
                     "field_counts": count_name(counts),
                     "matchete": "lambdaT3ExtraM" + "".join(str(x) for x in counts) + "I1",
                     "rgbeta_candidates": QUADRATIC_SYMBOLS[kind], "basis_rank": 1,
                     "correspondence": "CandidateOnly", "relative_sign": STATUS_UNCERTIFIED,
                     "normalization": STATUS_UNCERTIFIED,
                     "complex_conjugation": STATUS_UNCERTIFIED,
                     "group_tensor_map": STATUS_UNCERTIFIED})
    for y in ("y1", "y2") + ((general.expected_yukawa(model),) if general.expected_yukawa(model) else ()):
        rows.append({"model": model, "sector": "yukawa",
                     "field_counts": ("Bar[L] F^c S1" if y == "y1" else
                                      "Bar[L] F S2dag" if y == "y2" else
                                      "Bar[L] F^c S1dag" if y == "y1Conj" else
                                      "Bar[L] F S2"),
                     "matchete": y, "rgbeta_candidates": y, "basis_rank": 1,
                     "correspondence": "NamedVertexOnly", "relative_sign": STATUS_UNCERTIFIED,
                     "normalization": STATUS_UNCERTIFIED,
                     "complex_conjugation": STATUS_UNCERTIFIED,
                     "group_tensor_map": STATUS_UNCERTIFIED})
    return rows


def report_model(model: str, report_dir: Path, coeff: list[dict], ledger: list[dict],
                 rows: list[dict], structural: dict, data_run: Path, study: str) -> None:
    directory = report_dir / "models"
    directory.mkdir(parents=True, exist_ok=True)
    label = general.model_tag(model)
    group = [r for r in coeff if r["model"] == model]
    good = sum(r["status"] == STATUS_PASS for r in group)
    bad = [r for r in group if r["status"] != STATUS_PASS]
    dims, alpha = model_numbers(model)
    ref = f"../../GeneralSymbolic/models/{label}.md"
    lines = [f"# T3 {model}: unified general-model physics verification", "",
             f"S1 = ({dims[0]},{Fraction(alpha,2)}), "
             f"S2 = ({dims[1]},{Fraction(alpha+2,2)}), "
             f"F = ({dims[2]},{Fraction(alpha+1,2)}).", "",
             f"**Source-matched structural inventory:** {structural['structural_status']}.  "
             f"**Gauge/group coefficients verified within RGBeta:** {good}/{len(group)}.", "",
             "## Complete general BSM Lagrangian and UV/EFT1 RGE expressions", "",
             f"[Full model report]({ref}) displays the native Matchete BSM expression "
             "and all exported RGBeta UV/EFT1 renormalisable beta functions.", "",
             "These are candidate RGBeta RGEs in the extended interaction basis; "
             "the physical cross-package coupling dictionary is not yet certified.", "",
             "## Gauge group factors", "",
             "The tables independently project the coefficients of "
             "`g^3` (gauge), `g^2*lambda` (scalar), and "
             "`g^2*Matrix[y]` (Yukawa). These checks **do not** test all "
             "non-gauge scalar recoupling or matrix-product terms.", "",
             "| Stage | Beta | Gauge | Expected | Extracted | Status |",
             "|---|---|---|---:|---:|---|"]
    for r in group:
        lines.append(f"| {r['stage']} | `{r['beta']}` | {r['gauge']} | {r['expected']} | "
                     f"{r['extracted'] or '—'} | {r['status']} |")
    if bad:
        lines.extend(["", "## Unresolved or discrepant projections", ""])
        for r in bad:
            lines.append(f"- {r['stage']} `beta_{r['beta']}` / {r['gauge']}: "
                         f"{r['status']} {r['detail'] or r['residual']}")
    lines.extend(["", "## Physical coupling correspondence", "",
                  "The coupling ledger records the scalar-invariant basis block "
                  "and additional Yukawa vertices. The relative Matchete/RGBeta "
                  "action sign, conjugation, phase and normalisation remain **NotCertified**.", "",
                  "| Sector | Native Matchete coupling | RGBeta candidate(s) | Tensor rank |",
                  "|---|---|---|---:|"])
    for row in ledger:
        if row["model"] == model:
            lines.append(f"| {row['sector']} | `{row['matchete']}` | "
                         f"`{row['rgbeta_candidates']}` | {row['basis_rank']} |")
    lines.extend(["", "## EFT scope and open physics", "",
                  "- UV and EFT1 **renormalisable** RGEs are exported; EFT1 dimension-five Wilson-operator evolution is not rederived here.",
                  "- Full extended matching and final C5 are not certified. B(-1)/C(-1) require scalar mass diagonalisation.",
                  "- The separate legacy restricted-model matching must not be labelled as matching the general BSM Lagrangian.", ""])
    (directory / f"{label}.md").write_text("\n".join(lines), encoding="utf-8")


def write_output(destination: Path, models: list[str], entries: list[dict],
                 coefficient_rows: list[dict], scalar_rows: list[dict],
                 ledger_rows: list[dict], native_term_rows: list[dict], evidence: list[dict],
                 run_dir: Path, study: str) -> dict:
    destination.mkdir(parents=True, exist_ok=True)
    save_csv(destination / "group_factor_coefficients.csv", coefficient_rows,
             ["model", "stage", "sector", "beta", "gauge", "field_counts",
              "expected", "extracted", "residual", "status", "detail"])
    save_csv(destination / "scalar_gauge_structure.csv", scalar_rows,
             ["model", "stage", "beta", "sector", "fields", "g2_factor_expected",
              "gY_factor_expected", "physical_coupling_dictionary"])
    save_csv(destination / "physical_coupling_dictionary.csv", ledger_rows,
             ["model", "sector", "field_counts", "matchete", "rgbeta_candidates",
              "basis_rank", "correspondence", "relative_sign", "normalization",
              "complex_conjugation", "group_tensor_map"])
    save_csv(destination / "native_beta_term_index.csv", native_term_rows,
             ["model", "stage", "beta", "term", "native_expression",
              "source_format", "coefficient_validation"])
    save_csv(destination / "model_status.csv", entries,
             ["model", "structural_status", "gauge_projections", "passed", "mismatched",
              "unprojected", "extra_quartics", "extra_bilinears", "extra_yukawas", "error"])
    counts = Counter(x["status"] for x in coefficient_rows)
    all_pass = (len(entries) == len(models) and all("error" not in x for x in entries)
                and not counts[STATUS_MISMATCH] and not counts[STATUS_UNPROJECTABLE])
    summary = {
        "format": FORMAT, "model_count": len(models), "requested_models": models,
        "source_provenance_checked": all("error" not in x for x in entries),
        "structural_audit_reused": True,
        "projection_scope": "exact universal one-loop gauge/gauge-linear scalar and Yukawa coefficients only",
        "native_beta_additive_terms_indexed": len(native_term_rows),
        "gauge_factor_coefficients_total": len(coefficient_rows),
        "gauge_factor_coefficients_verified": counts[STATUS_PASS],
        "gauge_factor_coefficients_mismatched": counts[STATUS_MISMATCH],
        "gauge_factor_coefficients_unprojectable": counts[STATUS_UNPROJECTABLE],
        "all_projected_gauge_group_factors_equal": all_pass,
        "physical_matchete_rgbeta_dictionary_certified": False,
        "all_non_gauge_rge_coefficients_independently_verified": False,
        "extended_wilson_rge_and_c5_matching_certified": False,
        "model_status": entries,
        "prior_limited_certificates": evidence,
        "input_run_directory": str(run_dir),
    }
    (destination / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    for model in models:
        structural = next((r for r in entries if r["model"] == model), {})
        if "error" not in structural:
            report_model(model, destination, coefficient_rows, ledger_rows,
                         scalar_rows, structural, run_dir, study)
    lines = ["# T3: full general-interaction 16-model physics and group-factor batch", "",
             f"**Structural outputs checked:** {sum('error' not in r for r in entries)}/{len(models)}.",
             f"**Exact RGBeta gauge-group coefficients verified:** {counts[STATUS_PASS]}/{len(coefficient_rows)}. "
             f"Mismatched: {counts[STATUS_MISMATCH]}, unprojectable: {counts[STATUS_UNPROJECTABLE]}.", "",
             "This is a **read-only** comparison of the already saved and source-fingerprinted "
             "general BSM Lagrangians and native RGBeta UV/EFT1 beta functions.", "",
             "## Deliverables", "",
             "- `group_factor_coefficients.csv`: exact gauge and universal gauge-linear factors extracted from native beta functions, with independent predictions and residuals.",
             "- `scalar_gauge_structure.csv`: SU(2) and U(1) group-factor predictions for each tested scalar coupling.",
             "- `physical_coupling_dictionary.csv`: 50 extra quartic coupling slots, eight extra bilinears and 38 Yukawa slots, with the actual RGBeta basis block(s) and uncertified physical maps.",
             "- `native_beta_term_index.csv`: **all** native UV/EFT1 beta expressions split into additive terms for inspection. Non-gauge coefficients are indexed, **not independently certified**.",
             "- `model_status.csv`, `summary.json`: machine-readable status and source-qualified audit outcomes.",
             "- `models/`: consolidated per-model overview and links to complete native BSM Lagrangians and all UV/EFT1 RGEs.", "",
             "## Verified versus outstanding", "",
             "**Verified by this batch where status says VerifiedRGBetaGaugeCoefficient:** "
             "the native RGBeta one-loop coefficient of the indicated gauge monomial "
             "agrees exactly with its independent representation formula.", "",
             "**NOT certified:** the Matchete-to-RGBeta physical action dictionary, "
             "every non-gauge Yukawa/quartic group-factor coefficient, "
             "dimension-five Wilson RGEs, mass-eigenbasis threshold matching, and C5.", "",
             "The legacy `Reports/GroupFactors` output uses analytic formulas for the "
             "restricted operator basis and cannot automatically certify "
             "all newly allowed scalar invariants.", "",
             "## Per-model status", "",
             "| Model | Structural | Gauge checks | Matched | Unprojected | Details |",
             "|---|---|---:|---:|---:|---|"]
    for entry in entries:
        name = entry["model"]
        if "error" in entry:
            lines.append(f"| {name} | FAILED | — | — | — | `{entry['error']}` |")
        else:
            lines.append(f"| {name} | {entry['structural_status']} | {entry['gauge_projections']} | "
                         f"{entry['passed']} | {entry['unprojected']} | "
                         f"[Report](models/{general.model_tag(name)}.md) |")
    lines += ["", "## Next physics checks", "",
              "1. Reconstruct the physical scalar action for each independent RGBeta quartic tensor, including factorials, H.c., signs, and field-basis conventions, and compare to Matchete component amplitudes.",
              "2. Derive the full signed/complex Yukawa action dictionary from the saved 38 Yukawa CG witnesses, without identifying a CG norm ratio with a physical coupling map.",
              "3. After a certified dictionary exists, compare *all* beta monomial coefficients and repair genuine discrepancies (if found).", ""]
    (destination / "README.md").write_text("\n".join(lines), encoding="utf-8")
    return summary


def self_test() -> None:
    assert len(audit.general_models()) == 16
    assert expected_gauge_betas((1, 3, 2), -2, "UV") == {
        "gY": Fraction(55, 6), "g2": Fraction(-1, 2), "g3": -7}
    assert expected_gauge_betas((2, 2, 1), -1, "UV")["g2"] == Fraction(-17, 6)
    assert expected_gauge_betas((2, 2, 3), -1, "UV")["g2"] == Fraction(7, 6)
    assert expected_scalar_gauge_factors((1, 1, 1, 1, 0, 0), (1, 3, 2), -2) == {
        "g2": Fraction(-9, 2), "gY": Fraction(-15, 2)}
    assert project_native_beta("41*gY^3/6", gauge="gY") == Fraction(41, 6)
    assert project_native_beta("-3*g2^2*lambdaH1/2", gauge="g2", coupling="lambdaH1") == Fraction(-3, 2)
    assert project_native_beta("-9*g2^2*Matrix[y1][$i,$j]/4", gauge="g2", matrix_coupling="y1") == Fraction(-9, 4)
    print("PASS: catalogue, exact UV gauge factors, scalar factors, native-polynomial projections")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", help="One model (e.g. A_-2); default 16")
    parser.add_argument("--run-id", default="full16_general_v1")
    parser.add_argument("--study", default="full16_general_symbolic")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args(argv)
    if args.self_test:
        self_test()
        return 0
    for val in (args.run_id, args.study):
        require(re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]{0,63}", val) is not None,
                f"Invalid study or run ID: {val}")
    models = audit.general_models()
    if args.model:
        require(args.model in models, f"Unknown model: {args.model}")
        models = [args.model]
    registry = audit.read(scalar.REGISTRY)
    manifest = load_manifest(scalar.MANIFEST)
    basis = audit.quadratic_audit()
    scalar.check_sources(models, registry, manifest, basis)
    hashes = general.fingerprint()
    run_dir = general.DEFAULT_DATA_ROOT / args.run_id
    group_csv = ROOT / "Reports/output" / args.study / "GeneralSymbolic/group_invariants.csv"
    destination = ROOT / "Reports/output" / args.study / REPORT_NAME
    coefficient_rows: list[dict] = []
    scalar_rows: list[dict] = []
    ledger_rows: list[dict] = []
    native_term_rows: list[dict] = []
    model_rows: list[dict] = []
    for idx, model in enumerate(models, 1):
        try:
            entry = registry["models"][model]
            kind = quadratic_kind(manifest, model)
            structural, _, _ = audit.audit_model(model, registry, manifest, hashes,
                                                 run_dir, group_csv)
            coeff, sr, terms = project_model(model, entry, kind,
                                      run_dir / general.model_tag(model))
            dictionary = coupling_ledger(model, entry, manifest, kind)
            coefficient_rows.extend(coeff)
            scalar_rows.extend(sr)
            ledger_rows.extend(dictionary)
            native_term_rows.extend(terms)
            counts = Counter(r["status"] for r in coeff)
            result = {"model": model, "structural_status": structural["status"],
                      "gauge_projections": len(coeff), "passed": counts[STATUS_PASS],
                      "mismatched": counts[STATUS_MISMATCH],
                      "unprojected": counts[STATUS_UNPROJECTABLE],
                      "extra_quartics": structural["extra_quartics"],
                      "extra_bilinears": structural["extra_quadratics"],
                      "extra_yukawas": structural["extra_yukawas"]}
        except (OSError, KeyError, ValueError, TypeError, ProjectionError, json.JSONDecodeError) as exc:
            result = {"model": model, "error": str(exc),
                      "structural_status": "Failed"}
        model_rows.append(result)
        print(f"[{idx}/{len(models)}] {model}: "
              + (f"{result['passed']}/{result['gauge_projections']} gauge/group-factor checks"
                 if "error" not in result else f"FAILED: {result['error']}"), flush=True)
    summary = write_output(destination, models, model_rows, coefficient_rows,
                           scalar_rows, ledger_rows, native_term_rows,
                           audit.inspect_saved_bridge_reports(), run_dir, args.study)
    print(f"Gauge coefficient checks: {summary['gauge_factor_coefficients_verified']}/"
          f"{summary['gauge_factor_coefficients_total']}, "
          f"unprojectable={summary['gauge_factor_coefficients_unprojectable']}, "
          f"mismatched={summary['gauge_factor_coefficients_mismatched']}", flush=True)
    print(f"Report: {destination / 'README.md'}", flush=True)
    print("Physical Matchete/RGBeta dictionary + non-gauge beta coefficients: NOT CERTIFIED", flush=True)
    if any("error" in r for r in model_rows) or summary["gauge_factor_coefficients_mismatched"]:
        return 1
    return 2 if summary["gauge_factor_coefficients_unprojectable"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
