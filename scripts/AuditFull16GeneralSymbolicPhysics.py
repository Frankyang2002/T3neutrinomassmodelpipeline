"""Read-only physics/inventory audit of the saved 16-model general-symbolic study.

No Wolfram invocation and no production-code changes. This audits the SOURCE-MATCHED
Matchete/RGBeta exports, independent charge/Casimir arithmetic, declared SU(2)
operator multiplicities, and actual beta-function coverage. It does NOT equate
Matchete and RGBeta physical coupling normalizations or certify matching.

Run from repository root:
  python -X utf8 -m scripts.AuditFull16GeneralSymbolicPhysics --self-test
  python -X utf8 -m scripts.AuditFull16GeneralSymbolicPhysics --model A_-2
  python -X utf8 -m scripts.AuditFull16GeneralSymbolicPhysics

Audit output: Reports/output/<study>/GeneralPhysicsAudit/.
"""
from __future__ import annotations

import argparse
import csv
from fractions import Fraction
import json
from pathlib import Path
import re
from typing import Any

from scripts import RunT3GeneralSymbolicStudy as general
from scripts import RunT3CompleteGeneralScalarSymbolic as scalar
from scripts.AuditGeneralT3QuadraticBasis import audit as quadratic_audit
from Numerical.experimental.GeneralT3ScalarQuadraticMasses import REP, load_manifest, quadratic_kind

ROOT = Path(__file__).resolve().parents[1]
MODEL_CLASSES: dict[str, tuple[int, ...]] = {
    "A": (-4, -2, 0), "B": (-3, -1, 1), "C": (-3, -1, 1),
    "D": (-2, 0, 2), "E": (-4, -2, 0, 2),
}
EXPECTED_TOTALS = {"models": 16, "extra_quartics": 50,
                   "extra_quadratics": 8, "extra_yukawas": 6,
                   "rgbeta_new_quartics": 38, "rgbeta_reused_quartics": 12}
EXPECTED_RGE_BASE = {"gY", "g2", "g3", "yu", "yd", "ye", "mS1Sq",
                     "mS2Sq", "lambdaH", "lambdaS1", "lambdaS2", "lambdaH1",
                     "lambdaH2", "lambda12", "lambdaT3"}
FAILURE = "Failed"
PENDING = "NotCertified"
PASS = "VerifiedStructural"
FORMAT = "t3-general-symbolic-independent-physics-audit-v1"


def require(condition: bool, detail: str) -> None:
    if not condition:
        raise ValueError(detail)


def read(path: Path) -> dict[str, Any]:
    obj = json.loads(path.read_text(encoding="utf-8-sig"))
    require(isinstance(obj, dict), f"JSON object required: {path}")
    return obj


def tag(model: str) -> str:
    return model.replace("+", "p").replace("-", "m")


def exact_group_invariants(dim: int) -> tuple[Fraction, Fraction]:
    require(isinstance(dim, int) and dim >= 1, "Nonpositive SU(2) dimension")
    return (Fraction(dim * dim - 1, 4),
            Fraction(dim * (dim * dim - 1), 12))


def scalar_charges2(alpha: int) -> tuple[int, ...]:
    """Twice the U(1)_Y charge, order H,Hdag,S1,S1dag,S2,S2dag."""
    return (1, -1, alpha, -alpha, alpha + 2, -(alpha + 2))


def structure_key(counts: list[int]) -> str:
    require(len(counts) == 6 and all(type(x) is int and 0 <= x <= 4 for x in counts),
            f"Bad six-field structure {counts!r}")
    return "_".join(str(x) for x in counts)


def check_monomial(counts: list[int], alpha: int, degree: int) -> dict[str, str]:
    structure_key(counts)
    require(sum(counts) == degree, f"Not a degree-{degree} scalar operator")
    require((sum(counts[2:]) % 2) == 0, "Operator is Z2 odd")
    hypercharge2 = sum(n * y for n, y in zip(counts, scalar_charges2(alpha)))
    require(hypercharge2 == 0, f"Nonzero hypercharge 2Y={hypercharge2}")
    require(sum(counts[2:]) > 0, "SM-only operator in BSM registry")
    return {"counts": structure_key(counts), "hypercharge2": "0", "z2": "even"}


def check_yukawa(model: str) -> dict[str, Any]:
    cls, value = model.split("_")
    a = int(value)
    s1, s2, df = REP[cls]
    d = [s1, s2]
    yf = Fraction(a + 1, 2)
    y1, y2 = Fraction(a, 2), Fraction(a + 2, 2)
    # These are the actual Bar[L] F^c/F S1/Bar[S2] vertex orientations.
    ordinary = {"y1": Fraction(1, 2) - yf + y1,
                "y2": Fraction(1, 2) + yf - y2}
    require(all(q == 0 for q in ordinary.values()),
            f"{model}: ordinary Yukawa hypercharge mismatch")
    # jF and jS have to couple to SU(2) spin 1/2 (SM lepton doublet).
    require(all(abs(df - ds) == 1 for ds in d),
            f"{model}: SU(2) triple product lacks a singlet")
    extra = general.expected_yukawa(model)
    if a == 0:
        require(y1 == 0 and extra == "y1Conj", f"{model}: neutral S1 Yukawa missing")
        charge = Fraction(1, 2) - yf - y1
    elif a == -2 and cls in "ADE":
        require(y2 == 0 and extra == "y2Conj", f"{model}: neutral S2 Yukawa missing")
        charge = Fraction(1, 2) + yf + y2
    else:
        require(extra is None, f"{model}: spurious extra Yukawa")
        charge = None
    if charge is not None:
        require(charge == 0, f"{model}: extra Yukawa hypercharge mismatch")
    return {"extra_coupling": extra, "ordinary_hypercharge": "0",
            "extra_hypercharge": "0" if extra else None,
            "su2_triangle": "Pass"}


def declared_quartics(model: str, entry: dict, mentry: dict,
                      interactions: list[str]) -> tuple[list[dict], int, int]:
    """Verify exact independent invariant inventory, incl. reused B/C quartics."""
    a = int(model.split("_")[1])
    manifest = mentry["extra_quartic_structures"]
    new, reused = entry["new_rgbeta_symbols"], entry["preexisting_rgbeta_symbols"]
    wanted = {}
    for row in manifest:
        counts = row["counts"]
        check_monomial(counts, a, 4)
        key = structure_key(counts)
        require(key not in wanted, f"{model}: duplicate manifest quartic monomial")
        n = row["invariants"]
        require(type(n) is int and n >= 1, "Nonpositive quartic tensor rank")
        wanted[key] = n
    actual = {}
    for row in new:
        counts = row["counts"]
        check_monomial(counts, a, 4)
        key = structure_key(counts)
        require(re.fullmatch(r"lambdaT3GeneralQ\d{6}I\d+", row["name"]) is not None,
                f"{model}: malformed generated RGBeta coupling name")
        require(row.get("mapping_status") == "NotNormalised",
                f"{model}: unexpected physical quartic mapping claim")
        actual[key] = actual.get(key, 0) + 1
    for row in reused:
        counts = row["counts"]
        check_monomial(counts, a, 4)
        key = structure_key(counts)
        require(isinstance(row.get("existing_name"), str), "Missing reused quartic name")
        actual[key] = actual.get(key, 0) + 1
    require(actual == wanted,
            f"{model}: RGBeta invariant multiplicities differ from manifest: "
            f"declared={wanted}, RGBeta={actual}")
    require(sum(wanted.values()) == mentry["extra_quartic_complex_couplings"]
            == entry["total_additional_quartics"],
            f"{model}: scalar quartic counts disagree")
    expected_interactions = {
        "T3ExtraQ_" + key + "_" + str(i)
        for key, number in wanted.items() for i in range(1, number + 1)
    }
    observed = {x for x in interactions if x.startswith("T3ExtraQ_")}
    require(observed == expected_interactions,
            f"{model}: Matchete interaction inventory missing/extra quartic names: "
            f"missing={sorted(expected_interactions - observed)} "
            f"extra={sorted(observed - expected_interactions)}")
    rows = [{"model": model, "monomial": key, "invariants": number,
             "matchete_names": ";".join(f"T3ExtraQ_{key}_{i}" for i in range(1, number + 1)),
             "rank_source": "manifest-and-live-Matchete-operator-audit",
             "physical_coupling_map": PENDING}
            for key, number in sorted(wanted.items())]
    return rows, len(new), len(reused)


def check_quadratics(model: str, mentry: dict, interactions: list[str], kind: str) -> int:
    a = int(model.split("_")[1])
    manifest = mentry["quadratic_structures_not_enabled"]
    count = 0
    expected_names = set()
    for row in manifest:
        counts = row["counts"]
        check_monomial(counts, a, 2)
        n = row["invariants"]
        require(type(n) is int and n == 1, f"{model}: unexpected quadratic invariant rank")
        expected_names.add("T3ExtraM_" + structure_key(counts) + "_1")
        count += n
    require(count == mentry["quadratic_complex_couplings_not_enabled"] == int(kind != "none"),
            f"{model}: mass-operator inventory mismatch")
    actual_names = {x for x in interactions if x.startswith("T3ExtraM_")}
    require(actual_names == expected_names,
            f"{model}: scalar bilinear missing/spurious in Matchete Lagrangian")
    return count


def check_native_rges(model: str, folder: Path, entry: dict, extra: str | None) -> dict:
    results = {}
    quartics = {x["name"] for x in entry["new_rgbeta_symbols"]} | {
        x["existing_name"] for x in entry["preexisting_rgbeta_symbols"]}
    for stage in ("UV", "EFT1"):
        native = read(folder / f"{stage.lower()}_rgbeta.json")
        beta = native["report_betas"]
        keys = set(beta)
        require(EXPECTED_RGE_BASE <= keys,
                f"{model}/{stage}: missing baseline beta functions: {EXPECTED_RGE_BASE-keys}")
        require(quartics <= keys,
                f"{model}/{stage}: general quartic beta coverage incomplete: {quartics-keys}")
        present = {"y1Conj", "y2Conj"} & keys
        wanted = {extra} if stage == "UV" and extra else set()
        require(present == wanted, f"{model}/{stage}: wrong extra Yukawa beta set")
        if stage == "UV":
            require({"y1", "y2", "MF"} <= keys,
                    f"{model}: UV fermion beta functions absent")
        else:
            require(not ({"y1", "y2", "MF", "y1Conj", "y2Conj"} & keys),
                    f"{model}: integrated-out fermion still active in EFT1 RGEs")
        results[stage] = {"betas": len(beta), "general_quartic_betas": len(quartics),
                          "extra_yukawa_beta": extra if stage == "UV" else None,
                          "verified": PASS}
    return results


def check_rep_group_factors(model: str, report_path: Path | None) -> list[dict]:
    """Exact SU(2) Casimirs/Dynkin indices; not an RGE tensor-dictionary proof."""
    family, num = model.split("_")
    a = int(num)
    d1, d2, df = REP[family]
    vals = (("S1", d1, Fraction(a, 2)),
            ("S2", d2, Fraction(a + 2, 2)),
            ("F", df, Fraction(a + 1, 2)))
    rows = []
    for field, dim, y in vals:
        c2, dynkin = exact_group_invariants(dim)
        rows.append({"model": model, "field": field, "d_SU2": dim, "Y": str(y),
                     "C2_SU2": str(c2), "T_SU2": str(dynkin)})
    if report_path is not None and report_path.exists():
        with report_path.open(newline="", encoding="utf-8-sig") as stream:
            saved = {(r["model"], r["field"]): r for r in csv.DictReader(stream)}
        for r in rows:
            old = saved.get((model, r["field"]))
            require(old is not None, f"{model}: missing saved group invariant for {r['field']}")
            for key in ("d_SU2", "Y", "C2_SU2", "T_SU2"):
                require(str(old[key]) == str(r[key]),
                        f"{model}: wrong saved {r['field']} {key}: {old[key]}")
    return rows


def audit_model(model: str, registry: dict, manifest: dict, source_hashes: dict,
                run_dir: Path, group_csv: Path | None) -> tuple[dict, list[dict], list[dict]]:
    entry = registry["models"][model]
    mentry = manifest["models"][model]
    folder = run_dir / tag(model)
    kind = quadratic_kind(manifest, model)
    # This is the source-fingerprint, Matchete status, and native-RGBeta export gate.
    cached = general.validate_cached(model, folder, entry, kind, source_hashes)
    require(cached["status"] in ("SymbolicRGExportedMatchingUnverified",
                                 "SymbolicRGExportedMassBasisBlocked"),
            f"{model}: unexpected general-symbolic status")
    require(entry["alpha"] == int(model.split("_")[1]) and
            entry["dimensions"] == list(REP[model[0]]),
            f"{model}: registry representation changed")
    interactions = cached["lagrangian"]["interactions"]
    require(len(interactions) == len(set(interactions)), f"{model}: duplicate interaction names")
    require(any(x.startswith("Yukawa1_") for x in interactions) and
            any(x.startswith("Yukawa2_") for x in interactions),
            f"{model}: ordinary Yukawa interaction absent")
    yukawa = check_yukawa(model)
    quartics, n_new, n_reused = declared_quartics(model, entry, mentry, interactions)
    bilinears = check_quadratics(model, mentry, interactions, kind)
    rges = check_native_rges(model, folder, entry, yukawa["extra_coupling"])
    group = check_rep_group_factors(model, group_csv)
    row = {
        "model": model, "status": PASS,
        "extra_quartics": entry["total_additional_quartics"],
        "extra_quadratics": bilinears,
        "extra_yukawas": int(bool(yukawa["extra_coupling"])),
        "extra_yukawa_name": yukawa["extra_coupling"],
        "rgbeta_new_quartics": n_new, "rgbeta_reused_quartics": n_reused,
        "uv_betas": rges["UV"]["betas"], "eft1_betas": rges["EFT1"]["betas"],
        "su2_group_invariants": PASS,
        "charge_selection_rules": PASS,
        "matchete_detailed_audit": (
            "MassBasisBlocked" if model in general.MIXED_MASS_MODELS else "Passed"),
        "physical_quartic_coupling_dictionary": PENDING,
        "physical_yukawa_coupling_dictionary": PENDING,
        "physical_beta_coefficient_equality": PENDING,
        "general_threshold_matching": PENDING,
    }
    return row, quartics, group


def known_bridge_evidence() -> list[dict]:
    """Describe *actual* pre-existing audit locations; never upgrade their scope."""
    return [
        {"evidence": "Yukawa SU(2)/U(1) component Ward checks",
         "source": "Lagrangian/RunT3YukawaComponentWard.py",
         "saved_evidence": "Reports/output/workstream_a/yukawa_component_ward/yukawa_component_ward.json",
         "certifies_physical_dictionary": False},
        {"evidence": "Native Yukawa CG norms and internal RGBeta prefactors",
         "source": "Lagrangian/RunT3YukawaNormalizationBridge.py",
         "saved_evidence": "Reports/output/workstream_a/yukawa_normalization_bridge",
         "certifies_physical_dictionary": False},
        {"evidence": "B/C quartic invariant-tensor geometry",
         "source": "docs/T3_BC_RGBETA_NORMALISATION_NEXT_CHECK.md",
         "saved_evidence": "Source documents / previous geometric audit",
         "certifies_physical_dictionary": False},
        {"evidence": "A(-2) scalar-bilinear mass mapping",
         "source": "scripts/VerifyA2BilinearNormalization.py",
         "saved_evidence": "output/studies/general_scalar_extended/a2_bilinear_normalization",
         "certifies_physical_dictionary": False},
    ]


def inspect_saved_bridge_reports(root: Path = ROOT) -> list[dict]:
    """Read existing *partial* proof certificates, never relabel them physical."""
    sources = (
        ("Yukawa component Ward",
         "Reports/output/workstream_a/yukawa_component_ward/yukawa_component_ward.json",
         {"status": "ExactSU2U1ComponentWardPassedPhysicalCrossPackageDictionaryPending",
          "models_checked": 16, "baseline_yukawas_checked": 32,
          "new_native_extra_yukawas_checked": 6,
          "su2_ward_component_equations_checked": 1032,
          "physical_matchete_rgbeta_yukawa_normalization_certified": False}),
        ("Yukawa native-CG magnitudes / RGBeta internal coefficient",
         "Reports/output/workstream_a/yukawa_normalization_bridge/yukawa_normalization_bridge.json",
         {"models_checked": 16, "native_matchete_yukawas_checked": 38,
          "singlet_geometry_count": 14, "triplet_geometry_count": 24,
          "matchete_rgbeta_absolute_coupling_dictionary_certified": False,
          "matchete_rgbeta_physical_action_verified": False}),
    )
    evidence = []
    for label, rel, required in sources:
        path = root / rel
        if not path.is_file():
            evidence.append({"name": label, "path": rel, "status": "NotAvailable",
                             "physical_dictionary_certified": False})
            continue
        try:
            data = read(path)
            ok = all(data.get(k) == v and type(data.get(k)) is type(v)
                     for k, v in required.items())
            if "normalization_bridge" in rel:
                ok &= data.get("status") in (
                    "ExactCGMagnitudesAndNativeRGBetaInternalPrefactorsVerifiedActionDictionaryPending",
                    "ExactCGMagnitudesVerifiedRGBetaInternalProbePending")
            evidence.append({"name": label, "path": rel,
                             "status": "VerifiedLimitedScope" if ok else "UnexpectedOrStale",
                             "saved_status": data.get("status"),
                             "physical_dictionary_certified": False})
        except (OSError, ValueError, TypeError, json.JSONDecodeError) as exc:
            evidence.append({"name": label, "path": rel,
                             "status": "Unreadable", "error": str(exc),
                             "physical_dictionary_certified": False})
    return evidence


def save_report(destination: Path, rows: list[dict], quartic_rows: list[dict],
                group_rows: list[dict], requested: list[str], run_dir: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    totals = {"models": len([r for r in rows if r["status"] == PASS])}
    for key in ("extra_quartics", "extra_quadratics", "extra_yukawas",
                "rgbeta_new_quartics", "rgbeta_reused_quartics"):
        totals[key] = sum(int(r.get(key, 0)) for r in rows if r["status"] == PASS)
    total_gate = (len(requested) == 16 and set(requested) == set(general_models())
                  and totals == EXPECTED_TOTALS)
    all_checked = len(rows) == len(requested) and all(r["status"] == PASS for r in rows)
    data = {
        "format": FORMAT, "input_run": str(run_dir),
        "models_requested": requested, "model_count_requested": len(requested),
        "model_count_structurally_verified": totals["models"],
        "per_model_structural_checks_passed": all_checked,
        "full_16_expected_inventory_passed": total_gate,
        "totals_of_structurally_verified": totals,
        "physical_matchete_rgbeta_coupling_dictionary_certified": False,
        "one_loop_beta_coefficients_independently_equal_after_physical_mapping": False,
        "full_general_sequential_matching_certified": False,
        "overall": ("StructuralInventoryVerifiedPhysicalDictionaryPending" if all_checked
                    else "StructuralInventoryIncompleteOrFailed"),
        "models": rows,
        "existing_bridge_evidence": known_bridge_evidence(),
        "saved_bridge_report_checks": inspect_saved_bridge_reports(),
        "scope": "Source-matched inventory, exact charge/C2/T checks, RGBeta beta-key coverage; NOT coefficient identity",
    }
    (destination / "audit.json").write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    def table(name: str, entries: list[dict], keys: list[str]) -> None:
        with (destination / name).open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=keys, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(entries)
    table("audit.csv", rows, ["model", "status", "extra_quartics", "extra_quadratics",
                              "extra_yukawas", "rgbeta_new_quartics", "rgbeta_reused_quartics",
                              "uv_betas", "eft1_betas", "matchete_detailed_audit", "error"])
    table("quartic_invariants.csv", quartic_rows,
          ["model", "monomial", "invariants", "matchete_names", "rank_source", "physical_coupling_map"])
    table("group_invariants.csv", group_rows,
          ["model", "field", "d_SU2", "Y", "C2_SU2", "T_SU2"])
    lines = ["# General 16-model symbolic physics audit", "",
             f"Structurally verified: **{totals['models']}/{len(requested)}**.",
             f"Complete 16-model inventory totals gate: **{'PASS' if total_gate else 'not passed / not applicable'}**.",
             "", "## What has actually been checked", "",
             "- Current-source SHA-256 provenance and existing Matchete detailed-check statuses",
             "- Each independent additional scalar monomial's U(1) hypercharge and Z2 selection rules",
             "- Quartic multiplicities cross-checked among the manifest, Matchete interaction names and RGBeta registry",
             "- Additional scalar mass and Yukawa operator inventory (including six exceptional Yukawas)",
             "- Required UV and EFT1 beta-function *keys*, including all general quartics and massive operators",
             "- Independent exact SU(2) Casimirs and Dynkin indices compared with the symbolic-study report",
             "", "## What remains NOT CERTIFIED", "",
             "1. Physical action-level Matchete↔RGBeta coupling dictionary: tensor basis, phase, complex conjugation, signs, factorials and H.c.",
             "2. Equality of each group-factor/beta coefficient after applying that physical dictionary",
             "3. Scalar mass-basis matching for B(-1), C(-1) and complete extended sequential EFT/C5 matching",
             "", "The existing tensor-Ward and normalization-bridge tests are partial evidence; this audit never upgrades them to physical equivalence.",
             "Inspect `audit.json` -> `saved_bridge_report_checks` for the status of existing 38-Yukawa witness files.",
             "", "## Per-model structural gate", "",
             "| Model | Quartics | Mass terms | Extra Yukawa | UV betas | EFT1 betas | Status |",
             "|---|---:|---:|---|---:|---:|---|",]
    for r in rows:
        lines.append(f"| {r['model']} | {r.get('extra_quartics','—')} | "
                     f"{r.get('extra_quadratics','—')} | {r.get('extra_yukawa_name') or '—'} | "
                     f"{r.get('uv_betas','—')} | {r.get('eft1_betas','—')} | "
                     f"{r['status']} |")
        if r.get("error"):
            lines.append(f"\n**{r['model']} failure:** `{r['error']}`\n")
    lines += ["", "See `audit.json`, `quartic_invariants.csv`, and `group_invariants.csv` for machine-readable evidence.", ""]
    (destination / "README.md").write_text("\n".join(lines), encoding="utf-8")


def general_models() -> list[str]:
    return sorted(f"{cls}_{alpha:+d}" for cls, alphas in MODEL_CLASSES.items()
                  for alpha in alphas)


def self_test() -> None:
    require(len(general_models()) == 16, "wrong catalogue")
    require(exact_group_invariants(1) == (0, 0), "singlet group factors")
    require(exact_group_invariants(2) == (Fraction(3, 4), Fraction(1, 2)),
            "doublet group factors")
    require(exact_group_invariants(3) == (2, 2), "triplet group factors")
    for model in general_models():
        check_yukawa(model)
    require(sum(bool(general.expected_yukawa(m)) for m in general_models()) == 6,
            "extra Yukawa catalogue")
    check_monomial([0, 0, 0, 0, 0, 4], -2, 4)
    for args in (([1, 0, 0, 0, 0, 3], -2, 4),
                 ([0, 0, 1, 0, 0, 0], -2, 1)):
        try:
            check_monomial(*args)
        except ValueError:
            pass
        else:
            raise AssertionError("Invalid gauge/Z2 monomial was accepted")
    require(structure_key([0, 0, 0, 0, 0, 4]) == "0_0_0_0_0_4", "bad name")
    print("PASS: 16 representation/selection-rule witnesses, exact C2/T and negative controls")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", default=None, help="Single model e.g. A_-2 (pilot)")
    parser.add_argument("--run-id", default="full16_general_v1")
    parser.add_argument("--study", default="full16_general_symbolic")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args(argv)
    if args.self_test:
        self_test()
        return 0
    for name in ("run_id", "study"):
        require(re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]{0,63}", getattr(args, name)) is not None,
                f"Unsafe {name} directory name")
    names = general_models()
    if args.model is not None:
        require(args.model in names, f"Unknown model {args.model}")
        names = [args.model]
    # Cross-check the production-source and manifest inventory *without* running Wolfram.
    basis = quadratic_audit()
    require(set(basis["models"]) == set(general_models()), "Quadratic audit catalogue changed")
    registry = read(scalar.REGISTRY)
    manifest = load_manifest(scalar.MANIFEST)
    scalar.check_sources(names, registry, manifest, basis)
    current_sources = general.fingerprint()
    run_dir = general.DEFAULT_DATA_ROOT / args.run_id
    group_csv = ROOT / "Reports/output" / args.study / "GeneralSymbolic/group_invariants.csv"
    destination = ROOT / "Reports/output" / args.study / "GeneralPhysicsAudit"
    rows: list[dict] = []
    quartics: list[dict] = []
    groups: list[dict] = []
    for i, model in enumerate(names, 1):
        try:
            row, qrows, grows = audit_model(model, registry, manifest, current_sources,
                                            run_dir, group_csv)
            quartics.extend(qrows)
            groups.extend(grows)
        except (OSError, ValueError, KeyError, TypeError, json.JSONDecodeError) as exc:
            row = {"model": model, "status": FAILURE, "error": str(exc),
                   "physical_quartic_coupling_dictionary": PENDING,
                   "physical_yukawa_coupling_dictionary": PENDING}
        rows.append(row)
        print(f"[{i}/{len(names)}] {model}: {row['status']}"
              + (f" ({row['error']})" if "error" in row else ""), flush=True)
    save_report(destination, rows, quartics, groups, names, run_dir)
    structural_pass = all(r["status"] == PASS for r in rows)
    print(f"Structural inventory: {sum(r['status'] == PASS for r in rows)}/{len(names)}; "
          f"report: {destination/'README.md'}", flush=True)
    print("Physical Matchete/RGBeta action-level dictionary and RGE coefficient equality: NOT CERTIFIED", flush=True)
    return 0 if structural_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
