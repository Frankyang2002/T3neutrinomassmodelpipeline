r"""Verify the archived B(-1)/C(-1) native Y/X matching evidence; fail closed.

This is an evidence/provenance gate, NOT a calculation of physical C5.
The archive itself is read-only; no Matchete expression is evaluated, no pole
is silently deleted, and no source code/UV/RGE convention is modified.

Intended repository path: scripts/AuditT3BCPhysicalYMatchingEvidence.py
Example (from repository root):
  python scripts/AuditT3BCPhysicalYMatchingEvidence.py --archive .\xy_fresh_native_evidence_v3.zip --output output\audits\t3_bc_physical_y_matching_gate_v1\summary.json

This audit checks exact SHA256 metadata, raw sectors (including absence of
standalone X constants), Matchete check outcomes as reported in archived JSON,
and whether the X source inherited Y one-loop EOMEFT.  It deliberately cannot
establish algebraic equivalence or pole cancellation from an unpaired amplitude.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path, PurePosixPath
import zipfile

FORMAT = "t3-bc-physical-y-evidence-gate-v1"
ROOT = "output/general_t3_mass_eigen_matchete"
EXPECTED = {
    "Y/LLHH": {"B_-1": (132, 50, 50, 32), "C_-1": (132, 50, 50, 32)},
    "Y/LLXX": {"B_-1": (666, 240, 240, 186), "C_-1": (582, 198, 198, 186)},
    "X/LLHH": {"B_-1": (36, 18, 18, 0), "C_-1": (36, 18, 18, 0)},
}


class EvidenceAuditError(ValueError):
    pass


def require(ok: bool, message: str) -> None:
    if not ok:
        raise EvidenceAuditError(message)


def raw_terms(src: str) -> tuple[str, ...]:
    """Split additive Matchete InputForm at top level, not inside brackets.

    Delimiters in the archived exports are space-sign-space; this intentionally
    rejects unexpected formatting rather than attempting symbolic algebra.
    """
    src = src.strip()
    require(bool(src), "empty InputForm sector")
    out = []
    brackets = {"(": ")", "[": "]", "{": "}"}
    stack: list[str] = []
    last = 0
    for i, char in enumerate(src):
        if char in brackets:
            stack.append(brackets[char])
        elif char in ")]}":
            require(bool(stack) and stack.pop() == char,
                    f"unbalanced InputForm at position {i}")
        elif (char in "+-" and not stack and i > 0 and i + 1 < len(src)
              and src[i-1].isspace() and src[i+1].isspace()):
            out.append(src[last:i].strip())
            last = i
    require(not stack, "unbalanced trailing InputForm brackets")
    out.append(src[last:].strip())
    require(all(out), "empty additive term in raw Matchete expression")
    return tuple(out)


def audited_sector(src: str, *, sector: str, model: str) -> dict:
    ts = raw_terms(src)
    require(all("hbar" in t and "GammaCC" in t for t in ts),
            f"{model} {sector}: hbar/GammaCC inventory not established")
    pole = sum("ϵ" in t for t in ts)
    logs = sum("Log[" in t for t in ts)
    require(not any("ϵ" in t and "Log[" in t for t in ts),
            f"{model} {sector}: mixed pole/log term; split is ambiguous")
    constants = len(ts) - pole - logs
    observed = (len(ts), pole, logs, constants)
    require(observed == EXPECTED[sector][model],
            f"{model} {sector}: unexpected count {observed}, expected {EXPECTED[sector][model]}")
    return {"source_additive_terms": observed[0], "pole_terms": pole,
            "log_terms": logs, "mass_independent_additive_terms": constants,
            "pole_log_coefficient_equality_verified_here": False,
            "pole_subtraction_is_physical_matching": False}


def run(archive: Path) -> dict:
    require(archive.is_file(), f"missing archive: {archive}")
    archive_hash = sha256(archive.read_bytes()).hexdigest()
    with zipfile.ZipFile(archive) as z:
        entries = z.namelist()
        require(len(entries) == len(set(entries)), "duplicate ZIP paths")
        for p in entries:
            pp = PurePosixPath(p)
            require(not pp.is_absolute() and ".." not in pp.parts and "\\" not in p,
                    f"invalid archive member: {p}")
        require("MANIFEST.json" in entries, "missing MANIFEST.json")
        manifest = json.loads(z.read("MANIFEST.json"))
        require(manifest.get("format") == "t3-bc-partial-xy-native-evidence-v3",
                "unexpected native evidence format")
        require(manifest.get("models") == ["B_-1", "C_-1"], "unexpected model inventory")
        recorded = manifest.get("files")
        require(isinstance(recorded, dict), "missing file manifest")
        require(set(recorded) == set(entries) - {"MANIFEST.json"},
                "manifest files do not equal ZIP member inventory")
        for p, info in recorded.items():
            b = z.read(p)
            require(info.get("size_bytes") == len(b), f"size mismatch: {p}")
            require(info.get("sha256") == sha256(b).hexdigest(), f"SHA256 mismatch: {p}")
        require(all(manifest.get(k) is False for k in (
            "y_one_loop_inherited_in_x_source",
            "f_first_eomeft_equivalence_certified",
            "running_boundary_substituted",
            "physical_sequential_matching_certified")),
            "unexpected manifest certification status")
        def read(path: str) -> str:
            require(path in recorded, f"missing required evidence: {path}")
            return z.read(path).decode("utf-8")
        def js(path: str) -> dict:
            r = json.loads(read(path))
            require(isinstance(r, dict), f"invalid report: {path}")
            return r
        models = {}
        for prefix, model in (("B", "B_-1"), ("C", "C_-1")):
            saved_path = (f"{ROOT}/eft1_fermion_first_tree/"
                          f"{prefix}_m1_final_eft1_fermion_first_tree_eft1_transition_eom.inputform.txt")
            saved = read(saved_path)
            require("GammaCC" in saved and "Field[NewScalarX" in saved
                    and "Field[NewScalarY" in saved and "Field[NewFermion" not in saved,
                    f"{model}: invalid saved fermion-first EOMEFT field inventory")
            y_root = f"{ROOT}/eft1_scalar_Y_only_canonical_probe/{prefix}_m1_final_y_only_canonical_match"
            y_report = js(y_root + ".json")
            eq = js(y_root + "_reconstruction_audit.json")
            herm = js(y_root + "_hermiticity_diagnostics.json")
            require(y_report.get("model") == model and
                    y_report.get("status") == "YOnlyIndependentCanonicalRawEFT2NotCertified",
                    f"{model}: unexpected Y matching status")
            require(y_report.get("physical_matching_performed") is False and
                    y_report.get("hard_pole_counterterms_combined") is False and
                    y_report.get("eft1_running_performed") is False,
                    f"{model}: Y counterterm/running status changed")
            require(eq.get("algebraic_equality_proven") is False and
                    eq.get("accepted") is False and
                    eq.get("independent_matchete_check_passed") is True and
                    herm.get("full_saved_eft_hermiticity") is False,
                    f"{model}: EFT1 reconstruction status inconsistent")
            prefix_diff = eq.get("residual_nonzero_inputform_prefix")
            require(isinstance(prefix_diff, str) and
                    all(x in prefix_diff for x in ("Field[H, Scalar", "Field[NewScalarX", "Field[NewScalarY")),
                    f"{model}: residual prefix no longer HHXY; review explicitly")
            sectors = {}
            for key, suffix in (("Y/LLHH", "raw_llhh"), ("Y/LLXX", "raw_llxx")):
                path = y_root + "_" + suffix + ".inputform.txt"
                sectors[key] = {"sha256": recorded[path]["sha256"],
                                **audited_sector(read(path), sector=key, model=model)}
            x_root = f"{ROOT}/eft1_scalar_xy_fresh_kernel_v3/{prefix}_m1_final_x_from_saved_y_match"
            x_report = js(x_root + ".json")
            x_stage = js(x_root + "_x_followup.json")
            x_path = x_root + "_x_followup_x_raw_llhh.inputform.txt"
            sectors["X/LLHH"] = {"sha256": recorded[x_path]["sha256"],
                                  **audited_sector(read(x_path), sector="X/LLHH", model=model)}
            require(x_report.get("fresh_kernel_without_y_matching") is True and
                    x_report.get("y_one_loop_correction_inherited_in_x_source") is False and
                    x_stage.get("y_loop_inherited_in_x_matching") is False and
                    x_stage.get("source_is_y_tree_EOMEFT") is True and
                    x_stage.get("x_tree_matching_completed") is True,
                    f"{model}: X source/tree-stage provenance inconsistent")
            require("hbar" in read(x_root + "_x_followup_y_full_eomeft.inputform.txt"),
                    f"{model}: archived Y one-loop expression absent")
            models[model] = {
                "original_fermion_first_eomeft_sha256": recorded[saved_path]["sha256"],
                "original_fermion_first_eomeft_hermiticity_check_passed": False,
                "independent_canonical_matchete_checks_passed": True,
                "source_algebraic_equivalence_proven": False,
                "nonzero_reconstruction_residual_prefix": prefix_diff,
                "Y_native_tree_and_one_loop_completed": True,
                "Y_finite_and_pole_terms_are_only_formal": True,
                "Y_loop_EOMEFT_exists_but_inherited_into_X": False,
                "X_native_tree_completed_and_raw_loop_sector_exported": True,
                "X_absent_standalone_constant_implies_physical_zero": False,
                "sector_inventory": sectors,
                "physical_matching_certified": False,
                "physical_c5": None,
            }
    return {
        "format": FORMAT,
        "source_archive_sha256": archive_hash,
        "archive_manifest_hashes_verified": True,
        "source_used": "user-uploaded xy_fresh_native_evidence_v3.zip",
        "models": models,
        "first_blocker": "original fermion-first EOMEFT vs independent canonical scalar EFT algebraic equivalence unproven; source checks fail Hermiticity",
        "next_required_operation": "resolve full HHXY quartic residual with native CG/index identities and establish source equivalence before physical EFT1-EFT2 subtraction",
        "further_blockers": [
            "one-loop EFT1 minus EFT2 1PI amplitude subtraction with same infrared prescription not demonstrated",
            "counterterm and finite EOM/wavefunction/parameter shifts not combined",
            "Matchete GammaCC chirality and Majorana sign/factor map to ordered A5 not certified",
            "Y loop and LLXX inheritance at X not implemented at consistent fixed order",
            "EFT1 running and quartic/mixing physical normalization not certified",
        ],
        "renormalized_finite_y_c5_obtained": False,
        "physical_matching_certified": False,
        "physical_c5": None,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--archive", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    try:
        result = run(args.archive)
    except (EvidenceAuditError, ValueError, OSError, zipfile.BadZipFile) as exc:
        ap.exit(2, f"FAIL (audit not certified): {exc}\n")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    for model, model_info in result["models"].items():
        print(f"{model}: 3 raw sectors verified; source equivalence=UNPROVEN; physical Y C5=NOT CERTIFIED")
    print(f"PASS: archive hashes and B/C evidence gates; report={args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
