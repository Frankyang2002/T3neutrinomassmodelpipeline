"""Worker C: read-only gate for existing Matchete scalar-threshold exports.

No expression in a raw Matchete Match result is promoted automatically to a
renormalized, restricted-Yukawa, bubble-subtracted SMEFT C5 coefficient. This
module reads source-provenance files and returns *evidence*, not a Wilson matrix.

Runs independently of Wolfram, NumPy, Matchete and the numerical pipeline.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re
from typing import Any

FORMAT = "t3-worker-c-matchete-scalar-audit-v1"

# These files belong to a single *saved* diagnostic profile. Do not infer that
# the source still has the same state as a separate current local checkout.
FILES = {
    "match_result": "profile_full_scalar_one_loop_result.json",
    "mass_audit": "profile_full_scalar_one_loop_audit.json",
    "fixed_order": "profile_fixed_order_components_result.json",
    "flavour_candidate": "profile_c5_flavour_candidate_result.json",
    "pole_audit": "profile_c5_pole_cancellation_result.json",
    "normalization_audit": "profile_c5_normalization_reference_result.json",
    "candidate_inputform": "profile_full_scalar_one_loop_hard_c5_candidate.txt",
}
TOKEN_Y2 = re.compile(r"(?<![A-Za-z0-9_])y2(?![A-Za-z0-9_])")
TOKEN_Y2CONJ = re.compile(r"(?<![A-Za-z0-9_])y2Conj(?![A-Za-z0-9_])")
TOKEN_POLE = re.compile(r"(?:\bϵ\b|\bε\b|\\\[Epsilon\]|\bepsilon\b|\bEpsilon\b)")


def _hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _read(path: Path, *, text_only: bool) -> Any:
    value = path.read_text(encoding="utf-8")
    if text_only:
        return value
    doc = json.loads(value)
    if not isinstance(doc, dict):
        raise ValueError(f"Expected JSON object at {path}")
    return doc


def audit_matchete_scalar_profile(profile_dir: str | Path, *, model: str = "D_-2") -> dict:
    """Inspect a saved Matchete scalar-threshold profile without writing to it.

    Returns a portable JSON-serializable audit with evidence/blocks. Missing,
    malformed or contradictory files fail closed. In particular the candidate
    cannot be numerically added to the existing LLSS*HHSS bubble on this evidence.
    """
    profile = Path(profile_dir).resolve()
    evidence: dict[str, Any] = {}
    sources: dict[str, dict] = {}
    errors: list[str] = []
    for key, name in FILES.items():
        path = profile / name
        record: dict[str, Any] = {"filename": name, "exists": path.is_file()}
        if not path.is_file():
            errors.append(f"missing {name}")
            sources[key] = record
            continue
        try:
            content = _read(path, text_only=key == "candidate_inputform")
            record["sha256"] = _hash(path)
            record["size_bytes"] = path.stat().st_size
            evidence[key] = content
        except (OSError, ValueError, UnicodeError) as exc:
            errors.append(f"invalid {name}: {exc}")
        sources[key] = record

    match = evidence.get("match_result", {})
    mass = evidence.get("mass_audit", {})
    fixed = evidence.get("fixed_order", {})
    flavour = evidence.get("flavour_candidate", {})
    pole = evidence.get("pole_audit", {})
    norm = evidence.get("normalization_audit", {})
    candidate = evidence.get("candidate_inputform", "")

    hard_has_y2 = bool(TOKEN_Y2.search(candidate))
    hard_has_y2conj = bool(TOKEN_Y2CONJ.search(candidate))
    hard_has_pole = bool(TOKEN_POLE.search(candidate))
    source_y2conj_terms = mass.get("Y1Y2ConjWilsonTermCount")

    gates = {
        "source_files_complete": len(evidence) == len(FILES) and not errors,
        "exact_eft1_partition_passed": mass.get("Status") == "Passed"
            and mass.get("SMAlphaComparisonCertified") is True
            and mass.get("AllHeavyTermsPreservedExactly") is True,
        "all_scalars_heavy_same_threshold": mass.get("AllPhysicalScalarsRegisteredHeavy") is True
            and mass.get("PureHeavyQuadraticAlphaEqualToFreeLag") is True
            and mass.get("CheckLagrangianPassed") is True,
        "raw_one_loop_match_returned": match.get("OneLoopRawMatchReturned") is True
            and match.get("LoopOrder") == 1 and match.get("EFTOrder") == 5,
        "matchete_postprocessing_returned": match.get("SimplificationPerformed") is True
            and match.get("HardC5Extracted") is True
            and match.get("WeinbergExtractionStatus") == "Success",
        "one_loop_weinberg_terms_present": isinstance(match.get("WeinbergTermCount"), int)
            and match["WeinbergTermCount"] > 0,
        "no_independent_conjugate_scalar_yukawa_in_source": source_y2conj_terms == 0
            and not hard_has_y2conj and fixed.get("HardC5HasY2Conj") is False,
        "hard_candidate_has_no_unsubtracted_poles": bool(candidate) and not hard_has_pole,
        "pole_cancellation_audited_at_reported_nf": pole.get("ExactPoleCancellationVerifiedForNF1") is True,
        "full_flavour_hard_c5_verified": match.get("FullFlavorHardC5Verified") is True
            and fixed.get("FlavourAssemblyCertified") is True,
        "operator_normalization_certified": flavour.get("OperatorNormalizationIndependentlyCertified") is True
            and norm.get("FullPhysicalC5NormalizationCertified") is True,
        "finite_renormalization_certified": pole.get("RenormalizationCertified") is True,
        "physical_mass_relations_imposed": match.get("PhysicalMassRelationsImposed") is True,
        "independent_nonbubble_remainder_certified": match.get(
            "FiniteHardMinusBubbleRemainderCertified") is True,
        "physical_full_scalar_threshold_certified": match.get(
            "FullGeneralScalarThresholdCertified") is True,
    }

    # A matched one-loop hard coefficient may CONTAIN the already implemented
    # scalar bubble. A subtraction/projection proof is therefore mandatory.
    mandatory = (
        "source_files_complete", "exact_eft1_partition_passed",
        "all_scalars_heavy_same_threshold", "raw_one_loop_match_returned",
        "matchete_postprocessing_returned", "one_loop_weinberg_terms_present",
        "no_independent_conjugate_scalar_yukawa_in_source",
        "hard_candidate_has_no_unsubtracted_poles",
        "full_flavour_hard_c5_verified", "operator_normalization_certified",
        "finite_renormalization_certified", "physical_mass_relations_imposed",
        "independent_nonbubble_remainder_certified",
        "physical_full_scalar_threshold_certified",
    )
    source_ready = all(gates[k] for k in mandatory)
    # This module does NOT parse or certify a numerical finite C5 array. Even
    # perfect source flags would authorize *extraction work*, not addition to
    # the existing scalar bubble. Keep that final boundary fail-closed.
    safe = False
    reasons = errors + [f"unverified: {k}" for k in mandatory if not gates[k]]
    result = {
        "format": FORMAT,
        "model": model,
        "profile_dir": str(profile),
        "source_files": sources,
        "match_status": match.get("Status"),
        "candidate_inventory": {
            "total_weinberg_terms": match.get("WeinbergTermCount"),
            "holomorphic_hard_terms": fixed.get("HardWeinbergHolomorphicTermCount"),
            "conjugate_hard_terms": fixed.get("HardWeinbergConjugateTermCount"),
            "y1_y2_tree_wilson_terms": mass.get("Y1Y2WilsonTermCount"),
            "y1_y2conj_tree_wilson_terms": source_y2conj_terms,
            "raw_hard_has_y2": hard_has_y2,
            "raw_hard_has_y2conj": hard_has_y2conj,
            "raw_hard_has_explicit_epsilon_symbol": hard_has_pole,
            "pole_cancellation_nf": pole.get("NF"),
            "pole_cancellation_verified_for_nf1": gates[
                "pole_cancellation_audited_at_reported_nf"],
            "inherited_c5_already_added_to_match": match.get("InheritedOneLoopAdded"),
            "direct_running_already_added_to_match": match.get("DirectRunningAdded"),
        },
        "gates": gates,
        "source_ready_for_separate_numerical_extraction": source_ready,
        "finite_numeric_nonbubble_c5_supplied": False,
        "safe_to_add_as_independent_c5_correction": safe,
        "status": "SOURCE_READY_FOR_NUMERIC_EXTRACTION" if source_ready else "BLOCKED",
        "blocking_reasons": reasons,
        "physics_accounting": {
            "computed_by_this_audit": "none (read-only source check)",
            "existing_worker_c": "C5_EFT1(muS) + finite LLSS*HHSS joint scalar bubble",
            "candidate_isolation": "No verified UV/EFT1-to-SMEFT hard-minus-bubble remainder",
            "restriction": "Only Y1,Y2; independent y2Conj sector forbidden",
            "rule": "Do not add the Matchete candidate or its poles to c5_finite_scalar",
            "numeric_matrix": "not supplied or accepted by this read-only audit",
        },
    }
    return result


def save_audit(report: dict, output_dir: str | Path) -> Path:
    """Write a NEW report directory; never overwrite previous results."""
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=False)
    file = destination / "worker_c_scalar_audit.json"
    file.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return file
