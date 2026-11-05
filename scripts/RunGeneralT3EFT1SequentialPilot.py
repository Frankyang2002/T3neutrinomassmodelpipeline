"""General-interaction B(-3): F -> EFT1 running -> (S1,S2) -> SMEFT.

This is an opt-in pilot. It deliberately does NOT invoke the restricted
production pipeline or modify any existing model/matching source file.

Prerequisites (from repository root):

    python -X utf8 -m scripts.RunT3GeneralSymbolicStudy \
        --model B_-3 --study general_sequential_pilot \
        --run-id general_sequential_pilot

    $env:T3_EXPERIMENTAL_GENERAL_QUADRATICS = "1"
    $env:T3_EXPERIMENTAL_GENERAL_FERMIONS = "1"
    wolframscript -file Lagrangian/RunModel.wl \
        output/general_sequential_pilot/B_m3_matching 5 1 B m3 \
        'THRESHOLDS=F;S1,S2'

    python -X utf8 -m scripts.RunGeneralT3EFT1SequentialPilot

Options:
    --rge-only           Build the full one-loop Wilson tensor and direct C5
                         running insertion, without resuming scalar matching.
    --validate-threshold Enable the more expensive [A]/[B]/[C] provenance
                         diagnostics at the resumed scalar threshold.
    --matching-dir PATH  Use a different general two-stage Matchete output.
    --general-dir PATH   Use a different general-symbolic B(-3) output.

Fixed-order physics:
    C_LLSS(mu_s) = C_LLSS^(0) + hbar ln(mu_s/mu_f) beta_LLSS^(1)
    C_5(mu_s) = C_5^hard + hbar ln(mu_s/mu_f) beta_(5<-LLSS)^(1)

Only direct mixing into C5 is inserted into the scalar matching at O(hbar).
The heavy-operator O(hbar) self-running is computed for audit but is NOT
subsequently one-loop matched: that would be O(hbar^2).

This pilot does not certify general Matchete/RGBeta coupling normalisation or
an unrestricted all-16-model operator/flavour basis. It tests a specific
nontrivial general-Lagrangian benchmark and records explicit provenance.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MATCHING_DIR = ROOT / "output/general_sequential_pilot/B_m3_matching"
DEFAULT_GENERAL_DIR = (
    ROOT / "output/studies/general_scalar_extended/complete_scalar_symbolic"
    / "general_sequential_pilot/B_m3"
)
MODEL = "B_-3"
GENERAL_INTERACTION = "T3ExtraQ_0_2_0_0_0_2_1"
GENERAL_COUPLING = "lambdaT3ExtraQ020002I1"
GENERAL_RGBETA_COUPLING = "lambdaT3GeneralQ020002I1"
GENERAL_CG = "T3ExtraCGQ020002I1"
DIMENSIONS = (2, 2, 1)
ALPHA = -3


def _json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(f"Required input missing: {path}")
    value = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _written_json(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, default=str) + "\n", encoding="utf-8")


def validate_general_boundary(matching_dir: Path, general_dir: Path) -> dict[str, Any]:
    """Fail closed if matched UV is not the general B(-3) interaction basis."""
    general = _json(general_dir / "lagrangian.json")
    general_status = _json(general_dir / "status.json")
    general_eft1_rge = _json(general_dir / "eft1_rgbeta.json")
    summary = _json(matching_dir / "comparison_summary.json")
    stages = _json(matching_dir / "sequential_matching_debug.json")
    interval = _json(matching_dir / "data/eft1_after_F_rge_export.json")
    wilson = _json(matching_dir / "data/eft1_after_F_wilson_seed.json")
    quartic = _json(matching_dir / "data/eft1_after_F_scalar_quartic_seed.json")

    _require(general.get("model") == MODEL, "General symbolic model is not B_-3")
    _require(general.get("status") == "Success", "General UV Lagrangian did not validate")
    _require(general_status.get("status") == "SymbolicRGExportedMatchingUnverified"
             and not general_status.get("errors"), "General symbolic run is incomplete")
    _require(general.get("general_scalar_quadratics_requested") is True
             and general.get("extra_yukawas_requested") is True,
             "Complete general-interaction flags were not enabled")
    _require(general.get("quartic_expected") == 1
             and general.get("quartic_included") == 1,
             "B(-3) general quartic inventory incomplete")
    _require(general.get("quadratic_expected") == 0
             and general.get("quadratic_included") == 0
             and general.get("extra_yukawa_expected") == 0
             and general.get("extra_yukawa_included") == 0,
             "Unexpected B(-3) quadratic or extra Yukawa inventory")
    accepted = general.get("interactions", [])
    matched_accepted = summary.get("AcceptedInteractions", [])
    _require(isinstance(accepted, list) and len(accepted) == len(set(accepted)),
             "Invalid general interaction inventory")
    _require(set(accepted) == set(matched_accepted),
             "Matchete matching used a different UV interaction inventory")
    _require(GENERAL_INTERACTION in accepted,
             "The B(-3) extra allowed quartic is absent from matched UV")
    _require(summary.get("BuildStatus") == "Success"
             and summary.get("MatchingStatus") == "Success"
             and summary.get("WeinbergExtractionStatus") == "Success",
             "The general two-stage matching or Weinberg extraction failed")
    expected_plan = [["F"], ["S1", "S2"]]
    _require(stages.get("Status") == "Success"
             and stages.get("RequestedThresholdCount") == 2
             and stages.get("MatchedStageCount") == 2
             and summary.get("ThresholdPlan") == expected_plan,
             "Expected exact F -> (S1,S2) sequential threshold plan")
    _require(interval.get("IntegratedFields") == ["F"]
             and interval.get("ActiveHeavyFields") == ["S1", "S2"],
             "Incorrect intermediate EFT field-content metadata")
    first_fields = set(interval.get("TreeDiagnostics", {}).get("Fields", []))
    _require({"NewScalar1", "NewScalar2"} <= first_fields
             and "NewFermion" not in first_fields,
             "First EFT must retain both scalars and remove the fermion")
    _require(wilson.get("TreeWilsonTermCount", 0) > 0
             and any(t.get("Chirality") == "PL" for t in wilson.get("TreeWilsonTerms", [])),
             "Tree-level PL Wilson boundary is missing")
    _require(quartic.get("ScalarQuarticTermCount", 0) > 0,
             "Intermediate scalar quartic seed missing")

    # Couplings metadata may be empty in the current exporter: read the actual
    # InputForm expression and its explicit CG registry instead.
    general_terms = [
        item for item in quartic.get("ScalarQuarticTerms", [])
        if GENERAL_COUPLING in item.get("TermInputForm", "")
    ]
    _require(len(general_terms) == 2,
             "Expected both general quartic and its Hermitian conjugate in EFT1")
    _require(all(GENERAL_CG in item.get("CGNames", []) for item in general_terms),
             "Extra quartic CG not recorded in the EFT1 seed")
    registered = {
        item.get("Name") for item in quartic.get("CGRegistry", [])
    }
    _require(GENERAL_CG in registered,
             "General extra-quartic CG tensor was not exported")
    _require(all(item.get("TensorInputForm") and item.get("RepsInputForm")
                 for item in quartic.get("CGRegistry", []) if item.get("Name") == GENERAL_CG),
             "General CG registry has no tensor or representation")

    _require(general_eft1_rge.get("status") == "SuccessCandidateMassBeta"
             and general_eft1_rge.get("metadata", {}).get("GeneralScalarModel") == MODEL
             and general_eft1_rge.get("metadata", {}).get("GeneralScalarStage") == "EFT1",
             "General symbolic EFT1 renormalisable beta export is missing")
    _require(GENERAL_RGBETA_COUPLING in general_eft1_rge.get("report_betas", {}),
             "General EFT1 RGBeta omitted the additional quartic coupling")

    return {
        "model": MODEL,
        "general_uv_interaction_count": len(accepted),
        "general_quartic": GENERAL_INTERACTION,
        "matched_uv_agrees_with_general_symbolic_inventory": True,
        "fermion_absent_both_scalars_active_after_first_threshold": True,
        "tree_pl_wilson_term_count": sum(
            term.get("Chirality") == "PL" for term in wilson["TreeWilsonTerms"]
        ),
        "general_quartic_and_conjugate_in_eft1": len(general_terms),
        "general_quartic_cg_registry_present": True,
        "general_renormalisable_eft1_beta_includes_extra_quartic": True,
        "physical_matchete_rgbeta_mapping_certified": False,
    }


def make_tensor_metadata() -> dict[str, Any]:
    """Exact B(-3) field quantum numbers, NOT restricted RGBeta beta functions."""
    return {
        "status": "MetadataOnlyGeneralEFT1",
        "metadata": {
            "SharedScalar": False,
            "dS": None,
            "dS1": DIMENSIONS[0],
            "dS2": DIMENSIONS[1],
            "dF": DIMENSIONS[2],
            "alpha": ALPHA,
            "YS1": "-3/2",
            "YS2": "-1/2",
            "YS": "Null",
            "IntegratedField": "F",
            "ActiveBSMFields": ["S1", "S2"],
        },
        "note": (
            "Only field metadata is read by ScalarOnlyWilsonTensorRGE. "
            "The renormalisable beta functions are the separate GENERAL "
            "symbolic eft1_rgbeta.json, not the restricted RGBeta export."
        ),
    }


def _continuation(matching_dir: Path):
    from RGE.running.intermediate.ScalarThresholdMatching import ScalarThresholdContinuation

    base = matching_dir / "debug/fresh_kernel_stage_2"
    names = (
        "input_tree.wxf", "input_loop.wxf", "input_transition_tree.wxf",
        "input_transition_full.wxf", "input_cg_registry.wxf",
    )
    paths = [base / name for name in names]
    for path in paths:
        if not path.is_file():
            raise FileNotFoundError(f"Original stage-2 continuation missing: {path}")
    return ScalarThresholdContinuation(*paths)


def run_pilot(
    matching_dir: Path,
    general_dir: Path,
    *,
    mu_high: str = "MF",
    mu_low: str = "MS",
    rge_only: bool = False,
    validate_threshold: bool = False,
) -> dict[str, Any]:
    matching_dir = Path(matching_dir).resolve()
    general_dir = Path(general_dir).resolve()
    work = matching_dir / "data/general_eft1_pilot"
    work.mkdir(parents=True, exist_ok=True)
    summary_path = work / "pilot_summary.json"
    result: dict[str, Any] = {
        "status": "Started",
        "model": MODEL,
        "scope": "general B(-3) fixed-one-loop sequential EFT1 pilot",
        "matching_directory": str(matching_dir),
        "general_symbolic_directory": str(general_dir),
        "matching_source": "general symmetry-allowed UV Lagrangian",
        "renormalisable_running_mapping_certified": False,
        "physical_C5_numerics_required": False,
    }
    _written_json(summary_path, result)

    try:
        result["boundary_audit"] = validate_general_boundary(matching_dir, general_dir)
        print("[general B(-3)] General UV, EFT1 field, quartic/CG, and Wilson seeds validated.", flush=True)

        metadata_path = work / "general_eft1_tensor_metadata.json"
        _written_json(metadata_path, make_tensor_metadata())

        # Imports are intentionally lazy: source audits and lightweight tests
        # work even where optional symbolic physics dependencies are absent.
        import sympy as sp
        from RGE.running.intermediate.ScalarOnlyWilsonTensorRGE import (
            build_scalar_only_rge_context, run_scalar_only_wilson_rge,
        )
        from RGE.running.intermediate.FlavorWilsonBoundary import (
            export_full_flavor_wilson_boundary,
        )
        from RGE.running.intermediate.DirectWeinbergRunning import (
            build_direct_weinberg_running, export_direct_weinberg_insertion,
        )

        wilson_path = matching_dir / "data/eft1_after_F_wilson_seed.json"
        quartic_path = matching_dir / "data/eft1_after_F_scalar_quartic_seed.json"
        context = build_scalar_only_rge_context(wilson_path, quartic_path, metadata_path)
        general_quartic_entries = [
            (key, expr) for key, expr in context.quartic.entries.items()
            if sp.Symbol(GENERAL_COUPLING) in sp.sympify(expr).free_symbols
            or GENERAL_COUPLING in sp.sstr(expr)
        ]
        _require(bool(general_quartic_entries),
                 "General quartic disappeared during Matchete -> real-scalar RGE tensor conversion")
        result["general_quartic_tensor_component_count"] = len(general_quartic_entries)
        print(f"[general B(-3)] Extra quartic enters {len(general_quartic_entries)} real-tensor components.", flush=True)

        rge_path = work / "general_eft1_wilson_rge.json"
        rge = run_scalar_only_wilson_rge(
            wilson_path, quartic_path, metadata_path,
            output_path=rge_path, seed_support_only=False,
        )
        _require(rge.get("status") == "Success", "Wilson tensor RGE failed")
        check = rge.get("weinberg_subspace_validation") or {}
        _require(check.get("matches_weinberg_subspace") is True,
                 "Direct LLSS -> Weinberg beta failed the Weinberg-subspace test")
        result["rge"] = {
            "path": str(rge_path),
            "initial_component_count": rge.get("initial_independent_component_count"),
            "nonzero_beta_component_count": rge.get("nonzero_beta_component_count"),
            "generated_component_count": rge.get("generated_component_count"),
            "weinberg_subspace_validation": check,
            "includes_general_quartic_in_tensor": True,
        }
        print("[general B(-3)] Complete one-loop EFT1 Wilson RGE and Weinberg-subspace projection computed.", flush=True)

        flavor_seed_path = work / "general_eft1_flavor_seed.json"
        flavor = export_full_flavor_wilson_boundary(wilson_path, flavor_seed_path)
        _require(flavor.get("status") == "Success", "Full-flavor Wilson seed failed")

        transport_path = work / "general_eft1_direct_weinberg_running.json"
        transport = build_direct_weinberg_running(
            flavor_seed_path, rge_path, mu_high, mu_low,
            output_path=transport_path,
        )
        _require(transport.get("status") == "Success", "Full-flavor C5 transport failed")
        direct_coefficient = transport["direct_weinberg_beta"]["flavor_blind_prefactor"]
        # The existing Wolfram emitter has a known explicit allowed-coupling
        # dictionary. Never emit an unregistered extra quartic symbol silently.
        _require(not re.search(r"\blambdaT3Extra[QM][A-Za-z0-9_]*\b", direct_coefficient),
                 "Direct C5 beta needs new general-quartic Coupling[]/Bar[] emitter support")

        insertion_path = work / "general_eft1_direct_weinberg_insertion.wl"
        insertion = export_direct_weinberg_insertion(transport_path, insertion_path)
        _require(insertion.get("status") == "Success", "Wolfram C5 insertion export failed")
        result["direct_c5"] = {
            "transport_path": str(transport_path),
            "insertion_path": str(insertion_path),
            "full_flavor": True,
            "equal_scale_running_vanishes": transport.get("equal_scale_running_vanishes"),
            "only_direct_C5_insertion_at_one_loop": True,
            "heavy_LLSS_self_running_only_for_diagnostics": True,
        }
        print("[general B(-3)] Full-flavor direct C5 running insertion written.", flush=True)
        result["status"] = "RGEAndDirectC5InsertionReady"
        _written_json(summary_path, result)

        if rge_only:
            return result

        from RGE.running.intermediate.ScalarThresholdMatching import (
            resume_scalar_threshold_with_running,
        )
        from RGE.matching.FinalWeinbergCoefficient import (
            build_final_weinberg_coefficient,
        )

        continuation = _continuation(matching_dir)
        threshold_script = ROOT / "Lagrangian/RunThresholdStage.wl"
        threshold_path = work / "threshold_2_with_general_eft1_running.wxf"
        resumed = resume_scalar_threshold_with_running(
            output_dir=matching_dir,
            running_insertion=insertion_path,
            run_threshold_script=threshold_script,
            d_s1=DIMENSIONS[0], d_s2=DIMENSIONS[1],
            d_f=DIMENSIONS[2], alpha=ALPHA,
            eft_order=5, loop_order=1,
            continuation=continuation,
            result_path=threshold_path,
            validation_mode=validate_threshold,
        )
        _require(resumed.get("status") == "Success",
                 f"General scalar-threshold running resume failed: {resumed.get('stdout_tail')}")
        result["resumed_scalar_threshold"] = {
            "status": resumed["status"],
            "result_path": resumed["result_path"],
            "running_insertion_loaded": resumed.get("running_insertion_loaded"),
            "direct_weinberg_carried_separately": resumed.get("direct_weinberg_carried_separately"),
            "equal_scale_vanishes": resumed.get("direct_weinberg_equal_scale_vanishes"),
        }
        print("[general B(-3)] Resumed scalar matching with EFT1 running insertion.", flush=True)

        pole_path = work / "c5_pole_rge_consistency.json"
        _require(pole_path.is_file(), "Missing scalar-threshold pole/RGE diagnostic")
        pole = _json(pole_path)
        pole_equal = pole.get("DirectLogEqualsTwicePoleResidueOneGeneration")
        result["pole_rge_consistency"] = {
            "path": str(pole_path),
            "direct_log_equals_twice_pole_residue": pole_equal,
        }
        _written_json(summary_path, result)
        _require(pole_equal is True,
                 "General B(-3) hard-pole/direct-Weinberg-log consistency failed")
        final_path = work / "general_sequential_final_c5.json"
        final = build_final_weinberg_coefficient(
            threshold_c5_path=work / "c5_threshold_with_eft1_running.txt",
            flavor_transport_path=transport_path,
            pole_rge_consistency_path=pole_path,
            output_path=final_path,
        )
        _require(final.get("status") == "Success", "Final C5 assembly failed")
        combined = final.get("combined", {})
        result["final_c5"] = {
            "path": str(final_path),
            "status": final.get("status"),
            "structurally_full_flavor": combined.get("structurally_full_flavor"),
            "ready_for_physical_majorana_numerics": combined.get(
                "ready_for_physical_majorana_numerics"
            ),
            "not_claimed_as_general_16_model_certification": True,
        }
        result["status"] = "GeneralBMinus3SequentialEFT1Executed"
        _written_json(summary_path, result)
        print(f"[general B(-3)] Saved final audit: {summary_path}", flush=True)
        return result
    except Exception as exc:
        result["status"] = "Failed"
        result["error"] = f"{type(exc).__name__}: {exc}"
        _written_json(summary_path, result)
        raise


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--matching-dir", type=Path, default=DEFAULT_MATCHING_DIR)
    ap.add_argument("--general-dir", type=Path, default=DEFAULT_GENERAL_DIR)
    ap.add_argument("--mu-high", default="MF")
    ap.add_argument("--mu-low", default="MS")
    ap.add_argument("--rge-only", action="store_true")
    ap.add_argument("--validate-threshold", action="store_true")
    args = ap.parse_args()
    result = run_pilot(
        args.matching_dir, args.general_dir,
        mu_high=args.mu_high, mu_low=args.mu_low,
        rge_only=args.rge_only,
        validate_threshold=args.validate_threshold,
    )
    print(json.dumps({
        "model": result["model"], "status": result["status"],
        "rge": result.get("rge", {}),
        "direct_c5": result.get("direct_c5", {}),
        "resumed_scalar_threshold": result.get("resumed_scalar_threshold", {}),
        "final_c5": result.get("final_c5", {}),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
