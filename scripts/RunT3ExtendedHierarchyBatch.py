#!/usr/bin/env python3
"""Native Matchete hierarchical extra-Yukawa diagnostic (local-only).

The batch computes the F-only EFT1 with and without the independent Yukawa,
optionally follows the existing fresh-kernel scalar threshold, and compares
native expressions in Wolfram. It does NOT replace/approximate the missing
EFT1 Wilson-coefficient RGE or declare a finite physical hierarchical C5.

Usage (repo root):
  python scripts/RunT3ExtendedHierarchyBatch.py --self-test
  python scripts/RunT3ExtendedHierarchyBatch.py --preflight
  python scripts/RunT3ExtendedHierarchyBatch.py --pilot --stage1
  python scripts/RunT3ExtendedHierarchyBatch.py --all --stage1
  python scripts/RunT3ExtendedHierarchyBatch.py --pilot --hierarchy
  python scripts/RunT3ExtendedHierarchyBatch.py --all --hierarchy
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "Reports/output/workstream_a/extended_c5_hierarchy"
MODELS = ("A_p0", "A_m2", "D_p0", "D_m2", "E_p0", "E_m2")
EXTRAS = {m: ("y1Conj" if m.endswith("p0") else "y2Conj") for m in MODELS}
SUCCESS = {
    "stage1": "FermionEFT1NativeMatched",
    "hierarchy": "NativeHierarchicalC5ExtractedNotRenormalized",
}
NATIVE = ROOT / "Lagrangian/RunT3ExtendedHierarchyNative.wl"
COMPARE = ROOT / "Lagrangian/CompareT3ExtendedHierarchy.wl"


class GateError(RuntimeError):
    pass


def need(test: bool, msg: str) -> None:
    if not test:
        raise GateError(msg)


def json_read(path: Path) -> dict[str, Any]:
    need(path.is_file(), "Missing source evidence: " + str(path))
    obj = json.loads(path.read_text(encoding="utf-8"))
    need(isinstance(obj, dict), "Expected JSON object: " + str(path))
    return obj


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_gate(root: Path = ROOT, *, require_prior: bool = True) -> dict[str, Any]:
    files = {
        "F_index": ("Lagrangian/model/T3Fields.wl", "DefineFlavorIndex["),
        "F_mass": ("Lagrangian/model/T3Fields.wl", 'T3MassSpec["F"'),
        "builder": ("Lagrangian/model/LagrangianBuilder.wl", "BuildT3Lagrangian["),
        "extra_yukawa": ("Lagrangian/interactions/T3Topology.wl", "BuildT3ExtraYukawaCandidates["),
        "matching": ("Lagrangian/RunMatching.wl", "RunSequentialT3Matching["),
        "fresh_stage": ("Lagrangian/RunMatching.wl", "RunFreshThresholdStage["),
        "threshold": ("Lagrangian/RunThresholdStage.wl", "authoritativeThresholdC5"),
        "native_driver": ("Lagrangian/RunT3ExtendedHierarchyNative.wl", "SetT3HeavyFields[{\"F\"}]"),
        "compare_driver": ("Lagrangian/CompareT3ExtendedHierarchy.wl", "stage2_c5_ordered"),
    }
    hashes = {}
    for label, (rel, snippet) in files.items():
        path = root / rel
        need(path.is_file(), f"Required source missing: {rel}")
        text = path.read_text(encoding="utf-8-sig")
        need(snippet in text, f"Source contract changed: {label}: {rel}")
        hashes[rel] = digest(path)
    # Protect the actual two-threshold perturbative bookkeeping. Do not
    # use the stage1 expression as though it were a final Weinberg C5.
    matching = (root / "Lagrangian/RunMatching.wl").read_text(encoding="utf-8-sig")
    for needle in ("CompactExplicitEFT", "CanonicalTransitionEFT",
                   "InheritedLoopTreeMatchingFailed", "hbar Coefficient[",
                   '"TreeLagrangian"', '"OneLoopCorrection"'):
        need(needle in matching, "Native matching perturbative stage missing: " + needle)
    stage = (root / "Lagrangian/RunThresholdStage.wl").read_text(encoding="utf-8-sig")
    for needle in ('"AuthoritativeThresholdC5"', '"C5PoleRGEConsistency"',
                   '"EFT1WilsonRunningLoaded"'):
        need(needle in stage, "Threshold stage diagnostic missing: " + needle)
    native = (root / "Lagrangian/RunT3ExtendedHierarchyNative.wl").read_text(encoding="utf-8-sig")
    need(native.index('Needs["Matchete`"]') < native.index('RunSequentialT3Matching['),
         "Matchete load must precede the parsed matching body")
    need('RunSequentialT3Matching[uv,model,plan,5,1,hOut,hModel]' in native,
         "Native pipeline must call existing sequential Matchete matching")
    need('"eft1_Wilson_RGE_evolved"->False' in native,
         "Do not promote absent EFT1 Wilson running to completed")
    previous = {}
    if require_prior:
        base = root / "Reports/output/workstream_a/extended_c5_native/native"
        for model in MODELS:
            row = json_read(base / f"{model}_h3.json")
            need(row.get("model") == model and
                 row.get("status") == "NativeOneStepExtendedWeinbergExtracted" and
                 row.get("additional_yukawa") == EXTRAS[model] and
                 row.get("heavy_flavours") == 3 and
                 row.get("native_full_coefficient_contains_extra_name") is True and
                 row.get("absolute_matchete_rgbeta_dictionary_certified") is False,
                 "Prior native one-step extra Yukawa certificate absent/changed: " + model)
            previous[model] = digest(base / f"{model}_h3.json")
    return {
        "status": "ExistingSequentialSourceAndNativeOneStepEvidenceVerified",
        "source_sha256": hashes,
        "prior_one_step_sha256": previous,
        "physical_dictionary_certified": False,
        "hierarchical_extended_wilson_running_certified": False,
    }


def check_native(report: dict, model: str, mode: str, variant: str,
                 mass_mode: str = "indexed") -> dict:
    need(model in MODELS and mode in SUCCESS and variant in ("seeded", "zero"),
         "Invalid comparison identity")
    for label, expected in (("format", "t3-extra-yukawa-hierarchical-native-v1"),
                            ("model", model), ("mode", mode),
                            ("variant", variant), ("heavy_flavours", 3),
                            ("status", SUCCESS[mode])):
        need(report.get(label) == expected,
             f"Native Matchete report mismatch: {model}/{mode}/{variant}/{label}")
    need(report.get("mass_mode", "indexed") == mass_mode,
         "Native report belongs to a different heavy-mass mode")
    need(report.get("independent_heavy_masses_certified") is not True,
         "Independent heavy masses not certified by the diagnostic")
    if mass_mode == "common":
        need(report.get("mass_mode_common_is_diagnostic_only") is True,
             "Common mass result lacks diagnostic-only declaration")
    need(report.get("fermion_threshold_integrated_out") is True,
         "F threshold did not complete")
    need(report.get("scalar_threshold_integrated_out") is (mode == "hierarchy"),
         "Scalar threshold flag inconsistent with executed mode")
    need(report.get("eft1_Wilson_RGE_evolved") is False and
         report.get("absolute_matchete_rgbeta_dictionary_certified") is False and
         report.get("physical_low_energy_prediction_certified") is False and
         report.get("full_scalar_quadratic_basis_included") is False,
         "Native result contains unjustified physical certification")
    need(report.get("stage_count") == (1 if mode == "stage1" else 2),
         "Wrong number of Matchete stages")
    if mode == "hierarchy":
        need(report.get("stage2_c5_is_finite_renormalized") is False and
             report.get("EFT1_Wilson_running_insertion_added") is False,
             "EFT1 running or pole subtraction claimed without evidence")
    return report


def check_difference(report: dict, model: str, mode: str) -> dict:
    need(report.get("status") == "NativeEFTDifferenceExported" and
         report.get("model") == model and report.get("mode") == mode,
         "Invalid seeded-minus-zero comparison")
    expected = ["stage1_tree", "stage1_one_loop", "stage1_full"]
    if mode == "hierarchy":
        expected.extend(["stage2_full", "stage2_c5_ordered"])
    need(report.get("seeded_vs_uv_zero_expressions_compared") == expected,
         "Incomplete matched-EFT comparison")
    need(isinstance(report.get("difference_not_proven_zero_by_PossibleZeroQ"), dict) and
         set(report["difference_not_proven_zero_by_PossibleZeroQ"]) == set(expected),
         "Incomplete native difference witness")
    need(all(isinstance(v, bool) for v in report["difference_not_proven_zero_by_PossibleZeroQ"].values()),
         "Difference witness is not Boolean")
    for key in ("absolute_coupling_dictionary_certified",
                "full_renormalized_hierarchical_C5_certified", "EFT1_running_inserted",
                "cross_kernel_dummy_index_canonicalisation_certified"):
        need(report.get(key) is False, "Comparison claims unsupported status: " + key)
    return report


def invoke(command: list[str], log: Path, timeout: int, env: dict[str, str]) -> int:
    log.parent.mkdir(parents=True, exist_ok=True)
    try:
        process = subprocess.run(command, cwd=ROOT, env=env,
                                 text=True, stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT, errors="replace",
                                 timeout=timeout, check=False)
    except (OSError, subprocess.TimeoutExpired) as exc:
        log.write_text(str(exc) + "\n", encoding="utf-8")
        return 124
    log.write_text(process.stdout, encoding="utf-8")
    return process.returncode


def excerpt(log: Path, lines: int = 15) -> str:
    if not log.is_file():
        return "No native log saved"
    return "\n".join(log.read_text(encoding="utf-8", errors="replace").splitlines()[-lines:])


def run_model(model: str, mode: str, output: Path, args: argparse.Namespace) -> dict:
    need(model in MODELS, "Unsupported model")
    base = output / mode / model
    row: dict[str, Any] = {"model": model, "mode": mode,
                          "status": "Blocked", "variants": {},
                          "physical_hierarchical_prediction_certified": False}
    env = os.environ.copy()
    mass_mode = "common" if getattr(args, "common_mass_diagnostic", False) else "indexed"
    env.update({"T3_EXPERIMENTAL_GENERAL_FERMIONS": "1",
                "T3_EXPERIMENTAL_GENERAL_QUARTICS": "1",
                "T3_EXPERIMENTAL_FLAVOR_MASSES": "0" if mass_mode == "common" else "1",
                "T3_EXPERIMENTAL_GENERAL_QUADRATICS": "0"})
    row["mass_mode"] = mass_mode
    row["independent_heavy_masses_certified"] = False
    for variant in ("seeded", "zero"):
        dest = base / variant
        native = dest / "native.json"
        if args.use_existing and native.is_file():
            code = 0
            print(f"{model}/{variant}: loading previous native report", flush=True)
        else:
            dest.mkdir(parents=True, exist_ok=True)
            command = [args.wolframscript, "-file", str(NATIVE), str(dest),
                       model[0], model.split("_")[1], variant, mode, mass_mode]
            print(f"{model}/{variant}: Matchete {mode} matching...", flush=True)
            code = invoke(command, dest / "stdout.log", args.timeout_seconds, env)
        if code != 0 or not native.is_file():
            row["variants"][variant] = {"status": "Blocked",
                                        "native_exit": code,
                                        "log": str(dest / "stdout.log")}
            print(excerpt(dest / "stdout.log"), flush=True)
            return row
        try:
            obj = json_read(native)
            check_native(obj, model, mode, variant, mass_mode)
        except (GateError, KeyError, ValueError) as exc:
            row["variants"][variant] = {"status": "Blocked",
                                        "reason": str(exc),
                                        "native_report": obj if 'obj' in locals() else None}
            print(f"{model}/{variant}: {exc}", flush=True)
            return row
        row["variants"][variant] = {"status": SUCCESS[mode],
                                    "native_file": str(native),
                                    "sha256": digest(native)}
        print(f"{model}/{variant}: {SUCCESS[mode]}", flush=True)
    difference = base / "difference.json"
    if args.use_existing and difference.is_file():
        code = 0
    else:
        command = [args.wolframscript, "-file", str(COMPARE),
                   str(base / "seeded"), str(base / "zero"), str(base), mode]
        print(f"{model}: comparing native seeded/zero EFTs...", flush=True)
        code = invoke(command, base / "comparison.log", args.timeout_seconds, env)
    if code != 0 or not difference.is_file():
        row["comparison"] = {"status": "Blocked", "native_exit": code,
                             "log": str(base / "comparison.log")}
        print(excerpt(base / "comparison.log"), flush=True)
        return row
    try:
        diff = check_difference(json_read(difference), model, mode)
    except (GateError, ValueError, KeyError) as exc:
        row["comparison"] = {"status": "Blocked", "reason": str(exc)}
        return row
    row["comparison"] = {"status": diff["status"],
                         "difference_not_proven_zero_by_PossibleZeroQ": diff["difference_not_proven_zero_by_PossibleZeroQ"],
                         "native_difference_file": str(difference)}
    row["status"] = "NATIVE_STAGE1_AND_DIFFERENCE_COMPLETE" if mode == "stage1" else \
                    "NATIVE_SEQUENTIAL_DIFFERENCE_UNRENORMALIZED"
    return row


def self_test() -> None:
    need(len(MODELS) == 6 and len(set(MODELS)) == 6, "Incorrect catalogue")
    template = {
        "format": "t3-extra-yukawa-hierarchical-native-v1",
        "model": "A_p0", "mode": "stage1", "variant": "seeded",
        "heavy_flavours": 3, "status": SUCCESS["stage1"],
        "fermion_threshold_integrated_out": True,
        "scalar_threshold_integrated_out": False,
        "stage_count": 1, "eft1_Wilson_RGE_evolved": False,
        "absolute_matchete_rgbeta_dictionary_certified": False,
        "physical_low_energy_prediction_certified": False,
        "full_scalar_quadratic_basis_included": False,
    }
    check_native(template, "A_p0", "stage1", "seeded")
    for modify in ({"status": "Success"},
                   {"fermion_threshold_integrated_out": False},
                   {"eft1_Wilson_RGE_evolved": True},
                   {"absolute_matchete_rgbeta_dictionary_certified": True},
                   {"heavy_flavours": 1}):
        try:
            check_native({**template, **modify}, "A_p0", "stage1", "seeded")
        except GateError:
            pass
        else:
            raise AssertionError("Unsafe native report accepted: " + str(modify))
    hierarchy = {**template, "mode": "hierarchy", "stage_count": 2,
                 "scalar_threshold_integrated_out": True,
                 "status": SUCCESS["hierarchy"],
                 "stage2_c5_is_finite_renormalized": False,
                 "EFT1_Wilson_running_insertion_added": False}
    check_native(hierarchy, "A_p0", "hierarchy", "seeded")
    for key in ("stage2_c5_is_finite_renormalized", "EFT1_Wilson_running_insertion_added"):
        try:
            check_native({**hierarchy, key: True}, "A_p0", "hierarchy", "seeded")
        except GateError:
            pass
        else:
            raise AssertionError("Uncertified hierarchical claim accepted: " + key)
    print("PASS: six-model native threshold registry, F-only EFT1, fixed-order gate, seeded/zero negative controls")


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    action = p.add_mutually_exclusive_group(required=True)
    for opt in ("self-test", "preflight", "pilot", "all"):
        action.add_argument("--" + opt, action="store_true")
    mode = p.add_mutually_exclusive_group()
    mode.add_argument("--stage1", action="store_true")
    mode.add_argument("--hierarchy", action="store_true")
    p.add_argument("--use-existing", action="store_true")
    p.add_argument("--common-mass-diagnostic", action="store_true",
                   help="Three degenerate heavy flavours; diagnostic ONLY, not 3 independent masses")
    p.add_argument("--wolframscript", default="wolframscript")
    p.add_argument("--timeout-seconds", type=int, default=3600)
    p.add_argument("--output-dir", type=Path, default=OUT)
    args = p.parse_args()
    try:
        self_test()
        if args.self_test:
            return 0
        evidence = source_gate(ROOT, require_prior=True)
        if args.common_mass_diagnostic and args.output_dir == OUT:
            args.output_dir = OUT / "common_mass_diagnostic"
        args.output_dir.mkdir(parents=True, exist_ok=True)
        (args.output_dir / "source_gate.json").write_text(
            json.dumps(evidence, indent=2) + "\n", encoding="utf-8")
        print("PASS: local sequential matching source and six native one-step coefficients verified")
        if args.preflight:
            return 0
        need(args.timeout_seconds > 0, "Invalid Wolfram timeout")
        mode_name = "hierarchy" if args.hierarchy else "stage1"
        selected = MODELS if args.all else MODELS[:1]
        results = []
        for idx, model in enumerate(selected, 1):
            print(f"[{idx}/{len(selected)}] {model} / {mode_name}", flush=True)
            results.append(run_model(model, mode_name, args.output_dir, args))
        succeeded = sum(row["status"].startswith("NATIVE_") for row in results)
        summary = {
            "format": "t3-extra-yukawa-native-sequential-batch-v1",
            "status": "AllNativeDifferenceRunsCompleted" if succeeded == len(results)
                      else "NativeHierarchicalMatchingPartialOrBlocked",
            "mode": mode_name, "successful": succeeded, "requested": len(results),
            "mass_mode": "common" if args.common_mass_diagnostic else "indexed",
            "independent_heavy_masses_certified": False,
            "common_mass_diagnostic_only": bool(args.common_mass_diagnostic),
            "models": results,
            "absolute_matchete_rgbeta_dictionary_certified": False,
            "EFT1_Wilson_coefficient_running_implemented": False,
            "finite_renormalized_hierarchical_C5_certified": False,
            "low_energy_predictions_physical_certified": False,
        }
        out = args.output_dir / mode_name / "summary.json"
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
        print(f"SUMMARY {mode_name}: {succeeded}/{len(results)} native comparisons; {summary['status']}")
        print("Report:", out)
        print("Physical coupling dictionary, EFT1 Wilson RG running, finite hierarchical C5: NOT CERTIFIED")
        if args.common_mass_diagnostic:
            print("MASS MODE: three DEGENERATE flavours, diagnostic ONLY; independent heavy masses NOT represented.")
        return 0 if succeeded == len(results) else 4
    except (GateError, OSError, ValueError, KeyError, json.JSONDecodeError) as exc:
        print("BLOCKED: " + str(exc), file=sys.stderr)
        return 3


if __name__ == "__main__":
    sys.exit(main())
