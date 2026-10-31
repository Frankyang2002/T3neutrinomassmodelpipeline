#!/usr/bin/env python3
"""Experimental local Matchete F-threshold flavour-Wilson tensor transfer.

One F threshold integrates 3 diagonal heavy generations simultaneously.
The scalar-threshold experiment uses ONLY stage-1 tree EFT1 Wilson tensors.
The saved stage-1 hard one-loop EFT is deliberately not propagated, and no
EFT1 RGE is included; it cannot yield a finite hierarchical prediction.
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
MODELS = ("A_p0", "A_m2", "D_p0", "D_m2", "E_p0", "E_m2")
SUCCESS = {"prepare": "FermionThresholdWilsonTreeContracted",
           "tree": "ScalarThresholdWilsonSoftMatchedNotRenormalized",
           "soft": "ScalarThresholdWilsonSoftMatchedNotRenormalized"}
EXTRAS = {key: ("y1Conj" if key.endswith("p0") else "y2Conj") for key in MODELS}
OUT = ROOT / "Reports/output/workstream_a/eft1_wilson_tensor"
NATIVE = ROOT / "Lagrangian/RunT3EFT1WilsonHierarchyNative.wl"
EXPECTED_SOURCES = {
    "Lagrangian/RunMatching.wl": "2a3d7415695213b501e53023b19e8a4a7ee54275daa918ba2ce10064d0d87567",
    "Lagrangian/RunThresholdStage.wl": "9c81fc99982421c656f80050cd78dbdabf0e85d81c7062baf96f8cdd91af3d4c",
}

class Blocked(RuntimeError):
    pass


def require(condition: bool, msg: str) -> None:
    if not condition:
        raise Blocked(msg)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), "Missing file: " + str(path))
    obj = json.loads(path.read_text(encoding="utf-8-sig"))
    require(isinstance(obj, dict), "Not a JSON object: " + str(path))
    return obj


def preflight(root: Path = ROOT, *, verify_real: bool = True) -> dict[str, Any]:
    hashes = {}
    for rel, checksum in EXPECTED_SOURCES.items():
        path = root / rel
        require(path.is_file(), "Missing production source: " + rel)
        actual = sha(path)
        require(actual == checksum,
                "Production Matchete source differs from uploaded reference: " + rel
                + "; do not use experimental copied stage until it is reviewed")
        hashes[rel] = actual
    paths = (
        "Lagrangian/RunT3EFT1WilsonHierarchyNative.wl",
        "Lagrangian/RunT3EFT1WilsonStageLauncher.wl",
        "Lagrangian/RunThresholdStageEFT1Wilson.wl",
    )
    for rel in paths:
        require((root / rel).is_file(), "Missing experimental source: " + rel)
    driver = (root / paths[0]).read_text(encoding="utf-8-sig")
    launcher = (root / paths[1]).read_text(encoding="utf-8-sig")
    child = (root / paths[2]).read_text(encoding="utf-8-sig")
    checks = [
        ('Matchete loaded before function body', driver.index('Needs["Matchete`"]') < driver.index('wOutput =')),
        ('Three flavours explicit', 'T3ModelFromClass[wClass,wAlpha[wModel],3]' in driver),
        ('All F together', 'RunSequentialT3Matching[' in driver and '{{"F"}}' in driver),
        ('Both Wilson channels', 'T3EFT1W12' in driver and 'T3EFT1WExtra' in driver),
        ('Per-generation mass retained', 'Coupling[MF,{r:Index[_,NFlavor]},0]' in driver),
        ('Hard one-loop withheld', 'stage1_hard_one_loop_omitted_from_stage2.wxf' in driver),
        ('No fake EFT1 running', '"eft1_wilson_running_included"->False' in driver),
        ('No physical certification', '"renormalized_hierarchical_c5_certified"->False' in driver),
        ('Isolated fresh child', '"RunThresholdStageEFT1Wilson.wl"' in launcher),
        ('Wilson metadata in child', 'DefineCoupling[T3EFT1W12, Indices -> {Flavor, Flavor}' in child
             and 'DefineCoupling[T3EFT1WExtra, Indices -> {Flavor, Flavor}' in child),
        ('Preserved original stage2 tree matching', 'RunT3Matching[transitionTree, eftOrder, 0]' in child),
        ('Preserved original stage2 loop matching', 'RunT3Matching[transitionTree, eftOrder, 1]' in child),
    ]
    for desc, ok in checks:
        require(ok, "Experimental source contract failed: " + desc)
    # The experimental child must be a FULL copy of the uploaded production
    # threshold source apart from ONE inserted, local-only registration block.
    production = (root / "Lagrangian/RunThresholdStage.wl").read_text(encoding="utf-8-sig")
    before, after = child.split("(* Experimental contracted EFT1 Wilson tensors.", 1)
    require("Print[\"[fresh kernel] Restored contracted EFT1 Wilson tensor flavour indices.\"];" in after,
            "Missing experimental child registration end")
    after = after.split('Print["[fresh kernel] Restored contracted EFT1 Wilson tensor flavour indices."];', 1)[1]
    require(before + after.lstrip("\n") == production,
            "Experimental scalar stage differs from production outside Wilson registration")
    # And the launcher must match the original public fresh-stage function
    # exactly, apart from its function name and stage filename.
    matching = (root / "Lagrangian/RunMatching.wl").read_text(encoding="utf-8-sig")
    start = matching.index('RunFreshThresholdStage[\n')
    end = matching.index('RunSequentialT3Matching[\n', start)
    expected_launcher = matching[start:end].replace(
        'RunFreshThresholdStage[', 'RunT3EFT1WilsonFreshThresholdStage[', 1
    ).replace('"RunThresholdStage.wl"', '"RunThresholdStageEFT1Wilson.wl"')
    require(expected_launcher.strip() in launcher,
            "Experimental fresh-stage launcher no longer matches production source")
    native = {}
    if verify_real:
        d = root / "Reports/output/workstream_a/extended_c5_native/native"
        for model in MODELS:
            obj = read_json(d / f"{model}_h3.json")
            require(obj.get("model") == model
                    and obj.get("status") == "NativeOneStepExtendedWeinbergExtracted"
                    and obj.get("additional_yukawa") == EXTRAS[model]
                    and obj.get("heavy_flavours") == 3
                    and obj.get("native_full_coefficient_contains_extra_name") is True,
                    "Six-model native one-step input not validated: " + model)
            native[model] = sha(d / f"{model}_h3.json")
    return {"status": "EFT1WilsonTensorSourceGatePassed",
            "source_sha256": hashes, "one_step_native_sha256": native,
            "absolute_coupling_dictionary_certified": False,
            "hard_one_loop_and_eft1_wilson_running_certified": False}


def check_native(obj: dict, model: str, variant: str, mode: str) -> dict:
    for key, expected in (
        ("format", "t3-eft1-flavour-wilson-soft-v1"),
        ("model", model), ("variant", variant), ("mode", mode),
        ("status", SUCCESS[mode]),
        ("heavy_generations_integrated_together", 3),
        ("independent_heavy_masses_preserved_in_wilson_definition", True),
        ("tree_level_eft1_wilson_contraction_applied", True),
    ):
        require(obj.get(key) == expected,
                f"Invalid native output {model}/{variant}/{mode}: {key}")
    for flag in ("stage1_hard_one_loop_inherited_at_stage2",
                 "eft1_wilson_running_included",
                 "renormalized_hierarchical_c5_certified",
                 "absolute_matchete_rgbeta_dictionary_certified",
                 "physical_neutrino_prediction_certified"):
        require(obj.get(flag) is False, "Unjustified certification: " + flag)
    wc = obj.get("wilson_counts", {})
    require(wc.get("baseline") == 2 and
            wc.get("extra") == (2 if variant == "seeded" else 0) and
            wc.get("unmatched_mass") == 0,
            "Incomplete or invalid heavy-index Wilson contraction")
    if mode in ("tree", "soft"):
        require(obj.get("scalar_threshold_native_passed") is True and
                obj.get("stage2_loop_order") == (0 if mode == "tree" else 1) and
                obj.get("eft1_wilson_running_insertion_added") is False,
                "Scalar result flags not consistent with mode")
    return obj


def check_definitions(obj: dict, variant: str, model: str) -> None:
    require(obj.get("T3EFT1W12[p,q]") ==
            "Sum[y1[p,r] y2[q,r] / MF[r], {r,1,3}]",
            "Wrong baseline Wilson definition")
    name = "y1Conj" if model.endswith("p0") else "y2Conj"
    expected = ("Sum[y1Conj[p,r] y2[q,r] / MF[r], {r,1,3}]"
                if name == "y1Conj" else
                "Sum[y1[p,r] y2Conj[q,r] / MF[r], {r,1,3}]")
    require(obj.get("T3EFT1WExtra[p,q]") == expected,
            "Wrong extra Wilson definition")
    require("Three diagonal independent masses" in obj.get("mass_assumption", ""),
            "Mass basis / three-generation assumption missing")


def run_variant(model: str, variant: str, mode: str, args: argparse.Namespace,
                directory: Path) -> dict[str, Any]:
    directory.mkdir(parents=True, exist_ok=True)
    log = directory / "stdout.log"
    native_path = directory / "native.json"
    env = os.environ.copy()
    env.update({"T3_EXPERIMENTAL_GENERAL_FERMIONS": "1",
                "T3_EXPERIMENTAL_FLAVOR_MASSES": "1",
                "T3_EXPERIMENTAL_GENERAL_QUARTICS": "1",
                "T3_EXPERIMENTAL_GENERAL_QUADRATICS": "0"})
    if args.use_existing and native_path.is_file():
        ret = 0
        print(f"{model}/{variant}: reading saved native artifact", flush=True)
    else:
        # Prevent stale success from a prior run surviving an aborted child.
        if native_path.exists():
            native_path.unlink()
        cmd = [args.wolframscript, "-file", str(NATIVE), str(directory),
               model[0], model.split("_")[1], variant, mode]
        print(f"{model}/{variant}: F-threshold Wilson {mode}...", flush=True)
        try:
            proc = subprocess.run(cmd, cwd=ROOT, env=env, capture_output=True,
                                  text=True, errors="replace", timeout=args.timeout_seconds,
                                  check=False)
            ret = int(proc.returncode)
            log.write_text(proc.stdout + "\nSTDERR:\n" + proc.stderr,
                           encoding="utf-8")
        except (OSError, subprocess.TimeoutExpired) as exc:
            log.write_text(str(exc) + "\n", encoding="utf-8")
            ret = 124
    if ret != 0 or not native_path.is_file():
        note = ""
        if native_path.is_file():
            try:
                note = str(read_json(native_path).get("reason", ""))
            except (Blocked, ValueError):
                pass
        if log.is_file():
            note += "\n" + "\n".join(log.read_text(encoding="utf-8", errors="replace").splitlines()[-16:])
        raise Blocked(f"{model}/{variant} native exit={ret}: {note}")
    obj = check_native(read_json(native_path), model, variant, mode)
    check_definitions(read_json(directory / "wilson_definitions.json"), variant, model)
    for name in ("stage1_tree_wilson.wxf",
                 "stage1_hard_one_loop_omitted_from_stage2.wxf"):
        require((directory / name).is_file(), "Missing Matchete artifact: " + name)
    if mode != "prepare":
        require((directory / "stage2_weinberg_unrenormalized.wxf").is_file(),
                "Missing stage-2 Weinberg expression")
    return {"status": obj["status"], "native": str(native_path),
            "native_sha256": sha(native_path), "wilson_counts": obj["wilson_counts"]}


def run(args: argparse.Namespace) -> int:
    if args.self_test:
        preflight(ROOT, verify_real=False)
        print("PASS: exact local matching-source copies, six-model catalogue, experimental EFT1 Wilson gates")
        return 0
    result = preflight(ROOT, verify_real=True)
    out = args.output.resolve() if args.output else OUT
    out.mkdir(parents=True, exist_ok=True)
    if args.preflight:
        (out / "preflight.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
        print("PASS: six native reductions, exact source provenance, EFT1 Wilson tensor gates")
        print("Uncertified: full hierarchical C5, EFT1 running, physical Matchete-RGBeta dictionary")
        return 0
    require(args.pilot or args.all, "Specify --pilot or --all")
    models = MODELS if args.all else MODELS[:1]
    mode = args.mode
    summary = {"format": "t3-eft1-flavour-wilson-soft-batch-v1",
               "mode": mode, "status": "NativeMatchingPartialOrBlocked",
               "one_fermion_threshold_for_three_generations": True,
               "scalar_threshold_soft_matching_test_only": mode == "soft",
               "hard_one_loop_inherited": False, "eft1_wilson_running_evolved": False,
               "physical_matching_certified": False, "models": []}
    successes = 0
    for index, model in enumerate(models, 1):
        row: dict[str, Any] = {"model": model, "status": "Blocked", "variants": {}}
        print(f"[{index}/{len(models)}] {model}, {mode}", flush=True)
        for variant in ("seeded", "zero"):
            try:
                row["variants"][variant] = run_variant(
                    model, variant, mode, args, out / mode / model / variant)
            except Blocked as exc:
                row["variants"][variant] = {"status": "Blocked", "reason": str(exc)}
                print(str(exc), flush=True)
                break
        if len(row["variants"]) == 2 and all(
            v.get("status") == SUCCESS[mode] for v in row["variants"].values()
        ):
            row["status"] = "BothNativeEFT1WilsonVariantsCompleted"
            successes += 1
            print(f"PASS {model}: seeded/zero {mode}, independent heavy masses", flush=True)
        summary["models"].append(row)
        # Continue to next model only in the explicitly requested full batch.
        if row["status"] == "Blocked" and args.pilot:
            break
    summary["status"] = ("AllRequestedNativeEFT1WilsonRunsCompleted"
                         if successes == len(models) else "NativeMatchingPartialOrBlocked")
    target = out / mode / "summary.json"
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"SUMMARY {mode}: {successes}/{len(models)}; {summary['status']}")
    print("Report: " + str(target))
    print("Physical C5: NOT CERTIFIED (hard one-loop and EFT1 running not propagated)")
    return 0 if successes == len(models) else 4


def parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__)
    choices = ap.add_mutually_exclusive_group()
    choices.add_argument("--pilot", action="store_true")
    choices.add_argument("--all", action="store_true")
    choices.add_argument("--preflight", action="store_true")
    choices.add_argument("--self-test", action="store_true")
    ap.add_argument("--mode", choices=("prepare", "tree", "soft"), default="prepare")
    ap.add_argument("--wolframscript", default="wolframscript")
    ap.add_argument("--timeout-seconds", type=int, default=3600)
    ap.add_argument("--output", type=Path)
    ap.add_argument("--use-existing", action="store_true")
    return ap


def main() -> int:
    try:
        return run(parser().parse_args())
    except (Blocked, OSError, ValueError, IndexError, KeyError) as exc:
        print("BLOCKED: " + str(exc), file=sys.stderr)
        return 4


if __name__ == "__main__":
    raise SystemExit(main())
