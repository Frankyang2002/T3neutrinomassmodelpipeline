"""Run general D(-2) scalar matching on the physical real-S2 EFT1 boundary.

Existing stage-1 matching, one-generation RGE and full-flavour direct running
must already have succeeded.  This driver never reruns stage-1 matching.

The Wolfram stage keeps pure-SM direct LLSS->LLHH running separate from the
scalar hard Match and does not claim a fully assembled flavour-indexed C5.

Usage:
  python -X utf8 -m scripts.RunGeneralDMinus2ScalarThreshold
  python -X utf8 -m scripts.RunGeneralDMinus2ScalarThreshold --preflight-only
  python -X utf8 -m scripts.RunGeneralDMinus2ScalarThreshold --profile-tree
"""
from __future__ import annotations

import argparse
import json
import os
import signal
import time
from pathlib import Path
import subprocess
import sys
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MATCHING = ROOT / "output/general_sequential_pilot/D_m2_stage1"
STAGE_SUBDIR = Path("data/general_d_m2_scalar_threshold")
SCRIPT = ROOT / "Lagrangian/RunGeneralDMinus2ScalarThreshold.wl"
WXF_FILES = (
    "eft1_after_F_realbasis_tree.wxf",
    "eft1_after_F_realbasis_loop.wxf",
    "eft1_after_F_realbasis_transition_tree.wxf",
    "eft1_after_F_realbasis_transition_full.wxf",
)
EXPECTED_SCALARS = ["NewScalar1", "T3S2R", "T3S2I"]


def _load(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(f"Required D(-2) artifact missing: {path}")
    payload = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(payload, dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return payload


def _save(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def preflight(matching_dir: Path, script: Path = SCRIPT) -> dict[str, Any]:
    """Reject restricted or incompatible boundaries *before* calling Matchete."""
    matching_dir = Path(matching_dir)
    data = matching_dir / "data"
    first = _load(data / "general_d_m2_stage1_summary.json")
    if first.get("status") != "GeneralDMinus2RealMassBasisFermionThresholdVerified":
        raise ValueError("The general real-basis fermion threshold is not verified")
    uv = first.get("general_uv", {})
    boundary = first.get("first_matching", {})
    if (uv.get("extra_quartics") != 5 or uv.get("neutral_scalar_bilinears") != 1
            or uv.get("extra_yukawas") != 1 or boundary.get("tree_wilson_terms") != 8
            or boundary.get("remaining_physical_BSM_scalars") != EXPECTED_SCALARS):
        raise ValueError("The saved first-threshold inventory is not full general D(-2)")

    rge = _load(data / "general_d_m2_eft1_rge/pilot_summary.json")
    if rge.get("status") != "GeneralDMinus2RealBasisWeinbergBetaVerifiedOneGeneration":
        raise ValueError("The real-basis Weinberg beta has not passed its audit")
    full = _load(data / "general_d_m2_full_flavor/pilot_summary.json")
    if full.get("status") != "GeneralDMinus2FullFlavorDirectRunningPrepared":
        raise ValueError("The direct full-flavour running stage has not succeeded")
    direct = _load(data / "general_d_m2_full_flavor/full_flavor_direct_running.json")
    if (direct.get("status") != "Success"
        or not direct.get("one_generation_reduction_matches")
        or not direct.get("symmetric_under_external_flavor_exchange")
        or not direct.get("equal_scale_running_vanishes")
        or direct.get("real_scalar_hard_threshold_performed") is not False):
        raise ValueError("Direct running fails symmetry, reduction or power-counting checks")
    pre = _load(data / "general_d_m2_full_flavor/scalar_threshold_preflight.json")
    if (pre.get("status") != "FullFlavorDirectInsertionPrepared_ScalarHardThresholdNotRun"
            or pre.get("scalar_matching_integrate_out") != EXPECTED_SCALARS
            or pre.get("scalar_masses_matchete") != ["MS1", "T3S2RMass", "T3S2IMass"]
            or pre.get("physical_mass_relations_imposed") is not False
            or pre.get("insertion_is_pure_SM_and_is_not_to_be_one_loop_matched") is not True
            or pre.get("insertion_coefficient_omits_hbar") is not True):
        raise ValueError("Scalar matching preflight contains an incompatible mass/power convention")
    mass = _load(matching_dir / "scalar_mass_basis.json")
    if (mass.get("RealMassSquared") != "MS2^2-2*Abs[kappa]"
        or mass.get("ImagMassSquared") != "MS2^2+2*Abs[kappa]"
        or mass.get("RealFields") != ["T3S2R", "T3S2I"]):
        raise ValueError("D(-2) real-field mass convention is inconsistent")

    for name in WXF_FILES:
        if not (data / name).is_file():
            raise FileNotFoundError(f"Physical EFT1 WXF continuation missing: {name}")
    for name in ("eft1_after_F_realbasis_wilson_tensor_seed.json",
                 "eft1_after_F_realbasis_quartic_tensor_seed.json"):
        if not (data / name).is_file():
            raise FileNotFoundError(f"Physical EFT1 tensor/CG registry missing: {name}")
    insertion = data / "general_d_m2_full_flavor/full_flavor_weinberg_insertion.wl"
    if not insertion.is_file():
        raise FileNotFoundError(f"Direct running insertion missing: {insertion}")
    insertion_text = insertion.read_text(encoding="utf-8")
    for symbol in ("EFT1D2WilsonRunningHeavyInsertion = 0",
                   "EFT1D2WilsonRunningWeinbergCoefficient", "EFT1D2C12",
                   "EFT1D2C1Conj"):
        if symbol not in insertion_text:
            raise ValueError(f"Direct running insertion does not define {symbol!r}")
    if not Path(script).is_file():
        raise FileNotFoundError(f"Local Wolfram scalar-threshold driver missing: {script}")
    return {
        "status": "GeneralDMinus2ScalarThresholdPreflightPassed",
        "model": "D_-2",
        "physical_scalar_fields": EXPECTED_SCALARS,
        "physical_scalar_mass_relations_imposed": False,
        "common_fermion_mass_assumed": True,
        "one_loop_direct_running_separate_from_scalar_match": True,
        "full_flavour_hard_c5_certified": False,
        "matching_dir": str(matching_dir.resolve()),
        "insertion_file": str(insertion.resolve()),
        "wolfram_script": str(Path(script).resolve()),
    }


def make_command(matching_dir: Path, insertion_file: Path,
                 wolframscript: str = "wolframscript", script: Path = SCRIPT) -> list[str]:
    return [wolframscript, "-file", str(Path(script).resolve()),
            str(Path(matching_dir).resolve()), str(Path(insertion_file).resolve())]


# Each probe uses a fresh kernel: a timed-out Mathematica Match[] can leave
# private package state inconsistent, so do not reuse that kernel for a second
# diagnostic match. The full general matching is NEVER replaced by a probe.
PROFILE_SMOKE_PROBES = (
    "inspect", "sm-only", "s2r-only", "s1-only", "heavy-only", "free",
)
PROFILE_PROBES = (
    "inspect", "free", "scalar-standard", "scalar-general",
    "scalar-with-y2", "scalar-with-y2Conj",
)
# A separate diagnostic sequence that isolates the effect of sending the
# stand-alone Standard Model Lagrangian through Matchete's Match[] parser.
# This does not change any physics in the authoritative scalar threshold.
PROFILE_BSM_PROBES = (
    "inspect", "bsm-free", "bsm-standard", "bsm-general",
    "bsm-with-y2", "bsm-with-y2Conj", "bsm-with-both",
)
# Old probe names are retained for explicitly targeted comparisons only.
OPTIONAL_PROFILE_PROBES = ("full-tree",)
LEGACY_PROFILE_PROBES = ("renormalisable", "with_y2", "with_y2Conj")
PROFILE_BASELINE_PROBES = ("free", "scalar-standard", "scalar-general")
PROFILE_SUBDIR = STAGE_SUBDIR / "profile_tree"


def make_profile_command(
    matching_dir: Path, insertion_file: Path, probe: str,
    probe_timeout_seconds: int = 90,
    *, wolframscript: str = "wolframscript", script: Path = SCRIPT,
) -> list[str]:
    if probe not in PROFILE_PROBES + PROFILE_BSM_PROBES + PROFILE_SMOKE_PROBES + OPTIONAL_PROFILE_PROBES + LEGACY_PROFILE_PROBES:
        raise ValueError(f"Unknown D(-2) profile probe: {probe!r}")
    if not 5 <= probe_timeout_seconds <= 3600:
        raise ValueError("probe_timeout_seconds must be between 5 and 3600")
    return [
        wolframscript, "-file", str(Path(script).resolve()),
        str(Path(matching_dir).resolve()), str(Path(insertion_file).resolve()),
        f"profile:{probe}", str(probe_timeout_seconds),
    ]


def _stop_process_tree(process: subprocess.Popen) -> None:
    """Terminate a timed-out profile's *child Wolfram kernel*, not just its launcher."""
    if process.poll() is not None:
        return
    if os.name == "nt":
        # Windows wolframscript launches wolfram.exe as a child. Terminating
        # only wolframscript would leave that CPU-intensive kernel alive.
        try:
            subprocess.run(
                ["taskkill", "/PID", str(process.pid), "/T", "/F"],
                capture_output=True, text=True, timeout=20, check=False,
            )
        except (OSError, subprocess.TimeoutExpired):
            process.kill()
    else:
        try:
            os.killpg(process.pid, signal.SIGTERM)
        except ProcessLookupError:
            return
        try:
            process.wait(timeout=5)
        except subprocess.TimeoutExpired:
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
    try:
        process.wait(timeout=10)
    except subprocess.TimeoutExpired:
        process.kill()
        process.wait()


def _execute_profile_probe(
    matching_dir: Path, insertion_file: Path, probe: str,
    probe_timeout_seconds: int, startup_allowance_seconds: int,
    wolframscript: str, script: Path,
) -> dict[str, Any]:
    """Bounded diagnostic Match in a new kernel, with visible heartbeat."""
    profile_dir = matching_dir / PROFILE_SUBDIR
    profile_dir.mkdir(parents=True, exist_ok=True)
    log = profile_dir / f"profile_{probe}.log"
    status_path = profile_dir / f"wolfram_profile_status_{probe}.json"
    result_path = profile_dir / f"probe_result_{probe}.json"
    # A stale success from an earlier run must never be mistaken for this run.
    for target in (status_path, result_path):
        target.unlink(missing_ok=True)
    command = make_profile_command(
        matching_dir, insertion_file, probe, probe_timeout_seconds,
        wolframscript=wolframscript, script=script,
    )
    environment = os.environ.copy()
    environment["T3_EXPERIMENTAL_GENERAL_QUADRATICS"] = "1"
    environment["T3_EXPERIMENTAL_GENERAL_FERMIONS"] = "1"
    start = time.monotonic()
    # External safety cap: TimeConstrained[] may not interrupt a blocked
    # native/frontend operation, so use a separate OS process per probe.
    hard_limit = probe_timeout_seconds + startup_allowance_seconds
    timed_out = False
    print(f"[profile] {probe}: launching fresh Wolfram kernel; "
          f"Match limit {probe_timeout_seconds}s, total cap {hard_limit}s",
          flush=True)
    with log.open("wb") as stream:
        flags = getattr(subprocess, "CREATE_NEW_PROCESS_GROUP", 0) if os.name == "nt" else 0
        process = subprocess.Popen(
            command, cwd=ROOT, stdout=stream, stderr=subprocess.STDOUT,
            env=environment, creationflags=flags,
            start_new_session=(os.name != "nt"),
        )
        last_heartbeat = start
        while process.poll() is None:
            now = time.monotonic()
            if now - start >= hard_limit:
                timed_out = True
                print(f"[profile] {probe}: external limit reached; "
                      "terminating this probe's process tree", flush=True)
                _stop_process_tree(process)
                break
            if now - last_heartbeat >= 20:
                print(f"[profile] {probe}: active for {int(now-start)}s "
                      "(Match[] does not expose internal progress)", flush=True)
                last_heartbeat = now
            time.sleep(0.5)
        exit_code = process.poll()
    result: dict[str, Any] = {
        "probe": probe,
        "status": "ExternalTimeout" if timed_out else "Exited",
        "exit_code": exit_code,
        "elapsed_seconds": round(time.monotonic() - start, 2),
        "wolfram_log": str(log),
        "external_timeout_seconds": hard_limit,
        "diagnostic_only": True,
        "full_general_scalar_threshold_completed": False,
    }
    if status_path.is_file():
        try:
            result["wolfram_status"] = _load(status_path)
        except (OSError, ValueError) as exc:
            result["status_error"] = str(exc)
    if result_path.is_file():
        try:
            result["match_result"] = _load(result_path)
        except (OSError, ValueError) as exc:
            result["result_error"] = str(exc)
    if not timed_out:
        if exit_code != 0:
            result["status"] = "WolframFailed"
        elif probe == "inspect":
            result["status"] = (
                "Inspected" if result.get("wolfram_status", {}).get("Status")
                == "ProfileTreeInspected" else "OutputIncomplete"
            )
        elif result.get("match_result", {}).get("Status") in ("MatchReturned", "TimedOut"):
            result["status"] = result["match_result"]["Status"]
        else:
            result["status"] = "OutputIncomplete"
    print(f"[profile] {probe}: {result['status']} "
          f"(exit {exit_code}); log: {log}", flush=True)
    return result


def run_tree_profile(
    matching_dir: Path,
    *, wolframscript: str = "wolframscript", script: Path = SCRIPT,
    probe_timeout_seconds: int = 90,
    startup_allowance_seconds: int = 180,
    include_full_tree: bool = False,
    smoke_only: bool = False,
    bsm_only: bool = False,
) -> dict[str, Any]:
    """Profile scalar input sectors WITHOUT starting [B] or [C]."""
    matching_dir = Path(matching_dir).resolve()
    profile_dir = matching_dir / PROFILE_SUBDIR
    profile_dir.mkdir(parents=True, exist_ok=True)
    if smoke_only and bsm_only:
        raise ValueError("Choose only one diagnostic profile mode")
    summary_path = profile_dir / (
        "profile_smoke_summary.json" if smoke_only else
        "profile_bsm_summary.json" if bsm_only else "profile_tree_summary.json"
    )
    summary: dict[str, Any] = {
        "status": "ProfileStarted", "model": "D_-2",
        "scope": ("diagnostic registered-field free-sector Match[] smoke tests"
                  if smoke_only else
                  "diagnostic heavy-scalar Match[] without stand-alone L_SM"
                  if bsm_only else
                  "diagnostic physical real-scalar EFT1 tree Match[] decomposition"),
        "diagnostic_only": True,
        "full_general_threshold_matched": False,
        "general_interactions_removed_from_authoritative_matching": False,
        "probe_timeout_seconds": probe_timeout_seconds,
        "smoke_only": smoke_only,
        "bsm_only": bsm_only,
        "startup_allowance_seconds": startup_allowance_seconds,
        "probes": [],
    }
    _save(summary_path, summary)
    try:
        if not 5 <= probe_timeout_seconds <= 3600:
            raise ValueError("--probe-timeout-seconds must be between 5 and 3600")
        if not 30 <= startup_allowance_seconds <= 3600:
            raise ValueError("--startup-allowance-seconds must be between 30 and 3600")
        ready = preflight(matching_dir, script=script)
        _save(profile_dir / "profile_preflight.json", ready)
        probes = list(PROFILE_SMOKE_PROBES if smoke_only else
                      PROFILE_BSM_PROBES if bsm_only else PROFILE_PROBES)
        if include_full_tree and not (smoke_only or bsm_only):
            probes.extend(OPTIONAL_PROFILE_PROBES)
        for probe in probes:
            record = _execute_profile_probe(
                matching_dir, Path(ready["insertion_file"]), probe,
                probe_timeout_seconds, startup_allowance_seconds,
                wolframscript, script,
            )
            summary["probes"].append(record)
            _save(summary_path, summary)
            if probe == "inspect" and record["status"] != "Inspected":
                summary["status"] = "ProfileInspectionFailed"
                summary["skipped_probes"] = probes[1:]
                break
            # Smoke probes are intentionally independent: continue even when
            # an isolated sector fails or times out, to localise the problem.
            if smoke_only:
                continue
            # If a *smaller* renormalisable baseline already blocks Match[],
            # adding higher-order interactions cannot diagnose that bottleneck.
            # Record the finding and skip the unnecessary later probes.
            if probe in (PROFILE_BASELINE_PROBES + ("bsm-free", "bsm-standard", "bsm-general")
                          if bsm_only else PROFILE_BASELINE_PROBES) and record["status"] != "MatchReturned":
                summary["status"] = "TreeProfileRecorded"
                summary["baseline_bottleneck_probe"] = probe
                summary["baseline_bottleneck_status"] = record["status"]
                summary["skipped_probes"] = probes[probes.index(probe) + 1:]
                break
        else:
            summary["status"] = "TreeProfileRecorded"
        if smoke_only and summary["status"] == "TreeProfileRecorded":
            summary["status"] = "SmokeProfileRecorded"
        if bsm_only and summary["status"] == "TreeProfileRecorded":
            summary["status"] = "BSMProfileRecorded"
    except (OSError, RuntimeError, ValueError) as exc:
        summary["status"] = "Failed"
        summary["error"] = f"{type(exc).__name__}: {exc}"
    finally:
        _save(summary_path, summary)
    return summary


def run(matching_dir: Path, wolframscript: str = "wolframscript",
        *, preflight_only: bool = False, script: Path = SCRIPT) -> dict[str, Any]:
    matching_dir = Path(matching_dir).resolve()
    output_dir = matching_dir / STAGE_SUBDIR
    output_dir.mkdir(parents=True, exist_ok=True)
    summary_path = output_dir / "pilot_summary.json"
    state: dict[str, Any] = {
        "status": "Started", "model": "D_-2",
        "scope": "full-general real-mass-basis D(-2) scalar hard threshold at fixed one loop",
        "physical_mass_relations_imposed": False,
        "full_flavour_hard_c5_certified": False,
        "direct_running_double_matched": False,
        "full_wilson_self_rge_performed": False,
        "full_renormalisable_rgbeta_mapping_certified": False,
    }
    _save(summary_path, state)
    try:
        ready = preflight(matching_dir, script=script)
        _save(output_dir / "scalar_threshold_preflight_result.json", ready)
        if preflight_only:
            state["status"] = "GeneralDMinus2ScalarThresholdPreflightPassed"
            state["preflight_only"] = True
            return state

        command = make_command(matching_dir, Path(ready["insertion_file"]),
                               wolframscript, script=script)
        env = os.environ.copy()
        env["T3_EXPERIMENTAL_GENERAL_QUADRATICS"] = "1"
        env["T3_EXPERIMENTAL_GENERAL_FERMIONS"] = "1"
        log = output_dir / "wolfram_scalar_threshold.log"
        with log.open("w", encoding="utf-8") as stream:
            try:
                process = subprocess.Popen(command, cwd=ROOT,
                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                    env=env, text=True, encoding="utf-8", errors="replace")
            except OSError as exc:
                raise RuntimeError(f"Cannot start wolframscript: {exc}") from exc
            assert process.stdout is not None
            for line in process.stdout:
                stream.write(line)
                stream.flush()
                print(line, end="", flush=True)
            exit_code = process.wait()
        state["wolfram_exit_code"] = exit_code
        state["wolfram_log_path"] = str(log)
        wolfram_status = _load(output_dir / "wolfram_scalar_threshold_status.json")
        state["wolfram_status"] = wolfram_status
        if exit_code != 0:
            raise RuntimeError(
                f"Wolfram scalar threshold returned {exit_code}: "
                f"{wolfram_status.get('Failure', wolfram_status.get('Status'))}")
        if wolfram_status.get("Status") != "HardScalarThresholdMatched_FlavorAssemblyPending":
            raise RuntimeError("Wolfram did not certify completed scalar hard matching")
        if not all(wolfram_status.get(k) is True for k in
            ("TreeMatching", "OneLoopMatching", "InheritedLoopTreePropagation",
             "HardSMEFTExported", "AllScalarFieldsIntegratedOut")):
            raise RuntimeError("Wolfram scalar threshold is missing required stage checks")
        for filename in ("scalar_hard_smeft.wxf", "scalar_hard_c5_inputform.txt",
                         "direct_running_c5_inputform.txt", "smeft_c5_assembly.wl"):
            if not (output_dir / filename).is_file():
                raise FileNotFoundError(f"Completed scalar threshold missing {filename}")
        state["status"] = "GeneralDMinus2ScalarHardThresholdMatched_FlavorAssemblyPending"
        state["all_scalars_integrated_out"] = True
        state["full_flavour_hard_c5_certified"] = False
    except Exception as exc:
        state["status"] = "Failed"
        state["error"] = f"{type(exc).__name__}: {exc}"
    finally:
        _save(summary_path, state)
    return state


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matching-dir", type=Path, default=DEFAULT_MATCHING)
    parser.add_argument("--wolframscript", default="wolframscript")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--preflight-only", action="store_true")
    mode.add_argument("--profile-tree", action="store_true",
                      help="Incremental free / standard / general scalar / LLSS tree probes; no [B]/[C]")
    mode.add_argument("--profile-smoke", action="store_true",
                      help="Independent SM/free-heavy scalar Match probes; no quartics or [B]/[C]")
    mode.add_argument("--profile-bsm", action="store_true",
                      help="Incremental heavy-scalar Match[] without stand-alone L_SM; diagnostic only")
    parser.add_argument("--probe-timeout-seconds", type=int, default=90)
    parser.add_argument("--startup-allowance-seconds", type=int, default=180)
    parser.add_argument("--include-full-tree-probe", action="store_true",
                        help="Also profile the original full tree Match[] (bounded)")
    args = parser.parse_args(argv)
    if args.profile_tree or args.profile_smoke or args.profile_bsm:
        result = run_tree_profile(
            args.matching_dir, wolframscript=args.wolframscript,
            probe_timeout_seconds=args.probe_timeout_seconds,
            startup_allowance_seconds=args.startup_allowance_seconds,
            include_full_tree=args.include_full_tree_probe,
            smoke_only=args.profile_smoke,
            bsm_only=args.profile_bsm,
        )
        target = Path(args.matching_dir).resolve() / PROFILE_SUBDIR / (
            "profile_smoke_summary.json" if args.profile_smoke else
            "profile_bsm_summary.json" if args.profile_bsm else
            "profile_tree_summary.json"
        )
        print(f"{result['status']} -> {target}")
        if result["status"] not in ("TreeProfileRecorded", "SmokeProfileRecorded", "BSMProfileRecorded"):
            print("ERROR:", result.get("error", result["status"]), file=sys.stderr)
            return 1
        return 0
    result = run(args.matching_dir, wolframscript=args.wolframscript,
                 preflight_only=args.preflight_only)
    target = Path(args.matching_dir).resolve() / STAGE_SUBDIR / "pilot_summary.json"
    print(f"{result['status']} -> {target}")
    if result["status"] == "Failed":
        print("ERROR:", result["error"], file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
