"""End-to-end local general D(-2) fermion -> scalar sequential C5 run.

Default: rebuild physical UV -> EFT1, the real-basis EFT1 beta, full-flavour
running, the scalar hard one-loop EFT, the independent inherited term, and a
provisional flavour-symmetric C5 benchmark. Each Wolfram stage runs in its
own kernel. Existing legacy scalar-threshold files are not overwritten.

Use --reuse-boundary to skip only the UV/F/RGE producers and rerun the entire
second threshold starting from the already verified first-threshold WXF.
No GitHub writes. No claim of fully certified SMEFT operator basis or physical
model-specific Clebsch-Gordan normalization.
"""
from __future__ import annotations

import argparse
from collections import deque
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import signal
import subprocess
import sys
import time
from typing import Any, Callable

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "output/general_sequential_pilot/D_m2_stage1"
DEFAULT_CONFIG = ROOT / "benchmarks/general_d_m2_rough_benchmark_finitepart.json"
REL_PROFILE = Path("data/general_d_m2_scalar_threshold/profile_tree")
REL_RUN = Path("data/general_d_m2_scalar_threshold/sequential_run")
BOUNDARY_FILES = (
    "eft1_after_F_realbasis_tree.wxf",
    "eft1_after_F_realbasis_loop.wxf",
    "eft1_after_F_realbasis_transition_tree.wxf",
    "eft1_after_F_realbasis_transition_full.wxf",
)
STAGE_REPORTS = {
    "scalar_hard_one_loop": "profile_full_scalar_one_loop_result.json",
    "components": "profile_fixed_order_components_result.json",
    "operator_structure": "profile_c5_flavour_structure_result.json",
    "running_symmetry": "profile_c5_running_symmetry_audit_result.json",
    "candidate_assembly": "profile_c5_flavour_candidate_result.json",
    "warsaw_reference": "profile_c5_normalization_reference_result.json",
    "benchmark": "profile_c5_rough_benchmark_result.json",
    "exact_pole": "profile_c5_pole_cancellation_result.json",
}
EXPECTED_STATUSES = {
    "scalar_hard_one_loop": "FullEFT1OneLoopPostprocessed",
    "components": "ComponentsAudited",
    "operator_structure": "WeinbergStructureAudited",
    "running_symmetry": "RunningSymmetryAudited",
    "candidate_assembly": "C5FlavourCandidateExported",
    "warsaw_reference": "WarsawNormalizationReferenceAudited",
    "benchmark": "ProvisionalFinitePartBenchmarkEvaluated",
    "exact_pole": "ExactNF1PoleCancellationVerified",
}


class PipelineFailure(RuntimeError):
    pass


def read_object(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise PipelineFailure(f"Missing report: {path}")
    try:
        raw = path.read_text(encoding="utf-8-sig")
        data = json.loads(raw)
    except (OSError, ValueError, UnicodeError) as exc:
        raise PipelineFailure(f"Cannot decode JSON {path}: {exc}") from exc
    if not isinstance(data, dict):
        raise PipelineFailure(f"Expected JSON object: {path}")
    return data


def json_write(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    temporary.replace(path)


def boundary_hashes(out: Path) -> dict[str, str]:
    results: dict[str, str] = {}
    for name in BOUNDARY_FILES:
        file = out / "data" / name
        if not file.is_file():
            raise PipelineFailure(f"First-threshold WXF missing: {file}")
        h = hashlib.sha256()
        with file.open("rb") as fh:
            for block in iter(lambda: fh.read(1024 * 1024), b""):
                h.update(block)
        results[name] = h.hexdigest()
    return results


def check(condition: bool, message: str) -> None:
    if not condition:
        raise PipelineFailure(message)


def validate_report(stage: str, report: dict[str, Any]) -> None:
    expected = EXPECTED_STATUSES[stage]
    check(report.get("Status") == expected,
          f"{stage}: expected {expected!r}, got {report.get('Status')!r}: "
          f"{report.get('Reason', '')}")
    if stage == "scalar_hard_one_loop":
        check(report.get("OneLoopRawMatchReturned") is True and
              report.get("MatchResultContainsHbar") is True and
              report.get("ProcessedEFTHasHbar") is True and
              report.get("RemainingHeavyFields") == [] and
              report.get("WeinbergExtractionStatus") == "Success" and
              report.get("WeinbergTermCount") == 68 and
              len(report.get("ProcessingStages", [])) == 4 and
              all(x.get("Status") == "Returned" for x in report["ProcessingStages"]),
              "Scalar hard EFT failed output/one-loop/Weinberg validation")
    elif stage == "components":
        check(report.get("HardC5LinearInHbar") is True and
              report.get("DirectRunningHasExplicitHbar") is False and
              report.get("HardWeinbergTermCount") == 68 and
              report.get("InheritedPureSMWeinbergTermCount") == 12,
              "Fixed-order components changed")
    elif stage == "operator_structure":
        check(report.get("AllCandidateTermsHaveTwoLeptonsTwoHiggs") is True and
              report.get("AllCandidateTermsHaveUnambiguousChirality") is True and
              report.get("HardC5IsLinearInHbar") is True and
              report.get("InheritedC5IsLinearInHbar") is True,
              "Weinberg field-content or chirality audit failed")
    elif stage == "running_symmetry":
        check(report.get("SymmetricModuloSummationDummies") is True and
              report.get("RunningFlavourSymmetryCertifiedModuloDummyRenaming") is True and
              report.get("AlphaNormalizedResidualTermCount") == 0,
              "Running flavour symmetry not established modulo dummy variables")
    elif stage == "candidate_assembly":
        check(report.get("FixedOrderC5CandidateAssembled") is True and
              report.get("EachHolomorphicTermHasTwoDifferentLeptonFlavourLabels") is True and
              report.get("HardSymmetricByConstruction") is True and
              report.get("InheritedSymmetricByConstruction") is True and
              report.get("RunningSymmetricModuloSummationDummies") is True and
              report.get("TotalCandidateLinearInHbar") is True,
              "Candidate C5 flavour or fixed-order check failed")
    elif stage == "warsaw_reference":
        # The Wolfram reference exporter deliberately serializes its exact
        # InputForm factors as strings: "1" and "2" (not JSON numbers).
        # Accept exactly those canonical values or their numeric equivalents;
        # do not coerce arbitrary text or accept booleans as integers.
        unit = report.get("WarsawOrderedCoefficientWithReferenceCouplingsSetToOne")
        symmetric_unit = report.get("WarsawSymmetricUnitWilsonFactorAssumingEqualPQandQP")
        exact_unit = unit == "1" if isinstance(unit, str) else (
            type(unit) in (int, float) and unit == 1)
        exact_symmetric_unit = symmetric_unit == "2" if isinstance(symmetric_unit, str) else (
            type(symmetric_unit) in (int, float) and symmetric_unit == 2)
        check(report.get("ProjectExtractorRecognizesWarsawWeinbergOperator") is True and
              report.get("ProjectExtractorLinearOnWarsawReference") is True and
              report.get("WarsawOrderedCoefficientUnitWilsonFactorNumeric") is True and
              exact_unit and exact_symmetric_unit and
              report.get("WarsawOperatorNeutrinoVEVComponentIsPlusVSquaredOverTwo") is True,
              "Warsaw reference normalization check changed")
    elif stage == "benchmark":
        masses = report.get("NeutrinoMasses_eV")
        check(isinstance(masses, list) and len(masses) == 3 and
              all(isinstance(x, (float, int)) and math.isfinite(x) and x >= 0 for x in masses) and
              masses == sorted(masses) and
              report.get("FinitePartOnly") is True and
              report.get("RenormalizationCertified") is False and
              isinstance(report.get("C5RelativeSymmetryResidual"), (float, int)) and
              report["C5RelativeSymmetryResidual"] < 1e-8,
              "Benchmark masses, symmetry or finite-part metadata invalid")
    elif stage == "exact_pole":
        check(report.get("NF") == 1 and
              report.get("ExactPoleCancellationVerifiedForNF1") is True and
              report.get("ExactZeroByTogether") is True and
              report.get("ReducedResidueTermCount") == 0,
              "NF=1 exact pole cancellation not established")


def validate_mass_audit(profile: Path) -> None:
    a = read_object(profile / "profile_full_scalar_one_loop_audit.json")
    check(a.get("Status") == "Passed" and
          a.get("AllPhysicalScalarsRegisteredHeavy") is True and
          a.get("MassPowerCountingPromotionOnly") is True and
          a.get("PureHeavyQuadraticAlphaEqualToFreeLag") is True and
          a.get("ExactReversibleMassPromotion") is True and
          a.get("CheckLagrangianPassed") is True and
          a.get("SMAlphaComparisonCertified") is True and
          a.get("NoLinearHeavyScalarSources") is True and
          a.get("FullPhysicalSMIncluded") is True,
          "Full EFT1 heavy-scalar mass/registration audit failed")


def validate_resume_from_benchmark(out: Path, config: Path, profile: Path,
                                   manifest_path: Path) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Fail closed before modifying a prior failed run's manifest.

    This narrow recovery mode is intended only for the known failed Warsaw
    *report check* after a successful Wolfram Warsaw reference calculation.
    It rechecks every saved upstream report and the original EFT1 hashes.
    """
    previous = read_object(manifest_path)
    check(previous.get("status") == "Failed" and
          previous.get("failed_step") == "warsaw_reference" and
          previous.get("mode") == "ReuseVerifiedEFT1" and
          previous.get("model") == "D_-2",
          "Resume requires the preceding --reuse-boundary run stopped at Warsaw validation")
    check(previous.get("benchmark_config") == str(config),
          "Benchmark configuration differs from the failed Warsaw run")
    check(previous.get("eft1_boundary_sha256") == boundary_hashes(out),
          "First-threshold WXF files changed since the failed Warsaw run")

    stage_names = ["scalar_hard_one_loop", "components", "operator_structure",
                   "running_symmetry", "candidate_assembly", "warsaw_reference"]
    steps = previous.get("steps")
    check(isinstance(steps, list) and len(steps) == len(stage_names) and
          [x.get("stage") for x in steps if isinstance(x, dict)] == stage_names,
          "Failed-run stage sequence does not match the known Warsaw-check failure")
    for stage, record in zip(stage_names, steps):
        check(record.get("exit_code") == 0 and record.get("timed_out") is False,
              f"Prior stage {stage} did not finish successfully")
        report = read_object(profile / STAGE_REPORTS[stage])
        validate_report(stage, report)
        if stage != "warsaw_reference":
            check(record.get("reported_status") == EXPECTED_STATUSES[stage],
                  f"Prior stage {stage} lacks a validated completion record")
    validate_mass_audit(profile)
    for artifact in ("profile_full_scalar_one_loop_processed.wxf",
                     "profile_full_scalar_one_loop_hard_hbar_sector.wxf",
                     "profile_c5_flavour_candidate_total_provisional.wxf"):
        check((profile / artifact).is_file(), f"Missing prior validated output: {artifact}")

    validated_steps = [{**item, "reused_from_prior_attempt": True,
                        "reported_status": EXPECTED_STATUSES[item["stage"]]}
                       for item in steps]
    return previous, validated_steps


def kill_tree(proc: subprocess.Popen[bytes]) -> None:
    if proc.poll() is not None:
        return
    if os.name == "nt":
        try:
            subprocess.run(["taskkill", "/PID", str(proc.pid), "/T", "/F"],
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                           check=False, timeout=20)
        except (OSError, subprocess.TimeoutExpired):
            proc.kill()
    else:
        try:
            os.killpg(proc.pid, signal.SIGTERM)
        except ProcessLookupError:
            pass
    try:
        proc.wait(timeout=8)
    except subprocess.TimeoutExpired:
        if os.name != "nt":
            try:
                os.killpg(proc.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
        else:
            proc.kill()
        proc.wait()


def tail(path: Path, n: int = 25) -> str:
    if not path.is_file():
        return "(missing log)"
    with path.open("r", encoding="utf-8", errors="replace") as inp:
        return "".join(deque(inp, maxlen=n)).strip()


def run_command(stage: str, cmd: list[str], *, cap_seconds: int,
                log_dir: Path, on_pulse: Callable[[str], None] | None = None) -> dict[str, Any]:
    # wolframscript may exit 0 even when -file points at a missing source.
    # Validate the actual path before starting a kernel. This also catches
    # accidental variable shadowing in the stage launcher.
    if len(cmd) >= 3 and cmd[1] == "-file":
        wl_path = Path(cmd[2])
        if wl_path.suffix.lower() != ".wl" or not wl_path.is_file():
            raise PipelineFailure(
                f"{stage}: invalid Wolfram -file target: {wl_path}"
            )
    log = log_dir / f"{stage}.log"
    start = time.monotonic()
    print(f"[D(-2) sequential] {stage}: starting (cap {cap_seconds}s)", flush=True)
    flags = getattr(subprocess, "CREATE_NEW_PROCESS_GROUP", 0) if os.name == "nt" else 0
    with log.open("wb") as dst:
        try:
            proc = subprocess.Popen(cmd, cwd=ROOT, env={**os.environ,
                "T3_EXPERIMENTAL_GENERAL_QUADRATICS": "1",
                "T3_EXPERIMENTAL_GENERAL_FERMIONS": "1"},
                stdout=dst, stderr=subprocess.STDOUT,
                creationflags=flags, start_new_session=(os.name != "nt"))
        except OSError as exc:
            raise PipelineFailure(f"Cannot launch {stage}: {exc}") from exc
        timed_out = False
        last_pulse = start
        while proc.poll() is None:
            elapsed = time.monotonic() - start
            if elapsed > cap_seconds:
                timed_out = True
                kill_tree(proc)
                break
            if elapsed - (last_pulse - start) >= 20:
                print(f"[D(-2) sequential] {stage}: active for {int(elapsed)}s", flush=True)
                if on_pulse:
                    on_pulse(stage)
                last_pulse = time.monotonic()
            time.sleep(0.3)
        code = proc.poll()
    result = {"stage": stage, "exit_code": code, "timed_out": timed_out,
              "elapsed_seconds": round(time.monotonic() - start, 2),
              "log": str(log)}
    print(f"[D(-2) sequential] {stage}: exit {code}; {result['elapsed_seconds']}s", flush=True)
    if timed_out or code != 0:
        print(tail(log, 45), flush=True)
        raise PipelineFailure(f"Stage {stage} {'timed out' if timed_out else 'failed'} "
                              f"(exit {code}); see {log}")
    return result


def initial_preflight(out: Path) -> Path:
    # Existing audited preflight checks stage1, real-basis RGE, full flavour,
    # original tensor seeds and the direct-running insertion, but this time
    # verifies the new integrated scalar driver instead of the legacy driver.
    from scripts.RunGeneralDMinus2ScalarThreshold import preflight
    scalar_wl = ROOT / "Lagrangian/RunGeneralDMinus2SequentialScalarThreshold.wl"
    return Path(preflight(out, script=scalar_wl)["insertion_file"])


def validate_sources(config: Path) -> None:
    needed = [
        "scripts/RunGeneralDMinus2Stage1.py",
        "scripts/RunGeneralDMinus2RealEFT1RGE.py",
        "scripts/RunGeneralDMinus2FullFlavorRunning.py",
        "scripts/RunGeneralDMinus2ScalarThreshold.py",
        "Lagrangian/RunGeneralDMinus2Stage1.wl",
        "Lagrangian/RunGeneralDMinus2SequentialScalarThreshold.wl",
        "Lagrangian/RunGeneralDMinus2FixedOrderComponents.wl",
        "Lagrangian/RunGeneralDMinus2C5FlavourStructureAudit.wl",
        "Lagrangian/RunGeneralDMinus2RunningSymmetryAudit.wl",
        "Lagrangian/RunGeneralDMinus2C5FlavourAssemblyCandidate.wl",
        "Lagrangian/RunGeneralDMinus2C5NormalizationReference.wl",
        "Lagrangian/RunGeneralDMinus2C5RoughBenchmark.wl",
        "Lagrangian/RunGeneralDMinus2C5PoleCancellationAudit.wl",
    ]
    missing = [p for p in needed if not (ROOT / p).is_file()]
    check(not missing, "Missing project source files: " + ", ".join(missing))
    check(config.is_file(), f"Missing benchmark configuration: {config}")


def run_pipeline(out: Path, config: Path, *, reuse_boundary: bool,
                 wolframscript: str, stage1_timeout: int,
                 python_stage_timeout: int, scalar_stage_timeout: int,
                 wolfram_timeout: int, startup_allowance: int,
                 resume_from_benchmark: bool = False) -> dict[str, Any]:
    out = out.resolve()
    config = config.resolve()
    validate_sources(config)
    profile = out / REL_PROFILE
    log_dir = out / REL_RUN
    profile.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = log_dir / "sequential_matching_manifest.json"
    if resume_from_benchmark:
        check(reuse_boundary, "--resume-from-benchmark requires --reuse-boundary")
        previous, saved_steps = validate_resume_from_benchmark(
            out, config, profile, manifest_path)
    else:
        previous, saved_steps = None, []
    manifest: dict[str, Any] = {
        "status": "Started", "model": "D_-2",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "mode": ("ResumeValidatedWarsawReference" if resume_from_benchmark else
                 "ReuseVerifiedEFT1" if reuse_boundary else "FreshUVToC5Benchmark"),
        "matching_order": [["NewFermion"], ["NewScalar1", "T3S2R", "T3S2I"]],
        "hard_one_loop_input": "full saved EFT1 tree, with reversible heavy-mass EFT-order promotion",
        "inherited_one_loop_propagation": "pure-SM EFT1 O(hbar), given no linear heavy-scalar source",
        "direct_running_kept_out_of_match": True,
        "physical_c5_normalization_certified": False,
        "full_smeft_operator_basis_certified": False,
        "nf1_simple_pole_cancellation_required": True,
        "benchmark_finite_part_only": True,
        "benchmark_config": str(config),
        "steps": saved_steps,
        "resumed_from_previous_failed_run": resume_from_benchmark,
        "previous_failed_run_utc": previous.get("generated_utc") if previous else None,
        "scalar_matching_executed_this_run": False if resume_from_benchmark else True,
    }
    json_write(manifest_path, manifest)
    python = sys.executable
    step = "initialization"
    boundary_before: dict[str, str] | None = None
    try:
        def execute(name: str, cmd: list[str], seconds: int) -> None:
            nonlocal step
            step = name
            result = run_command(name, cmd, cap_seconds=seconds, log_dir=log_dir)
            manifest["steps"].append(result)
            json_write(manifest_path, manifest)

        if not reuse_boundary:
            execute("fermion_threshold", [python, "-X", "utf8", "-m",
                     "scripts.RunGeneralDMinus2Stage1", "--output", str(out),
                     "--wolframscript", wolframscript], stage1_timeout)
            execute("real_basis_eft1_rge", [python, "-X", "utf8", "-m",
                     "scripts.RunGeneralDMinus2RealEFT1RGE", "--output", str(out)],
                    python_stage_timeout)
            execute("full_flavour_running", [python, "-X", "utf8", "-m",
                     "scripts.RunGeneralDMinus2FullFlavorRunning", "--matching-dir", str(out)],
                    python_stage_timeout)
        step = "preflight"
        insertion = initial_preflight(out)
        boundary_before = boundary_hashes(out)
        manifest["eft1_boundary_sha256"] = boundary_before
        json_write(manifest_path, manifest)

        def wolfram(name: str, script_filename: str, args: list[str],
                    seconds: int, report: str) -> dict[str, Any]:
            target = profile / report
            target.unlink(missing_ok=True)  # prevent stale previous success
            if name == "scalar_hard_one_loop":
                (profile / "profile_full_scalar_one_loop_audit.json").unlink(missing_ok=True)
                for artifact_filename in (
                    "profile_full_scalar_one_loop_raw.wxf",
                    "profile_full_scalar_one_loop_green.wxf",
                    "profile_full_scalar_one_loop_eom.wxf",
                    "profile_full_scalar_one_loop_evaluated.wxf",
                    "profile_full_scalar_one_loop_processed.wxf",
                    "profile_full_scalar_one_loop_hard_hbar_sector.wxf",
                    "profile_full_scalar_one_loop_hard_c5_candidate.txt",
                ):
                    (profile / artifact_filename).unlink(missing_ok=True)
            if name == "candidate_assembly":
                (profile / "profile_c5_flavour_candidate_total_provisional.wxf").unlink(missing_ok=True)
            script_path = ROOT / "Lagrangian" / script_filename
            check(script_path.is_file() and script_path.suffix.lower() == ".wl",
                  f"{name}: missing Wolfram source: {script_path}")
            execute(name, [wolframscript, "-file", str(script_path), *args], seconds)
            data = read_object(target)
            validate_report(name, data)
            if name == "scalar_hard_one_loop":
                validate_mass_audit(profile)
                check((profile / "profile_full_scalar_one_loop_processed.wxf").is_file() and
                      (profile / "profile_full_scalar_one_loop_hard_hbar_sector.wxf").is_file(),
                      "Hard SMEFT WXF outputs missing")
            if name == "candidate_assembly":
                check((profile / "profile_c5_flavour_candidate_total_provisional.wxf").is_file(),
                      "Assembled C5 WXF missing")
            manifest["steps"][-1]["reported_status"] = data["Status"]
            json_write(manifest_path, manifest)
            return data

        wcap = wolfram_timeout + startup_allowance
        if resume_from_benchmark:
            # Earlier stages were already revalidated against the failed
            # manifest and unchanged EFT1 boundary before any files were
            # overwritten. Do not repeat the matching or the Warsaw test.
            print("[D(-2) sequential] Reusing six validated scalar/C5 stages; "
                  "continuing with benchmark and exact pole", flush=True)
        else:
            wolfram("scalar_hard_one_loop", "RunGeneralDMinus2SequentialScalarThreshold.wl",
                    [str(out), str(insertion), str(scalar_stage_timeout)],
                    7 * scalar_stage_timeout + startup_allowance,
                    STAGE_REPORTS["scalar_hard_one_loop"])
            wolfram("components", "RunGeneralDMinus2FixedOrderComponents.wl",
                    [str(out), str(insertion)], wcap, STAGE_REPORTS["components"])
            wolfram("operator_structure", "RunGeneralDMinus2C5FlavourStructureAudit.wl",
                    [str(out)], wcap, STAGE_REPORTS["operator_structure"])
            wolfram("running_symmetry", "RunGeneralDMinus2RunningSymmetryAudit.wl",
                    [str(out), str(insertion)], wcap, STAGE_REPORTS["running_symmetry"])
            wolfram("candidate_assembly", "RunGeneralDMinus2C5FlavourAssemblyCandidate.wl",
                    [str(out), str(insertion)], wcap, STAGE_REPORTS["candidate_assembly"])
            wolfram("warsaw_reference", "RunGeneralDMinus2C5NormalizationReference.wl",
                    [str(out)], wcap, STAGE_REPORTS["warsaw_reference"])
        numerical = wolfram("benchmark", "RunGeneralDMinus2C5RoughBenchmark.wl",
                            [str(out), str(config)], wcap, STAGE_REPORTS["benchmark"])
        pole = wolfram("exact_pole", "RunGeneralDMinus2C5PoleCancellationAudit.wl",
                       [str(out), "1"], max(wcap, 420), STAGE_REPORTS["exact_pole"])

        step = "boundary_integrity"
        check(boundary_hashes(out) == boundary_before,
              "One of the four authoritative first-threshold WXF files changed during scalar processing")
        manifest.update(status="SequentialC5BenchmarkReproduced_ProvisionalNormalization",
                        fermion_threshold_executed_this_run=not reuse_boundary,
                        scalar_one_loop_matched=True,
                        fixed_order_c5_candidate_assembled=True,
                        nf1_simple_pole_cancellation_verified=True,
                        boundary_hashes_unchanged=True,
                        neutrino_masses_eV=numerical["NeutrinoMasses_eV"],
                        mass_squared_differences_eV2=numerical.get("MassSquaredDifferences_eV2"),
                        c5_relative_symmetry_residual=numerical["C5RelativeSymmetryResidual"],
                        exact_pole_residual_terms=pole["ReducedResidueTermCount"],
                        full_smeft_operator_basis_certified=False,
                        physical_c5_normalization_certified=False,
                        match_result_file=str(profile / STAGE_REPORTS["scalar_hard_one_loop"]),
                        candidate_file=str(profile / "profile_c5_flavour_candidate_total_provisional.wxf"),
                        benchmark_report=str(profile / STAGE_REPORTS["benchmark"]))
        json_write(manifest_path, manifest)
        print("[D(-2) sequential] SUCCESS: " + manifest["status"], flush=True)
        print("[D(-2) sequential] m_nu Takagi masses (eV): " +
              str(manifest["neutrino_masses_eV"]), flush=True)
        print("[D(-2) sequential] Manifest: " + str(manifest_path), flush=True)
        return manifest
    except Exception as exc:
        manifest.update(status="Failed", failed_step=step,
                        error=f"{type(exc).__name__}: {exc}")
        json_write(manifest_path, manifest)
        print(f"[D(-2) sequential] FAILED at {step}: {exc}", flush=True)
        print("[D(-2) sequential] Manifest: " + str(manifest_path), flush=True)
        raise


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matching-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--reuse-boundary", action="store_true",
                        help="Reuse verified UV->F EFT1 and RGE inputs; execute entire scalar and C5 stack")
    parser.add_argument("--resume-from-benchmark", action="store_true",
                        help="After the specific failed Warsaw check, revalidate all saved stages and "
                             "run only benchmark and exact pole; requires --reuse-boundary")
    parser.add_argument("--wolframscript", default="wolframscript")
    parser.add_argument("--stage1-timeout-seconds", type=int, default=3600)
    parser.add_argument("--python-stage-timeout-seconds", type=int, default=1800)
    parser.add_argument("--scalar-stage-timeout-seconds", type=int, default=90)
    parser.add_argument("--wolfram-stage-timeout-seconds", type=int, default=360)
    parser.add_argument("--startup-allowance-seconds", type=int, default=180)
    args = parser.parse_args(argv)
    if any(x < 30 for x in (args.stage1_timeout_seconds,
                            args.python_stage_timeout_seconds,
                            args.wolfram_stage_timeout_seconds,
                            args.startup_allowance_seconds)) or not 5 <= args.scalar_stage_timeout_seconds <= 3600:
        parser.error("Invalid stage timeout")
    try:
        run_pipeline(args.matching_dir, args.config, reuse_boundary=args.reuse_boundary,
                     wolframscript=args.wolframscript,
                     stage1_timeout=args.stage1_timeout_seconds,
                     python_stage_timeout=args.python_stage_timeout_seconds,
                     scalar_stage_timeout=args.scalar_stage_timeout_seconds,
                     wolfram_timeout=args.wolfram_stage_timeout_seconds,
                     startup_allowance=args.startup_allowance_seconds,
                     resume_from_benchmark=args.resume_from_benchmark)
        return 0
    except (OSError, PipelineFailure, ValueError) as exc:
        print(f"[D(-2) sequential] ERROR: {exc}", file=sys.stderr, flush=True)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
