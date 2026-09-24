from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

from common.RunRecords import RunRecord
from common.Thresholds import (
    default_threshold_plan,
    threshold_plan_for_wolfram,
)
from common.T3Model import (
    T3_CLASSES,
    encode_alpha,
    identify_t3_class,
    shared_scalar_formal_dimensions,
    t3_has_neutral_bsm_component,
    t3_neutral_component_fields,
    valid_shared_scalar_dimensions,
    valid_shared_scalar_topology_dimensions,
    valid_t3_dimensions,
    valid_t3_topology_dimensions,
)


# Repository paths and matching order used by the Wolfram model runner.
PROJECT_ROOT = Path(__file__).resolve().parent.parent
OUTPUT_DIR = PROJECT_ROOT / "output"
RUN_MODEL_SCRIPT = PROJECT_ROOT / "Lagrangian" / "RunModel.wl"

# Weinberg dimension-five matching at one loop.
EFT_ORDER = 5
LOOP_ORDER = 1


def _physicalize_shared_scalar_summary(
    summary: dict,
    *,
    d_s: int,
    d_f: int,
) -> dict:
    """Return a shared-scalar summary expressed in physical field language.

    Wolfram matching always works with the formal T3 roles S1 and S2. For the
    shared-scalar branch those two roles represent one physical scalar S, so
    Python-facing metadata must collapse simultaneous S1/S2 occurrences back
    to S.

    Only summary metadata is rewritten here. Matching expressions and files
    produced by Wolfram are untouched.
    """
    physical_summary = dict(summary)
    physical_summary["SharedScalar"] = True
    physical_summary["PhysicalScalarField"] = "S"
    physical_summary["FormalScalarIdentification"] = "S1=C*S*, S2=S"
    physical_summary["PhysicalDimensions"] = {"dS": d_s, "dF": d_f}

    physical_stages: list[dict] = []

    for raw_stage in summary.get("EFTStages", []):
        stage = dict(raw_stage)
        integrated = list(stage.get("IntegratedFields", []))
        active = list(stage.get("ActiveHeavyFields", []))

        if "S1" in integrated and "S2" in integrated:
            integrated = [
                field
                for field in integrated
                if field not in {"S1", "S2"}
            ] + ["S"]

        if "S1" in active and "S2" in active:
            active = [
                field
                for field in active
                if field not in {"S1", "S2"}
            ] + ["S"]

        stage["IntegratedFields"] = integrated
        stage["ActiveHeavyFields"] = active

        for key in ("Label", "label"):
            if key in stage:
                stage[key] = str(stage[key]).replace("S1_S2", "S")

        physical_stages.append(stage)

    if "EFTStages" in summary:
        physical_summary["EFTStages"] = physical_stages

    return physical_summary


def run_model(
    name: str,
    alpha: int,
    d_s1: int,
    d_s2: int,
    d_f: int,
    output_dir: Path,
    model_args: list[str],
    debug_reports: bool = False,
    export_rge_tensors: bool = False,
    threshold_plan: tuple[tuple[str, ...], ...] | None = None,
    shared_scalar: bool = False,
) -> RunRecord:
    """What this does is 
    1. Delete previous output directory and recreate for new results
    2. Run Runmodel.wl with out inputs
    3. Get its output and errors into a file
    4. Get debug reports and summaries
    5. Return a RunRecord object with all the data."""

    # Delete the previous output directory and recreate it.
    # This prevents an old successful result being mistaken for a new result
    # if the current Wolfram run fails before producing its summary.
    shutil.rmtree(output_dir, ignore_errors=True)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Run the model with inputs:
    #   output directory
    #   EFT order
    #   loop order
    #   model arguments
    # We use RunModel.wl+
    threshold_plan = threshold_plan or default_threshold_plan(
        shared_scalar=shared_scalar
    )
    wolfram_threshold_plan = threshold_plan_for_wolfram(
        threshold_plan, shared_scalar=shared_scalar
    )
    threshold_token = "THRESHOLDS=" + ";".join(
        ",".join(group) for group in wolfram_threshold_plan
    )

    command = [
        "wolframscript",
        "-file",
        str(RUN_MODEL_SCRIPT),
        str(output_dir),
        str(EFT_ORDER),
        str(LOOP_ORDER),
        *model_args,
        threshold_token,
        *(["DEBUG"] if debug_reports else []),
        *(["RGETENSORS"] if export_rge_tensors else []),
    ]

    # Launch Wolfram and stream its output live.
    #
    # Sequential matching can take substantially longer than the old single
    # common-threshold Match.  Using subprocess.run(..., capture_output=True)
    # hid every Wolfram message until the entire model finished, which made a
    # long stage-2 Match look like Python had frozen at "Running T3-B".
    #
    # We merge stderr into stdout, print each line immediately, and also retain
    # the complete transcript for the existing debug/error logs.
    process_handle = subprocess.Popen(
        command,
        cwd=PROJECT_ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )

    output_lines: list[str] = []

    assert process_handle.stdout is not None

    try:
        for line in process_handle.stdout:
            output_lines.append(line)
            print(line, end="", flush=True)
    except KeyboardInterrupt:
        # Make Ctrl+C stop wolframscript as well instead of leaving an orphan
        # Mathematica kernel running in the background.
        process_handle.terminate()
        try:
            process_handle.wait(timeout=5)
        except subprocess.TimeoutExpired:
            process_handle.kill()
            process_handle.wait()
        raise
    finally:
        process_handle.stdout.close()

    return_code = process_handle.wait()
    wolfram_output = "".join(output_lines)

    # Small compatibility object for the rest of this function.
    class _ProcessResult:
        pass

    process = _ProcessResult()
    process.returncode = return_code
    process.stdout = wolfram_output
    process.stderr = ""

    # Save Wolfram output for debugging.
    if debug_reports:
        debug_dir = output_dir / "debug"
        debug_dir.mkdir(parents=True, exist_ok=True)
        (debug_dir / "wolfram_stdout.log").write_text(
            process.stdout,
            encoding="utf-8",
        )
        (debug_dir / "wolfram_stderr.log").write_text(
            process.stderr,
            encoding="utf-8",
        )
    elif process.returncode != 0:
        (output_dir / "run_error.log").write_text(
            process.stdout + "\n" + process.stderr,
            encoding="utf-8",
        )

    # Wolfram output has already been streamed live above.  On failure the
    # complete transcript is still saved to run_error.log / debug logs.

    # Read the summary produced by the Wolfram side.
    summary_path = output_dir / "comparison_summary.json"

    summary = (
        json.loads(summary_path.read_text(encoding="utf-8"))
        if summary_path.exists()
        else {
            "BuildStatus": "ProcessFailed",
            "MatchingStatus": "NotRun",
        }
    )

    if shared_scalar:
        summary = _physicalize_shared_scalar_summary(
            summary,
            d_s=d_s1,
            d_f=d_f,
        )

    # Sequential matching must never silently degrade to the historical
    # common-threshold result.  Surface the actual Wolfram stage count here.
    requested_stage_count = len(threshold_plan)
    wolfram_stages = summary.get("EFTStages", [])
    sequential_status = summary.get("SequentialMatchingStatus", "Missing")

    print(
        f"  {name}: sequential matching={sequential_status}; "
        f"EFT stages={len(wolfram_stages)}/{requested_stage_count}",
        flush=True,
    )

    if requested_stage_count > 1 and len(wolfram_stages) != requested_stage_count:
        debug_dir = output_dir / "debug"
        debug_dir.mkdir(parents=True, exist_ok=True)

        (debug_dir / "wolfram_stdout.log").write_text(
            process.stdout,
            encoding="utf-8",
        )
        (debug_dir / "wolfram_stderr.log").write_text(
            process.stderr,
            encoding="utf-8",
        )

        print(
            f"ERROR: {name} requested {requested_stage_count} threshold stages "
            f"but Wolfram exported {len(wolfram_stages)}. "
            f"See {debug_dir / 'wolfram_stdout.log'}",
            flush=True,
        )

        # Keep the summary for inspection, but mark the run as invalid so
        # reports cannot silently omit the requested EFT levels.
        summary["MatchingStatus"] = "SequentialStageExportFailed"

    # Put everything associated with this run into one object.
    return RunRecord(
        name=name,
        alpha=alpha,
        d_s1=d_s1,
        d_s2=d_s2,
        d_f=d_f,
        return_code=process.returncode,
        summary=summary,
        output_dir=output_dir,
        shared_scalar=shared_scalar,
    )


def validate_dimensions(
    d_s1: int,
    d_s2: int,
    d_f: int,
    alpha: int,
    debug_reports: bool = False,
    export_rge_tensors: bool = False,
    threshold_plan: tuple[tuple[str, ...], ...] | None = None,
    output_root: Path | None = None,
    force: bool = False,
) -> RunRecord:
    """All it does is 
    1. Check if dimensions are correct, if not then return error
    2. Identify if its an T3-A..E model and name the output folder after it
    3. Use run_model"""

    # ``--force`` bypasses only production-scope and neutrality restrictions.
    # It never bypasses the representation-theory conditions that define T3.
    topology_ok = valid_t3_topology_dimensions(d_s1, d_s2, d_f)
    if not topology_ok:
        raise ValueError(
            f"({d_s1}, {d_s2}, {d_f}) does not form the required T3 topology. "
            "Dimensions must be positive, each scalar must satisfy dS=dF±1, "
            "and S1⊗S2 must contain the triplet."
        )

    if not force and not valid_t3_dimensions(d_s1, d_s2, d_f):
        raise ValueError(
            f"({d_s1}, {d_s2}, {d_f}) is outside the current production support. "
            "Normal mode supports only SU(2) dimensions 1, 2, and 3. "
            "Use --force to attempt a larger representation that still satisfies "
            "the T3 topology conditions."
        )

    if not force and not t3_has_neutral_bsm_component(d_s1, d_s2, d_f, alpha):
        raise ValueError(
            f"({d_s1}, {d_s2}, {d_f}), alpha={alpha} has no electrically neutral "
            "BSM component. Normal mode requires at least one neutral state. "
            "Use --force to run this charged-only point explicitly."
        )

    # Check whether these dimensions correspond to one of the known
    # T3-A ... T3-E models from the original classification. (Its for output report names)
    output_root = output_root or OUTPUT_DIR

    model_class = identify_t3_class(d_s1, d_s2, d_f)

    if model_class is not None:
        # Known model: keep its familiar A-E name.
        name = f"T3-{model_class}"
        output_dir = output_root / (
            f"T3_{model_class}_alpha_{encode_alpha(alpha)}"
        )
    else:
        # New/generalised representation: identify it directly by dimensions.
        name = f"T3-d{d_s1}-d{d_s2}-F{d_f}"
        output_dir = output_root / (
            f"T3_d{d_s1}_d{d_s2}_F{d_f}_alpha_{encode_alpha(alpha)}"
        )

    # Give arguments for the alpha and the dimensions
    model_args = [
        "DIMS",
        str(d_s1),
        str(d_s2),
        str(d_f),
        encode_alpha(alpha),
    ]

    print(
        f"Running {name}, "
        f"dims=({d_s1}, {d_s2}, {d_f}), "
        f"alpha={alpha} ...",
        flush=True,
    )

    return run_model(
        name,
        alpha,
        d_s1,
        d_s2,
        d_f,
        output_dir,
        model_args,
        debug_reports,
        export_rge_tensors,
        threshold_plan,
    )


def validate_shared_dimensions(
    d_s: int,
    d_f: int,
    *,
    debug_reports: bool = False,
    export_rge_tensors: bool = False,
    threshold_plan: tuple[tuple[str, ...], ...] | None = None,
    output_root: Path | None = None,
    force: bool = False,
) -> RunRecord:
    """Run the one-physical-scalar scotogenic branch."""
    if not valid_shared_scalar_topology_dimensions(d_s, d_f):
        raise ValueError(
            f"({d_s}, {d_f}) does not form the shared-scalar T3 topology. "
            "Dimensions must be positive, dS=dF±1, and S⊗S must contain the triplet."
        )

    if not force and not valid_shared_scalar_dimensions(d_s, d_f):
        raise ValueError(
            f"({d_s}, {d_f}) is outside supported shared-scalar mode. "
            "Current production support is dS=2 with dF=1 or 3. "
            "Use --force to attempt a larger topology-compatible representation."
        )

    d_s1, d_s2, d_f = shared_scalar_formal_dimensions(d_s, d_f, force=force)
    alpha = -1
    output_root = output_root or OUTPUT_DIR
    fermion_label = "N" if d_f == 1 else f"F{d_f}"
    name = f"Scotogenic-dS{d_s}-{fermion_label}"
    output_dir = output_root / f"Scotogenic_dS{d_s}_F{d_f}"
    model_args = ["DIMS", str(d_s1), str(d_s2), str(d_f), encode_alpha(alpha)]
    print(f"Running {name}, shared scalar dims=({d_s}, {d_f}), alpha=-1 ...", flush=True)
    return run_model(
        name, alpha, d_s1, d_s2, d_f, output_dir, model_args,
        debug_reports, export_rge_tensors, threshold_plan, shared_scalar=True,
    )


def obtain_class_dimensions(
    model_class: str,
    alpha: int,
    debug_reports: bool = False,
    export_rge_tensors: bool = False,
    threshold_plan: tuple[tuple[str, ...], ...] | None = None,
    output_root: Path | None = None,
    force: bool = False,
) -> RunRecord:
    """Convert a known A-E class into dimensions and then run normally."""

    if model_class not in T3_CLASSES:
        raise ValueError(f"Unknown T3 model class: {model_class}")

    d_s1, d_s2, d_f = T3_CLASSES[model_class]

    return validate_dimensions(
        d_s1,
        d_s2,
        d_f,
        alpha,
        debug_reports,
        export_rge_tensors,
        threshold_plan,
        output_root,
        force,
    )
