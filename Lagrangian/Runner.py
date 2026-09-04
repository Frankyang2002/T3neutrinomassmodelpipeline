from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

from common.Paths import (
    EFT_ORDER,
    LOOP_ORDER,
    OUTPUT_DIR,
    PROJECT_ROOT,
    RUN_MODEL_SCRIPT,
)
from common.Records import RunRecord
from common.T3Model import (
    T3_CLASSES,
    encode_alpha,
    identify_t3_class,
    valid_t3_dimensions,
)


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
    command = [
        "wolframscript",
        "-file",
        str(RUN_MODEL_SCRIPT),
        str(output_dir),
        str(EFT_ORDER),
        str(LOOP_ORDER),
        *model_args,
        *(["DEBUG"] if debug_reports else []),
        *(["RGETENSORS"] if export_rge_tensors else []),
    ]

    # Launch Wolfram and capture both normal output and errors as text.
    process = subprocess.run(
        command,
        cwd=PROJECT_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )

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

    # If Wolfram fails, show only the final part of its output so the
    # terminal remains readable while still giving useful debug information.
    if process.returncode != 0:
        if process.stdout:
            print("\n".join(process.stdout.splitlines()[-40:]))

        if process.stderr:
            print(
                "\n".join(process.stderr.splitlines()[-20:]),
                file=sys.stderr,
            )

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
    )

def run_dimensions(
    d_s1: int,
    d_s2: int,
    d_f: int,
    alpha: int,
    debug_reports: bool = False,
    export_rge_tensors: bool = False,
) -> RunRecord:
    """All it does is 
    1. Check if dimensions are correct, if not then return error
    2. Identify if its an T3-A..E model and name the output folder after it
    3. Use run_model"""

    # First check whether the dimensions can form the required T3
    # Yukawa and scalar interactions.
    if not valid_t3_dimensions(d_s1, d_s2, d_f):
        raise ValueError(
            f"({d_s1}, {d_s2}, {d_f}) is not supported by the current T3 pipeline. "
            "For now only SU(2) dimensions 1, 2, and 3 are supported, and the "
            "assignment must satisfy dS=dF±1 with S1⊗S2 containing the triplet."
        )

    # Check whether these dimensions correspond to one of the known
    # T3-A ... T3-E models from the original classification. (Its for output report names)
    model_class = identify_t3_class(d_s1, d_s2, d_f)

    if model_class is not None:
        # Known model: keep its familiar A-E name.
        name = f"T3-{model_class}"
        output_dir = OUTPUT_DIR / (
            f"T3_{model_class}_alpha_{encode_alpha(alpha)}"
        )
    else:
        # New/generalised representation: identify it directly by dimensions.
        name = f"T3-d{d_s1}-d{d_s2}-F{d_f}"
        output_dir = OUTPUT_DIR / (
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
    )

def run_known_class(
    model_class: str,
    alpha: int,
    debug_reports: bool = False,
    export_rge_tensors: bool = False,
) -> RunRecord:
    """Convert a known A-E benchmark into dimensions and run normally."""

    if model_class not in T3_CLASSES:
        raise ValueError(f"Unknown T3 model class: {model_class}")

    d_s1, d_s2, d_f = T3_CLASSES[model_class]

    return run_dimensions(
        d_s1,
        d_s2,
        d_f,
        alpha,
        debug_reports,
        export_rge_tensors,
    )

