"""Wolfram process execution for T3 Lagrangian construction and matching.

This module owns only the external-process boundary:

- prepare a clean model output directory;
- translate the physical threshold plan to formal Wolfram roles;
- construct and execute the ``wolframscript`` command;
- stream stdout live;
- preserve the existing debug/error log contract.

It does not validate T3 representation choices and does not interpret matching
summary JSON.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import shutil
import subprocess

from common.Thresholds import (
    ThresholdPlan,
    default_threshold_plan,
    threshold_plan_for_wolfram,
)


PROJECT_ROOT = Path(__file__).resolve().parent.parent
OUTPUT_DIR = PROJECT_ROOT / "output"
RUN_MODEL_SCRIPT = PROJECT_ROOT / "Lagrangian" / "RunModel.wl"

EFT_ORDER = 5
LOOP_ORDER = 1


@dataclass(frozen=True, slots=True)
class WolframProcessResult:
    """Captured result of one streamed ``wolframscript`` process."""

    returncode: int
    stdout: str
    stderr: str = ""


def _prepare_output_directory(output_dir: Path) -> None:
    """Remove stale model outputs and create a fresh run directory."""

    shutil.rmtree(output_dir, ignore_errors=True)
    output_dir.mkdir(parents=True, exist_ok=True)


def _wolfram_command(
    output_dir: Path,
    model_args: list[str] | tuple[str, ...],
    *,
    threshold_plan: ThresholdPlan,
    shared_scalar: bool,
    debug_reports: bool,
    export_rge_tensors: bool,
) -> list[str]:
    """Build the existing ``RunModel.wl`` command line."""

    wolfram_threshold_plan = threshold_plan_for_wolfram(
        threshold_plan,
        shared_scalar=shared_scalar,
    )
    threshold_token = "THRESHOLDS=" + ";".join(
        ",".join(group) for group in wolfram_threshold_plan
    )

    return [
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


def _stream_process(command: list[str]) -> WolframProcessResult:
    """Run Wolfram, stream merged output, and make Ctrl+C terminate it."""

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

    return WolframProcessResult(
        returncode=return_code,
        stdout="".join(output_lines),
        stderr="",
    )


def _write_process_logs(
    output_dir: Path,
    process: WolframProcessResult,
    *,
    debug_reports: bool,
) -> None:
    """Preserve the existing Wolfram debug/error logging contract."""

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


def run_wolfram_model(
    output_dir: Path,
    model_args: list[str] | tuple[str, ...],
    *,
    debug_reports: bool = False,
    export_rge_tensors: bool = False,
    threshold_plan: ThresholdPlan | None = None,
    shared_scalar: bool = False,
) -> tuple[WolframProcessResult, ThresholdPlan]:
    """Execute ``RunModel.wl`` and return its process result and resolved plan."""

    _prepare_output_directory(output_dir)

    resolved_plan = threshold_plan or default_threshold_plan(
        shared_scalar=shared_scalar
    )
    command = _wolfram_command(
        output_dir,
        model_args,
        threshold_plan=resolved_plan,
        shared_scalar=shared_scalar,
        debug_reports=debug_reports,
        export_rge_tensors=export_rge_tensors,
    )

    process = _stream_process(command)
    _write_process_logs(
        output_dir,
        process,
        debug_reports=debug_reports,
    )

    return process, resolved_plan
