from __future__ import annotations

"""Re-run threshold 2 after EFT1 running has been constructed.

This is the first orchestration bridge between the existing monolithic
Wolfram matching run and the new Python EFT1-running chain.

The current pipeline still performs threshold 2 too early.  However, the
stage-1 run already leaves the exact fresh-kernel continuation payload on
disk.  This helper reuses that payload and launches RunThresholdStage.wl a
second time with the full-flavor running insertion as its optional 15th
argument and a validation-mode token as its optional 16th argument.

Ordinary result runs use validation_mode=False.  The authoritative [A]/[B]/[C]
physics is still performed, but provenance-only Matchete re-matches are
skipped.  Final/debug runs use validation_mode=True to restore those expensive
cross-checks.

Physics ordering of THIS re-run is therefore

    stage-1 continuation package
        + EFT1 running insertion
        -> fresh threshold-2 kernel
        -> [A] tree match L0
        -> [B] one-loop match L0
        -> [C] tree propagate (L1_threshold + L1_running).

The earlier premature stage-2 result is not used by this helper.

This is intentionally a small bridge before changing the whole pipeline
driver.  Once regression passes, pipeline.py can call this function and use
its result as the authoritative stage-2 output.
"""

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
import subprocess
from typing import Iterable


@dataclass(frozen=True)
class Threshold2Continuation:
    tree: Path
    loop: Path
    transition_tree: Path
    transition_full: Path
    cg_registry: Path


def _all_wxf(output_dir: Path) -> list[Path]:
    return sorted(
        p for p in output_dir.rglob("*.wxf")
        if p.is_file()
    )


def _score(path: Path, required: tuple[str, ...], forbidden: tuple[str, ...]) -> int:
    name = path.name.lower()
    full = str(path).lower()
    if any(word not in full for word in required):
        return -10_000
    if any(word in full for word in forbidden):
        return -10_000

    score = 0
    if "fresh_kernel_stage_2" in full:
        score += 20
    if "stage_2" in full or "stage2" in full:
        score += 5
    score -= len(name) // 20
    return score


def _pick(
    files: Iterable[Path],
    *,
    required: tuple[str, ...],
    forbidden: tuple[str, ...] = (),
    label: str,
) -> Path:
    ranked = sorted(
        ((_score(path, required, forbidden), path) for path in files),
        key=lambda item: item[0],
        reverse=True,
    )
    ranked = [item for item in ranked if item[0] > -10_000]

    if not ranked:
        candidates = "\n".join(f"  {p}" for p in files)
        raise FileNotFoundError(
            f"Could not auto-discover {label}. WXF files found:\n{candidates}"
        )

    best_score = ranked[0][0]
    best = [path for score, path in ranked if score == best_score]

    if len(best) != 1:
        candidates = "\n".join(f"  {p}" for p in best)
        raise RuntimeError(
            f"Ambiguous auto-discovery for {label}:\n{candidates}\n"
            "Pass the path explicitly."
        )

    return best[0]


def discover_continuation(output_dir: Path) -> Threshold2Continuation:
    files = _all_wxf(output_dir)
    if not files:
        raise FileNotFoundError(
            f"No WXF continuation files were found below {output_dir}."
        )

    # RunMatching has used slightly different filenames during development,
    # so discovery is semantic rather than tied to one exact spelling.
    transition_tree = _pick(
        files,
        required=("transition", "tree"),
        forbidden=("result",),
        label="canonical transition tree",
    )
    transition_full = _pick(
        files,
        required=("transition",),
        forbidden=("tree", "result", "cg"),
        label="canonical transition full/O(hbar) EFT",
    )
    cg_registry = _pick(
        files,
        required=("cg",),
        forbidden=("result",),
        label="CG registry",
    )

    # The stage input tree/loop files can have generic names.  Prefer files
    # in the fresh-kernel directory and avoid transition/result/CG payloads.
    non_transition = [
        p for p in files
        if "transition" not in p.name.lower()
        and "result" not in p.name.lower()
        and "cg" not in p.name.lower()
    ]

    tree = _pick(
        non_transition,
        required=("tree",),
        forbidden=("transition",),
        label="stage-2 input tree EFT",
    )
    loop = _pick(
        non_transition,
        required=("loop",),
        forbidden=("transition",),
        label="stage-2 inherited one-loop EFT",
    )

    return Threshold2Continuation(
        tree=tree,
        loop=loop,
        transition_tree=transition_tree,
        transition_full=transition_full,
        cg_registry=cg_registry,
    )


def _alpha_token(alpha: int) -> str:
    if alpha < 0:
        return f"m{abs(alpha)}"
    return f"p{alpha}"


def rerun_threshold2_with_running(
    *,
    output_dir: Path,
    running_insertion: Path,
    run_threshold_script: Path,
    d_s1: int,
    d_s2: int,
    d_f: int,
    alpha: int,
    eft_order: int = 5,
    loop_order: int = 1,
    continuation: Threshold2Continuation | None = None,
    result_path: Path | None = None,
    shared_scalar: bool = False,
    validation_mode: bool = False,
) -> dict:
    output_dir = Path(output_dir).resolve()
    running_insertion = Path(running_insertion).resolve()
    run_threshold_script = Path(run_threshold_script).resolve()

    if not running_insertion.exists():
        raise FileNotFoundError(running_insertion)
    if not run_threshold_script.exists():
        raise FileNotFoundError(run_threshold_script)

    if continuation is None:
        continuation = discover_continuation(output_dir)

    if result_path is None:
        result_path = (
            output_dir
            / "data"
            / "threshold_2_with_eft1_running.wxf"
        )
    result_path = Path(result_path).resolve()
    result_path.parent.mkdir(parents=True, exist_ok=True)

    debug_dir = output_dir / "debug" / "threshold_2_with_eft1_running"
    debug_dir.mkdir(parents=True, exist_ok=True)

    command = [
        "wolframscript",
        "-file",
        str(run_threshold_script),
        str(result_path),
        str(eft_order),
        str(loop_order),
        str(d_s1),
        str(d_s2),
        str(d_f),
        _alpha_token(alpha),
        "S1,S2",
        "S1,S2",
        str(continuation.tree),
        str(continuation.loop),
        str(continuation.transition_tree),
        str(continuation.transition_full),
        str(continuation.cg_registry),
        str(running_insertion),
        "validation" if validation_mode else "results",
    ]

    # Stream Wolfram output live. The old capture_output=True implementation
    # made long Matchete calculations look frozen until the subprocess ended.
    stdout_log = debug_dir / "stdout.log"
    stderr_log = debug_dir / "stderr.log"

    (debug_dir / "command.json").write_text(
        json.dumps(command, indent=2),
        encoding="utf-8",
    )

    output_lines: list[str] = []
    print("[resume] Launching threshold-2 Wolfram kernel...")
    print(f"[resume] Live log: {stdout_log}")

    with stdout_log.open("w", encoding="utf-8") as log_handle:
        process = subprocess.Popen(
            command,
            cwd=run_threshold_script.parent.parent,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )

        assert process.stdout is not None
        try:
            for line in process.stdout:
                clean = line.rstrip("\n")
                output_lines.append(clean)
                print(clean, flush=True)
                log_handle.write(line)
                log_handle.flush()
        except KeyboardInterrupt:
            print(
                "\n[resume] Interrupted; terminating Wolfram process...",
                flush=True,
            )
            process.terminate()
            try:
                process.wait(timeout=10)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait()
            raise

        return_code = process.wait()

    stderr_log.write_text(
        "stderr was merged into stdout.log for live streaming.\n",
        encoding="utf-8",
    )

    combined_output = "\n".join(output_lines)

    loaded_marker = (
        "[fresh kernel] EFT1 running heavy insertion loaded."
        in combined_output
    )
    c_marker = (
        "[fresh kernel] Added EFT1 leading-log heavy insertion to [C] only."
        in combined_output
    )
    direct_weinberg_marker = (
        "[fresh kernel] Direct EFT1 Weinberg running coefficient carried separately to final C5."
        in combined_output
    )
    equal_scale_marker = (
        "[fresh kernel] Direct Weinberg equal-scale check: True"
        in combined_output
    )

    status = (
        "Success"
        if return_code == 0
        and result_path.exists()
        and loaded_marker
        and c_marker
        and direct_weinberg_marker
        and equal_scale_marker
        else "Failed"
    )

    result = {
        "status": status,
        "return_code": return_code,
        "stdout_tail": output_lines[-40:],
        "stderr_tail": [],
        "result_path": str(result_path),
        "running_insertion": str(running_insertion),
        "running_insertion_loaded": loaded_marker,
        "running_inserted_in_C_only": c_marker,
        "direct_weinberg_carried_separately": direct_weinberg_marker,
        "direct_weinberg_equal_scale_vanishes": equal_scale_marker,
        "validation_mode": bool(validation_mode),
        "shared_scalar": bool(shared_scalar),
        "continuation": {
            "tree": str(continuation.tree),
            "loop": str(continuation.loop),
            "transition_tree": str(continuation.transition_tree),
            "transition_full": str(continuation.transition_full),
            "cg_registry": str(continuation.cg_registry),
        },
        "stdout_log": str(stdout_log),
        "stderr_log": str(stderr_log),
    }

    (output_dir / "data" / "threshold_2_running_resume_summary.json").write_text(
        json.dumps(result, indent=2),
        encoding="utf-8",
    )

    return result


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Re-run the S1,S2 threshold using the saved stage-1 continuation "
            "package and a full-flavor EFT1 running insertion."
        )
    )
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("running_insertion", type=Path)
    parser.add_argument("--dims", nargs="+", type=int, required=True, metavar="D")
    parser.add_argument("--alpha", type=int, default=None)
    parser.add_argument(
        "--script",
        type=Path,
        default=Path("Lagrangian") / "RunThresholdStage.wl",
    )
    parser.add_argument("--eft-order", type=int, default=5)
    parser.add_argument("--loop-order", type=int, default=1)
    parser.add_argument(
        "--validation",
        action="store_true",
        help=(
            "enable expensive threshold-2 provenance/consistency re-matches; "
            "ordinary result runs leave this disabled"
        ),
    )

    # Optional explicit continuation paths if auto-discovery is ambiguous.
    parser.add_argument("--tree", type=Path)
    parser.add_argument("--loop", type=Path)
    parser.add_argument("--transition-tree", type=Path)
    parser.add_argument("--transition-full", type=Path)
    parser.add_argument("--cg-registry", type=Path)

    args = parser.parse_args()

    explicit = [
        args.tree,
        args.loop,
        args.transition_tree,
        args.transition_full,
        args.cg_registry,
    ]
    if any(path is not None for path in explicit):
        if not all(path is not None for path in explicit):
            parser.error(
                "If one continuation path is supplied, all five must be supplied."
            )
        continuation = Threshold2Continuation(
            tree=args.tree.resolve(),
            loop=args.loop.resolve(),
            transition_tree=args.transition_tree.resolve(),
            transition_full=args.transition_full.resolve(),
            cg_registry=args.cg_registry.resolve(),
        )
    else:
        continuation = None

    if len(args.dims) == 2:
        d_s1, d_f = args.dims
        d_s2 = d_s1
        alpha = -1 if args.alpha is None else args.alpha
        if alpha != -1:
            parser.error("Two-number shared-scalar mode requires alpha=-1.")
        shared_scalar = True
    elif len(args.dims) == 3:
        d_s1, d_s2, d_f = args.dims
        alpha = 0 if args.alpha is None else args.alpha
        shared_scalar = False
    else:
        parser.error("--dims requires DS DF or DS1 DS2 DF.")

    result = rerun_threshold2_with_running(
        output_dir=args.output_dir,
        running_insertion=args.running_insertion,
        run_threshold_script=args.script,
        d_s1=d_s1, d_s2=d_s2, d_f=d_f, alpha=alpha,
        eft_order=args.eft_order, loop_order=args.loop_order,
        continuation=continuation, shared_scalar=shared_scalar,
        validation_mode=args.validation,
    )

    print(json.dumps(result, indent=2))
    return 0 if result["status"] == "Success" else 1


if __name__ == "__main__":
    raise SystemExit(main())
