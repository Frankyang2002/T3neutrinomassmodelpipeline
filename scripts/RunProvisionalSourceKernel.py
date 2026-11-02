"""Read-only all-16 bridge: Chat-1 provisional LLSS -> symbolic quartic C5 pivot.

No existing certified outputs are written; no physical C5 is claimed.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from RGE.running.intermediate.ProvisionalSourceKernel import (
    MODELS, SourceKernelError, file_digest, kernel_for_boundary, load_json, tag,
)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release-root", required=True, type=Path,
                        help="Root of Chat 1's immutable source_<hash> release")
    parser.add_argument("--output-root", required=True, type=Path,
                        help="New subdirectory in output/general_observable/lo_eft1_mixing")
    parser.add_argument("--model", choices=MODELS, default=None,
                        help="One model only; omit for all 16")
    args = parser.parse_args(argv)
    models = (args.model,) if args.model else MODELS
    release = args.release_root.resolve()
    output = args.output_root.resolve()
    owned_output_root = (ROOT / "output/general_observable/lo_eft1_mixing").resolve()
    if not output.is_relative_to(owned_output_root) or output == owned_output_root:
        raise SourceKernelError(
            f"Output must be a dedicated subdirectory of {owned_output_root}"
        )
    # Check that no result is already present before writing ANY result.
    for model in models:
        target = output / tag(model) / "quartic_pivot_kernel.json"
        if target.exists():
            raise FileExistsError(f"Refusing to overwrite existing source-kernel result: {target}")
    summary_target = output / "summary.json"
    if summary_target.exists():
        raise FileExistsError(f"Refusing to overwrite existing summary: {summary_target}")
    source_summary = release / "summary.json"
    if not source_summary.is_file():
        raise SourceKernelError("Chat 1 release summary.json is required")
    summary = load_json(source_summary)
    if args.model is None:
        if (summary.get("n_models_verified") != 16 or
                {x.get("model") for x in summary.get("models", [])} != set(MODELS)):
            raise SourceKernelError("Chat 1 summary does not verify exactly 16 models")
    if summary.get("absolute_spinor_normalization_certified") is not False or \
       summary.get("physical_eft1_workers_must_reject_this_as_certified_input") is not True:
        raise SourceKernelError("Source normalization gate changed; inspect before transport")
    inventory = {row.get("model"): row for row in summary.get("models", [])}
    results = []
    for model in models:
        row = inventory.get(model)
        if not isinstance(row, dict) or row.get("status") != "NumericalLLSSBySourceVerified_SpinorNormalizationOpen":
            raise SourceKernelError(f"Chat 1 release does not certify source structure for {model}")
        source_path = release / tag(model) / "boundary.json"
        if row.get("path") != f"{tag(model)}/boundary.json" or not source_path.is_file():
            raise SourceKernelError(f"Missing/mismatched Chat1 source file for {model}")
        raw = load_json(source_path)
        payload = kernel_for_boundary(raw, source_sha256=file_digest(source_path))
        if payload["source_count"] != row.get("PL_sources"):
            raise SourceKernelError(f"{model}: Chat1 source inventory count mismatch")
        results.append(payload)
    # Only enter the write phase after all input validation and numerical checks passed.
    for result in results:
        destination = output / tag(result["model"]) / "quartic_pivot_kernel.json"
        destination.parent.mkdir(parents=True, exist_ok=True)
        with destination.open("x", encoding="utf-8") as handle:
            json.dump(result, handle, indent=2, allow_nan=False)
            handle.write("\n")
    summary_target.parent.mkdir(parents=True, exist_ok=True)
    with summary_target.open("x", encoding="utf-8") as handle:
        json.dump({"format": "t3-eft1-conditional-quartic-kernel-batch-v1",
                   "status": "SourceResolvedQuarticKernelsOnly_NotPhysicalC5",
                   "source_release_summary_sha256": file_digest(source_summary),
                   "model_count": len(results),
                   "model_rows": [{"model": r["model"], "sources": r["source_count"],
                                   "nonzero_scalar_pairs": r["distinct_nonzero_quartic_pairs"],
                                   "source_sha256": r["source_boundary_sha256"]}
                                  for r in results],
                   "gamma_cc_spinor_normalization_certified": False,
                   "physical_C5_computed": False,
                   "joint_scalar_matching_included": False}, handle, indent=2)
        handle.write("\n")
    print(f"Extracted {len(results)} provisional source-resolved quartic kernels: {output}")
    print("Physical C5 results: 0 (normalization and quartic dictionary still uncertified)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
