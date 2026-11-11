"""Read-only inspection of saved Matchete scalar-threshold result files.

Run from project root: python -m scripts.RunT3WorkerCScalarAudit --model D_-2
"""
from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path
from uuid import uuid4

from Numerical.joint_scalar_matching.worker_c_matchete_audit import (
    audit_matchete_scalar_profile, save_audit,
)

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_D_M2 = Path("output/general_sequential_pilot/D_m2_stage1/data/"
                    "general_d_m2_scalar_threshold/profile_tree")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", default="D_-2", help="Model name recorded in the audit")
    parser.add_argument("--profile-dir", type=Path, default=None,
                        help="Existing Matchete one-loop profile directory (read-only)")
    parser.add_argument("--output-dir", type=Path, default=None,
                        help="Fresh audit directory; never overwrite")
    parser.add_argument("--require-ready", action="store_true",
                        help="Exit 2 if the source is not ready for an additional C5 term")
    args = parser.parse_args()
    if args.profile_dir is None and args.model != "D_-2":
        parser.error("--profile-dir required unless --model D_-2")
    source = args.profile_dir or DEFAULT_D_M2
    source = source if source.is_absolute() else ROOT / source
    report = audit_matchete_scalar_profile(source, model=args.model)
    if args.output_dir is None:
        run_id = datetime.now().strftime("%Y%m%d_%H%M%S") + "_" + uuid4().hex[:8]
        target = ROOT / "output" / "audits" / f"worker_c_scalar_{args.model}_{run_id}"
    else:
        target = args.output_dir if args.output_dir.is_absolute() else ROOT / args.output_dir
    written = save_audit(report, target)
    print(f"Worker C scalar audit: {report['status']}")
    print(f"Result: {written}")
    print(f"One-loop terms: {report['candidate_inventory']['total_weinberg_terms']}")
    print("Independent finite scalar correction safe to add:",
          report["safe_to_add_as_independent_c5_correction"])
    for reason in report["blocking_reasons"]:
        print("  -", reason)
    return 2 if args.require_ready and not report["safe_to_add_as_independent_c5_correction"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
