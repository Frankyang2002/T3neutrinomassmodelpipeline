"""Read-only Batch 3 comparison of independent UV and conditional EFT outputs."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from Numerical.lo_certification.benchmark.uv_eft_compare import (  # noqa: E402
    ComparisonError, compare_all, emit,
)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source-mode", choices=("github-reference", "local"), default="github-reference")
    p.add_argument("--project-root", type=Path, default=ROOT)
    p.add_argument("--check-only", action="store_true")
    p.add_argument("--run-id", default="source01")
    args = p.parse_args(argv)
    try:
        result = compare_all(args.project_root, args.source_mode)
        if not args.check_only:
            result["output"] = str(emit(args.project_root, result, args.run_id))
        print(json.dumps({"status": result["status"], "source_mode": args.source_mode,
                          "output": result.get("output"),
                          "models": [{"model": x["model"],
                                      "UV_masses_ev": x["uv_neutrino_mass_singular_values_ev_desc"],
                                      "EFT_masses_ev_per_abs_N": x["eft_mass_singular_values_ev_per_abs_N_desc"],
                                      "residual_after_best_complex_scale": x["residual_after_best_complex_scale"]}
                                     for x in result["models"]]}, indent=2))
        return 0
    except (ComparisonError, KeyError, TypeError, ValueError, OSError) as exc:
        print(f"BATCH3 COMPARISON BLOCKED: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
