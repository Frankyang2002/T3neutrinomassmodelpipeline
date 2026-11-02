"""Run certified/synthetic master-Weyl EFT1 C5 source mixing for 1 or 16 T3 models.

Reads inputs without modifying them. All generated output is under an
explicit --output-root and uses exclusive file creation (never overwrites).
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from RGE.running.intermediate.GeneralEFT1C5Mixing import MODELS, solve_boundary


def _load(path: Path) -> tuple[dict, str]:
    data = path.read_bytes()
    return json.loads(data.decode("utf-8-sig")), sha256(data).hexdigest()


def _write_exclusive(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8") as f:
        json.dump(value, f, indent=2, allow_nan=False)
        f.write("\n")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--boundary-root", type=Path, required=True,
                        help="Directory of <model_tag>/boundary.json files")
    parser.add_argument("--output-root", type=Path,
                        default=Path("output/general_observable/lo_eft1_mixing"))
    parser.add_argument("--model", choices=MODELS,
                        help="Otherwise process all 16 models")
    parser.add_argument("--allow-synthetic", action="store_true")
    args = parser.parse_args(argv)
    models = (args.model,) if args.model else MODELS
    results = []
    # Preflight EVERYTHING before writing any file. Existing destinations block.
    for model in models:
        tag = model.replace("+", "p").replace("-", "m")
        boundary = args.boundary_root / tag / "boundary.json"
        dest = args.output_root / tag / "c5_eft1_log.json"
        if dest.exists():
            raise FileExistsError(f"Refusing to overwrite existing output: {dest}")
        raw, fingerprint = _load(boundary)
        if raw.get("model") != model:
            raise ValueError(f"Model tag/boundary mismatch: {model}")
        result = solve_boundary(raw, allow_synthetic=args.allow_synthetic)
        result["source_boundary_sha256"] = fingerprint
        result["source_boundary_path"] = str(boundary)
        results.append((model, dest, result))
    summary_path = args.output_root / "summary.json"
    if summary_path.exists():
        raise FileExistsError(f"Refusing to overwrite existing summary: {summary_path}")
    for _, dest, result in results:
        _write_exclusive(dest, result)
    _write_exclusive(summary_path, {
        "format": "t3-eft1-c5-source-batch-summary-v1",
        "count": len(results), "all16": len(results) == 16,
        "models": [{"model": model, "status": result["status"],
                    "output": str(path), "physical": result["physics_completeness"]["physical_eft1_source_verified"]}
                   for model, path, result in results],
        "finite_joint_scalar_matching_included": False,
        "all16_physical_predictions": all(r["physics_completeness"]["physical_eft1_source_verified"]
                                      for _, _, r in results),
    })
    print(f"Wrote {len(results)} model results: {args.output_root}")
    print("Physical sources:", sum(r["physics_completeness"]["physical_eft1_source_verified"]
                                    for _, _, r in results))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
