"""Standalone fitted-Yukawa grouping and optional candidate configuration writer.

Usage: python -m Numerical.orchestration.RunYukawaGrouping [--write-candidates]
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from Numerical.fitting.YukawaSimilarity import (group_collection, load_collection,
                                               make_candidate)

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_MANIFEST = PROJECT_ROOT / "Numerical" / "benchmarks" / "fitted_reference" / "benchmark_manifest.json"
DEFAULT_OUTPUT = PROJECT_ROOT / "Reports" / "output" / "full" / "yukawa_grouping"


def run_grouping(*, manifest: Path = DEFAULT_MANIFEST, output: Path = DEFAULT_OUTPUT,
                 threshold: float = 0.1, write_candidates: bool = False) -> dict:
    collection = load_collection(manifest, PROJECT_ROOT)
    report = group_collection(collection, threshold=threshold)
    output.mkdir(parents=True, exist_ok=True)
    (output / "grouping_summary.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    with (output / "group_assignments.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["group", "size", "reference_model", "member_model", "distance_to_reference", "maximum_pairwise_distance"])
        for number, group in enumerate(report["groups"], start=1):
            for model in group["members"]:
                writer.writerow([number, group["size"], group["reference_model"], model,
                                 group["reference_distances"][model], group["maximum_pairwise_distance"]])
    if write_candidates:
        by_name = {x["model"]: x for x in collection}
        group = report["groups"][0]
        ref = by_name[group["reference_model"]]
        destination = output / "UNVALIDATED_reference_candidates"
        destination.mkdir(parents=True, exist_ok=True)
        for name in group["members"]:
            payload = make_candidate(by_name[name], ref)
            (destination / f"{name}.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--threshold", type=float, default=0.1)
    parser.add_argument("--write-candidates", action="store_true",
                        help="Write unvalidated Yukawa-replacement configs for largest group only")
    args = parser.parse_args()
    report = run_grouping(manifest=args.manifest.resolve(), output=args.output.resolve(),
                          threshold=args.threshold, write_candidates=args.write_candidates)
    for number, group in enumerate(report["groups"], 1):
        print(f"Group {number}: {group['size']} models; reference={group['reference_model']}")
    print("Results:", args.output.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
