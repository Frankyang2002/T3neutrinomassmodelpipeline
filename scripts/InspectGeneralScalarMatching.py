"""Inspect baseline vs extended symbolic matching, without numerical claims.

Reads full output text exported by RunModel.wl. In particular checks whether
lambdaT3Extra* appears in extracted Weinberg coefficient or holomorphic terms.
Presence shows symbolic sensitivity, not that the UV coupling has a nonzero
numerical boundary value. Absence in a failed/incomplete extraction proves
nothing; the result explicitly distinguishes missing artifacts.

  python scripts/InspectGeneralScalarMatching.py --root output/general_t3_trial
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[1]
MANIFEST = PROJECT / "docs/general_scalar_trial_operator_manifest.json"

COEFFICIENT_FILES = ("c5_coefficient.txt", "weinberg_holomorphic_terms.txt")
COUPLING = re.compile(r"lambdaT3Extra[A-Za-z0-9]+")


def inspect_model(model_path: Path) -> dict:
    expected_extra = None
    try:
        classes, alpha = model_path.name.split('_', 1)
        model_key = f"{classes}_{'-' if alpha.startswith('m') else '+'}{alpha[1:]}"
        expected_extra = json.loads(MANIFEST.read_text(encoding='utf-8'))['models'][model_key]['extra_quartic_complex_couplings']
    except (KeyError, OSError, ValueError):
        pass
    accepted_extra = None
    extended_log = model_path / 'extended_matching.log'
    if extended_log.is_file():
        log = extended_log.read_text(encoding='utf-8', errors='replace')
        accepted_match = re.search(r'Accepted interactions:\s*\{([^}]*)\}', log)
        if accepted_match:
            accepted_extra = len(re.findall(r'\bT3ExtraQ_[0-9_]+\b', accepted_match.group(1)))
    variants: dict[str, dict] = {}
    for variant in ("baseline", "extended"):
        folder = model_path / variant
        artifacts = {}
        detected: set[str] = set()
        coefficient_text = None
        for filename in COEFFICIENT_FILES:
            path = folder / filename
            if path.is_file():
                text = path.read_text(encoding="utf-8", errors="replace")
                if filename == 'c5_coefficient.txt':
                    coefficient_text = re.sub(r'\s+', '', text)
                matches = sorted(set(COUPLING.findall(text)))
                detected.update(matches)
                artifacts[filename] = {"exists": True, "couplings": matches,
                                       "size_bytes": path.stat().st_size}
            else:
                artifacts[filename] = {"exists": False}
        variants[variant] = {
            "matched_c5_artifact_available": any(
                item["exists"] for item in artifacts.values()),
            "extra_couplings_observed": sorted(detected),
            "artifacts": artifacts,
            "coefficient_text_normalized": coefficient_text,
        }
    baseline = variants["baseline"]
    extended = variants["extended"]
    text_equal = (baseline['coefficient_text_normalized'] == extended['coefficient_text_normalized']
                  if baseline['coefficient_text_normalized'] is not None
                  and extended['coefficient_text_normalized'] is not None else None)
    status = (
        "MissingBaselineArtifacts" if not baseline["matched_c5_artifact_available"]
        else "MissingExtendedArtifacts" if not extended["matched_c5_artifact_available"]
        else "IncompleteExtendedInteractions" if (
            accepted_extra is not None and expected_extra is not None
            and accepted_extra != expected_extra)
        else "ExtendedC5DependsOnNewCouplings" if extended["extra_couplings_observed"]
        else "NoNewCouplingsSeenInAvailableC5Artifacts"
    )
    for variant in variants.values():
        # No need to copy large symbolic expression text into a JSON report.
        variant.pop('coefficient_text_normalized')
    return {
        "model": model_path.name,
        "expected_extra_quartics": expected_extra,
        "accepted_extra_quartics_in_extended_log": accepted_extra,
        "extended_log_operator_count_matches_manifest": (
            accepted_extra == expected_extra if accepted_extra is not None
            and expected_extra is not None else None
        ),
        "baseline_extended_c5_text_equal_ignoring_whitespace": text_equal,
        "variants": variants,
        "status": status,
        "warning": "Textual C5 and interaction-log inspection only. Textual equality does not prove algebraic equality or general RGE/numerical consistency.",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("output/general_t3_trial"))
    parser.add_argument("--output", type=Path)
    parser.add_argument('--require-complete', action='store_true',
                        help='fail if a model has missing C5 artifacts or an incomplete extended interaction log')
    args = parser.parse_args()
    project = Path(__file__).resolve().parents[1]
    root = args.root if args.root.is_absolute() else project / args.root
    rows = [inspect_model(folder) for folder in sorted(root.iterdir())
            if folder.is_dir()] if root.is_dir() else []
    for row in rows:
        extra = row["variants"]["extended"]["extra_couplings_observed"]
        print(f"{row['model']:10} {row['status']:46} new_C5_couplings={len(extra)} "
              f"quartics={row['accepted_extra_quartics_in_extended_log']}/{row['expected_extra_quartics']} "
              f"C5_text_equal={row['baseline_extended_c5_text_equal_ignoring_whitespace']}")
        for item in extra:
            print(f"  {item}")
    if args.output:
        dest = args.output if args.output.is_absolute() else project / args.output
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_text(json.dumps({"results": rows}, indent=2) + "\n", encoding="utf-8")
        print(f"JSON report: {dest}")
    if not rows:
        return 1
    if args.require_complete and any(row['status'].startswith('Missing')
           or row['status'] == 'IncompleteExtendedInteractions' for row in rows):
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
