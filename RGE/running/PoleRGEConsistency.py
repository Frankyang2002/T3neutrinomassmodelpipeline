from __future__ import annotations

"""Strict C5 pole/RGE consistency validation.

The authoritative fixed-one-loop regression is

    direct logarithmic coefficient = 2 * hard UV-pole residue

in the same one-generation normalization.  No flavor-orientation or
representation-dependent fallback is allowed here: a failure should expose an
upstream normalization error rather than be normalized away diagnostically.
"""

import json
from pathlib import Path
from typing import Any


def normalize_pole_rge_consistency(path: Path) -> dict[str, Any]:
    path = Path(path)
    payload = json.loads(path.read_text(encoding="utf-8"))

    if payload.get("Status") != "Success":
        raise ValueError("Pole/RGE diagnostic itself is not successful.")

    validated = bool(
        payload.get("DirectLogEqualsTwicePoleResidueOneGeneration", False)
    )

    payload["ConsistencyConvention"] = "strict_one_generation_direct_vs_hard"
    payload["ConsistencyValidated"] = validated
    payload["ConsistencyReason"] = (
        "Validated: direct running log coefficient = 2 x hard pole residue"
        if validated
        else
        "Failed strict regression: direct running log coefficient != "
        "2 x hard pole residue. No normalization fallback was applied."
    )

    # Remove fields left by the temporary orientation-resolved workaround so
    # old cached JSON cannot make a failed strict comparison look successful.
    for key in (
        "OriginalDirectLogEqualsTwicePoleResidueOneGeneration",
        "DirectFlavorOrientationCount",
        "HardPoleFlavorOrientationCount",
        "FlavorOrientationMultiplicity",
        "OrientationResolvedDirectLogToPoleRatio",
        "OrientationResolvedDirectLogEqualsTwicePoleResidue",
    ):
        payload.pop(key, None)

    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return payload


def main() -> int:
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("diagnostic_json", type=Path)
    args = parser.parse_args()

    result = normalize_pole_rge_consistency(args.diagnostic_json)
    print(json.dumps({
        "ConsistencyValidated": result.get("ConsistencyValidated"),
        "ConsistencyConvention": result.get("ConsistencyConvention"),
        "ConsistencyReason": result.get("ConsistencyReason"),
        "DirectLogToPoleRatioOneGenerationInputForm":
            result.get("DirectLogToPoleRatioOneGenerationInputForm"),
    }, indent=2))
    return 0 if result.get("ConsistencyValidated") else 1


if __name__ == "__main__":
    raise SystemExit(main())
