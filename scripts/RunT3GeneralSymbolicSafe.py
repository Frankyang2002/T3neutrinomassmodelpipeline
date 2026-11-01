"""Run the existing canonical general-T3 symbolic workflow with collision recovery.

Only a prematurely created EFT1 export is quarantined. The exporter must still
produce a fresh EFT1 JSON during the actual EFT1 subprocess; no cached beta
functions are substituted. Use as an alternative CLI entry point, without
modifying the canonical symbolic or numerical modules.

Run from the project root:
    python -X utf8 -m scripts.RunT3GeneralSymbolicSafe --model A_-2 \
        --run-id a_m2_general_recovery_v1 --study a_m2_general
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

from scripts import RunT3CompleteGeneralScalarSymbolic as scalar
from scripts import RunT3GeneralSymbolicStudy as canonical


_original_execute = scalar.execute


def _execute_with_early_eft1_recovery(
    name: str, command: list[str], path: Path, log: Path, *, env=None
) -> dict:
    path = Path(path)
    if name == "EFT1" and path.is_file():
        # The normal flow must not create this path before the EFT1 command.
        # Preserve the unexpected file rather than deleting or reusing it.
        destination = None
        for index in range(1, 10000):
            candidate = path.with_name(f"{path.stem}.premature_{index:03d}{path.suffix}")
            if not candidate.exists():
                destination = candidate
                break
        if destination is None:
            raise RuntimeError(f"No free quarantine filename for {path}")
        path.rename(destination)
        print(
            f"EFT1 pre-existing export preserved at {destination}; "
            "requiring a fresh EFT1 subprocess export",
            flush=True,
        )
    result = _original_execute(name, command, path, log, env=env)
    if name in ("UV", "EFT1"):
        # Check that the subprocess really wrote the requested stage/model.
        if not isinstance(result, dict):
            raise ValueError(f"{name}: exporter did not return a JSON object")
        metadata = result.get("metadata")
        if not isinstance(metadata, dict) or metadata.get("GeneralScalarStage") != name:
            raise ValueError(f"{name}: export metadata stage mismatch")
        model = command[-2]
        if metadata.get("GeneralScalarModel") != model:
            raise ValueError(f"{name}: export metadata model mismatch: {model}")
        if not isinstance(result.get("report_betas"), dict) or not result["report_betas"]:
            raise ValueError(f"{name}: no beta functions exported")
    return result


def main(argv: list[str] | None = None) -> int:
    scalar.execute = _execute_with_early_eft1_recovery
    try:
        return canonical.main(argv)
    finally:
        scalar.execute = _original_execute


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
