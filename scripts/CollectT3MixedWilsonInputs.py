"""Collect exact local B_-1/C_-1 synchronized Matchete sources for Wilson conversion.

Read-only inspection/collection. Never recomputes matching or asserts Wilson
normalization. Run from the repository using ``python -m scripts...``.

The resulting ZIP preserves repository paths, carries SHA256 provenance, and
contains only the two requested models plus specific source files used to
interpret the numerical mass basis and native Matchete LLSS operators.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import sys
from typing import Any
from zipfile import ZIP_DEFLATED, ZipFile

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SYNC = ROOT / "output/general_observable/mixed_f_synchronised_v1"
DEFAULT_ARCHIVE = ROOT / "mixed_wilson_inputs.zip"
MODELS = ("B_-1", "C_-1")
STATUS = "CurrentRunningMassBasisFThresholdTreeLLSSExported"
SOURCE_FILES = (
    "scripts/RunT3GeneralOrderedWilsonTensor.py",
    "Numerical/running/GeneralT3OrderedTensor.py",
    "tests/test_t3_general_ordered_wilson_tensor.py",
    "scripts/RunT3GeneralMixedFSynchronised.py",
    "scripts/RunT3GeneralWilsonBoundary.py",
    "Numerical/experimental/GeneralT3MixedScalarMassBasis.py",
    "scripts/RunGeneralT3MassEigenMatchete.py",
    "Lagrangian/RunGeneralT3FermionFirstEFT1Tree.wl",
    "Lagrangian/RunMatching.wl",
    "RGE/running/intermediate/ScalarOnlyTensorAdapters.py",
    "RGE/running/intermediate/IntermediateMatcheteParsing.py",
)


def _file_hash(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _resolve(value: str | Path) -> Path:
    p = Path(value)
    return p.resolve() if p.is_absolute() else (ROOT / p).resolve()


def _json(path: Path) -> dict[str, Any]:
    obj = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(obj, dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return obj


def _arcname(path: Path) -> str:
    try:
        return "project/" + path.relative_to(ROOT).as_posix()
    except ValueError:
        # External source paths must not leak the user's full filesystem path.
        short_hash = hashlib.sha256(str(path).encode("utf-8")).hexdigest()[:12]
        return f"external/{short_hash}_{path.name}"


def build_archive(sync_root: Path, archive: Path, *, max_file_mb: int = 100) -> dict:
    """Collect present files; failures are surfaced explicitly in a manifest."""
    if max_file_mb < 1:
        raise ValueError("--max-file-mb must be >= 1")
    archive = archive.resolve()
    sync_root = sync_root.resolve()
    limit = max_file_mb * 1024 * 1024
    report: dict[str, Any] = {
        "format": "t3-mixed-wilson-local-source-collection-v1",
        "scope": list(MODELS),
        "physical_wilson_tensor_computed": False,
        "physical_eft1_rge_computed": False,
        "source_summary": {}, "files": [], "missing_required": [],
        "skipped_optional": [], "too_large": [], "errors": [],
    }
    chosen: dict[str, tuple[Path, str, bool]] = {}

    def add(path: Path, role: str, required: bool) -> None:
        path = path.resolve()
        key = _arcname(path)
        if path == archive:
            return
        if not path.is_file():
            (report["missing_required"] if required else report["skipped_optional"]).append(
                {"role": role, "path": key}
            )
            return
        if path.stat().st_size > limit:
            report["too_large"].append({"role": role, "path": key, "bytes": path.stat().st_size})
            if required:
                report["missing_required"].append({"role": role, "path": key, "reason": "too_large"})
            return
        if key in chosen:
            old_path, old_role, old_required = chosen[key]
            chosen[key] = (old_path, old_role, old_required or required)
        else:
            chosen[key] = (path, role, required)

    for model in MODELS:
        tag = model.replace("-", "m").replace("+", "p")
        status_path = sync_root / tag / "status.json"
        add(status_path, f"{model}: synchronized matching status", True)
        if not status_path.is_file():
            continue
        try:
            status = _json(status_path)
            if status.get("model") != model or status.get("status") != STATUS:
                report["errors"].append(f"{model}: unexpected status/model in {status_path}")
            if status.get("source_operator_matching_at_running_mu_F") is not True:
                report["errors"].append(f"{model}: current mu_F source not certified")
            if status.get("historical_mu_S_basis_reused") is not False:
                report["errors"].append(f"{model}: historical mu_S basis may have been reused")
            report["source_summary"][model] = {
                "status": status.get("status"),
                "mu_F_gev": status.get("mu_F_gev"),
                "llss_term_count": status.get("llss_term_count"),
                "scalar_eigenmasses_gev": status.get("scalar_eigenmasses_gev"),
            }
            for key in ("matched_llss_file", "mass_basis_file", "running_file"):
                value = status.get(key)
                if not isinstance(value, str) or not value.strip():
                    report["missing_required"].append({"role": f"{model}: {key}", "reason": "not_in_status"})
                else:
                    add(_resolve(value), f"{model}: {key}", True)
            # Include the untruncated matched EFT / handoff when available, not
            # only the filtered LLSS text: it may contain phase conventions.
            for sibling in status_path.parent.glob("*.json"):
                add(sibling, f"{model}: sibling certificate", False)
            for sibling in status_path.parent.glob("*.inputform.txt"):
                add(sibling, f"{model}: sibling Matchete InputForm", False)
            for sibling in status_path.parent.glob("*.wl"):
                add(sibling, f"{model}: sibling Wolfram source", False)
        except (OSError, ValueError, KeyError) as exc:
            report["errors"].append(f"{model}: status parsing failed: {type(exc).__name__}: {exc}")

    for relative in SOURCE_FILES:
        add(ROOT / relative, "source: " + relative, relative in SOURCE_FILES[:2])

    archive.parent.mkdir(parents=True, exist_ok=True)
    with ZipFile(archive, "w", compression=ZIP_DEFLATED, compresslevel=6) as zf:
        for key, (path, role, required) in sorted(chosen.items()):
            try:
                digest = _file_hash(path)
                zf.write(path, arcname=key)
                report["files"].append({"archive_path": key, "role": role,
                                        "bytes": path.stat().st_size, "sha256": digest,
                                        "required": required})
            except OSError as exc:
                report["errors"].append(f"Unable to pack {key}: {exc}")
        report["ready_for_source_review"] = (
            not report["missing_required"] and not report["errors"]
        )
        zf.writestr("manifest.json", json.dumps(report, indent=2, allow_nan=False) + "\n")
    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--synchronised-root", type=Path, default=DEFAULT_SYNC)
    parser.add_argument("--output", type=Path, default=DEFAULT_ARCHIVE)
    parser.add_argument("--max-file-mb", type=int, default=100)
    args = parser.parse_args(argv)
    report = build_archive(_resolve(args.synchronised_root), _resolve(args.output),
                           max_file_mb=args.max_file_mb)
    print(f"Source ZIP: {_resolve(args.output)}")
    print(f"Included files: {len(report['files'])}")
    for model in MODELS:
        s = report["source_summary"].get(model, {})
        print(f"{model}: {s.get('status', 'missing')}; LLSS terms: {s.get('llss_term_count', '?')}")
    print(f"Missing required: {len(report['missing_required'])}; errors: {len(report['errors'])}")
    if report["missing_required"] or report["errors"]:
        print("See manifest.json in the ZIP for the missing sources or invalid status.")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
