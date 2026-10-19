"""Temporary, reversible reorganisation of Numerical/orchestration.

Run from the repository root:
  python OrganizeT3Orchestration.py --plan
  python OrganizeT3Orchestration.py --apply
  python OrganizeT3Orchestration.py --restore

Original import paths remain as forwarding compatibility modules. No permanent
file deletion is performed. Archived modules remain importable via launchers.
"""
from __future__ import annotations

import argparse
import ast
import hashlib
import json
from pathlib import Path
import shutil
import sys
from datetime import datetime, timezone

CORE = {
    "ModelNames.py", "NumericalConfig.py", "PipelineNumericalResults.py",
    "RunYukawaGrouping.py", "__init__.py",
}
GROUPS = {
    "studies/common_inputs": (
        "OptimizeCommonYukawas.py", "RunCommonInputCouplingScan.py",
        "RunCommonInputEnvelopes.py", "RunCommonInputFermionScale.py",
        "RunCommonInputRGE.py", "RunCommonParameterComparisonV2.py",
        "RunFourBenchmarkCorrelations.py", "RunScalarCouplingCorrelations.py",
    ),
    "studies/scalar_rge": (
        "RunCommonScalarBetaDominance.py", "RunCommonScalarComparison.py",
        "RunCommonScalarRGControls.py",
    ),
    "studies/fitted_rge": (
        "RunFittedRGEDisplay.py", "RunFittedSensitivityComparison.py",
        "RunFittedStagedRGComparison.py", "RunFittedWeinbergRunningComparison.py",
    ),
    "studies/fermion_scales": (
        "RunFittedFermionScaleScan.py", "RunHeavyFermionScaleComparison.py",
    ),
    "presentation": (
        "BuildThreePanelInteractive.py", "CollectChapter4Interactive.py",
        "ExportChapter4Figures.py",
    ),
    "support": (
        "ResumePipelineNumericalResults.py",  # heavily imported: NOT obsolete as a library
    ),
    "archive_candidates": (
        "BuildComparisonReports.py", "CollectThesisFigures.py",
        "ModelReportLayout.py", "RunCommonParameterComparison.py",
    ),
}
MANIFEST = ".orchestration_organizer_manifest.json"
STUB_MARKER = "# AUTO-GENERATED T3 ORCHESTRATION FORWARDER"


def project_root() -> Path:
    here = Path.cwd().resolve()
    if (here / "Numerical" / "orchestration" / "PipelineNumericalResults.py").is_file():
        return here
    raise SystemExit("Run this script from the project root containing Numerical/orchestration/.")


def checksum(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def relative_parent_fix(text: str, depth: int) -> str:
    """Retain semantics of __file__-based project-root references after relocation."""
    import re
    # Only static Path(__file__).resolve().parents[N] expressions are remapped.
    return re.sub(
        r"(Path\(\s*__file__\s*\)\.resolve\(\)\.parents\[)(\d+)(\])",
        lambda m: m.group(1) + str(int(m.group(2)) + depth) + m.group(3),
        text,
    )


def source_for_moving(path: Path, depth: int) -> str:
    text = path.read_text(encoding="utf-8-sig")
    ast.parse(text, filename=str(path))
    return relative_parent_fix(text, depth)


def forwarding_module(target_module: str) -> str:
    return (f"{STUB_MARKER}\n"
            '"""Backward-compatible import and `python -m` entry point."""\n'
            'from importlib import import_module as _import_module\n'
            f'_target = _import_module("{target_module}")\n'
            'for _name, _value in vars(_target).items():\n'
            '    if _name not in {"__name__", "__package__", "__file__", "__spec__", "__loader__", "__cached__", "__builtins__"}:\n'
            '        globals()[_name] = _value\n'
            'if __name__ == "__main__":\n'
            '    if not hasattr(_target, "main"):\n'
            '        raise SystemExit("This utility has no command-line main(); import it instead.")\n'
            '    raise SystemExit(_target.main())\n')


def mapping():
    return {filename: group for group, files in GROUPS.items() for filename in files}


def plan(root: Path) -> list[dict]:
    base = root / "Numerical" / "orchestration"
    names = mapping()
    inventory = {p.name for p in base.glob("*.py")}
    expected = CORE | set(names)
    missing = sorted(expected - inventory)
    additional = sorted(inventory - expected)
    if missing or additional:
        raise SystemExit(
            "Source inventory differs from the reviewed 30-file snapshot. "
            f"Missing: {missing}; new/unclassified: {additional}. No changes made."
        )
    result = []
    for name, group in names.items():
        source = base / name
        destination = base / group / name
        if destination.exists():
            raise SystemExit(f"Destination exists; refusing to overwrite: {destination}")
        new_source = source_for_moving(source, len(Path(group).parts))
        ast.parse(new_source, filename=str(destination))
        result.append({"name": name, "group": group,
                       "old": str(source.relative_to(root)),
                       "new": str(destination.relative_to(root)),
                       "original_sha256": checksum(source.read_bytes())})
    return result


def apply(root: Path):
    base = root / "Numerical" / "orchestration"
    manifest_path = root / MANIFEST
    if manifest_path.exists():
        raise SystemExit("An organizer manifest already exists. Restore first or inspect it manually.")
    operations = plan(root)
    # Stage original files in a backup before touching any originals.
    backup = root / ".orchestration_organizer_backup"
    if backup.exists():
        raise SystemExit(f"Backup already exists; refusing to overwrite: {backup}")
    backup.mkdir()
    try:
        for op in operations:
            src = root / op["old"]
            shutil.copy2(src, backup / op["name"])
        manifest = {
            "schema": 1, "created_utc": datetime.now(timezone.utc).isoformat(),
            "operations": operations, "core_untouched": sorted(CORE),
            "note": "Old module paths are launchers. archive_candidates are NOT deleted.",
        }
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        for op in operations:
            src = root / op["old"]
            dst = root / op["new"]
            dst.parent.mkdir(parents=True, exist_ok=True)
            # Every new Python package folder gets an __init__.py.
            folder = dst.parent
            while folder != base:
                (folder / "__init__.py").touch(exist_ok=True)
                folder = folder.parent
            dst.write_text(source_for_moving(src, len(Path(op["group"]).parts)), encoding="utf-8")
            module = ".".join(Path(op["new"]).with_suffix("").parts)
            src.write_text(forwarding_module(module), encoding="utf-8")
        print(f"Organized {len(operations)} utilities. Manifest: {manifest_path}")
        print("Core files unchanged. Legacy candidates moved to archive_candidates (NOT deleted).")
        print("Original import/module entry points remain available through forwarding files.")
        print("Run Python import smoke tests and your usual project tests before keeping changes.")
        print("To undo: python OrganizeT3Orchestration.py --restore")
    except Exception:
        print("Error during apply. Attempting automatic restoration.", file=sys.stderr)
        restore(root)
        raise


def restore(root: Path):
    base = root / "Numerical" / "orchestration"
    manifest_path = root / MANIFEST
    backup = root / ".orchestration_organizer_backup"
    if not manifest_path.is_file() or not backup.is_dir():
        raise SystemExit("Missing manifest or backup; cannot safely restore.")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    for op in manifest["operations"]:
        src = root / op["old"]
        if not (backup / op["name"]).is_file():
            raise SystemExit(f"Original backup is missing: {op['name']}")
        if src.is_file() and STUB_MARKER not in src.read_text(encoding="utf-8-sig"):
            raise SystemExit(f"Old location was edited since apply; refusing to overwrite: {src}")
    for op in manifest["operations"]:
        shutil.copy2(backup / op["name"], root / op["old"])
        new = root / op["new"]
        if new.is_file():
            new.unlink()
    # Remove only empty new group directories and our empty init markers.
    for op in manifest["operations"]:
        directory = (root / op["new"]).parent
        while directory != base:
            children = list(directory.iterdir())
            if len(children) == 1 and children[0].name == "__init__.py" and not children[0].read_text(encoding="utf-8").strip():
                children[0].unlink()
            try:
                directory.rmdir()
            except OSError:
                break
            directory = directory.parent
    shutil.rmtree(backup)
    manifest_path.unlink()
    print("Restored original orchestration files. Original code is unchanged.")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    options = ap.add_mutually_exclusive_group(required=True)
    options.add_argument("--plan", action="store_true", help="Read-only inventory and move plan")
    options.add_argument("--apply", action="store_true", help="Back up and reorganize")
    options.add_argument("--restore", action="store_true", help="Restore backed-up originals")
    args = ap.parse_args()
    root = project_root()
    if args.restore:
        restore(root)
    elif args.apply:
        apply(root)
    else:
        ops = plan(root)
        for op in ops:
            print(f"{op['old']}  ->  {op['new']}")
        print(f"\n{len(ops)} moves planned; {len(CORE)} core files untouched. No changes made.")


if __name__ == "__main__":
    main()
