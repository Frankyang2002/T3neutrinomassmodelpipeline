"""Reversible, conservative removal of legacy packaging material from the project view.

This utility never renames physics code, touches generated numerical results, or
writes to remote Git. Preview and review the plan before using --apply.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import shutil
import subprocess
import sys

FORMAT = "t3-repository-presentation-v1"
CONFIRM = "ARCHIVE_PRESENTATION_ONLY"
PLAN_REL = Path("output/maintenance/repository_presentation_plan.json")
ARCHIVE_PARENT = "T3_external_archives"

# These are *exact* tracked paths, not filename patterns.
CANDIDATES = {
    "CLEANUP_README.md": "obsolete root-level maintenance instructions",
    "CONSOLIDATION_README.md": "obsolete package installation instructions",
    "README_REPLACEMENT.md": "earlier proposed front-page replacement",
    "REGENERATE_README.md": "earlier package-specific regeneration instructions",
    "consolidation_preview.json": "historical presentation preview",
    "generated_archive_preview.json": "historical generated-output preview",
    "t -q --tb=line": "accidental command-output filename",
    "PACKAGE_SHA256.txt": "historical distribution checksum file",
    "SHA256_MANIFEST.json": "historical distribution manifest",
    "SOURCE_MANIFEST.json": "historical distribution manifest",
    "backup_normalization_20261031_215003": "dated source backup snapshot",
    "backup_probe_20261031_215644": "dated source probe snapshot",
    "report_cleanup_archive": "saved TeX compilation artifacts",
    "scripts/maintenance/FinishRepositoryConsolidation.py": "superseded front-page installation utility",
    "docs/WORKSTREAM_A_RGBETA_YUKAWA_BRIDGE.md.backup-20261031_205437": "dated document backup",
}

SOURCE_ROOTS = frozenset({
    "common", "model", "Lagrangian", "RGE", "Numerical", "physics",
    "studies", "validation", "tests", "scripts", "configs", "docs", "tools",
})
TEXT_EXTENSIONS = frozenset({
    ".py", ".wl", ".wls", ".m", ".ps1", ".md", ".tex",
    ".txt", ".json", ".toml", ".ini", ".yaml", ".yml",
    ".cfg", ".sh", ".bat",
})
# The source files defining this very inventory naturally mention its paths.
SCAN_EXCLUSIONS = {
    "scripts/maintenance/OrganizeRepositoryPresentation.py",
    "tests/python/test_repository_presentation.py",
}
# Some older maintenance code lists candidate names or names the *output* of
# a preview. These are not input dependencies, so excluding those specific
# references is safe while retaining all other source-reference protection.
SOFT_REFERENCES = {
    "CLEANUP_README.md": {"scripts/maintenance/ExtremeRepositoryCleanup.py"},
    "CONSOLIDATION_README.md": {"scripts/maintenance/ExtremeRepositoryCleanup.py"},
    "REGENERATE_README.md": {"scripts/maintenance/ExtremeRepositoryCleanup.py"},
    "generated_archive_preview.json": {"scripts/maintenance/ArchiveGeneratedResults.py"},
}


class PresentationError(RuntimeError):
    pass


def _run_git(root: Path, *args: str) -> subprocess.CompletedProcess[bytes]:
    p = subprocess.run(["git", *args], cwd=root, capture_output=True, check=False)
    if p.returncode not in (0, 1):
        raise PresentationError(f"Git command failed: git {' '.join(args)}; {p.stderr.decode(errors='replace')}")
    return p


def repo_head(root: Path) -> str:
    root = root.resolve()
    p = _run_git(root, "rev-parse", "--show-toplevel")
    if p.returncode:
        raise PresentationError("Not a Git checkout")
    checkout = Path(os.fsdecode(p.stdout.strip())).resolve()
    if os.path.normcase(str(checkout)) != os.path.normcase(str(root)):
        raise PresentationError(f"Expected the Git repository root, got {root}")
    p = _run_git(root, "rev-parse", "HEAD")
    if p.returncode:
        raise PresentationError("Cannot resolve Git HEAD")
    return p.stdout.decode("ascii").strip()


def tracked_paths(root: Path) -> set[str]:
    p = _run_git(root, "ls-files", "-z", "--cached")
    if p.returncode:
        raise PresentationError("Cannot list tracked files")
    return {os.fsdecode(x).replace("\\", "/") for x in p.stdout.split(b"\x00") if x}


def _digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _safe_rel(rel: str) -> Path:
    posix = PurePosixPath(rel)
    if posix.is_absolute() or not rel or posix.as_posix() != rel or any(
        part in ("", ".", "..") for part in posix.parts
    ) or ":" in rel or "\\" in rel:
        raise PresentationError(f"Unsafe repository path: {rel!r}")
    return Path(*posix.parts)


def _candidate_files(root: Path, rel: str) -> list[str]:
    target = root / _safe_rel(rel)
    if not target.exists() and not target.is_symlink():
        return []
    if target.is_symlink():
        raise PresentationError("Symlink candidate")
    if target.is_file():
        return [rel]
    if not target.is_dir():
        raise PresentationError("Unsupported candidate type")
    files: list[str] = []
    for p in sorted(target.rglob("*")):
        if p.is_symlink():
            raise PresentationError("Symlink inside candidate")
        if p.is_file():
            files.append(p.relative_to(root).as_posix())
        elif not p.is_dir():
            raise PresentationError("Unsupported entry in candidate")
    if not files:
        raise PresentationError("Empty directory: do not archive blindly")
    return files


def _candidate_owner(rel: str) -> bool:
    return any(rel == p or rel.startswith(p + "/") for p in CANDIDATES)


def _reference_hits(root: Path, tracked: set[str]) -> dict[str, list[str]]:
    """Block candidates referred to in retained source, config or documentation."""
    tokens = {rel: (rel.casefold(), PurePosixPath(rel).name.casefold()) for rel in CANDIDATES}
    hits: dict[str, list[str]] = {rel: [] for rel in CANDIDATES}
    for rel in sorted(tracked):
        if rel in SCAN_EXCLUSIONS or _candidate_owner(rel):
            continue
        relative = PurePosixPath(rel)
        if relative.suffix.lower() not in TEXT_EXTENSIONS:
            continue
        if len(relative.parts) > 1 and relative.parts[0] not in SOURCE_ROOTS:
            continue
        path = root / _safe_rel(rel)
        if not path.is_file() or path.is_symlink():
            continue
        # Very large data files should not be parsed as executable references.
        try:
            if path.stat().st_size > 1_000_000:
                continue
            content = path.read_text(encoding="utf-8-sig", errors="replace").casefold()
        except OSError:
            continue
        for candidate, parts in tokens.items():
            if rel not in SOFT_REFERENCES.get(candidate, set()) and any(token in content for token in parts):
                hits[candidate].append(rel)
    return hits


def _unchanged_in_git(root: Path, rel: str) -> bool:
    # Git's normal clean filters handle Windows CRLF correctly. This checks
    # both staged and unstaged changes against HEAD, without raw-byte guesses.
    p = _run_git(root, "diff", "--quiet", "HEAD", "--", rel)
    return p.returncode == 0


def build_plan(root: Path) -> dict:
    root = root.resolve()
    head = repo_head(root)
    tracked = tracked_paths(root)
    refs = _reference_hits(root, tracked)
    approved = []
    held = []
    absent = []
    for rel, reason in sorted(CANDIDATES.items()):
        try:
            files = _candidate_files(root, rel)
        except (PresentationError, OSError) as exc:
            held.append({"path": rel, "reason": str(exc)})
            continue
        if not files:
            absent.append(rel)
            continue
        if refs[rel]:
            held.append({"path": rel, "reason": "referenced by retained files", "referenced_by": refs[rel][:25]})
            continue
        untracked = [p for p in files if p not in tracked]
        if untracked:
            held.append({"path": rel, "reason": "contains untracked files", "examples": untracked[:10]})
            continue
        dirty = [p for p in files if not _unchanged_in_git(root, p)]
        if dirty:
            held.append({"path": rel, "reason": "tracked file differs from Git HEAD", "examples": dirty[:10]})
            continue
        entries = [
            {"path": p, "sha256": _digest(root / _safe_rel(p)), "bytes": (root / _safe_rel(p)).stat().st_size}
            for p in files
        ]
        approved.append({"path": rel, "reason": reason, "files": entries})
    return {
        "format": FORMAT,
        "git_head": head,
        "approved": approved,
        "held": held,
        "absent": absent,
        "approved_files": sum(len(x["files"]) for x in approved),
        "approved_bytes": sum(f["bytes"] for x in approved for f in x["files"]),
    }


def write_plan(root: Path) -> dict:
    plan = build_plan(root)
    destination = root / PLAN_REL
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(plan, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return plan


def apply_plan(root: Path, *, confirm: str) -> Path:
    root = root.resolve()
    if confirm != CONFIRM:
        raise PresentationError(f"Explicit confirmation required: {CONFIRM}")
    plan_path = root / PLAN_REL
    if not plan_path.is_file():
        raise PresentationError("First run preview and review the generated plan")
    saved = json.loads(plan_path.read_text(encoding="utf-8"))
    live = build_plan(root)
    if saved != live:
        raise PresentationError("The reviewed plan is stale. Preview again before archiving")
    if not live["approved"]:
        raise PresentationError("There are no reviewed eligible files")
    archive_parent = root.parent / ARCHIVE_PARENT
    archive_parent.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    archive = archive_parent / f"presentation_{stamp}"
    archive.mkdir(exist_ok=False)
    items = [f for group in live["approved"] for f in group["files"]]
    # Stage and verify *all* archived copies before deleting any original.
    try:
        for item in items:
            rel = item["path"]
            src = root / _safe_rel(rel)
            dst = archive / _safe_rel(rel)
            if _digest(src) != item["sha256"] or src.stat().st_size != item["bytes"]:
                raise PresentationError(f"Changed since review: {rel}")
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dst)
            if _digest(dst) != item["sha256"]:
                raise PresentationError(f"Archive checksum failed: {rel}")
        manifest = {
            "format": FORMAT,
            "repository_root": str(root),
            "git_head": live["git_head"],
            "files": items,
            "groups": [{"path": g["path"], "reason": g["reason"]} for g in live["approved"]],
        }
        (archive / "manifest.json").write_text(
            json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
        )
        for item in items:
            src = root / _safe_rel(item["path"])
            if _digest(src) != item["sha256"]:
                raise PresentationError(f"Source changed during copy: {item['path']}")
        for item in items:
            (root / _safe_rel(item["path"])).unlink()
        for group in live["approved"]:
            directory = root / _safe_rel(group["path"])
            if directory.is_dir():
                for nested in sorted(directory.rglob("*"), key=lambda p: len(p.parts), reverse=True):
                    if nested.is_dir() and not nested.is_symlink():
                        try:
                            nested.rmdir()
                        except OSError:
                            pass
                try:
                    directory.rmdir()
                except OSError:
                    pass
    except Exception:
        # Copies are left intact for recovery; the sources may also still exist.
        raise
    return archive


def restore(root: Path, archive: Path) -> int:
    root = root.resolve()
    archive = archive.resolve()
    permitted = (root.parent / ARCHIVE_PARENT).resolve()
    if permitted not in archive.parents:
        raise PresentationError("Restore source must be in the sibling external archive directory")
    manifest_path = archive / "manifest.json"
    if not manifest_path.is_file():
        raise PresentationError("Archive manifest not found")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("format") != FORMAT or manifest.get("repository_root") != str(root):
        raise PresentationError("Archive does not belong to this repository")
    items = manifest.get("files", [])
    if not isinstance(items, list):
        raise PresentationError("Malformed manifest")
    for item in items:
        rel = item["path"]
        target = root / _safe_rel(rel)
        archived = archive / _safe_rel(rel)
        if not archived.is_file() or archived.is_symlink() or _digest(archived) != item["sha256"]:
            raise PresentationError(f"Archive checksum or file missing: {rel}")
        if target.exists() or target.is_symlink():
            raise PresentationError(f"Refusing to overwrite existing file: {rel}")
    for item in items:
        rel = item["path"]
        target = root / _safe_rel(rel)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(archive / _safe_rel(rel), target)
        if _digest(target) != item["sha256"]:
            raise PresentationError(f"Restore verification failed: {rel}")
    return len(items)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--apply", action="store_true")
    parser.add_argument("--confirm", default="")
    parser.add_argument("--restore", type=Path)
    args = parser.parse_args(argv)
    try:
        if args.restore is not None:
            count = restore(args.root, args.restore)
            print(f"Restored {count} files. No existing files were overwritten.")
        elif args.apply:
            archive = apply_plan(args.root, confirm=args.confirm)
            print(f"Verified external archive: {archive}")
            print("Only approved, unchanged historical packaging files were removed from the working tree.")
        else:
            plan = write_plan(args.root)
            print(f"Git HEAD: {plan['git_head']}")
            print(f"Approved: {len(plan['approved'])} groups, {plan['approved_files']} files, {plan['approved_bytes']} bytes")
            print(f"Reference/modified/untracked holds: {len(plan['held'])}; absent: {len(plan['absent'])}")
            for group in plan["approved"]:
                print(f"  APPROVED {group['path']}")
            for group in plan["held"]:
                print(f"  HELD {group['path']}: {group['reason']}")
            print(f"Review: {args.root.resolve() / PLAN_REL}")
            print("Preview only. No files moved.")
    except (PresentationError, OSError, ValueError, KeyError, json.JSONDecodeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
