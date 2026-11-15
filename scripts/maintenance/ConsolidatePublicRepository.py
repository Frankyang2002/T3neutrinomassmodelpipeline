#!/usr/bin/env python3
"""Audit and optionally archive root-level research packaging notes.

This is a local-only, external-archive-first utility. It does not rename any
executable, touch numerical data, delete research outputs or contact GitHub.

Usage, from the project root:
    python -X utf8 scripts/maintenance/ConsolidatePublicRepository.py plan
    python -X utf8 scripts/maintenance/ConsolidatePublicRepository.py apply --plan <plan.json> --confirm ARCHIVE_UNREFERENCED_ROOT_NOTES
    python -X utf8 scripts/maintenance/ConsolidatePublicRepository.py restore --manifest <manifest.json>

The explicit selection is deliberately narrow. A file with an incoming literal
reference from a retained tracked source is not automatically archived.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import subprocess
import sys
from typing import Any

FORMAT = "t3-reader-facing-root-cleanup-v1"
CONFIRM = "ARCHIVE_UNREFERENCED_ROOT_NOTES"
CORE_ROOT = frozenset({"README.md", ".gitignore", "pipeline.py", "requirements.txt", "full16_symbolic.py"})
ROOT_SPECIAL = frozenset({"README.txt", "INSTALL_AND_RUN.txt", "T3_COMMAND_REFERENCE.md"})
ROOT_NOTE = re.compile(r"^(?:BATCH\d+[_-]|README[_-]|EFT1[_-]|WORKER[_-]|STREAM\d+[_-]).*\.(?:md|txt)$", re.I)
TEXT_SUFFIXES = frozenset({
    ".py", ".md", ".txt", ".wl", ".wls", ".m", ".json", ".jsonl", ".ps1", ".psm1",
    ".csv", ".tsv", ".rst", ".tex", ".toml", ".yaml", ".yml", ".ini", ".cfg",
    ".sh", ".bat", ".cmd", ".html", ".ipynb", ".xml",
})
MAX_SCAN_BYTES = 32 * 1024 * 1024
SCRIPT_REL = "scripts/maintenance/ConsolidatePublicRepository.py"


def run_git(root: Path, *args: str) -> bytes:
    proc = subprocess.run(["git", "-C", str(root), *args], capture_output=True, check=False)
    if proc.returncode:
        raise RuntimeError("git " + " ".join(args) + ": " + proc.stderr.decode("utf-8", "replace").strip())
    return proc.stdout


def git_root(root: Path) -> Path:
    actual = Path(os.fsdecode(run_git(root, "rev-parse", "--show-toplevel").strip())).resolve()
    if actual != root.resolve():
        raise RuntimeError(f"Expected the repository root {root}, not {actual}")
    return actual


def tracked_paths(root: Path) -> list[str]:
    paths = [os.fsdecode(item).replace("\\", "/")
             for item in run_git(root, "ls-files", "-z", "--cached").split(b"\x00") if item]
    return sorted(paths)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def is_candidate(relative: str) -> bool:
    if "/" in relative or relative in CORE_ROOT:
        return False
    return relative in ROOT_SPECIAL or bool(ROOT_NOTE.fullmatch(relative))


def safe_path(root: Path, relative: str) -> Path:
    p = PurePosixPath(relative)
    if (not relative or p.is_absolute() or p.as_posix() != relative
            or any(segment in ("", ".", "..") for segment in p.parts)
            or ":" in relative or "\\" in relative):
        raise RuntimeError(f"Unsafe repository path: {relative!r}")
    file = root.joinpath(*p.parts)
    if file.is_symlink() or not file.is_file() or not file.resolve().is_relative_to(root.resolve()):
        raise RuntimeError(f"Missing, symlinked or unsafe candidate: {relative}")
    return file


def external_base(root: Path, supplied: str | None) -> Path:
    target = Path(supplied).expanduser().resolve() if supplied else (
        root.parent / "T3_external_archives" / root.name / "reader_facing").resolve()
    if target == root or target.is_relative_to(root):
        raise RuntimeError("The archive must be outside the Git repository")
    return target


def check_candidates_unmodified(root: Path, candidates: list[str]) -> None:
    # A note that someone has edited locally must never be selected for silent removal.
    if not candidates:
        return
    for args in (("diff", "--quiet", "--", *candidates),
                 ("diff", "--cached", "--quiet", "--", *candidates)):
        p = subprocess.run(["git", "-C", str(root), *args], capture_output=True)
        if p.returncode == 1:
            raise RuntimeError("A candidate has local tracked modifications; review them before archiving")
        if p.returncode:
            raise RuntimeError(p.stderr.decode("utf-8", "replace"))


def inventory(root: Path) -> dict[str, Any]:
    paths = tracked_paths(root)
    untracked = [os.fsdecode(item).replace("\\", "/")
                 for item in run_git(root, "ls-files", "--others", "--exclude-standard", "-z").split(b"\x00") if item]
    candidates = [p for p in paths if is_candidate(p)]
    references: dict[str, list[str]] = {p: [] for p in candidates}
    # Filename matching only: avoids treating a similar-sounding experimental source
    # as a dependency. All non-candidate tracked textual paths are scanned, including
    # small saved text outputs, tests, configuration files and documentation.
    keys = {p: p.casefold() for p in candidates}
    skipped: list[str] = []
    errors: list[str] = []
    scanned = 0
    for relative in sorted(set(paths + untracked)):
        if relative in keys or relative == SCRIPT_REL:
            continue
        suffix = Path(relative).suffix.lower()
        if suffix not in TEXT_SUFFIXES:
            continue
        file = root / Path(*PurePosixPath(relative).parts)
        if file.is_symlink() or not file.is_file():
            errors.append(relative)
            continue
        try:
            if file.stat().st_size > MAX_SCAN_BYTES:
                skipped.append(relative)
                continue
            content = file.read_text(encoding="utf-8-sig", errors="replace").casefold()
            scanned += 1
        except OSError:
            errors.append(relative)
            continue
        for candidate, token in keys.items():
            if token in content:
                references[candidate].append(relative)
    safe_coverage = not skipped and not errors
    candidates_out: list[dict[str, Any]] = []
    for relative in candidates:
        file = safe_path(root, relative)
        hits = references[relative]
        candidates_out.append({
            "path": relative, "bytes": file.stat().st_size, "sha256": sha256(file),
            "incoming_references": hits[:50], "incoming_reference_count": len(hits),
            "eligible": bool(safe_coverage and not hits),
        })
    groups = Counter(p.split("/", 1)[0] if "/" in p else "(repository root)" for p in paths)
    # Filename inventory for a second, manual audit. NEVER auto-archive these.
    research_name = re.compile(r"(?:^|[_./-])(?:batch\d*|stage\d*|worker\w*|legacy|experimental|trial|probe)(?:$|[_./-])", re.I)
    review = [p for p in paths if "/" in p and research_name.search(Path(p).name)]
    return {
        "format": FORMAT + "-plan", "repo_root": str(root.resolve()),
        "git_head": run_git(root, "rev-parse", "HEAD").decode("ascii").strip(),
        "tracked_count": len(paths), "untracked_checked_count": len(untracked),
        "tracked_by_top_level": dict(sorted(groups.items())),
        "reference_scan_count": scanned, "reference_scan_skipped": skipped,
        "reference_scan_errors": errors, "reference_scan_complete": safe_coverage,
        "root_note_candidates": candidates_out,
        "non_root_research_name_review_count": len(review),
        "non_root_research_name_review_sample": review[:100],
        "note": ("This plan only auto-selects unreferenced root text notes. "
                 "No modules, tests, data or generated assets are selected for removal."),
    }


def build_plan(root: Path, archive_base: Path) -> Path:
    data = inventory(root)
    archive_base.mkdir(parents=True, exist_ok=True)
    path = archive_base / "root_notes_plan.json"
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    rows = data["root_note_candidates"]
    eligible = sum(r["eligible"] for r in rows)
    blocked = sum(not r["eligible"] for r in rows)
    print(f"Read-only: {len(rows)} root notes reviewed: {eligible} eligible, {blocked} blocked.")
    print(f"Tracked files: {data['tracked_count']}; text files scanned: {data['reference_scan_count']}.")
    if not data["reference_scan_complete"]:
        print("Reference scan incomplete: no automatic archival is eligible.")
    print(f"Review additional research-named files manually: {data['non_root_research_name_review_count']}.")
    print(f"Plan: {path}")
    return path


def apply(root: Path, archive_base: Path, plan_file: Path, confirmation: str) -> Path | None:
    if confirmation != CONFIRM:
        raise RuntimeError(f"Explicit confirmation is required: --confirm {CONFIRM}")
    data = json.loads(plan_file.read_text(encoding="utf-8-sig"))
    if data.get("format") != FORMAT + "-plan" or Path(data.get("repo_root", "")).resolve() != root.resolve():
        raise RuntimeError("Plan has a different format or repository root")
    actual_head = run_git(root, "rev-parse", "HEAD").decode("ascii").strip()
    if actual_head != data.get("git_head"):
        raise RuntimeError("Git HEAD changed after the plan; regenerate the plan")
    if not data.get("reference_scan_complete"):
        raise RuntimeError("Plan was produced with an incomplete reference scan")
    fresh = inventory(root)
    if fresh["root_note_candidates"] != data["root_note_candidates"] or not fresh["reference_scan_complete"]:
        raise RuntimeError("Candidate files or incoming references changed; regenerate the plan")
    rows = [r for r in data["root_note_candidates"] if r["eligible"]]
    if not rows:
        print("No unreferenced root notes found; nothing was removed.")
        return None
    paths = [r["path"] for r in rows]
    check_candidates_unmodified(root, paths)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    archive = archive_base / f"root_notes_{stamp}_{actual_head[:10]}"
    if archive.exists():
        raise RuntimeError(f"Archive already exists: {archive}")
    original_root = archive / "original"
    original_root.mkdir(parents=True)
    manifest = {
        "format": FORMAT + "-archive", "repo_root": str(root), "git_head": actual_head,
        "archived_at_utc": stamp, "files": rows,
        "restoration": "Use restore --manifest with the path to this manifest.json",
    }
    # Two-phase operation: copy and verify every source before any removal.
    for row in rows:
        relative = row["path"]
        source = safe_path(root, relative)
        dest = original_root / relative
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, dest)
        if sha256(dest) != row["sha256"]:
            raise RuntimeError(f"Failed archive verification for {relative}; original intact")
    manifest_file = archive / "manifest.json"
    manifest_file.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    for row in rows:
        relative = row["path"]
        if sha256(safe_path(root, relative)) != row["sha256"]:
            raise RuntimeError(f"Source changed while archiving: {relative}")
    for row in rows:
        safe_path(root, row["path"]).unlink()
    print(f"Archived {len(rows)} reviewed notes outside the repository: {archive}")
    print(f"Restore manifest: {manifest_file}")
    print("Local Git tracks the resulting deletions; no Git commit or remote action was performed.")
    return manifest_file


def restore(root: Path, manifest_file: Path) -> None:
    data = json.loads(manifest_file.read_text(encoding="utf-8-sig"))
    if data.get("format") != FORMAT + "-archive" or Path(data.get("repo_root", "")).resolve() != root.resolve():
        raise RuntimeError("Archive manifest does not belong to this repository")
    entries = data.get("files")
    if not isinstance(entries, list):
        raise RuntimeError("Invalid archive entries")
    # Validate every archive and destination before restoring any file.
    for row in entries:
        relative = row["path"]
        if not is_candidate(relative):
            raise RuntimeError(f"Unexpected restore path: {relative}")
        target = root / relative
        if target.exists() or target.is_symlink():
            raise RuntimeError(f"Refusing to overwrite existing path: {relative}")
        saved = manifest_file.parent / "original" / relative
        if saved.is_symlink() or not saved.is_file() or sha256(saved) != row["sha256"]:
            raise RuntimeError(f"Archive missing or checksum mismatch: {relative}")
    for row in entries:
        relative = row["path"]
        shutil.copy2(manifest_file.parent / "original" / relative, root / relative)
        if sha256(root / relative) != row["sha256"]:
            raise RuntimeError(f"Restore checksum mismatch: {relative}")
    print(f"Restored {len(entries)} notes; no tracked changes were staged or committed.")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("plan", "apply", "restore"), nargs="?", default="plan")
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--archive-base", help="External directory for the plan and archives")
    parser.add_argument("--plan", type=Path, help="Previously reviewed JSON plan (for apply)")
    parser.add_argument("--manifest", type=Path, help="Archive manifest (for restore)")
    parser.add_argument("--confirm", help="Explicit opt-in token (for apply)")
    args = parser.parse_args(argv)
    try:
        root = git_root(args.root.resolve())
        if args.action == "restore":
            if not args.manifest:
                parser.error("restore requires --manifest")
            restore(root, args.manifest.resolve())
        else:
            archive = external_base(root, args.archive_base)
            if args.action == "plan":
                build_plan(root, archive)
            else:
                if not args.plan:
                    parser.error("apply requires --plan")
                apply(root, archive, args.plan.resolve(), args.confirm or "")
        return 0
    except (OSError, ValueError, RuntimeError, json.JSONDecodeError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
