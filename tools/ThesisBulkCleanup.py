"""Conservative, local-only archival of reproducible T3 LaTeX build byproducts.

Only selects tracked Reports/output/ .aux, .out, .toc, .log and .compile.log
when a tracked same-stem .tex file exists. No PDFs, TeX sources, numerical
outputs, source programs, matching, RGEs, fitted data or saved physics JSONs
are selected. SHA-256 backup manifest and full preflight before removals.

From repository root:
  python -X utf8 tools/ThesisBulkCleanup.py preview
  python -X utf8 tools/ThesisBulkCleanup.py apply --plan PATH --confirm ARCHIVE_TEX_BUILD_ARTIFACTS
  python -X utf8 tools/ThesisBulkCleanup.py restore --archive PATH
  python -X utf8 tools/ThesisBulkCleanup.py duplicates
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

FORMAT = "t3-thesis-bulk-tex-build-v1"
CONFIRM = "ARCHIVE_TEX_BUILD_ARTIFACTS"
SUFFIXES = (".compile.log", ".aux", ".out", ".toc", ".log")


def git(root: Path, *args: str, input_text: str | None = None) -> str:
    cp = subprocess.run(["git", "-C", str(root), *args], input=input_text,
                        capture_output=True, text=True, encoding="utf-8",
                        errors="surrogateescape", check=False)
    if cp.returncode:
        raise RuntimeError(f"git {' '.join(args)} failed: {cp.stderr.strip()}")
    return cp.stdout


def head(root: Path) -> str:
    return git(root, "rev-parse", "HEAD").strip()


def ensure_clean(root: Path) -> None:
    if git(root, "status", "--porcelain", "--untracked-files=no").strip():
        raise RuntimeError("Tracked working tree is not clean: commit or restore prior changes first")


def tracked_blobs(root: Path) -> dict[str, tuple[str, str]]:
    """Path -> (Git blob id, file mode); works with spaces in filenames."""
    raw = git(root, "ls-files", "--stage", "-z")
    entries: dict[str, tuple[str, str]] = {}
    for record in raw.split("\0"):
        if not record:
            continue
        metadata, path = record.split("\t", 1)
        mode, blob, stage = metadata.split(" ")
        if stage != "0":
            raise RuntimeError(f"Unmerged index entry: {path}")
        entries[path.replace("\\", "/")] = (blob, mode)
    return entries


def safe_file(root: Path, relative: str) -> Path:
    if not relative or relative.startswith("/") or "\\" in relative:
        raise RuntimeError(f"Unsafe project-relative path: {relative!r}")
    parts = relative.split("/")
    if any(p in ("", ".", "..") for p in parts):
        raise RuntimeError(f"Unsafe path components: {relative!r}")
    path = root.joinpath(*parts)
    if not path.resolve().is_relative_to(root.resolve()):
        raise RuntimeError(f"Path leaves repository: {relative}")
    if path.is_symlink() or not path.is_file():
        raise RuntimeError(f"Not a regular existing file: {relative}")
    return path


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def tex_source(relative: str) -> str | None:
    if not relative.startswith("Reports/output/"):
        return None
    for suffix in SUFFIXES:
        if relative.endswith(suffix):
            return relative[:-len(suffix)] + ".tex"
    return None


def find_candidates(root: Path, tracked: dict[str, tuple[str, str]] | None = None) -> list[dict]:
    tracked = tracked if tracked is not None else tracked_blobs(root)
    rows: list[dict] = []
    for rel, (blob, mode) in sorted(tracked.items()):
        source = tex_source(rel)
        if source is None or source not in tracked or mode != "100644":
            continue
        path = safe_file(root, rel)
        if not root.joinpath(*source.split("/")).is_file():
            continue
        rows.append({"path": rel, "tex_source": source, "sha256": sha256(path),
                     "bytes": path.stat().st_size, "git_blob": blob})
    return rows


def external_root(root: Path, override: Path | None) -> Path:
    return override.resolve() if override else (root.parent / "T3_external_archives" / "thesis_bulk_cleanup").resolve()


def preview(root: Path, external: Path) -> Path:
    ensure_clean(root)
    rows = find_candidates(root)
    report = external / "preview"
    report.mkdir(parents=True, exist_ok=True)
    plan = {"format": FORMAT + "-plan", "repo_root": str(root.resolve()),
            "git_head": head(root), "kind": "tracked_tex_build_byproducts",
            "count": len(rows), "bytes": sum(r["bytes"] for r in rows), "entries": rows}
    target = report / "plan.json"
    target.write_text(json.dumps(plan, indent=2) + "\n", encoding="utf-8")
    with (report / "candidates.csv").open("w", newline="", encoding="utf-8-sig") as f:
        fields = ("path", "tex_source", "bytes", "sha256", "git_blob")
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"Git HEAD: {plan['git_head']}")
    print(f"READ-ONLY: {len(rows)} tracked LaTeX build byproducts, {plan['bytes']:,} bytes")
    print(f"Plan: {target}")
    return target


def _verify_repo_binding(root: Path, plan: dict, kind: str) -> None:
    if plan.get("format") != FORMAT + kind:
        raise RuntimeError("Unexpected cleanup plan/archive format")
    if Path(plan.get("repo_root", "")).resolve() != root.resolve():
        raise RuntimeError("Plan/archive belongs to another local checkout")


def apply(root: Path, external: Path, plan_file: Path, confirm: str) -> Path:
    if confirm != CONFIRM:
        raise RuntimeError(f"Apply requires --confirm {CONFIRM}")
    plan = json.loads(plan_file.read_text(encoding="utf-8"))
    _verify_repo_binding(root, plan, "-plan")
    ensure_clean(root)
    if plan.get("git_head") != head(root):
        raise RuntimeError("Git HEAD differs from preview; rerun preview")
    current = find_candidates(root)
    if current != plan.get("entries"):
        raise RuntimeError("Candidate paths/content differ from preview; rerun preview")
    if not current:
        raise RuntimeError("No matching compiler artifacts to archive")
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    destination = external / ("tex_build_" + stamp)
    destination.mkdir(parents=True, exist_ok=False)
    # Copy and checksum ALL files before unlinking ANY source.
    for row in current:
        src = safe_file(root, row["path"])
        dst = destination / "files" / row["path"]
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
        if sha256(dst) != row["sha256"]:
            raise RuntimeError(f"Backup verification failed: {row['path']}; no originals removed")
    manifest = {"format": FORMAT + "-archive", "repo_root": str(root.resolve()),
                "git_head": head(root), "created_utc": datetime.now(timezone.utc).isoformat(),
                "count": len(current), "bytes": sum(r["bytes"] for r in current),
                "entries": current}
    (destination / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    # All copies and the restore manifest are ready. Abort on a changed source.
    for row in current:
        src = safe_file(root, row["path"])
        if sha256(src) != row["sha256"]:
            raise RuntimeError(f"Original changed during archive: {row['path']}; archive retained")
        src.unlink()
    print(f"Archived {len(current)} compiled LaTeX artifacts: {destination}")
    print("Run core tests and review git status before local commit. No GitHub write occurred.")
    return destination


def restore(root: Path, archive_dir: Path) -> int:
    manifest = json.loads((archive_dir / "manifest.json").read_text(encoding="utf-8"))
    _verify_repo_binding(root, manifest, "-archive")
    entries = manifest.get("entries", [])
    if not isinstance(entries, list) or len(entries) != manifest.get("count"):
        raise RuntimeError("Invalid archive entries/count")
    seen: set[str] = set()
    # Preflight entire restore first; refuse overwrite.
    for row in entries:
        rel = row["path"]
        if rel in seen:
            raise RuntimeError(f"Duplicate archive path: {rel}")
        seen.add(rel)
        out = root.joinpath(*rel.split("/"))
        if not out.resolve().is_relative_to(root.resolve()) or out.exists() or out.is_symlink():
            raise RuntimeError(f"Restore destination exists or is unsafe: {rel}")
        source = safe_file(archive_dir / "files", rel)
        if sha256(source) != row["sha256"]:
            raise RuntimeError(f"Archived file checksum mismatch: {rel}")
    for row in entries:
        rel = row["path"]
        dest = root.joinpath(*rel.split("/"))
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(safe_file(archive_dir / "files", rel), dest)
        if sha256(dest) != row["sha256"]:
            raise RuntimeError(f"Restored file checksum mismatch: {rel}")
    print(f"Restored {len(entries)} files from {archive_dir}")
    return len(entries)


def duplicates(root: Path, external: Path) -> Path:
    """Report identical tracked Git blobs; NEVER remove or link these files."""
    items = tracked_blobs(root)
    groups: dict[str, list[str]] = {}
    for path, (sha, mode) in items.items():
        if mode == "100644":
            groups.setdefault(sha, []).append(path)
    repeated = {sha: sorted(paths) for sha, paths in groups.items() if len(paths) > 1}
    if not repeated:
        raise RuntimeError("No identical tracked blobs found")
    cp = git(root, "cat-file", "--batch-check=%(objectname) %(objectsize)",
             input_text="\n".join(sorted(repeated)) + "\n")
    sizes = {line.split()[0]: int(line.split()[1]) for line in cp.splitlines()}
    report = external / "duplicates_review"
    report.mkdir(parents=True, exist_ok=True)
    results = []
    for sha, paths in repeated.items():
        bytes_each = sizes[sha]
        all_output = all(p.startswith(("output/", "Reports/output/")) for p in paths)
        results.append({"blob_sha": sha, "identical_path_count": len(paths),
                        "bytes_each": bytes_each,
                        "extra_worktree_bytes": bytes_each * (len(paths)-1),
                        "all_generated_output_paths": all_output,
                        "paths": json.dumps(paths, ensure_ascii=False)})
    results.sort(key=lambda r: r["extra_worktree_bytes"], reverse=True)
    csv_path = report / "identical_git_blobs.csv"
    with csv_path.open("w", newline="", encoding="utf-8-sig") as f:
        fields = ("blob_sha", "identical_path_count", "bytes_each", "extra_worktree_bytes",
                  "all_generated_output_paths", "paths")
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(results)
    summary = {"git_head": head(root), "duplicate_groups": len(results),
               "duplicate_paths_beyond_one": sum(r["identical_path_count"]-1 for r in results),
               "extra_worktree_bytes": sum(r["extra_worktree_bytes"] for r in results),
               "generated_output_groups": sum(r["all_generated_output_paths"] for r in results),
               "note": "No files removed. Git already shares identical blobs in its object store. "
                       "Dynamic readers and thesis provenance may depend on every path."}
    (report / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"READ-ONLY: {summary['duplicate_groups']} identical groups, "
          f"{summary['duplicate_paths_beyond_one']} paths beyond first, "
          f"{summary['extra_worktree_bytes']:,} apparent duplicate working-tree bytes")
    print(f"Review: {csv_path}. NO AUTOMATIC DUPLICATE DELETIONS.")
    return csv_path


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("command", choices=("preview", "apply", "restore", "duplicates"), nargs="?", default="preview")
    ap.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    ap.add_argument("--external-root", type=Path)
    ap.add_argument("--plan", type=Path)
    ap.add_argument("--archive", type=Path)
    ap.add_argument("--confirm", default="")
    args = ap.parse_args(argv)
    root = args.root.resolve()
    if not (root / ".git").exists():
        ap.error(f"Not a Git repository checkout root: {root}")
    outside = external_root(root, args.external_root)
    if outside.is_relative_to(root):
        ap.error("External archive location must be outside the repository")
    try:
        if args.command == "preview":
            preview(root, outside)
        elif args.command == "apply":
            if args.plan is None:
                ap.error("--plan is required for apply")
            apply(root, outside, args.plan.resolve(), args.confirm)
        elif args.command == "restore":
            if args.archive is None:
                ap.error("--archive is required for restore")
            restore(root, args.archive.resolve())
        else:
            duplicates(root, outside)
    except (RuntimeError, OSError, ValueError, KeyError, json.JSONDecodeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
