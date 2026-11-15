#!/usr/bin/env python3
"""Conservative, local-only T3 thesis output archival utility.

Only explicitly listed, unreferenced, unchanged generated-output directories
can be archived. Expected file counts and bytes refer to the committed Git
blobs, not working-tree byte lengths (which differ under Windows autocrlf).
Every local file is individually compared with its committed blob using Git's
path-aware clean filters. Untracked files, changed content, unreadable files,
protected paths and manually reviewed cases always remain in place.

Usage from the repository root:
    python -X utf8 scripts/maintenance/ArchiveThesis45UnusedOutputs.py
    python -X utf8 scripts/maintenance/ArchiveThesis45UnusedOutputs.py --archive-copy
    python -X utf8 scripts/maintenance/ArchiveThesis45UnusedOutputs.py --apply
    python -X utf8 scripts/maintenance/ArchiveThesis45UnusedOutputs.py --restore-zip PATH

No remote operations, commits, recalculations or physics-data edits.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import shutil
import subprocess
import sys
import zipfile

ROOT = Path(__file__).resolve().parents[2]
PLAN = ROOT / "audits/thesis45_archive_candidates.csv"
REPORT = ROOT / "audits/thesis45_archive_preflight.csv"
OUTSIDE_ARCHIVES = ROOT.parent / "T3_thesis45_local_archives"

# Protected scientific source, figure provenance and numerical dependencies.
# This is intentionally the original project's exclusion list.
PROTECTED = (
    "Lagrangian/", "RGE/", "Numerical/", "physics/", "configs/",
    "tests/", "scripts/", "docs/", "input/", "archives/",
    "Reports/output/thesis_general_5_1/", "Reports/output/thesis_general_5_2/",
    "Reports/output/thesis_general_5_3/", "Reports/output/thesis_qualitative/",
    "output/studies/general_scalar_extended/thesis_general_all16/",
    "output/studies/general_scalar_extended/optimal_yukawa_thesis_refresh/",
    "output/studies/general_scalar_extended/experimental_simultaneous_eft1_",
    "output/general_observable/", "output/full16_symbolic/",
    "output/parallel_t3/thesis_two_benchmarks_all16_stagefreeze_v2/",
    "output/parallel_t3/eft1_llss_thesis_consistent_v2/",
    "output/parallel_t3/eft1_llss_thesis_consistent_v1/",
    "output/parallel_t3/thesis_generalised_tables_v1/",
    "output/parallel_t3/rough_inputs/", "output/parallel_t3/rough_wilson/",
    "output/parallel_t3/rough_results/",
)
SOURCE_SUFFIXES = {".py", ".wl", ".wls", ".m", ".ps1", ".sh", ".bat", ".cmd",
                   ".tex", ".bib", ".json", ".toml", ".yaml", ".yml", ".ini", ".cfg"}
OUTPUT_METADATA_NAMES = ("manifest", "provenance", "summary", "settings", "status", "index")
EXCLUDED_SCAN_PREFIXES = (".git/", "archives/", ".cleanup_quarantine/",
                          "output/cleanup_audit/", "audits/", "docs/")


def clean_rel(value: str) -> str:
    p = PurePosixPath(value.strip().replace("\\", "/"))
    if p.is_absolute() or ".." in p.parts or not p.parts or p.parts[0] not in ("output", "Reports"):
        raise ValueError(f"Non-output or unsafe path in candidate manifest: {value!r}")
    result = p.as_posix().rstrip("/")
    if result in {"output", "Reports", "Reports/output"}:
        raise ValueError(f"Unsafe top-level candidate: {result}")
    if not result.startswith(("output/", "Reports/output/")):
        raise ValueError(f"Not a generated output path: {result}")
    if any(result == q.rstrip("/") or result.startswith(q) for q in PROTECTED):
        raise ValueError(f"Protected thesis/source path in candidate manifest: {result}")
    return result


def load_candidates() -> list[dict]:
    with PLAN.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    required = {"path", "classification", "expected_files", "expected_bytes", "reason"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError("Missing columns or entries in thesis45_archive_candidates.csv")
    seen: set[str] = set()
    for row in rows:
        row["path"] = clean_rel(row["path"])
        if row["path"] in seen:
            raise ValueError(f"Duplicate candidate: {row['path']}")
        seen.add(row["path"])
        if row["classification"] not in {"archive_if_unreferenced", "manual_review"}:
            raise ValueError(f"Unknown classification for {row['path']}")
        row["expected_files"] = int(row["expected_files"])
        row["expected_bytes"] = int(row["expected_bytes"])
        if row["expected_files"] < 1 or row["expected_bytes"] < 0:
            raise ValueError(f"Invalid expected inventory: {row['path']}")
    if any(x != y and (x.startswith(y + "/") or y.startswith(x + "/"))
           for x in seen for y in seen):
        raise ValueError("Overlapping candidate directories")
    return rows


def git(*args: str) -> bytes:
    result = subprocess.run(["git", "-C", str(ROOT), *args], capture_output=True, check=False)
    if result.returncode:
        raise RuntimeError(f"git {' '.join(args)}: {result.stderr.decode('utf-8', 'replace')}")
    return result.stdout


def git_root() -> str | None:
    try:
        return os.fsdecode(git("rev-parse", "--show-toplevel")).strip()
    except (OSError, RuntimeError):
        return None


def git_status(path: str) -> str:
    try:
        return os.fsdecode(git("status", "--porcelain", "--", path)).strip()
    except (OSError, RuntimeError):
        return "git-status-unavailable"


def git_inventory(folder: str) -> dict[str, tuple[str, int]]:
    """Read exact tracked paths, SHA-1/256 object IDs and canonical blob sizes.

    The binary -z output safely supports spaces and Unicode paths. This
    compares the archived manifest against HEAD, not Windows CRLF bytes.
    """
    records = git("ls-tree", "-r", "-l", "-z", "HEAD", "--", folder)
    files: dict[str, tuple[str, int]] = {}
    for entry in records.split(b"\0"):
        if not entry:
            continue
        metadata, raw_name = entry.split(b"\t", 1)
        fields = metadata.split()
        if len(fields) != 4 or fields[1] != b"blob" or fields[3] == b"-":
            raise RuntimeError(f"Not a Git blob in candidate: {folder}")
        name = raw_name.decode("utf-8", errors="surrogateescape")
        if not name.startswith(folder + "/"):
            raise RuntimeError(f"Unexpected Git tree path for {folder}: {name}")
        files[name] = (fields[2].decode("ascii"), int(fields[3]))
    return files


def git_clean_blob_hash(rel: str) -> str:
    """Hash worktree file with the same attributes/EOL filters Git uses.

    `git hash-object --path=...` computes the canonical Git blob hash even
    if the Windows worktree contains CRLF rather than tracked LF. Unlike
    merely removing CR characters, this also respects .gitattributes and
    detects changes to numbers, CSV cells and binary files.
    """
    return git("hash-object", "--path=" + rel, "--", rel).decode("ascii").strip()


def file_inventory(folder: Path) -> tuple[list[Path], int, str | None]:
    if not folder.is_dir() or folder.is_symlink():
        return [], 0, "missing-or-symlink"
    files: list[Path] = []
    total_bytes = 0
    for base, dirs, names in os.walk(folder, followlinks=False):
        for name in dirs:
            if (Path(base) / name).is_symlink():
                return [], 0, "symlink-directory"
        for name in names:
            f = Path(base) / name
            if not f.is_file() or f.is_symlink():
                return [], 0, "symlink-or-nonregular-file"
            files.append(f)
            total_bytes += f.stat().st_size
    return sorted(files), total_bytes, None


def should_scan(rel: str, filename: str, suffix: str) -> bool:
    if rel.startswith(EXCLUDED_SCAN_PREFIXES):
        return False
    if rel.startswith(("output/", "Reports/output/")):
        stem = Path(filename).stem.lower()
        return suffix in {".json", ".csv", ".txt"} and any(
            stem.startswith(name) for name in OUTPUT_METADATA_NAMES
        )
    if filename.lower().startswith("readme") or suffix == ".md":
        return False
    return suffix in SOURCE_SUFFIXES


def references(rows: list[dict]) -> dict[str, list[str]]:
    """Conservative literal-reference scan, preserving original semantics."""
    lookups = {
        row["path"]: (row["path"].encode().lower(), row["path"].split("/")[-1].encode().lower())
        for row in rows
    }
    hits: dict[str, list[str]] = {row["path"]: [] for row in rows}
    for root, dirs, filenames in os.walk(ROOT):
        base = Path(root)
        dirs[:] = [name for name in dirs if not name.startswith(".git")
                   and not name.startswith(".cleanup_quarantine")
                   and not (base == ROOT and name in {"archives", "audits", "docs", "node_modules"})]
        for filename in filenames:
            path = base / filename
            try:
                rel = path.relative_to(ROOT).as_posix()
            except ValueError:
                continue
            if path.is_symlink() or not should_scan(rel, filename, path.suffix.lower()):
                continue
            if not path.is_file():
                continue
            try:
                with path.open("rb") as handle:
                    previous = b""
                    while chunk := handle.read(1024 * 1024):
                        data = (previous + chunk).lower().replace(b"\\\\", b"/").replace(b"\\", b"/")
                        for name, (whole, basename) in lookups.items():
                            if rel == name or rel.startswith(name + "/") or len(hits[name]) >= 8:
                                continue
                            if whole in data or (len(basename) >= 9 and basename in data):
                                hits[name].append(rel)
                        previous = chunk[-512:]
            except (OSError, PermissionError):
                for name in hits:
                    if len(hits[name]) < 8:
                        hits[name].append("UNREADABLE_SCAN_INPUT:" + rel)
    return hits


def evaluate_one(row: dict, source_refs: dict[str, list[str]]) -> dict:
    rel = row["path"]
    files, local_bytes, inventory_error = file_inventory(ROOT / rel)
    reasons: list[str] = []
    status = "READY"
    if row["classification"] != "archive_if_unreferenced":
        status = "MANUAL_REVIEW"
        reasons.append("not-auto-approved")
    if inventory_error:
        status = "HOLD"
        reasons.append(inventory_error)

    committed_files: dict[str, tuple[str, int]] = {}
    try:
        committed_files = git_inventory(rel)
        committed_bytes = sum(size for _, size in committed_files.values())
        if (len(committed_files) != row["expected_files"] or
                committed_bytes != row["expected_bytes"]):
            status = "HOLD"
            reasons.append("committed-inventory-differs-from-audited-main")
        local_names = {file.relative_to(ROOT).as_posix() for file in files}
        if set(committed_files) != local_names:
            status = "HOLD"
            reasons.append("local-file-set-differs-from-HEAD")
        elif not inventory_error:
            for name, (expected_hash, _) in committed_files.items():
                if git_clean_blob_hash(name) != expected_hash:
                    status = "HOLD"
                    reasons.append("local-content-differs-from-HEAD")
                    break
    except (OSError, ValueError, RuntimeError, UnicodeError) as exc:
        status = "HOLD"
        reasons.append("git-verification-failed:" + str(exc)[:150])
        committed_bytes = -1

    if source_refs[rel]:
        status = "HOLD"
        reasons.append("referenced-outside-own-folder")
    gs = git_status(rel)
    if gs:
        status = "HOLD"
        reasons.append("local-git-changes-or-unavailable")
    return {
        **row,
        "found_files": len(files),
        "found_bytes": local_bytes,
        "committed_files": len(committed_files),
        "committed_bytes": committed_bytes,
        "status": status,
        "reason_codes": ";".join(reasons),
        "referencing_paths": ";".join(source_refs[rel]),
        "git_status": gs,
    }


def evaluate(rows: list[dict]) -> list[dict]:
    source_refs = references(rows)
    return [evaluate_one(row, source_refs) for row in rows]


def write_report(rows: list[dict]) -> None:
    fields = ["path", "classification", "status", "expected_files", "found_files",
              "expected_bytes", "found_bytes", "committed_files", "committed_bytes",
              "reason", "reason_codes", "referencing_paths", "git_status"]
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    with REPORT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for data in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(data)
    return digest.hexdigest()


def verify_archive_selection(rows: list[dict]) -> list[tuple[str, Path]]:
    """Revalidate Git objects and local content before making an archive."""
    items: list[tuple[str, Path]] = []
    for row in rows:
        if row["classification"] != "archive_if_unreferenced":
            raise RuntimeError("Manual-review result cannot be archived automatically: " + row["path"])
        rel = clean_rel(row["path"])
        fs, _, error = file_inventory(ROOT / rel)
        committed = git_inventory(rel)
        if error or len(committed) != row["expected_files"] or sum(
            size for _, size in committed.values()
        ) != row["expected_bytes"]:
            raise RuntimeError("Changed committed/working tree inventory: " + rel)
        names = {p.relative_to(ROOT).as_posix() for p in fs}
        if names != set(committed):
            raise RuntimeError("Changed local file set: " + rel)
        if git_status(rel):
            raise RuntimeError("Git status changed since audit: " + rel)
        for name, (blob_hash, _) in committed.items():
            if git_clean_blob_hash(name) != blob_hash:
                raise RuntimeError("Changed content since audit: " + name)
        items.extend((p.relative_to(ROOT).as_posix(), p) for p in fs)
    return items


def make_verified_zip(rows: list[dict], zip_path: Path) -> tuple[dict, list[str]]:
    files = verify_archive_selection(rows)
    if not files:
        raise RuntimeError("No qualifying candidate files")
    entries = [{"path": rel, "bytes": src.stat().st_size, "sha256": sha256_file(src)}
               for rel, src in files]
    record = {
        "schema": "t3-thesis45-conservative-local-archive-v1",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "project_head": git("rev-parse", "HEAD").decode("ascii").strip(),
        "physics_source_moved": False,
        "candidate_directories": [row["path"] for row in rows],
        "files": entries,
        "note": "Git blob verification accounts for platform EOL. Verified ZIP is reversible. No remote operations.",
    }
    zip_path.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(zip_path, "x", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as archive:
        for rel, src in files:
            archive.write(src, arcname=rel)
        archive.writestr("ARCHIVE_MANIFEST.json", json.dumps(record, indent=2) + "\n")
    with zipfile.ZipFile(zip_path, "r") as archive:
        actual = set(archive.namelist())
        if actual != {entry["path"] for entry in entries} | {"ARCHIVE_MANIFEST.json"}:
            raise RuntimeError("ZIP contents differ from manifest")
        for entry in entries:
            digest = hashlib.sha256()
            with archive.open(entry["path"]) as source:
                for data in iter(lambda: source.read(1024 * 1024), b""):
                    digest.update(data)
            if digest.hexdigest() != entry["sha256"]:
                raise RuntimeError("ZIP SHA-256 mismatch: " + entry["path"])
    return record, [entry["path"] for entry in entries]


def restore_zip(path: Path) -> None:
    if not path.is_file():
        raise FileNotFoundError(path)
    with zipfile.ZipFile(path) as archive:
        record = json.loads(archive.read("ARCHIVE_MANIFEST.json"))
        if record.get("schema") != "t3-thesis45-conservative-local-archive-v1":
            raise ValueError("Unknown archive manifest format")
        entries = record["files"]
        expected = {item["path"] for item in entries} | {"ARCHIVE_MANIFEST.json"}
        if set(archive.namelist()) != expected or len(entries) != len(expected) - 1:
            raise ValueError("Unexpected members in archive")
        for item in entries:
            rel = clean_rel(item["path"])
            dest = ROOT / rel
            if dest.exists() or any(parent.is_symlink() for parent in dest.parents if parent.is_relative_to(ROOT)):
                raise FileExistsError("Refusing existing or symlink destination: " + rel)
        for item in entries:
            rel = clean_rel(item["path"])
            dest = ROOT / rel
            dest.parent.mkdir(parents=True, exist_ok=True)
            digest = hashlib.sha256()
            with archive.open(rel) as source, dest.open("xb") as target:
                for chunk in iter(lambda: source.read(1024 * 1024), b""):
                    target.write(chunk)
                    digest.update(chunk)
            if digest.hexdigest() != item["sha256"]:
                raise RuntimeError("Restored SHA-256 mismatch: " + rel)
    print("RESTORED", len(entries), "files to", ROOT)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--archive-copy", action="store_true", help="ZIP qualifying folders without deleting them")
    modes.add_argument("--apply", action="store_true", help="ZIP then move qualifying folders to adjacent quarantine")
    modes.add_argument("--restore-zip", type=Path, help="Restore archive members to their original paths")
    args = parser.parse_args()
    if git_root() is None or Path(git_root()).resolve() != ROOT.resolve():
        parser.error("Run from a local checkout with this script at scripts/maintenance")
    if args.restore_zip:
        restore_zip(args.restore_zip.resolve())
        return 0
    rows = evaluate(load_candidates())
    write_report(rows)
    for row in rows:
        print(f"{row['status']:14} {row['found_files']:4} files  {row['found_bytes']/1e6:6.2f} MB  {row['path']}")
        if row["reason_codes"]:
            print("               reasons:", row["reason_codes"])
    ready = [row for row in rows if row["status"] == "READY"]
    print(f"\nREADY: {len(ready)} directories, {sum(row['found_files'] for row in ready)} files, "
          f"{sum(row['found_bytes'] for row in ready)/1e6:.2f} MB")
    print("Full preflight:", REPORT)
    if not (args.archive_copy or args.apply):
        print("AUDIT ONLY: no original files moved. --archive-copy creates a backup; --apply also moves files.")
        return 0
    if not ready:
        print("No directories passed conservative audit. Nothing archived.")
        return 0
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    zip_path = OUTSIDE_ARCHIVES / f"t3_thesis45_archive_{stamp}.zip"
    record, paths = make_verified_zip(ready, zip_path)
    print("VERIFIED ARCHIVE:", zip_path, "files:", len(paths), "SHA256:", sha256_file(zip_path))
    if not args.apply:
        print("COPY ONLY: original files were not moved.")
        return 0
    # Recheck immediately before moving. Refuse any changed data since ZIP creation.
    for row in ready:
        verify_archive_selection([row])
    quarantine = OUTSIDE_ARCHIVES / f"quarantine_{stamp}"
    quarantine.mkdir(parents=True, exist_ok=False)
    moved: list[dict] = []
    try:
        for row in ready:
            src = ROOT / row["path"]
            dst = quarantine / row["path"]
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(src), str(dst))
            moved.append({"original": row["path"], "quarantined": str(dst)})
        (quarantine / "MOVE_MANIFEST.json").write_text(json.dumps(
            {"zip": str(zip_path), "moved": moved, "record": record["schema"]}, indent=2
        ), encoding="utf-8")
    except Exception:
        for entry in reversed(moved):
            src = Path(entry["quarantined"])
            dst = ROOT / entry["original"]
            if src.exists() and not dst.exists():
                dst.parent.mkdir(parents=True, exist_ok=True)
                shutil.move(str(src), str(dst))
        raise
    print("MOVED to quarantine:", len(moved), "directories. Restore from ZIP if needed.")
    print("GitHub not modified; inspect local git status before any commit.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
