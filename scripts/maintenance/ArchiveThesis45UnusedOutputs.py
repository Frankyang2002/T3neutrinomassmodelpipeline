#!/usr/bin/env python3
"""Conservative, local-only T3 thesis artifact archive tool.

The GitHub repository is NEVER modified. Only listed generated-output
DIRECTORIES are eligible. No source, tests, thesis data, configs or symbolic
Lagrangians are archived by this tool.

Default mode: read-only audit (also writes a preflight CSV).
--archive-copy: create a verified ZIP; leave original paths untouched.
--apply: create a verified ZIP, then MOVE eligible directories to an adjacent
         quarantine folder, preserving their project relative paths.
--restore-zip PATH: restore ZIP contents, refusing all overwrites.

Safety checks: exact tracked inventory count/bytes, no local modifications,
no textual references in source/tests/configs or other live result metadata,
no symlinks, and protected-path guardrails.

NOTE: Absence of a literal reference is not proof of absence of dynamic usage.
Use the CSV classifications and preflight report before --apply.
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

# Anything touching these subtrees must NOT be archived via this utility.
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
# Source code and data-bearing configs; stale Markdown notes intentionally
# excluded as evidence of active runtime dependency.
SOURCE_SUFFIXES = {".py", ".wl", ".wls", ".m", ".ps1", ".sh", ".bat", ".cmd",
                   ".tex", ".bib", ".json", ".toml", ".yaml", ".yml", ".ini", ".cfg"}
# Only metadata among output trees; historical outputs may contain huge text
# dumps that are not active imports. These keys are still used for lineage.
OUTPUT_METADATA_NAMES = ("manifest", "provenance", "summary", "settings", "status", "index")
EXCLUDED_SCAN_PREFIXES = (".git/", "archives/", ".cleanup_quarantine/",
                          "output/cleanup_audit/", "audits/", "docs/")


def clean_rel(value: str) -> str:
    p = PurePosixPath(value.strip().replace("\\", "/"))
    if p.is_absolute() or ".." in p.parts or not p.parts or p.parts[0] not in ("output", "Reports"):
        raise ValueError(f"Non-output or unsafe path in candidate manifest: {value!r}")
    result = p.as_posix().rstrip("/")
    if any(result == q.rstrip("/") or result.startswith(q) for q in PROTECTED):
        raise ValueError(f"Protected thesis/source path in candidate manifest: {result}")
    return result


def load_candidates() -> list[dict]:
    with PLAN.open("r", encoding="utf-8-sig", newline="") as f:
        rows = list(csv.DictReader(f))
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
    if any(x != y and (x.startswith(y + "/") or y.startswith(x + "/"))
           for x in seen for y in seen):
        raise ValueError("Overlapping candidate directories")
    return rows


def git_root() -> str | None:
    try:
        p = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "--show-toplevel"],
                           capture_output=True, text=True, check=True)
        return p.stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def git_status(path: str) -> str:
    try:
        p = subprocess.run(["git", "-C", str(ROOT), "status", "--porcelain", "--", path],
                           capture_output=True, text=True, check=True)
        return p.stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "git-status-unavailable"


def file_inventory(folder: Path) -> tuple[list[Path], int, str | None]:
    if not folder.is_dir() or folder.is_symlink():
        return [], 0, "missing-or-symlink"
    files: list[Path] = []
    bytes_total = 0
    for base, dirs, names in os.walk(folder, followlinks=False):
        for name in dirs:
            if (Path(base) / name).is_symlink():
                return [], 0, "symlink-directory"
        for name in names:
            f = Path(base) / name
            if not f.is_file() or f.is_symlink():
                return [], 0, "symlink-or-nonregular-file"
            files.append(f)
            bytes_total += f.stat().st_size
    return sorted(files), bytes_total, None


def should_scan(rel: str, filename: str, suffix: str) -> bool:
    if rel.startswith(EXCLUDED_SCAN_PREFIXES):
        return False
    if rel.startswith("output/") or rel.startswith("Reports/output/"):
        stem = Path(filename).stem.lower()
        return suffix in {".json", ".csv", ".txt"} and any(stem.startswith(n) for n in OUTPUT_METADATA_NAMES)
    if filename.lower().startswith("readme") or suffix == ".md":
        return False
    return suffix in SOURCE_SUFFIXES


def references(rows: list[dict]) -> dict[str, list[str]]:
    """Conservative literal-reference scan for all potential candidate basenames.

    Excludes candidate self-references but scans other candidate provenance.
    Unicode-agnostic byte search handles large reports without loading them
    into memory. Case folded matches in UTF-8/ASCII source files.
    """
    lookups = {r["path"]: (r["path"].encode().lower(), r["path"].split("/")[-1].encode().lower())
               for r in rows}
    hits: dict[str, list[str]] = {r["path"]: [] for r in rows}
    for root, dirs, filenames in os.walk(ROOT):
        base = Path(root)
        dirs[:] = [d for d in dirs if not d.startswith(".git") and not d.startswith(".cleanup_quarantine")
                   and not (base == ROOT and d in {"archives", "audits", "docs", "node_modules"})]
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
                # Search rolling overlapping chunks for UTF-8 filename tokens.
                with path.open("rb") as f:
                    prev = b""
                    while chunk := f.read(1024 * 1024):
                        data = (prev + chunk).lower().replace(b"\\\\", b"/").replace(b"\\", b"/")
                        for name, (pbytes, bbytes) in lookups.items():
                            if rel == name or rel.startswith(name + "/"):
                                continue
                            if len(hits[name]) >= 8:
                                continue
                            if pbytes in data or (len(bbytes) >= 9 and bbytes in data):
                                hits[name].append(rel)
                        prev = chunk[-512:]
            except (OSError, PermissionError):
                # Inaccessible live metadata means never archive automatically.
                for name in hits:
                    if len(hits[name]) < 8:
                        hits[name].append("UNREADABLE_SCAN_INPUT:" + rel)
    return hits


def evaluate(rows: list[dict]) -> list[dict]:
    source_refs = references(rows)
    results = []
    for row in rows:
        rel = row["path"]
        fs, total_bytes, inventory_error = file_inventory(ROOT / rel)
        status = "READY"
        reasons = []
        if row["classification"] != "archive_if_unreferenced":
            status = "MANUAL_REVIEW"
            reasons.append("not-auto-approved")
        if inventory_error:
            status = "HOLD"
            reasons.append(inventory_error)
        if len(fs) != row["expected_files"] or total_bytes != row["expected_bytes"]:
            status = "HOLD"
            reasons.append("local-inventory-differs-from-audited-main")
        if source_refs[rel]:
            status = "HOLD"
            reasons.append("referenced-outside-own-folder")
        gs = git_status(rel)
        if gs:
            status = "HOLD"
            reasons.append("local-git-changes-or-unavailable")
        results.append({**row, "found_files": len(fs), "found_bytes": total_bytes,
                        "status": status, "reason_codes": ";".join(reasons),
                        "referencing_paths": ";".join(source_refs[rel]), "git_status": gs})
    return results


def write_report(rows: list[dict]) -> None:
    fields = ["path", "classification", "status", "expected_files", "found_files",
              "expected_bytes", "found_bytes", "reason", "reason_codes", "referencing_paths", "git_status"]
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    with REPORT.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for data in iter(lambda: f.read(1024 * 1024), b""):
            h.update(data)
    return h.hexdigest()


def make_verified_zip(rows: list[dict], zip_path: Path) -> tuple[dict, list[str]]:
    files: list[tuple[str, Path]] = []
    for row in rows:
        fs, total, error = file_inventory(ROOT / row["path"])
        if error or len(fs) != row["expected_files"] or total != row["expected_bytes"]:
            raise RuntimeError("Candidate changed since audit: " + row["path"])
        files.extend((f.relative_to(ROOT).as_posix(), f) for f in fs)
    if not files:
        raise RuntimeError("No qualifying candidate files")
    entries = [{"path": rel, "bytes": src.stat().st_size, "sha256": sha256_file(src)}
               for rel, src in files]
    record = {"schema": "t3-thesis45-conservative-local-archive-v1",
              "timestamp_utc": datetime.now(timezone.utc).isoformat(),
              "project_head": subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"],
                                              capture_output=True, text=True).stdout.strip(),
              "physics_source_moved": False,
              "candidate_directories": [x["path"] for x in rows],
              "files": entries,
              "note": "Conservative source-reference scan. ZIP is reversible; no remote GitHub operations."}
    zip_path.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(zip_path, "x", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as z:
        for rel, src in files:
            z.write(src, arcname=rel)
        z.writestr("ARCHIVE_MANIFEST.json", json.dumps(record, indent=2) + "\n")
    # Independently verify every member's SHA256 after ZIP closure.
    with zipfile.ZipFile(zip_path, "r") as z:
        listed = set(z.namelist())
        if listed != {x["path"] for x in entries} | {"ARCHIVE_MANIFEST.json"}:
            raise RuntimeError("ZIP contents not identical to manifest")
        for item in entries:
            h = hashlib.sha256()
            with z.open(item["path"]) as f:
                for data in iter(lambda: f.read(1024 * 1024), b""):
                    h.update(data)
            if h.hexdigest() != item["sha256"]:
                raise RuntimeError("ZIP hash verification failed for " + item["path"])
    return record, [x["path"] for x in entries]


def restore_zip(path: Path) -> None:
    if not path.is_file():
        raise FileNotFoundError(path)
    with zipfile.ZipFile(path) as z:
        record = json.loads(z.read("ARCHIVE_MANIFEST.json"))
        if record.get("schema") != "t3-thesis45-conservative-local-archive-v1":
            raise ValueError("Unknown archive manifest format")
        entries = record["files"]
        for item in entries:
            rel = clean_rel(item["path"])
            if not rel.startswith(("output/", "Reports/output/")):
                raise ValueError("Unsafe archive content: " + rel)
            if (ROOT / rel).exists():
                raise FileExistsError("Refusing to overwrite existing path: " + rel)
        for item in entries:
            rel = clean_rel(item["path"])
            dest = ROOT / rel
            dest.parent.mkdir(parents=True, exist_ok=True)
            h = hashlib.sha256()
            with z.open(rel) as source, dest.open("xb") as target:
                for chunk in iter(lambda: source.read(1024 * 1024), b""):
                    target.write(chunk)
                    h.update(chunk)
            if h.hexdigest() != item["sha256"]:
                raise RuntimeError("Restored file hash mismatch: " + rel)
    print("RESTORED", len(entries), "files to", ROOT)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    modes = p.add_mutually_exclusive_group()
    modes.add_argument("--archive-copy", action="store_true", help="ZIP qualifying folders; leave originals in place")
    modes.add_argument("--apply", action="store_true", help="ZIP and move qualifying folders to adjacent quarantine")
    modes.add_argument("--restore-zip", type=Path, help="Restore archived files to their original paths")
    args = p.parse_args()
    if git_root() is None or Path(git_root()).resolve() != ROOT.resolve():
        p.error("Run from a local Git checkout with this script at scripts/maintenance; GitHub is not modified")
    if args.restore_zip:
        restore_zip(args.restore_zip.resolve())
        return 0
    rows = evaluate(load_candidates())
    write_report(rows)
    for row in rows:
        print(f"{row['status']:14} {row['found_files']:4} files  {row['found_bytes']/1e6:6.2f} MB  {row['path']}")
        if row["reason_codes"]:
            print("               reasons:", row["reason_codes"])
    ready = [r for r in rows if r["status"] == "READY"]
    print(f"\nREADY: {len(ready)} directories, {sum(r['found_files'] for r in ready)} files, "
          f"{sum(r['found_bytes'] for r in ready)/1e6:.2f} MB")
    print("Full preflight:", REPORT)
    if not (args.archive_copy or args.apply):
        print("AUDIT ONLY: no original files moved. To create a verified archive, use --archive-copy or --apply.")
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
    quarantine = OUTSIDE_ARCHIVES / f"quarantine_{stamp}"
    quarantine.mkdir(parents=True, exist_ok=False)
    moved = []
    try:
        for row in ready:
            src = ROOT / row["path"]
            dst = quarantine / row["path"]
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(src), str(dst))
            moved.append({"original": row["path"], "quarantined": str(dst)})
        (quarantine / "MOVE_MANIFEST.json").write_text(json.dumps({"zip": str(zip_path), "moved": moved,
                                     "record": record["schema"]}, indent=2), encoding="utf-8")
    except Exception:
        for item in reversed(moved):
            src = Path(item["quarantined"])
            dst = ROOT / item["original"]
            if src.exists() and not dst.exists():
                dst.parent.mkdir(parents=True, exist_ok=True)
                shutil.move(str(src), str(dst))
        raise
    print("MOVED to quarantine:", len(moved), "directories. Restore from ZIP via --restore-zip if needed.")
    print("GitHub was NOT modified; review git status before any local commit.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
