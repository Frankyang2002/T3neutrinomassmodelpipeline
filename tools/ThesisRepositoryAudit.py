#!/usr/bin/env python3
"""Inventory the local T3 checkout, review source utility roles, and archive reviewed files.

Default: READ-ONLY scientific/source audit, writing reports outside the checkout.
--apply: requires a reviewed selection, a matching Git HEAD and a confirmation token.
--restore: restores SHA-256 verified archived bytes without overwriting existing files.

This tool never writes Git history, contacts GitHub, or executes scientific scripts.
"""
from __future__ import annotations

import argparse
import ast
from collections import Counter, defaultdict
import csv
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

FORMAT = "t3-thesis-repository-review-v1"
CONFIRM = "ARCHIVE_REVIEWED_THESIS_FILES"
REPORT_DIRNAME = "thesis_repository_review"
SELF = "tools/ThesisRepositoryAudit.py"
SELF_TEST = "tests/python/test_thesis_repository_audit.py"
IGNORED_INDEXES = frozenset({
    "scripts/INDEX.md", "docs/research/INDEX.md", "tools/README.md",
    "output/maintenance/research_organization_plan.json",
    "output/maintenance/research_organization_manifest.json",
    "tools/AuditScriptDependencies.py",  # its static preserve-list is not a call
    "tools/RefineScriptDependencyAudit.py",  # its static inventory is not a call
})
TEXT_SUFFIXES = frozenset({
    ".py", ".wl", ".wls", ".m", ".ps1", ".md", ".rst", ".txt",
    ".tex", ".bib", ".json", ".jsonl", ".csv", ".tsv", ".yaml",
    ".yml", ".toml", ".ini", ".cfg", ".sh", ".bat", ".cmd", ".html",
})
SOURCE_SUFFIXES = frozenset({".py", ".wl", ".wls", ".m", ".ps1", ".sh", ".bat", ".cmd"})

# Checked against read-only GitHub main at 806cdf08608dbe3f5e0dcd8af785cfad16f26af6.
# These are one-time *repository operations*, not numerical/physics implementations.
# Even these paths may be archived ONLY if the current local reference scan clears them.
ONE_TIME_TOOL_REVIEW = frozenset({
    "tools/AuditDocumentationRetirement.py",
    "tools/AuditLegacyNotes.py",
    "tools/AuditOneTimeArtifacts.py",
    "tools/ConsolidateHistoricalDocs.py",
    "tools/FinishHistoricalDocsConsolidation.py",
    "tools/RenameNumberedResearchFiles.py",
    "tools/MoveStudyFigures.py",
})

# Standalone installers, literal-name migrations and report-retirement scripts.
# These were inspected individually, not selected by a generic Fix/Old/Batch pattern.
ONE_TIME_SCRIPT_REVIEW = frozenset({
    "scripts/InstallT3RGBetaDoubletMassFix.py",
    "scripts/RenameT3ModelFolders.py",
    "scripts/RemoveOldGeneralScalarDisplay.ps1",
    "scripts/RetireOldStudyDisplays.ps1",
    "scripts/RetireSupersededT3HTMLGenerators.ps1",
    "scripts/maintenance/Remove-LegacyReadmes.ps1",
})

# Preserve facilities needed for scientific reproducibility or archive restore.
PROTECTED = frozenset({
    "tools/AuditScriptDependencies.py", "tools/RefineScriptDependencyAudit.py",
    "tools/RebuildCurrentReports.py", "tools/RebuildReports.py",
    "tools/RebuildSavedStudyPlots.py", "tools/RedirectStudyFigures.py",
    "tools/LabelSensitivityDirections.py", "tools/PlotDisplayWrapper.py",
    "tools/ConsolidateFittedBenchmarks.py",
    SELF,
    "scripts/RunCoreTests.ps1", "scripts/RunBaselineTests.ps1",
    "scripts/RunT3GeneralSymbolicStudy.py",
    "scripts/RunT3ReducedScalarSobolComparison.py",
    "scripts/RunT3GeneralScalar48Campaign.py",
    "scripts/RunT3GeneralFermionThresholdScan.py",
    "scripts/maintenance/ApplyCuratedArchive.py",
    "scripts/maintenance/ArchiveGeneratedResults.py",  # may restore previous output archives
    "scripts/maintenance/CleanRepository.py",  # may restore previous archives
    "scripts/maintenance/ExtremeRepositoryCleanup.py",  # may restore previous archives
    "scripts/maintenance/ReduceHistoricalTests.py",  # restore previous test archive
    "scripts/maintenance/OrganizeRepositoryPresentation.py",  # restore previous archive
    "scripts/maintenance/OrganizeResearchFiles.py",  # restore previous reorganisations
    "scripts/maintenance/SmokeTestFullPipeline.py",
})
ROOT_JUNK = "t -q --tb=line"
# These inspected files only list the pager artifact as a retirement target or
# construct it in temporary test fixtures; none consumes the real root file.
PAGER_CLEANUP_MENTIONS = frozenset({
    "tools/AuditOneTimeArtifacts.py",
    "tests/python/test_one_time_artifact_cleanup.py",
    "scripts/maintenance/OrganizeRepositoryPresentation.py",
    "scripts/maintenance/ExtremeRepositoryCleanup.py",
})


def checked_path(root: Path, value: str) -> Path:
    if not isinstance(value, str) or not value or "\\" in value:
        raise ValueError(f"Invalid project-relative path: {value!r}")
    posix = PurePosixPath(value)
    if posix.is_absolute() or any(part in ("..", ".") for part in posix.parts):
        raise ValueError(f"Unsafe project-relative path: {value!r}")
    result = root.joinpath(*posix.parts)
    if not result.resolve().is_relative_to(root.resolve()):
        raise ValueError(f"Path escapes repository: {value!r}")
    return result


def git(root: Path, *args: str, binary: bool = False) -> str | bytes:
    p = subprocess.run(["git", *args], cwd=root, capture_output=True, check=False)
    if p.returncode:
        raise RuntimeError(f"git {' '.join(args)} failed: {p.stderr.decode('utf-8', 'replace')[:350]}")
    return p.stdout if binary else p.stdout.decode("utf-8", "surrogateescape").strip()


def git_names(root: Path, *args: str) -> list[str]:
    raw = git(root, *args, binary=True)
    assert isinstance(raw, bytes)
    return [x.decode("utf-8", "surrogateescape").replace("\\", "/")
            for x in raw.split(b"\x00") if x]


def verify_repo(root: Path) -> str:
    if Path(git(root, "rev-parse", "--show-toplevel")).resolve() != root.resolve():
        raise RuntimeError("--root must point to the Git working-tree root")
    return str(git(root, "rev-parse", "HEAD"))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def file_role(path: str) -> str:
    parts = path.split("/")
    top = parts[0]
    if top == "scripts":
        if len(parts) > 2 and parts[1] == "maintenance":
            return "REPOSITORY_MAINTENANCE"
        if len(parts) > 2 and parts[1] == "research":
            return "EXPERIMENTAL_RESEARCH_SCRIPT"
        if path.endswith(".experimental.py"):
            return "EXPERIMENTAL_RESEARCH_SCRIPT"
        return "RESEARCH_SCRIPT"
    if top == "tools":
        return "REPORT_OR_MAINTENANCE_TOOL"
    if top in {"Lagrangian", "RGE", "Numerical", "physics", "common", "model", "studies"}:
        return "SCIENTIFIC_IMPLEMENTATION"
    if top in {"tests", "validation"}:
        return "VERIFICATION"
    if top in {"configs", "input", "benchmarks"}:
        return "SCIENTIFIC_INPUT"
    if top in {"output", "Reports"}:
        if path.startswith(("output/", "Reports/output/")):
            return "GENERATED_OR_SAVED_RESULT"
        return "REPORTING_IMPLEMENTATION"
    if top == "docs":
        return "RESEARCH_DOCUMENT" if path.startswith("docs/research/") else "TECHNICAL_DOCUMENT"
    if top == "report_cleanup_archive":
        return "HISTORICAL_REPORT_ARCHIVE"
    return "ROOT_OR_MISCELLANEOUS"


def script_topic(path: str) -> str:
    name = Path(path).name.lower()
    if path.startswith("tools/"):
        if path in ONE_TIME_TOOL_REVIEW:
            return "one-time repository migration or retrospective audit"
        if name.startswith("rebuild") or "plot" in name or "figure" in name or "sensitivity" in name:
            return "saved-result reporting and figures"
        if "dependencies" in name or "audit" in name:
            return "dependency and reproducibility review"
        if "restore" in name:
            return "historical scientific report reconstruction"
        return "study/repository utility"
    if "/maintenance/" in path:
        return "repository maintenance and archival"
    if name.startswith(("test", "check", "verify")):
        return "validation or diagnostic executable"
    if name.startswith(("plot", "build", "summarize", "summarise", "report")):
        return "reports, derived results or figure generation"
    if name.startswith(("audit", "probe", "inspect", "review", "diagnose")):
        return "physics investigation and source audit"
    if name.startswith(("export", "generate")):
        return "symbolic or numerical source generation"
    if name.startswith(("fit", "search", "refit", "sweep", "continue")):
        return "parameter fitting and numerical exploration"
    if name.startswith(("run", "evaluate", "integrate", "compute")):
        return "study execution and physics calculation"
    return "other research driver or helper"


def description(path: Path) -> str:
    if path.suffix.lower() == ".py":
        try:
            txt = path.read_text(encoding="utf-8-sig")
            module = ast.parse(txt)
            doc = ast.get_docstring(module) or ""
            if doc:
                return " ".join(doc.split())[:240]
        except (UnicodeError, SyntaxError, OSError):
            pass
    try:
        txt = path.read_text(encoding="utf-8-sig", errors="replace")[:1600]
    except OSError:
        return ""
    for line in txt.splitlines()[:12]:
        s = line.strip(" #(*);\t")
        if len(s) > 20 and not s.startswith(("from ", "import ", "<", "#!")):
            return s[:240]
    return ""


def targets_from_inventory(files: list[str]) -> list[str]:
    out = []
    for p in files:
        if p == ROOT_JUNK:
            out.append(p)
        elif p.startswith(("tools/", "scripts/")) and Path(p).suffix.lower() in SOURCE_SUFFIXES:
            if Path(p).name != "__init__.py":
                out.append(p)
    return sorted(out)


def incoming_class(source: str) -> str:
    if source.startswith("tests/"):
        return "test"
    if source.startswith(("output/", "Reports/output/")):
        return "saved"
    if source.startswith(("configs/", "input/")):
        return "config"
    if source.startswith("docs/") or source in {"README.md", "T3_COMMAND_REFERENCE.md"}:
        return "documentation"
    if Path(source).suffix.lower() in SOURCE_SUFFIXES:
        return "executable"
    if Path(source).suffix.lower() in {".json", ".yaml", ".yml", ".toml", ".ini", ".cfg"}:
        return "config"
    return "documentation"


def reference_scan(root: Path, tracked: list[str], extra: list[str], targets: list[str],
                   max_text_bytes: int = 8 * 1024 * 1024) -> tuple[dict, list, list]:
    # A filename, module name, or explicit project-relative path is treated as
    # an incoming mention; these are conservative witnesses, not execution proofs.
    tokens: dict[str, set[str]] = defaultdict(set)
    for p in targets:
        if p == ROOT_JUNK:
            tokens[p].add(p)
            continue
        name = Path(p).name
        tokens[name].add(p)
        tokens[Path(p).stem].add(p)
        tokens[p].add(p)
        if p.startswith("scripts/") and len(PurePosixPath(p).parts) == 2:
            tokens["scripts." + Path(p).stem].add(p)
        if p.startswith("tools/"):
            tokens["tools." + Path(p).stem].add(p)
    words = sorted(tokens, key=lambda x: (-len(x), x))
    pattern = re.compile(r"(?<![A-Za-z0-9_])(?:" + "|".join(re.escape(x) for x in words) + r")(?![A-Za-z0-9_])")
    counts: dict[str, Counter] = defaultdict(Counter)
    witnesses: list[dict] = []
    seen: set[tuple[str, str]] = set()
    issues: list[dict] = []
    source_list = sorted(set(tracked) | set(extra))
    for rel in source_list:
        if rel in {SELF, SELF_TEST} or rel in IGNORED_INDEXES:
            continue
        if rel.startswith(("output/maintenance/", "output/cleanup_audit/")):
            continue
        path = checked_path(root, rel)
        if not path.is_file() or path.is_symlink() or path.suffix.lower() not in TEXT_SUFFIXES:
            continue
        try:
            size = path.stat().st_size
            if size > max_text_bytes:
                issues.append({"path": rel, "problem": "oversized_text", "bytes": size})
                continue
            with path.open("r", encoding="utf-8-sig", errors="replace") as stream:
                for lineno, line in enumerate(stream, 1):
                    # Prevent repeated mention of the same target from one file
                    # from overwhelming the audit, without hiding that reference.
                    for match in pattern.finditer(line):
                        term = match.group(0)
                        for target in tokens[term]:
                            if target == rel or (target == ROOT_JUNK and rel in PAGER_CLEANUP_MENTIONS):
                                continue
                            key = (target, rel)
                            if key in seen:
                                continue
                            seen.add(key)
                            kind = incoming_class(rel)
                            counts[target][kind] += 1
                            witnesses.append({"target": target, "source": rel,
                                              "kind": kind, "line": lineno,
                                              "excerpt": line.strip()[:160]})
        except (OSError, UnicodeError) as exc:
            issues.append({"path": rel, "problem": type(exc).__name__, "detail": str(exc)[:160]})
    return counts, witnesses, issues


def is_git_blob_unchanged(root: Path, path: str) -> bool:
    try:
        expected = git(root, "rev-parse", "HEAD:" + path)
        actual = git(root, "hash-object", "--path=" + path, "--", path)
        return expected == actual
    except RuntimeError:
        return False


def evaluate(root: Path) -> dict:
    head = verify_repo(root)
    tracked = sorted(set(git_names(root, "ls-files", "-z")))
    modified = set(git_names(root, "diff", "--name-only", "-z", "HEAD", "--"))
    deleted = set(git_names(root, "ls-files", "-d", "-z"))
    untracked = sorted(set(git_names(root, "ls-files", "--others", "--exclude-standard", "-z")))
    targets = targets_from_inventory(tracked)
    counts, witnesses, coverage = reference_scan(root, tracked, untracked, targets)
    # If the scan is incomplete, no file gets an automatic archive recommendation.
    coverage_complete = not coverage
    bypath = []
    for p in tracked:
        f = checked_path(root, p)
        size = f.stat().st_size if f.is_file() else 0
        bypath.append({"path": p, "role": file_role(p), "bytes": size,
                       "modified": p in modified or p in deleted or not f.is_file(),
                       "source_kind": "tracked"})
    for p in untracked:
        if p.startswith(("output/maintenance/", "output/cleanup_audit/")):
            continue
        f = checked_path(root, p)
        bypath.append({"path": p, "role": file_role(p), "bytes": f.stat().st_size if f.is_file() else 0,
                       "modified": True, "source_kind": "untracked"})
    bypath.sort(key=lambda x: x["path"])
    rows = []
    for p in targets:
        f = checked_path(root, p)
        c = counts[p]
        clean = p not in modified and p not in deleted and f.is_file() and not f.is_symlink()
        role = file_role(p)
        locked = p in PROTECTED
        hard = c["executable"] + c["test"] + c["config"]
        saved = c["saved"]
        docs = c["documentation"]
        # This file is a verified non-project less(1) command-help capture.
        junk = False
        if p == ROOT_JUNK and f.is_file():
            try:
                txt = f.read_text(encoding="utf-8", errors="replace")
                junk = "Commands marked with * may be preceded" in txt and "Exit." in txt
            except OSError:
                pass
        if locked:
            status, reason = "KEEP", "active/restore/dependency infrastructure protected"
        elif not clean:
            status, reason = "HOLD_LOCAL_CHANGE", "changed, missing or symlinked file"
        elif hard:
            status, reason = "KEEP_REFERENCED", "mentioned by executable, tests or input/config"
        elif saved:
            status, reason = "KEEP_PROVENANCE", "mentioned by saved numerical/result artifacts"
        elif not coverage_complete:
            status, reason = "HOLD_INCOMPLETE_SCAN", "reference scan has unreadable/oversized text"
        elif junk:
            status, reason = "SUGGEST_ARCHIVE", "verified accidental terminal pager capture"
        elif p in ONE_TIME_TOOL_REVIEW or p in ONE_TIME_SCRIPT_REVIEW:
            status, reason = "SUGGEST_ARCHIVE", "inspected one-time install, rename, report retirement or migration tool"
        elif docs:
            status, reason = "REVIEW_DOCUMENTED", "historical/current documentation mention; investigate external use"
        else:
            status, reason = "REVIEW_STANDALONE", "no scanned incoming mention; external invocation remains possible"
        rows.append({"path": p, "role": role, "topic": script_topic(p),
                     "description": description(f) if f.is_file() else "",
                     "bytes": f.stat().st_size if f.is_file() else 0,
                     "sha256": sha256(f) if f.is_file() and not f.is_symlink() else None,
                     "clean": clean, "executable_refs": c["executable"],
                     "test_refs": c["test"], "config_refs": c["config"],
                     "saved_refs": saved, "doc_refs": docs,
                     "status": status, "reason": reason})
    suggested = [r["path"] for r in rows if r["status"] == "SUGGEST_ARCHIVE" and r["doc_refs"] == 0]
    return {"format": FORMAT, "head": head,
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "tracked_files": len(tracked), "untracked_files": len(untracked),
            "audited_targets": len(rows), "coverage_complete": coverage_complete,
            "coverage_issues": coverage, "files": bypath,
            "rows": rows, "witnesses": witnesses,
            "suggested_selection": suggested,
            "caution": "Static references cannot prove absence of dynamically generated calls or external user commands."}


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        w = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def report_directory(root: Path) -> Path:
    return root.parent / "T3_external_archives" / REPORT_DIRNAME


def emit(root: Path, report: dict) -> Path:
    out = report_directory(root)
    out.mkdir(parents=True, exist_ok=True)
    doc = {k: v for k, v in report.items() if k not in ("files", "witnesses")}
    (out / "plan.json").write_text(json.dumps(doc, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    write_csv(out / "all_files.csv", report["files"],
              ["path", "role", "bytes", "modified", "source_kind"])
    fields = ["path", "role", "topic", "description", "bytes", "sha256", "clean",
              "executable_refs", "test_refs", "config_refs", "saved_refs", "doc_refs", "status", "reason"]
    write_csv(out / "scripts_and_tools.csv", report["rows"], fields)
    write_csv(out / "reference_witnesses.csv", report["witnesses"],
              ["target", "source", "kind", "line", "excerpt"])
    selection = out / "reviewed_selection.txt"
    # Never overwrite a user's manually curated selection on later previews.
    if not selection.exists():
        header = ["# Exact relative paths approved for external archival.",
                  "# An empty selection archives nothing. Do not uncomment a research script without checking its scientific role.",
                  "# Automatic low-risk suggestions from initial local audit follow:"]
        selection.write_text("\n".join(header + report["suggested_selection"]) + "\n", encoding="utf-8")
    groups = Counter(x["role"] for x in report["files"] if x["source_kind"] == "tracked")
    statuses = Counter(x["status"] for x in report["rows"])
    lines = ["# T3 thesis repository: file-level inventory", "",
             f"- Git HEAD: `{report['head']}`",
             f"- Tracked paths inventoried: **{report['tracked_files']}**",
             f"- Untracked paths inventoried: **{report['untracked_files']}**",
             f"- Scripts/tools individually classified: **{report['audited_targets']}**",
             f"- Full text reference coverage: **{'complete' if report['coverage_complete'] else 'INCOMPLETE — no automatic approval'}**", "",
             "## Tracked-file roles", "", "| Role | Files |", "|---|---:|"]
    lines.extend(f"| {k} | {v} |" for k, v in sorted(groups.items()))
    lines.extend(["", "## Source review statuses", "", "| Status | Files |", "|---|---:|"])
    lines.extend(f"| {k} | {v} |" for k, v in sorted(statuses.items()))
    lines += ["", "## Suggested candidates (not automatically removed)", ""]
    for r in report["rows"]:
        if r["status"] == "SUGGEST_ARCHIVE":
            lines.append(f"- `{r['path']}` — {r['reason']}; documentation references: {r['doc_refs']}")
    if not any(r["status"] == "SUGGEST_ARCHIVE" for r in report["rows"]):
        lines.append("No file cleared the strict suggested-candidate checks.")
    lines += ["", "**Important:** Empty static references are not evidence that a research program has no thesis value.",
              "The archive only acts on exact filenames in `reviewed_selection.txt` after explicit confirmation.",
              "Read `scripts_and_tools.csv` and `reference_witnesses.csv` before selecting anything."]
    if report["coverage_issues"]:
        lines += ["", "## Coverage limitations", ""]
        lines.extend(f"- `{x['path']}`: {x['problem']}" for x in report["coverage_issues"][:80])
    (out / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return out


def load_selection(path: Path) -> list[str]:
    if not path.is_file():
        raise ValueError(f"Selection file not found: {path}")
    selected = []
    for line in path.read_text(encoding="utf-8-sig").splitlines():
        name = line.strip()
        if not name or name.startswith("#"):
            continue
        if name in selected:
            raise ValueError(f"Duplicate selection entry: {name}")
        selected.append(name)
    return selected


def archive(root: Path, old: dict, selected: list[str], *, allow_research: bool,
            allow_documented: bool) -> Path | None:
    if not selected:
        print("No files selected; no changes made.")
        return None
    fresh = evaluate(root)
    if fresh["head"] != old["head"]:
        raise RuntimeError("Git HEAD changed since preview; rerun the audit")
    if not fresh["coverage_complete"]:
        raise RuntimeError("Reference scan incomplete; nothing may be archived")
    planned = {r["path"]: r for r in old["rows"]}
    current = {r["path"]: r for r in fresh["rows"]}
    eligible = {"SUGGEST_ARCHIVE", "REVIEW_STANDALONE", "REVIEW_DOCUMENTED"}
    for p in selected:
        checked_path(root, p)
        if p not in planned or p not in current:
            raise RuntimeError(f"Selected file not in reviewed inventory: {p}")
        prev, now = planned[p], current[p]
        if prev["status"] not in eligible or now["status"] not in eligible:
            raise RuntimeError(f"Selected file is not archive-eligible: {p} ({now['status']})")
        if prev["sha256"] != now["sha256"] or not now["clean"]:
            raise RuntimeError(f"Selected source changed since the reviewed plan: {p}")
        for key in ("executable_refs", "test_refs", "config_refs", "saved_refs", "doc_refs"):
            if prev[key] != now[key]:
                raise RuntimeError(f"Incoming references changed for {p}; rerun review")
        if now["doc_refs"] and not allow_documented:
            raise RuntimeError(f"Documented path {p}: use --allow-documented only after reviewing references")
        if (p.startswith("scripts/") and not p.startswith("scripts/maintenance/")
                and p not in ONE_TIME_SCRIPT_REVIEW and not allow_research):
            raise RuntimeError(f"Research script {p} needs --allow-research-scripts and explicit review")
        if not is_git_blob_unchanged(root, p):
            raise RuntimeError(f"Git blob differs from HEAD: {p}")
    external = root.parent / "T3_external_archives"
    external.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    target = external / ("thesis_review_" + stamp)
    if target.exists():
        raise RuntimeError("Archive directory already exists")
    target.mkdir(parents=False)
    rows = []
    try:
        for p in selected:
            source = checked_path(root, p)
            dest = target / "files" / p
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, dest)
            expected = current[p]["sha256"]
            if sha256(dest) != expected:
                raise RuntimeError(f"Archive copy checksum mismatch: {p}")
            rows.append({"path": p, "sha256": expected, "bytes": dest.stat().st_size,
                         "status_before": current[p]["status"]})
        manifest = {"format": FORMAT + "-archive", "git_head": fresh["head"],
                    "root": str(root.resolve()), "created_utc": datetime.now(timezone.utc).isoformat(),
                    "entries": rows}
        (target / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        # All source removals happen after every file is copied and checked.
        for row in rows:
            source = checked_path(root, row["path"])
            if sha256(source) != row["sha256"]:
                raise RuntimeError("Original changed during copy; archive retained: " + row["path"])
            source.unlink()
    except Exception:
        # Preserve the archive as soon as copies exist, even after partial unlink.
        raise
    print(f"Archived {len(rows)} reviewed files to: {target}")
    print("Run the core tests, and inspect Git status before committing locally.")
    return target


def restore(root: Path, folder: Path) -> None:
    manifest_file = folder / "manifest.json"
    if not manifest_file.is_file():
        raise RuntimeError("Missing archive manifest")
    manifest = json.loads(manifest_file.read_text(encoding="utf-8"))
    if manifest.get("format") != FORMAT + "-archive":
        raise RuntimeError("Unsupported archive manifest format")
    if Path(manifest.get("root", "")).resolve() != root.resolve():
        raise RuntimeError("Archive belongs to another checkout")
    entries = manifest.get("entries")
    if not isinstance(entries, list) or not entries:
        raise RuntimeError("Archive contains no entries")
    seen = set()
    # Preflight the ENTIRE restore: fail before copying anything on a conflict.
    for row in entries:
        p = row["path"]
        if p in seen:
            raise RuntimeError("Repeated entry in archive: " + p)
        seen.add(p)
        dest = checked_path(root, p)
        source = checked_path(folder / "files", p)
        if not source.is_file() or sha256(source) != row["sha256"]:
            raise RuntimeError("Archive checksum invalid: " + p)
        if dest.exists() or dest.is_symlink():
            raise RuntimeError("Refusing to overwrite existing path: " + p)
    for row in entries:
        p = row["path"]
        dest = checked_path(root, p)
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(checked_path(folder / "files", p), dest)
        if sha256(dest) != row["sha256"]:
            raise RuntimeError("Restored checksum mismatch: " + p)
    print(f"Restored {len(entries)} files from: {folder}")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    ap.add_argument("--apply", action="store_true", help="Archive explicitly selected files")
    ap.add_argument("--selection", type=Path, help="Text file containing exact project-relative filenames")
    ap.add_argument("--confirm", default="", help="Required exact archival confirmation token")
    ap.add_argument("--allow-research-scripts", action="store_true",
                    help="Allow individually selected unreferenced research scripts")
    ap.add_argument("--allow-documented", action="store_true",
                    help="Allow documented-only paths after checking affected links")
    ap.add_argument("--restore", type=Path, help="Restore external archive, refusing overwrites")
    args = ap.parse_args(argv)
    root = args.root.resolve()
    if args.restore is not None:
        if args.apply:
            ap.error("Cannot apply and restore in one invocation")
        verify_repo(root)
        restore(root, args.restore.resolve())
        return 0
    if args.apply:
        if args.confirm != CONFIRM:
            ap.error("--apply requires --confirm " + CONFIRM)
        planned = report_directory(root) / "plan.json"
        if not planned.is_file():
            ap.error("Run the default preview before applying")
        old = json.loads(planned.read_text(encoding="utf-8"))
        if old.get("format") != FORMAT:
            ap.error("Unsupported/incomplete review plan")
        selection = args.selection or (report_directory(root) / "reviewed_selection.txt")
        if not selection.is_absolute():
            selection = root / selection
        selected = load_selection(selection)
        archive(root, old, selected, allow_research=args.allow_research_scripts,
                allow_documented=args.allow_documented)
        return 0
    result = evaluate(root)
    out = emit(root, result)
    print("Git HEAD:", result["head"])
    print("Tracked files:", result["tracked_files"])
    print("Research scripts and tools reviewed:", result["audited_targets"])
    print("Reference coverage issues:", len(result["coverage_issues"]))
    for status, count in sorted(Counter(r["status"] for r in result["rows"]).items()):
        print(f"  {status}: {count}")
    print("READ-ONLY preview; no scientific files moved.")
    print("Reports:", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
