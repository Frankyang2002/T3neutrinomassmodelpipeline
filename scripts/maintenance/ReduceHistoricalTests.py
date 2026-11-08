#!/usr/bin/env python3
"""Reversible, Git-pinned retirement of narrow historical pytest modules.

Does not modify physics code, tests selected for retention, or GitHub.
Run without flags for read-only planning, then review the generated JSON.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
from typing import Any

CONFIRM = "ARCHIVE_HISTORICAL_TESTS"
ARCHIVE_SUBDIR = "T3_external_archives"

# These tests are central to the general UV -> all fermions -> all scalars -> SMEFT
# numerical programme, or safeguard core conventions. Pattern guards are kept
# deliberately broader than individual names so no such test is retired by
# a coincidental name match.
PROTECT = tuple(re.compile(x) for x in (
    r"test_pipeline_plan\.py$",
    r"test_scalar_first_safety\.py$",
    r"test_eft_transitions\.py$",
    r"test_t3_field_content\.py$",
    r"test_full16_general_symbolic_physics_audit\.py$",
    r"test_general_coupling_(inventory|layout)\.py$",
    r"test_general_t3_fermion_first_eft1_(tree|one_loop)\.py$",
    r"test_general_d_m2_",
    r"test_general_t3_",
    r"test_batch34_",  # contains physical covariance and gauge tests
    r"test_t3_mixed_batch",
    r"test_t3_scalar_(vertex|relevance)_batch",
    r"test_t3_quadratic_(impact|rg_closure)_batch",
    r"test_lo_completeness_",
    r"test_lo_couplings_",
    r"test_lo_direct_uv_",
    r"test_lo_spinor_",
    r"test_t3_complete_general_scalar_symbolic\.py$",
    r"test_numerical_state\.py$",
    r"test_parameter_scan\.py$",
    r"test_uv_runner\.py$",
    r"test_pipeline_backbone\.py$",
    r"test_pipeline_cli\.py$",
    r"test_low_energy_neutrino_architecture\.py$",
    r"test_weinberg_flavor_matching_architecture\.py$",
    r"test_final_c5_bridge\.py$",
    r"test_sobol_benchmark_search\.py$",
    r"test_thesis_self_mass_basis\.py$",
    r"test_one_generation_weinberg_benchmark_architecture\.py$",
))

# Name-based triage, NOT a claim that each test is redundant or scientifically
# wrong. The archival operation is reversible and intentionally narrow.
RETIRE_RULES = tuple((label, re.compile(pattern)) for label, pattern in (
    ("missing-pdf-writer", r"test_(general_latex_report_layout|general_report_structure_tables|t3_general_pdf_fallback|t3_general_pdf_functionality)\.py$"),
    ("presentation-regression", r"test_.*(pdf|report|figure|display|plot|html|latex|visual|layout).*\.py$"),
    ("historical-benchmark", r"test_.*(frozen|zero_extra|zero_boundary|reduced_scalar|uniform|sobol|thesis|qualitative|refit|benchmark|overnight).*\.py$"),
    ("batch-or-pilot-regression", r"test_.*(batch\d|batch_\d|batch_v\d|stage[234]|continuation|pilot|prototype|repair|closeout).*\.py$"),
    ("historical-docs-test", r"test_.*(historical_docs_consolidation|descriptive_research_file_names|legacy_notes_cleanup|one_time_artifact_cleanup|analytical_housekeeping_auto).*\.py$"),
    ("lo-analysis-snapshot", r"test_lo_.*(pattern|synthesis|factorization|texture|splitting|preliminary|analytic_explanation).*\.py$"),
    ("historical-exploration", r"test_.*(best_fit_diagnostics|running_diagnostics|t3b_exploratory|general_t3_zee_sensitivity|t3_general_scalar_diagnostics|general_t3_scalar_quartic_response_audit).*\.py$"),
))

CODE_SOURCE_SUFFIXES = frozenset({".py", ".ps1", ".wl", ".toml", ".ini", ".cfg", ".yml", ".yaml"})


def run_git(root: Path, *args: str) -> bytes:
    result = subprocess.run(["git", *args], cwd=root, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"git {' '.join(args)} failed: {result.stderr.decode(errors='replace').strip()}")
    return result.stdout


def checked_repository(root: Path, expected_head: str | None = None) -> str:
    root = root.resolve()
    git_root = Path(run_git(root, "rev-parse", "--show-toplevel").decode().strip()).resolve()
    if git_root != root:
        raise RuntimeError(f"Expected repository at {root}; Git root is {git_root}")
    head = run_git(root, "rev-parse", "HEAD").decode().strip()
    if expected_head is not None and head != expected_head:
        raise RuntimeError(f"Git HEAD changed ({head}); requires fresh review of test inventory")
    return head


def tracked_head_entries(root: Path) -> dict[str, tuple[str, str]]:
    result: dict[str, tuple[str, str]] = {}
    for raw in run_git(root, "ls-tree", "-r", "-z", "HEAD", "--", "tests").split(b"\x00"):
        if not raw:
            continue
        meta, path = raw.split(b"\t", 1)
        mode, obj_type, sha = meta.decode().split(" ")
        if obj_type == "blob":
            result[path.decode("utf-8", errors="surrogateescape")] = (mode, sha)
    return result


def test_group(path: str) -> str | None:
    if not re.fullmatch(r"tests/(?:python/)?test_[^/]+\.py", path):
        return None
    if any(pattern.search(path) for pattern in PROTECT):
        return None
    for name, pattern in RETIRE_RULES:
        if pattern.search(path):
            return name
    return None


def git_clean_blob_sha1(root: Path, rel: str) -> str:
    # `hash-object --path` normally applies the same clean filters as Git.
    return run_git(root, "hash-object", f"--path={rel}", rel).decode().strip()


def git_sha1_bytes(root: Path, data: bytes) -> str:
    process = subprocess.run(["git", "hash-object", "--stdin"], cwd=root,
                             input=data, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, check=False)
    if process.returncode:
        raise RuntimeError(process.stderr.decode(errors="replace"))
    return process.stdout.decode().strip()


def same_as_head_blob(root: Path, rel: str, head_blob: str) -> bool:
    if git_clean_blob_sha1(root, rel) == head_blob:
        return True
    # Some Windows Git configurations do not apply the expected clean filter
    # to hash-object --path if core.autocrlf was changed after checkout.
    # Accept *only* an exact LF-normalised match with HEAD, never arbitrary
    # content changes. This is valid for Git-managed text with CRLF conversion.
    attributes = run_git(root, "check-attr", "text", "eol", "--", rel).decode(
        "utf-8", errors="replace"
    )
    if "text: unset" in attributes:
        return False  # explicitly binary via -text
    config = subprocess.run(["git", "config", "--get", "core.autocrlf"],
                            cwd=root, capture_output=True, text=True, check=False)
    autocrlf = config.stdout.strip().lower() if config.returncode == 0 else ""
    text_managed = (autocrlf in ("true", "input") or
                    "text: set" in attributes or
                    "text: auto" in attributes or
                    "eol: lf" in attributes or
                    "eol: crlf" in attributes)
    if not text_managed:
        return False
    data = (root / rel).read_bytes()
    if b"\r\n" not in data or b"\x00" in data:
        return False
    return git_sha1_bytes(root, data.replace(b"\r\n", b"\n")) == head_blob


def safe_repo_path(root: Path, rel: str) -> Path:
    if not rel.startswith("tests/") or "\\" in rel or "\x00" in rel or ".." in Path(rel).parts:
        raise ValueError(f"Invalid test path: {rel!r}")
    path = root / rel
    if path.is_symlink() or root.resolve() not in path.resolve().parents:
        raise ValueError(f"Unsafe test path: {rel!r}")
    return path


def eligible_sources(root: Path, selected: set[str]) -> list[Path]:
    # Scan maintained code/config and retained test modules, not historical
    # markdown prose (a reference in an old article is not an execution edge).
    names = run_git(root, "ls-files", "-z").split(b"\x00")
    result = []
    for raw in names:
        if not raw:
            continue
        rel = raw.decode("utf-8", errors="surrogateescape")
        if rel in selected or Path(rel).suffix.lower() not in CODE_SOURCE_SUFFIXES:
            continue
        file_path = root / rel
        if file_path.is_file() and not file_path.is_symlink():
            result.append(file_path)
    # Newly installed retained files may not be tracked yet; scan only
    # the current maintenance/test-runner files that can reference tests.
    for rel in ("scripts/RunCoreTests.ps1", "scripts/maintenance/ReduceHistoricalTests.py"):
        p = root / rel
        if p.is_file() and p not in result:
            result.append(p)
    return result


def build_plan(root: Path, expected_head: str | None = None) -> dict[str, Any]:
    root = root.resolve()
    head = checked_repository(root, expected_head)
    entries = tracked_head_entries(root)
    preliminary = {p: test_group(p) for p in entries}
    prelim = {p: reason for p, reason in preliminary.items() if reason}
    candidate_names = set(prelim)
    # To avoid self-blocking on the test cleanup tool, exclude its source:
    # its category expressions do not contain entire old filenames.
    sources = eligible_sources(root, candidate_names)
    references: dict[str, list[str]] = {p: [] for p in candidate_names}
    aliases: dict[str, list[str]] = {}
    for p in candidate_names:
        aliases.setdefault(Path(p).stem, []).append(p)
    for file_path in sources:
        try:
            data = file_path.read_text(encoding="utf-8-sig", errors="replace")
        except OSError:
            continue
        for stem, paths in aliases.items():
            # Detect both literal paths and Python module imports such as
            # `from tests.test_old_campaign import ...` (no `.py` suffix).
            if re.search(r"(?<![\w])" + re.escape(stem) + r"(?![\w])", data):
                src = file_path.relative_to(root).as_posix()
                for p in paths:
                    references[p].append(src)
    approved: list[dict[str, Any]] = []
    held: list[dict[str, Any]] = []
    missing: list[str] = []
    for rel, reason in sorted(prelim.items()):
        path = safe_repo_path(root, rel)
        if not path.is_file():
            missing.append(rel)
            continue
        content = path.read_bytes()
        mode, head_blob = entries[rel]
        row = {"path": rel, "reason": reason, "bytes": len(content),
               "sha256": hashlib.sha256(content).hexdigest(),
               "head_blob_sha1": head_blob}
        issues = []
        if mode not in ("100644", "100755"):
            issues.append(f"unsupported Git mode {mode}")
        if not same_as_head_blob(root, rel, head_blob):
            issues.append("different from Git HEAD (modified or staged content)")
        if references[rel]:
            issues.append("explicitly referenced by retained source")
            row["references"] = references[rel]
        if issues:
            row["hold_reason"] = "; ".join(issues)
            held.append(row)
        else:
            approved.append(row)
    total = sum(bool(re.fullmatch(r"tests/(?:python/)?test_[^/]+\.py", p)) for p in entries)
    return {
        "schema": "t3-curated-test-retirement-v1", "git_head": head,
        "repository": str(root), "original_python_test_modules": total,
        "selected_by_name": len(prelim), "approved": approved, "held": held,
        "already_missing": missing, "still_retained": total - len(approved),
        "approved_bytes": sum(a["bytes"] for a in approved),
        "note": "Names are a triage heuristic; archive is reversible. The remaining tests still require actual execution.",
    }


def plan_path(root: Path) -> Path:
    return root / "output/cleanup_audit/historical_test_retirement_plan.json"


def load_reviewed_plan(root: Path) -> dict[str, Any]:
    path = plan_path(root)
    if not path.is_file():
        raise RuntimeError("No reviewed plan exists. Run without --apply, inspect the JSON, then apply.")
    plan = json.loads(path.read_text(encoding="utf-8"))
    if plan.get("schema") != "t3-curated-test-retirement-v1":
        raise RuntimeError("Unexpected test retirement plan schema")
    if Path(plan.get("repository", "")).resolve() != root.resolve():
        raise RuntimeError("Test retirement plan belongs to another repository")
    if not isinstance(plan.get("approved"), list) or not isinstance(plan.get("held"), list):
        raise RuntimeError("Invalid or incomplete reviewed plan")
    return plan


def save_plan(root: Path, plan: dict[str, Any]) -> Path:
    dest = plan_path(root)
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(plan, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return dest


def apply_plan(root: Path, plan: dict[str, Any]) -> Path | None:
    # Recheck the entire plan, not merely files, immediately before applying.
    root = root.resolve()
    fresh = build_plan(root, plan["git_head"])
    compared = ("approved", "held", "already_missing", "selected_by_name",
                "original_python_test_modules", "still_retained")
    if any(fresh.get(key) != plan.get(key) for key in compared):
        raise RuntimeError("Test inventory or dependencies changed after preview; rerun dry run")
    items = fresh["approved"]
    if not items:
        return None
    external = root.parent / ARCHIVE_SUBDIR
    external.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("test_cleanup_%Y%m%dT%H%M%S%fZ")
    archive = external / stamp
    archive.mkdir(parents=False, exist_ok=False)
    manifest = {"schema": "t3-test-retirement-archive-v1",
                "source_repository": str(root), "git_head": fresh["git_head"],
                "files": items}
    try:
        # Copy and verify ALL files before removing any from the checkout.
        for item in items:
            source = safe_repo_path(root, item["path"])
            dest = archive / item["path"]
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, dest)
            if hashlib.sha256(dest.read_bytes()).hexdigest() != item["sha256"]:
                raise RuntimeError(f"Archive hash mismatch: {item['path']}")
        (archive / "manifest.json").write_text(
            json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        for item in items:
            source = safe_repo_path(root, item["path"])
            if hashlib.sha256(source.read_bytes()).hexdigest() != item["sha256"]:
                raise RuntimeError("Source changed during archive copy; no further files will be deleted")
            source.unlink()
    except Exception:
        # Leave the external archive (including its checked files) available
        # for recovery if a failure happens during the move stage.
        raise
    return archive


def restore_archive(root: Path, archive: Path) -> int:
    root = root.resolve()
    checked_repository(root, expected_head=None)
    archive = archive.resolve()
    # Never take an arbitrary path under repository as the archive.
    expected_parent = (root.parent / ARCHIVE_SUBDIR).resolve()
    if archive.parent != expected_parent or archive == root:
        raise ValueError(f"Archive must be a direct child of {expected_parent}")
    manifest = json.loads((archive / "manifest.json").read_text(encoding="utf-8"))
    if manifest.get("schema") != "t3-test-retirement-archive-v1":
        raise RuntimeError("Unrecognized archive manifest")
    if Path(manifest["source_repository"]).resolve() != root:
        raise RuntimeError("Archive belongs to a different checkout")
    for item in manifest["files"]:
        target = safe_repo_path(root, item["path"])
        if target.exists():
            raise FileExistsError(f"Restore would overwrite existing file: {target}")
        origin = archive / item["path"]
        if origin.is_symlink() or not origin.is_file():
            raise RuntimeError(f"Missing or unsafe archived file: {origin}")
        if hashlib.sha256(origin.read_bytes()).hexdigest() != item["sha256"]:
            raise RuntimeError(f"Archive content changed: {origin}")
    for item in manifest["files"]:
        target = safe_repo_path(root, item["path"])
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(archive / item["path"], target)
    return len(manifest["files"])


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--apply", action="store_true", help="Move reviewed tests to external verified archive")
    parser.add_argument("--confirm", default="", help=f"Required for --apply: {CONFIRM}")
    parser.add_argument("--restore", type=Path, help="Restore a verified test archive into original paths")
    args = parser.parse_args(argv)
    root = args.repo.resolve()
    if args.restore:
        if args.apply:
            parser.error("--restore and --apply are mutually exclusive")
        print(f"Restored {restore_archive(root, args.restore)} tests from {args.restore}")
        return 0
    if args.apply:
        if args.confirm != CONFIRM:
            parser.error(f"--apply requires --confirm {CONFIRM}")
        # Apply the plan previously written and inspected by the user.
        # Never implicitly approve a changed HEAD or recomputed candidate list.
        plan = load_reviewed_plan(root)
        where = plan_path(root)
        archive = apply_plan(root, plan)
        print(f"Reviewed plan: {where}")
        print(f"Git HEAD: {plan['git_head']}")
        print(f"Approved for archive: {len(plan['approved'])}")
    else:
        plan = build_plan(root, expected_head=None)
        where = save_plan(root, plan)
        print(f"Git HEAD: {plan['git_head']}")
        print(f"Python test modules at HEAD: {plan['original_python_test_modules']}")
        print(f"Selected by rules: {plan['selected_by_name']}")
        print(f"Approved for archive: {len(plan['approved'])}")
        print(f"Held for safety: {len(plan['held'])}; already missing: {len(plan['already_missing'])}")
        print(f"Remaining tracked test modules after archive: {plan['still_retained']}")
        print(f"Plan: {where}")
        print("Preview only. No files changed.")
        return 0
    print(f"Archived to: {archive}" if archive else "No eligible tests to archive")
    print("Run the core test suite and then full pytest. No GitHub changes were made.")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValueError, RuntimeError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)
