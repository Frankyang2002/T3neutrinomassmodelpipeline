"""Local, reversible audit/archive of obsolete documentation and disposable artifacts.

No GitHub writes and no physics-code edits. Python 3.10+.

Default: create output/cleanup_audit/extreme_repo_plan.json only.
Apply: archive candidates outside checkout (hash-verified) before removing them.
Restore: verify archive hashes, refuse to overwrite existing paths.

The reference scan is intentionally conservative but cannot prove absence of a
runtime dependency (dynamic paths may not appear as literal strings).
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
from typing import Iterable

ROOT = Path(__file__).resolve().parents[2]
PLAN_PATH = Path("output/cleanup_audit/extreme_repo_plan.json")
KEEP_DOCS = {
    "README.md",
    "docs/MODEL_AND_PIPELINE.md",
    "docs/REPOSITORY_CLEANUP_POLICY.md",
}
SOURCE_SUFFIXES = {".py", ".wl", ".wls", ".m", ".ps1", ".json", ".yaml", ".yml", ".toml", ".md", ".txt", ".tex"}
BINARY_PRESENTATIONS = {".png", ".pdf", ".html", ".svg"}
BACKUP_NAME = re.compile(r"(?i)(?:\.backup[-_]|\.bak$|\.orig$|\.tmp$|\.before_[a-z0-9_]+$|~$)")
NOTE_MARKERS = {
    "restricted_or_reduced_model": re.compile(r"\breduced\b|\brestricted\b", re.I),
    "progress_or_temporary_instruction": re.compile(r"\b(?:workstream|milestone|stage [0-9]+|next step|currently|status update|batch [0-9]+)\b", re.I),
    "strong_generality_claim": re.compile(r"\b(?:fully general|complete general|complete physical|all allowed|certified)\b", re.I),
    "old_install_instruction": re.compile(r"Expand-Archive|(?:\.zip\b)|README_REPLACEMENT", re.I),
}
GENERALITY_MARKERS = {
    "opt_in_general": re.compile(r"T3_EXPERIMENTAL_GENERAL|--general-symbolic", re.I),
    "restricted_route": re.compile(r"--restricted-production|restricted[_ -](?:matching|production|lagrangian|scalar)|reduced[_ -](?:potential|scalar|lagrangian)", re.I),
    "uncertified_basis": re.compile(r"FullSymmetryCompletenessCertified|FermionFlavorMassGeneralityCertified|RequiresMixedScalarMassDiagonalization", re.I),
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def git_output(root: Path, *args: str) -> bytes | None:
    try:
        p = subprocess.run(["git", "-C", str(root), *args], capture_output=True, check=True)
        return p.stdout
    except (FileNotFoundError, subprocess.CalledProcessError):
        return None


def is_git_checkout_root(root: Path) -> bool:
    """Check Git's *actual* worktree root, not just whether an ancestor has .git.

    Git discovery walks upwards. A synthetic temporary directory under some
    other Git working tree is *not* a checkout of its own. In that situation
    `git ls-files` can succeed but return zero entries, hiding every fixture.
    """
    raw = git_output(root, "rev-parse", "--show-toplevel")
    if not raw:
        return False
    try:
        reported = Path(raw.decode("utf-8", errors="surrogateescape").strip()).resolve()
        expected = root.resolve()
    except (OSError, ValueError):
        return False
    return os.path.normcase(str(reported)) == os.path.normcase(str(expected))


def _git_paths(raw: bytes) -> list[str]:
    return sorted(set(s.decode("utf-8", errors="surrogateescape").replace("\\", "/")
                      for s in raw.split(b"\0") if s))


def tracked_paths(root: Path) -> list[str]:
    if is_git_checkout_root(root):
        raw = git_output(root, "ls-files", "-z")
        if raw is None:
            raise RuntimeError("Cannot enumerate Git-tracked files; refusing cleanup plan")
        return _git_paths(raw)
    if (root / ".git").exists():
        # A real checkout whose Git metadata is unreadable must not fall back
        # to scanning untracked files and accidentally offering them for deletion.
        raise RuntimeError("Git metadata exists but checkout root could not be validated")
    # Synthetic/non-git fixture. Git may detect an *ancestor* repository;
    # only the fixture's own files are relevant in this fallback.
    return sorted(p.relative_to(root).as_posix() for p in root.rglob("*")
                  if p.is_file() and ".git" not in p.relative_to(root).parts)


def untracked_reference_paths(root: Path) -> list[str]:
    """Inspect new local code/docs for references, without treating them as deletable."""
    if not is_git_checkout_root(root):
        return []
    raw = git_output(root, "ls-files", "--others", "--exclude-standard", "-z")
    if raw is None:
        raise RuntimeError("Cannot enumerate untracked reference files")
    return _git_paths(raw)


def modified_paths(root: Path) -> set[str]:
    """Return tracked paths with changes. Untracked files are never removed."""
    if not is_git_checkout_root(root):
        return set()
    raw = git_output(root, "status", "--porcelain=v1", "-z", "--untracked-files=no")
    if raw is None:
        raise RuntimeError("Cannot read Git working-tree status; refusing cleanup plan")
    changed: set[str] = set()
    chunks = raw.split(b"\0")
    index = 0
    while index < len(chunks):
        chunk = chunks[index]
        index += 1
        if not chunk or len(chunk) < 4:
            continue
        xy = chunk[:2]
        changed.add(chunk[3:].decode("utf-8", errors="surrogateescape").replace("\\", "/"))
        if b"R" in xy or b"C" in xy:
            if index < len(chunks):
                changed.add(chunks[index].decode("utf-8", errors="surrogateescape").replace("\\", "/"))
                index += 1
    return changed


def category(relative: str, *, nested_notes: bool, presentations: bool) -> str | None:
    name = Path(relative).name
    suffix = Path(relative).suffix.lower()
    if relative in KEEP_DOCS:
        return None
    if "/" not in relative:
        if relative in {"CLEANUP_README.md", "CONSOLIDATION_README.md", "REGENERATE_README.md", "T3_COMMAND_REFERENCE.md"}:
            return "root-historical-document"
        if re.fullmatch(r"README_.*\.md", name, re.I) or re.fullmatch(r"INSTALL.*\.txt", name, re.I):
            return "root-historical-document"
        if name.endswith(".zip"):
            return "root-export-archive"
        if name in {"t -q --tb=line"} or name.startswith("olidate historical"):
            return "root-stray-file"
    if relative.startswith("docs/") and suffix == ".md":
        return "historical-document"
    if nested_notes and suffix == ".md" and any(relative.startswith(x) for x in
            ("Numerical/", "RGE/", "Reports/notes/")):
        return "nested-research-document"
    if BACKUP_NAME.search(name):
        return "obsolete-backup"
    if presentations and suffix in BINARY_PRESENTATIONS and (
            relative.startswith("Reports/output/") or relative.startswith("output/")):
        return "generated-presentation"
    return None


def source_texts(root: Path, paths: Iterable[str], excluded: set[str]) -> list[tuple[str, str]]:
    sources = []
    for rel in paths:
        if rel in excluded or rel.startswith(("output/", "Reports/output/", ".git/")):
            continue
        f = root / rel
        if not f.is_file() or f.suffix.lower() not in SOURCE_SUFFIXES or f.stat().st_size > 2_000_000:
            continue
        try:
            sources.append((rel, f.read_text(encoding="utf-8-sig", errors="replace")))
        except OSError:
            pass
    return sources


def referenced_by(relative: str, sources: list[tuple[str, str]]) -> list[str]:
    base = Path(relative).name
    # Long filenames can be dynamically opened by basename rather than full path.
    needles = [relative, relative.replace("/", "\\")]
    if len(base) >= 18:
        needles.append(base)
    results = []
    for source, content in sources:
        if source == relative:
            continue
        if any(n in content for n in needles):
            results.append(source)
            if len(results) >= 8:
                break
    return results


def general_source_flags(sources: list[tuple[str, str]]) -> list[dict]:
    results = []
    for rel, content in sources:
        if Path(rel).suffix.lower() not in {".py", ".wl", ".wls", ".m", ".ps1"}:
            continue
        matched = {}
        for kind, regex in GENERALITY_MARKERS.items():
            lines = [n for n, line in enumerate(content.splitlines(), 1) if regex.search(line)]
            if lines:
                matched[kind] = lines[:12]
        if matched:
            results.append({"path": rel, "markers": matched})
    return results


def note_review(root: Path, paths: list[str]) -> list[dict]:
    """Triages every tracked non-generated Markdown note; NOT physics verification."""
    review = []
    for rel in paths:
        if Path(rel).suffix.lower() != ".md" or rel.startswith(("output/", "Reports/output/")):
            continue
        f = root / rel
        if not f.is_file() or f.stat().st_size > 500_000:
            continue
        try:
            lines = f.read_text(encoding="utf-8-sig", errors="replace").splitlines()
        except OSError:
            continue
        flags = {}
        for name, regex in NOTE_MARKERS.items():
            hits = [idx for idx, line in enumerate(lines, 1) if regex.search(line)]
            if hits:
                flags[name] = hits[:25]
        review.append({"path": rel, "lines": len(lines), "review_flags": flags})
    return review


def manual_code_review(paths: list[str], sources: list[tuple[str, str]]) -> list[dict]:
    """Only list probable debugging source; never classify it as safe to delete."""
    entries = []
    for path in paths:
        if not path.endswith((".py", ".wl", ".wls")) or path.startswith(("output/", "Reports/output/", "tests/")):
            continue
        name = Path(path).name.lower()
        if not any(word in name for word in ("probe", "experiment", "legacy", "debug", "backup")):
            continue
        entries.append({"path": path, "literal_references": referenced_by(path, sources),
                        "action": "manual-review-only; no automatic source deletion"})
    return entries


def build_plan(root: Path, *, nested_notes: bool = False, presentations: bool = False) -> dict:
    if not (root / "pipeline.py").is_file() or not (root / "Lagrangian").is_dir():
        raise ValueError(f"Not a T3 repository checkout: {root}")
    tracked = tracked_paths(root)
    modified = modified_paths(root)
    possible = {rel: reason for rel in tracked
                if (reason := category(rel, nested_notes=nested_notes, presentations=presentations))
                and (root / rel).is_file()}
    # Exclude candidates from the reference scan: links between retired records
    # do not keep the entire obsolete-document network alive.
    sources = source_texts(
        root,
        sorted(set(tracked) | set(untracked_reference_paths(root)) | KEEP_DOCS),
        set(possible),
    )
    notes = note_review(root, tracked)
    code_review = manual_code_review(tracked, sources)
    moves, held = [], []
    for rel, kind in sorted(possible.items()):
        p = root / rel
        item = {"path": rel, "category": kind, "bytes": p.stat().st_size,
                "sha256": sha256(p)}
        if rel in modified:
            held.append({**item, "reason": "locally modified tracked file"})
            continue
        hits = referenced_by(rel, sources)
        if hits:
            held.append({**item, "reason": "named by retained source/config/document",
                         "references": hits})
        else:
            moves.append(item)
    return {
        "schema": "t3-extreme-cleanup-plan-v1",
        "repository": str(root),
        "git_head": ((git_output(root, "rev-parse", "HEAD") or b"")
                     if is_git_checkout_root(root) else b"").decode("ascii", "replace").strip(),
        "options": {"include_nested_notes": nested_notes, "include_presentations": presentations},
        "archive_candidates": moves,
        "held_for_references": held,
        "general_source_flags": general_source_flags(sources),
        "note_review": notes,
        "manual_code_review": code_review,
        "counts": {"tracked_paths": len(tracked), "candidate_paths": len(possible),
                   "archivable": len(moves), "held": len(held),
                   "archivable_bytes": sum(item["bytes"] for item in moves),
                   "notes_reviewed": len(notes),
                   "notes_flagged": sum(bool(item["review_flags"]) for item in notes),
                   "manual_code_candidates": len(code_review)},
        "warnings": [
            "Name/path search cannot detect all dynamic dependencies.",
            "No physics source, configuration, fitted input, JSON/CSV data, or live matching/RGE export is removed.",
            "The general Lagrangian is not certified by documentation cleanup.",
            "Inspect before using --apply; test local numerical source dependencies afterward.",
        ],
    }


def write_plan(root: Path, plan: dict) -> Path:
    target = root / PLAN_PATH
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(plan, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return target


def apply_plan(root: Path, plan: dict) -> Path | None:
    moves = plan["archive_candidates"]
    if not moves:
        return None
    archive_parent = root.parent / "T3_external_archives"
    archive_parent.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    archive = archive_parent / ("cleanup_" + stamp)
    archive.mkdir(exist_ok=False)
    copied = []
    try:
        # Archive and verify *everything* before removing *anything*.
        for item in moves:
            rel = item["path"]
            src = root / rel
            dst = archive / rel
            if not src.is_file() or sha256(src) != item["sha256"]:
                raise RuntimeError(f"File changed during planning: {rel}")
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dst)
            if sha256(dst) != item["sha256"]:
                raise RuntimeError(f"Archive verification failed: {rel}")
            copied.append(item)
        manifest = {"schema": "t3-extreme-cleanup-archive-v1", "source_root": str(root),
                    "git_head": plan["git_head"], "files": copied,
                    "created_utc": datetime.now(timezone.utc).isoformat()}
        (archive / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        for item in copied:
            rel = item["path"]
            src = root / rel
            if sha256(src) != item["sha256"]:
                raise RuntimeError(f"File changed before removal: {rel}")
            src.unlink()
        return archive
    except Exception:
        # An interrupted failure should leave original files accessible.
        for item in copied:
            rel = item["path"]
            src, dst = root / rel, archive / rel
            if not src.exists() and dst.is_file() and sha256(dst) == item["sha256"]:
                src.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(dst, src)
        raise


def restore(root: Path, archive: Path) -> int:
    archive = archive.resolve()
    if archive.parent != (root.parent / "T3_external_archives").resolve():
        raise ValueError("Archive must be inside the checkout's sibling T3_external_archives folder")
    manifest = json.loads((archive / "manifest.json").read_text(encoding="utf-8"))
    if manifest.get("schema") != "t3-extreme-cleanup-archive-v1" or Path(manifest["source_root"]).resolve() != root:
        raise ValueError("Archive manifest does not belong to this checkout")
    items = manifest["files"]
    for item in items:
        dst = root / item["path"]
        saved = archive / item["path"]
        if dst.exists():
            raise FileExistsError(f"Restore would overwrite: {dst}")
        if not saved.is_file() or sha256(saved) != item["sha256"]:
            raise RuntimeError(f"Archive checksum failure: {item['path']}")
    for item in items:
        dst = root / item["path"]
        saved = archive / item["path"]
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(saved, dst)
        if sha256(dst) != item["sha256"]:
            raise RuntimeError(f"Restore verification failure: {item['path']}")
    return len(items)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--include-nested-notes", action="store_true",
                        help="Opt in to notes under Numerical/, RGE/ and Reports/notes/")
    parser.add_argument("--include-presentations", action="store_true",
                        help="Opt in to potentially used generated PNG/PDF/HTML/SVG; review carefully")
    parser.add_argument("--apply", action="store_true")
    parser.add_argument("--confirm", default="")
    parser.add_argument("--restore", type=Path, metavar="EXTERNAL_ARCHIVE")
    args = parser.parse_args()
    root = args.root.resolve()
    if args.restore is not None:
        print(f"Restored {restore(root, args.restore)} files")
        return 0
    if args.apply and args.confirm != "ARCHIVE_REVIEWED_FILES":
        parser.error("--apply requires --confirm ARCHIVE_REVIEWED_FILES")
    plan = build_plan(root, nested_notes=args.include_nested_notes,
                      presentations=args.include_presentations)
    plan_file = write_plan(root, plan)
    print("Plan:", plan_file)
    print("Candidates:", plan["counts"]["candidate_paths"],
          "archivable:", plan["counts"]["archivable"],
          "held:", plan["counts"]["held"],
          "bytes:", plan["counts"]["archivable_bytes"])
    print("Files requiring general-physics source review:", len(plan["general_source_flags"]))
    print("Notes triaged:", plan["counts"]["notes_reviewed"], "flagged:", plan["counts"]["notes_flagged"],
          "code candidates for manual review:", plan["counts"]["manual_code_candidates"])
    if args.apply:
        destination = apply_plan(root, plan)
        print("Archived outside repository:", destination or "(no candidates)")
        print("Run repository tests and review git status; nothing has been committed or pushed.")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, RuntimeError) as exc:
        print("ERROR:", exc, file=sys.stderr)
        raise SystemExit(1) from exc
