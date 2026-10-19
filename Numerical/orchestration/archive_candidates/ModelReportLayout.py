"""Consolidate symbolic PDF/TeX reports by model, without physics reruns.

Scenario-specific numerical outputs and raw calculation artifacts are untouched.
If contents differ for a nominally shared report, all versions are preserved:
the first scenario becomes the canonical model report, and differing copies go
into _scenario_variants/<scenario>/. No differing content is discarded.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from Numerical.orchestration.ModelNames import is_model, shorten
import shutil

SCENARIOS = ("smallY_smallL", "smallY_largeL", "largeY_smallL", "largeY_largeL")
SYMBOLIC_DIRS = ("Lagrangian", "RGE", "GroupFactors", "Matching")
SYMBOLIC_SUFFIXES = {".pdf", ".tex"}


@dataclass(frozen=True)
class Candidate:
    scenario: str
    model: str
    relative: Path
    source: Path
    digest: str


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1048576), b""):
            h.update(block)
    return h.hexdigest()


def candidates(root: Path) -> list[Candidate]:
    found = []
    for scenario in SCENARIOS:
        scenario_root = root / "full" / "comparison" / scenario
        if not scenario_root.is_dir():
            continue
        for model_dir in sorted(scenario_root.iterdir()):
            if not model_dir.is_dir() or not is_model(model_dir.name):
                continue
            for category in SYMBOLIC_DIRS:
                folder = model_dir / category
                if not folder.is_dir():
                    continue
                for path in sorted(folder.rglob("*")):
                    if path.is_file() and path.suffix.lower() in SYMBOLIC_SUFFIXES:
                        found.append(Candidate(scenario, shorten(model_dir.name),
                                               path.relative_to(model_dir), path,
                                               sha256(path)))
    return found


def migrate(report_root: Path, *, execute: bool = False, prune: bool = True) -> dict:
    """Plan or execute a collision-safe consolidation. Dry-run by default.

    prune=True makes actual moves (copy, verify, remove). Set prune=False to
    retain the originals. Existing destination conflicts are preserved as
    variants rather than overwritten.
    """
    root = Path(report_root).resolve()
    dest_root = root / "models"
    srcs = candidates(root)
    groups: dict[tuple[str, Path], list[Candidate]] = {}
    for candidate in srcs:
        groups.setdefault((candidate.model, candidate.relative), []).append(candidate)
    plan: list[tuple[Candidate, Path]] = []
    variants = []
    for (model, relative), rows in sorted(groups.items(), key=lambda pair: str(pair[0])):
        canonical = dest_root / model / relative
        canonical_hash = sha256(canonical) if canonical.is_file() else rows[0].digest
        for row in rows:
            if row.digest == canonical_hash:
                destination = canonical
            else:
                destination = dest_root / model / "_scenario_variants" / row.scenario / relative
                variants.append({"source": str(row.source), "destination": str(destination),
                                 "reason": "content differs from canonical report"})
            plan.append((row, destination))

    # Abort before writes on any pre-existing destination collision.
    conflicts = []
    by_dest: dict[Path, set[str]] = {}
    for row, destination in plan:
        by_dest.setdefault(destination, set()).add(row.digest)
    for destination, hashes in by_dest.items():
        if destination.exists():
            if not destination.is_file():
                conflicts.append(str(destination) + " is not a file")
            else:
                hashes.add(sha256(destination))
        if len(hashes) > 1:
            conflicts.append(str(destination) + " has incompatible contents")
    result = {"status": "Conflict" if conflicts else "Ready",
              "mode": "apply" if execute else "dry-run",
              "report_root": str(root), "source_files": len(srcs),
              "unique_destinations": len(by_dest),
              "models": sorted({r.model for r in srcs}),
              "different_content_variants": variants,
              "conflicts": conflicts, "copied": 0, "already_present": 0,
              "removed_source_files": 0,
              "note": "Only symbolic PDF/TeX files are relocated; raw outputs and numerical reports are untouched."}
    if not execute or conflicts:
        return result
    for destination, rows in sorted(((d, [(r, target) for r, target in plan if target == d])
                                     for d in by_dest), key=lambda x: str(x[0])):
        if destination.exists():
            result["already_present"] += 1
            continue
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(rows[0][0].source, destination)
        if sha256(destination) != rows[0][0].digest:
            raise IOError(f"Hash verification failed: {destination}")
        result["copied"] += 1
    if prune:
        for row, destination in plan:
            if sha256(destination) != row.digest:
                raise IOError(f"Cannot remove unverified source: {row.source}")
            row.source.unlink()
            result["removed_source_files"] += 1
        for scenario in SCENARIOS:
            base = root / "full" / "comparison" / scenario
            if not base.is_dir():
                continue
            for model_dir in base.iterdir():
                if not model_dir.is_dir() or not is_model(model_dir.name):
                    continue
                for category in SYMBOLIC_DIRS:
                    folder = model_dir / category
                    if folder.is_dir():
                        for directory in sorted((p for p in folder.rglob("*") if p.is_dir()),
                                                key=lambda p: len(p.parts), reverse=True):
                            if not any(directory.iterdir()):
                                directory.rmdir()
                        if not any(folder.iterdir()):
                            folder.rmdir()
    result["status"] = "Complete"
    dest_root.mkdir(parents=True, exist_ok=True)
    (dest_root / "report_migration_manifest.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result
