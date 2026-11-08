#!/usr/bin/env python3
"""Organise unreferenced research documents and optional experimental scripts.

Default is preview only. Physical files are moved within the *local* checkout
only after an exact, Git-HEAD-bound plan is reviewed. No network or GitHub write.
All moved bytes are hash-verified, and a restoration manifest is retained.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import subprocess
import sys
from typing import Any

FORMAT = "t3-research-organization-v1"
CONFIRM = "ORGANIZE_RESEARCH_FILES"
PLAN = Path("output/maintenance/research_organization_plan.json")
MANIFEST = Path("output/maintenance/research_organization_manifest.json")
DOC_GROUP_ROOT = "docs/research"
SCRIPT_GROUP_ROOT = "scripts/research/experimental"

# Documents in the main documentation directory are the reader-facing guides
# or consolidated research volumes. Do not move them through name heuristics.
PROTECTED_DOCS = frozenset({
    "README.md", "MODEL_AND_PIPELINE.md", "PHYSICS.md", "SCALAR_SECTOR.md",
    "NUMERICAL.md", "RESULTS.md", "PROJECT_GUIDE.md",
    "REPOSITORY_STRUCTURE.md", "TEST_STRATEGY.md",
    "REGENERATE_CURRENT_REPORTS.md", "GENERAL_SYMBOLIC_PHYSICS_AUDIT.md",
    "RESEARCH_HISTORY_INDEX.md", "RESEARCH_GENERAL_T3_TENSORS_AND_QUARTICS.md",
    "RESEARCH_GENERAL_T3_RUNNING_AND_MATCHING.md",
    "RESEARCH_MATCHETE_AND_MASS_EIGEN.md",
    "RESEARCH_EFT1_THRESHOLDS_AND_WEINBERG.md",
    "RESEARCH_NUMERICAL_FITS_AND_BENCHMARKS.md",
    "RESEARCH_GENERAL_SCALAR_AND_PILOTS.md",
})

TEXT_EXTS = frozenset({
    ".py", ".wl", ".wls", ".m", ".ps1", ".md", ".tex", ".txt",
    ".json", ".toml", ".ini", ".cfg", ".yaml", ".yml", ".csv",
    ".sh", ".bat",
})
SOURCE_DIRS = frozenset({
    "Lagrangian", "Numerical", "RGE", "Reports", "scripts", "tests",
    "docs", "configs", "common", "physics", "model", "studies",
    "input", "tools", "validation", "benchmarks",
})
# The organizer and its tests necessarily contain candidate filename examples.
SCAN_EXCLUSIONS = frozenset({
    "scripts/maintenance/OrganizeResearchFiles.py",
    "tests/python/test_research_organization.py",
})
DOC_LINK = re.compile(r"!?\[[^\]]*\]\(\s*<?([^\s)>]+)", re.IGNORECASE)
HTML_LINK = re.compile(r"(?:href|src)\s*=\s*['\"]([^'\"]+)['\"]", re.IGNORECASE)
REFERENCE_DEF = re.compile(r"^\s*\[[^\]]+\]:\s*<?([^\s>]+)", re.MULTILINE)
BARE_MARKDOWN_NAME = re.compile(r"(?<![A-Za-z0-9_])(?:[A-Za-z0-9_.-]+/)*[A-Za-z0-9_.-]+\.md\b", re.IGNORECASE)


class OrganizationError(RuntimeError):
    """An unsafe or stale organization operation was refused."""


def _git(root: Path, *args: str, allow_diff: bool = False) -> subprocess.CompletedProcess[bytes]:
    p = subprocess.run(("git", *args), cwd=root, capture_output=True, check=False)
    good = (0, 1) if allow_diff else (0,)
    if p.returncode not in good:
        raise OrganizationError(f"Git failed: {' '.join(args)}: {p.stderr.decode(errors='replace')}")
    return p


def _root(root: Path) -> Path:
    root = root.resolve()
    top = Path(os.fsdecode(_git(root, "rev-parse", "--show-toplevel").stdout.strip())).resolve()
    if os.path.normcase(str(top)) != os.path.normcase(str(root)):
        raise OrganizationError("Run from the root of this local Git checkout")
    return root


def _head(root: Path) -> str:
    return _git(root, "rev-parse", "HEAD").stdout.decode("ascii").strip()


def _tracked(root: Path) -> list[str]:
    out = _git(root, "ls-files", "-z", "--cached").stdout
    return sorted({os.fsdecode(x).replace("\\", "/") for x in out.split(b"\0") if x})


def _safe(root: Path, rel: str) -> Path:
    p = PurePosixPath(rel)
    if (not rel or not p.parts or p.is_absolute() or rel != p.as_posix()
            or "\\" in rel or ":" in rel or any(x in ("", ".", "..") for x in p.parts)):
        raise OrganizationError(f"Unsafe relative path: {rel!r}")
    path = root.joinpath(*p.parts)
    # Reject symlink parents; neither the inputs nor destinations may escape.
    part = root
    for name in p.parts:
        part = part / name
        if part.is_symlink():
            raise OrganizationError(f"Symlink in path: {rel}")
    if not path.resolve().is_relative_to(root):
        raise OrganizationError(f"Path escapes checkout: {rel}")
    return path


def _sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def _doc_group(name: str) -> str | None:
    if name in PROTECTED_DOCS or not name.endswith(".md"):
        return None
    u = name.upper()
    if u.startswith(("WORKSTREAM_A_", "A_MINUS4_", "T3_ALL16_", "T3_A_MINUS2_")):
        return "model-and-symmetry"
    if u.startswith(("WORKSTREAM_B_", "WORKSTREAMC_", "T3_BC_", "T3_DIRECT_",
                     "T3_EXTENDED_", "T3_SELF_BILINEAR_", "T3_SIX_SELF_BILINEAR_",
                     "GENERAL_D_M2_", "T3_SCALAR_BOUNDARY_", "T3B_YUKAWA_")):
        return "eft-matching"
    if u.startswith(("WORKSTREAMD_", "T3_FROZEN_", "T3_GENERAL_FIT_", "T3_GENERAL_ZEE_",
                     "T3_REDUCED_", "T3_REFITTED_", "T3_ZERO_", "T3_FULL_SCALAR_",
                     "T3_FULL_STACK_", "T3_GENERAL_DUAL_", "T3_GENERAL_SCALAR_",
                     "T3_GENERAL_SHARED_", "T3_JOINT_", "T3_FROZEN64_",
                     "T3_FULLUV_", "T3_FULL_UV_", "T3_GENERAL_FERMION_",
                     "FERMION_START_", "HIGH_SCALAR_", "GENERAL_SCALAR_",
                     "SIX_GENERAL_")):
        return "numerical-studies"
    if u.startswith(("THESIS_", "T3_THESIS_", "T3_GENERAL_SCALAR_INTERACTIVE_",
                     "T3_FROZEN_DISPLAY_", "T3_REFITTED_FROZEN_SCALAR_DISPLAY_")):
        return "scientific-reports"
    if u.startswith(("T3_GENERAL_GROUP_FACTOR_", "T3_COMPLETE_GENERAL_SCALAR_",
                     "T3_GENERAL_SCALAR_STUDY_", "T3_GENERAL_SCALAR_MATCHING_",
                     "T3_BC_110101_")):
        return "model-and-symmetry"
    if u.startswith(("REPOSITORY_CLEANUP_", "REPOSITORY_MAINTENANCE_",
                     "HISTORICAL_DOCS_", "TEST_BASELINE_", "TEST_FAILURE_")):
        return "project-records"
    if u.startswith(("T3_GENERAL_", "T3_SCALAR_", "T3_FROZEN", "T3_", "WORKSTREAM_")):
        return "research-notes"
    if u in ("GENERAL_WILSON_NORMALIZATION_STAGE.MD", "T3_UV_AUDIT.MD"):
        return "research-notes"
    return None


def _candidates(tracked: list[str], include_experiments: bool) -> dict[str, tuple[str, str]]:
    result: dict[str, tuple[str, str]] = {}
    for rel in tracked:
        parts = PurePosixPath(rel).parts
        if len(parts) == 2 and parts[0] == "docs":
            group = _doc_group(parts[1])
            if group:
                result[rel] = (f"{DOC_GROUP_ROOT}/{group}/{parts[1]}", "historical research documentation")
        if (include_experiments and len(parts) == 2 and parts[0] == "scripts"
                and parts[1].endswith(".experimental.py")):
            clean_name = parts[1].removesuffix(".experimental.py") + ".py"
            result[rel] = (f"{SCRIPT_GROUP_ROOT}/{clean_name}", "standalone experimental study")
        if rel == "mixed_wilson_inputs.zip":
            result[rel] = ("input/reference/mixed_wilson_inputs.zip", "saved reference input bundle")
    return result


def _has_relative_doc_links(text: str) -> bool:
    # Even prose like "see FILENAME.md" becomes misleading after relocation.
    if BARE_MARKDOWN_NAME.search(text):
        return True
    for pattern in (DOC_LINK, HTML_LINK, REFERENCE_DEF):
        for match in pattern.finditer(text):
            url = match.group(1).strip()
            if not url or url.startswith(("#", "/", "//")):
                continue
            if re.match(r"^[a-z][a-z0-9+.-]*:", url, flags=re.I):
                continue
            return True
    return False


def _references(root: Path, tracked: list[str], sources: dict[str, tuple[str, str]]) -> dict[str, list[str]]:
    """Conservative literal scan of retained source, configuration, and docs.

    Dynamic Python/Wolfram references cannot be proven absent by this scan;
    human review and the original Git objects remain important.
    """
    names: dict[str, tuple[str, ...]] = {}
    for rel in sources:
        filename = PurePosixPath(rel).name.casefold()
        tokens = (rel.casefold(), filename)
        if rel.startswith("scripts/"):
            tokens += (filename.removesuffix(".experimental.py"),)
        names[rel] = tokens
    hits: dict[str, list[str]] = {key: [] for key in sources}
    for rel in tracked:
        if rel in sources or rel in SCAN_EXCLUSIONS:
            continue  # separate candidate-to-candidate scan follows below
        parts = PurePosixPath(rel).parts
        if len(parts) > 1 and parts[0] not in SOURCE_DIRS:
            continue
        if len(parts) >= 2 and parts[:2] == ("Reports", "output"):
            continue
        if PurePosixPath(rel).suffix.lower() not in TEXT_EXTS:
            continue
        p = root / rel
        if not p.is_file() or p.is_symlink():
            continue
        try:
            if p.stat().st_size > 2_000_000:
                continue
            text = p.read_text(encoding="utf-8-sig", errors="replace").replace("\\", "/").casefold()
        except OSError:
            continue
        for source, tokens in names.items():
            if any(token in text for token in tokens):
                hits[source].append(rel)
    # A candidate must also not be referenced by another candidate: the
    # second candidate might be held for unrelated reasons, or moved elsewhere.
    for referrer in sources:
        p = root / referrer
        if not p.is_file():
            continue
        try:
            text = p.read_text(encoding="utf-8-sig", errors="replace").replace("\\", "/").casefold()
        except OSError:
            continue
        for target, tokens in names.items():
            if target != referrer and any(token in text for token in tokens):
                hits[target].append(referrer)
    return hits


def build_plan(root: Path, *, include_experiments: bool = False) -> dict[str, Any]:
    root = _root(root)
    tracked = _tracked(root)
    candidates = _candidates(tracked, include_experiments)
    hits = _references(root, tracked, candidates)
    approved: list[dict[str, Any]] = []
    held: list[dict[str, Any]] = []
    for source, (destination, reason) in sorted(candidates.items()):
        try:
            p = _safe(root, source)
            q = _safe(root, destination)
            if not p.is_file():
                raise OrganizationError("original tracked file is missing")
            if q.exists():
                raise OrganizationError("destination already exists")
            if _git(root, "diff", "--quiet", "HEAD", "--", source, allow_diff=True).returncode:
                raise OrganizationError("original differs from Git HEAD (staged or unstaged)")
            if hits[source]:
                raise OrganizationError("referenced by retained project files")
            if source.startswith("docs/") and _has_relative_doc_links(p.read_text(encoding="utf-8-sig", errors="replace")):
                raise OrganizationError("document contains relative hyperlinks requiring migration")
            if source.startswith("scripts/"):
                content = p.read_text(encoding="utf-8-sig", errors="replace")
                if ("__file__" in content or "sys.path" in content
                        or re.search(r"\bfrom\s+\.\w*\s+import\b", content)
                        or "importlib" in content):
                    raise OrganizationError("script has location-dependent or dynamic-import semantics")
            approved.append({
                "source": source, "destination": destination, "reason": reason,
                "sha256": _sha(p), "bytes": p.stat().st_size,
            })
        except (OrganizationError, OSError) as exc:
            record: dict[str, Any] = {"source": source, "destination": destination, "reason": str(exc)}
            if hits.get(source):
                record["referenced_by"] = hits[source][:12]
            held.append(record)
    return {
        "format": FORMAT, "git_head": _head(root),
        "include_experiments": include_experiments,
        "approved": approved, "held": held,
        "approved_count": len(approved),
        "approved_bytes": sum(x["bytes"] for x in approved),
    }


def _json_write(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    tmp.replace(path)


def preview(root: Path, *, include_experiments: bool = False) -> dict[str, Any]:
    root = _root(root)
    result = build_plan(root, include_experiments=include_experiments)
    _json_write(root / PLAN, result)
    return result


def apply(root: Path, *, confirm: str) -> dict[str, Any]:
    root = _root(root)
    if confirm != CONFIRM:
        raise OrganizationError(f"Confirmation required: {CONFIRM}")
    plan_file = root / PLAN
    if not plan_file.is_file():
        raise OrganizationError("Preview and review the plan before applying")
    stored = json.loads(plan_file.read_text(encoding="utf-8"))
    if stored.get("format") != FORMAT:
        raise OrganizationError("Unsupported organization plan format")
    live = build_plan(root, include_experiments=bool(stored["include_experiments"]))
    if stored != live:
        raise OrganizationError("Reviewed plan is stale: preview again")
    if not live["approved"]:
        raise OrganizationError("No safe moves approved in this plan")
    manifest_path = root / MANIFEST
    if manifest_path.exists():
        status = json.loads(manifest_path.read_text(encoding="utf-8")).get("status")
        if status != "restored":
            raise OrganizationError("An applied manifest already exists; restore or preserve it first")
    for item in live["approved"]:
        if _sha(_safe(root, item["source"])) != item["sha256"]:
            raise OrganizationError(f"Source hash changed: {item['source']}")
        if _safe(root, item["destination"]).exists():
            raise OrganizationError(f"Destination appeared: {item['destination']}")
    moved: list[dict[str, Any]] = []
    try:
        for item in live["approved"]:
            src = _safe(root, item["source"])
            dst = _safe(root, item["destination"])
            dst.parent.mkdir(parents=True, exist_ok=True)
            os.replace(src, dst)
            moved.append(item)
            if _sha(dst) != item["sha256"]:
                raise OrganizationError(f"Hash mismatch after moving: {item['source']}")
        manifest = {
            "format": FORMAT, "status": "applied",
            "git_head": live["git_head"],
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "moves": moved,
        }
        _json_write(manifest_path, manifest)
    except Exception:
        for item in reversed(moved):
            dst = root / item["destination"]
            src = root / item["source"]
            if dst.exists() and not src.exists():
                src.parent.mkdir(parents=True, exist_ok=True)
                os.replace(dst, src)
        raise
    return manifest


def restore(root: Path) -> dict[str, Any]:
    root = _root(root)
    path = root / MANIFEST
    if not path.is_file():
        raise OrganizationError(f"Restoration manifest not found: {MANIFEST}")
    manifest = json.loads(path.read_text(encoding="utf-8"))
    if manifest.get("format") != FORMAT or manifest.get("status") != "applied":
        raise OrganizationError("Manifest is not an applied organization record")
    # Fail before any change if any file has been changed or overwritten.
    for item in manifest["moves"]:
        src = _safe(root, item["source"])
        dst = _safe(root, item["destination"])
        if src.exists() or not dst.is_file() or _sha(dst) != item["sha256"]:
            raise OrganizationError(f"Restore blocked for changed/missing destination: {item['destination']}")
    reversed_items: list[dict[str, Any]] = []
    try:
        for item in reversed(manifest["moves"]):
            src = root / item["source"]
            dst = root / item["destination"]
            src.parent.mkdir(parents=True, exist_ok=True)
            os.replace(dst, src)
            reversed_items.append(item)
            if _sha(src) != item["sha256"]:
                raise OrganizationError(f"Restore hash mismatch: {item['source']}")
        manifest["status"] = "restored"
        manifest["restored_utc"] = datetime.now(timezone.utc).isoformat()
        _json_write(path, manifest)
    except Exception:
        for item in reversed(reversed_items):
            src = root / item["source"]
            dst = root / item["destination"]
            if src.exists() and not dst.exists():
                dst.parent.mkdir(parents=True, exist_ok=True)
                os.replace(src, dst)
        raise
    return manifest


def _script_category(name: str) -> str:
    n = name.casefold()
    if n.endswith(".experimental.py"):
        return "Experimental research programs"
    if n.startswith(("export", "generate", "analyzegeneralt3rgbeta", "runphysicalb2rgbeta")) or n.endswith(".wl"):
        return "Symbolic model and RGE exports"
    if n.startswith(("rungeneral", "runt3bc", "runaminus", "runt3direct", "runlo", "runphysical")):
        return "Matching, EFT and model-specific studies"
    if any(x in n for x in ("fit", "scan", "benchmark", "sobol", "neutrino", "observable")):
        return "Numerical studies and parameter fitting"
    if n.startswith(("build", "plot", "report", "polish", "summar", "readable")):
        return "Reports and visualisations"
    if n.startswith(("audit", "check", "test", "verify", "probe", "inspect", "diagnos", "review", "validate")):
        return "Scientific validation and diagnostics"
    if n.startswith(("clean", "archive", "retire", "remove", "rename", "install")):
        return "Maintenance and historical utilities"
    return "Specialised research commands"


def _write_indexes(root: Path) -> tuple[Path, Path]:
    """Create navigational lists; do not change or certify calculation modules."""
    root = _root(root)
    scripts = root / "scripts"
    categorized: dict[str, list[str]] = {}
    for p in sorted(scripts.iterdir()):
        if p.is_file() and p.suffix.lower() in (".py", ".wl", ".ps1"):
            group = _script_category(p.name)
            categorized.setdefault(group, []).append(p.name)
    out = [
        "# Research command index", "",
        "This index groups existing scripts by **purpose**, not by completion status.",
        "The root Python entry points are `pipeline.py` and `full16_symbolic.py`.",
        "Scientific routines remain at their existing paths when moving them",
        "would invalidate imports, tests or saved provenance.", "",
    ]
    for group, files in sorted(categorized.items()):
        out.extend((f"## {group}", ""))
        out.extend(f"- [`{name}`]({name})" for name in files)
        out.append("")
    out.extend((
        "## Experimental research programs", "",
        "Standalone experimental scripts relocated after dependency review are",
        "listed under [`research/experimental/`](research/experimental/).",
        "From the repository root, run a relocated Python program via",
        "`python -X utf8 -m scripts.research.experimental.ModuleName`.", "",
    ))
    script_index = scripts / "INDEX.md"
    script_index.write_text("\n".join(out), encoding="utf-8")
    directory = root / DOC_GROUP_ROOT
    docs_out = [
        "# Research document index", "",
        "Historical technical notes are arranged by the physics question they",
        "address. They are research records, not independent certification.", "",
    ]
    names = {
        "model-and-symmetry": "Model representations, group invariants and UV completeness",
        "eft-matching": "Heavy-particle matching, scalar thresholds and Wilson coefficients",
        "numerical-studies": "RGE evolution, benchmarks and exploratory numerical scans",
        "scientific-reports": "Figures, thesis material and reporting studies",
        "project-records": "Maintenance and verification history",
        "research-notes": "Additional specialised studies",
    }
    for group, description in names.items():
        docs_out.extend((f"## {description}", ""))
        folder = directory / group
        for p in sorted(folder.glob("*.md")) if folder.is_dir() else []:
            docs_out.append(f"- [{p.stem.replace('_', ' ').title()}]({group}/{p.name})")
        if not folder.is_dir() or not any(folder.glob("*.md")):
            docs_out.append("No documents moved into this category yet.")
        docs_out.append("")
    doc_index = directory / "INDEX.md"
    doc_index.parent.mkdir(parents=True, exist_ok=True)
    doc_index.write_text("\n".join(docs_out), encoding="utf-8")
    return script_index, doc_index


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2])
    action = ap.add_mutually_exclusive_group()
    action.add_argument("--apply", action="store_true")
    action.add_argument("--restore", action="store_true")
    action.add_argument("--index", action="store_true")
    ap.add_argument("--include-experiments", action="store_true",
                    help="Also consider standalone .experimental.py scripts; never the default")
    ap.add_argument("--confirm", default="")
    args = ap.parse_args(argv)
    try:
        if args.index:
            a, b = _write_indexes(args.root)
            print(f"Updated navigation: {a} and {b}")
        elif args.restore:
            manifest = restore(args.root)
            _write_indexes(args.root)
            print(f"Restored {len(manifest['moves'])} files. Verify the core suite.")
        elif args.apply:
            manifest = apply(args.root, confirm=args.confirm)
            _write_indexes(args.root)
            print(f"Moved and SHA-256 verified {len(manifest['moves'])} files.")
            print(f"Restore using: python -X utf8 scripts/maintenance/OrganizeResearchFiles.py --restore")
        else:
            plan = preview(args.root, include_experiments=args.include_experiments)
            print(f"Git HEAD: {plan['git_head']}")
            print(f"Approved: {plan['approved_count']} files ({plan['approved_bytes']:,} bytes); held: {len(plan['held'])}")
            print(f"Review exact paths in {PLAN.as_posix()}")
            print("Preview only: nothing moved.")
    except (OrganizationError, OSError, ValueError, KeyError, json.JSONDecodeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
