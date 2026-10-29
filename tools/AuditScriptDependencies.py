"""Read-only inventory of scripts/ dependencies and provenance references.

The scanner deliberately NEVER archives or deletes files. A script with no
static references is not automatically safe to remove: string-generated imports,
subprocesses, and external invocations may not be statically visible.

Run from any directory: python -X utf8 -m tools.AuditScriptDependencies
"""
from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import json
import os
import re
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = Path("output/cleanup_audit/script_dependencies_v1")
SCANNED_EXTENSIONS = frozenset({
    ".py", ".md", ".rst", ".txt", ".json", ".jsonl", ".csv",
    ".yaml", ".yml", ".toml", ".ini", ".cfg", ".ps1", ".sh",
    ".wl", ".wls", ".tex", ".bib", ".bat", ".cmd",
})
SKIP_DIRS = frozenset({
    ".git", ".venv", "venv", ".pytest_cache", "__pycache__",
    "node_modules", ".mypy_cache", ".ruff_cache", ".tox",
})
# Declared preserved workflow entrypoints, not a claim of dynamic verification.
PRESERVE_ENTRYPOINTS = frozenset({
    "RunT3ReducedScalarSobolComparison",
    "BuildT3ReducedScalarFinalResults",
    "BuildReducedScalarExistingInteractiveReports",
    "BuildReducedScalarLambdaStyleReports",
    "RunT3GeneralScalar48Campaign",
    "BuildGeneralScalar48FinalReports",
    "RunT3ZeroExtraBoundaryBenchmark",
    "RunT3ZeroBoundaryLowEnergyC5",
    "RunGeneralT3FrozenScalarC5",
    "RunGeneralT3NumericalC5",
    "Stage36FinalWeinbergRunning",
})
REF_CLASSES = (
    "STATIC_IMPORT", "TEST_REFERENCE", "CODE_REFERENCE",
    "SHELL_REFERENCE", "DOC_REFERENCE", "CONFIG_REFERENCE",
    "SAVED_RESULT_REFERENCE", "TEXT_REFERENCE",
)


def rel(path: Path, root: Path) -> str:
    return path.relative_to(root).as_posix()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def head(root: Path) -> str:
    try:
        p = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=root,
            check=False, capture_output=True, text=True, timeout=10,
        )
        return p.stdout.strip() if p.returncode == 0 else "unavailable"
    except (OSError, subprocess.TimeoutExpired):
        return "unavailable"


def classify_source(path: str) -> str:
    suffix = Path(path).suffix.lower()
    if path.startswith("tests/") or "/tests/" in path:
        return "TEST_REFERENCE"
    if suffix == ".py":
        return "CODE_REFERENCE"
    if suffix in {".ps1", ".sh", ".bat", ".cmd"}:
        return "SHELL_REFERENCE"
    if suffix in {".md", ".rst", ".tex", ".txt"}:
        return "DOC_REFERENCE"
    if path.startswith("output/") or path.startswith("Reports/output/"):
        return "SAVED_RESULT_REFERENCE"
    if suffix in {".json", ".jsonl", ".yaml", ".yml", ".toml", ".ini", ".cfg"}:
        return "CONFIG_REFERENCE"
    return "TEXT_REFERENCE"


def ast_imports(source: str, filename: str, inventory: set[str]) -> tuple[list[tuple[str, int]], str | None]:
    """Return statically resolvable imports from the scripts package."""
    try:
        tree = ast.parse(source, filename=filename)
    except (SyntaxError, ValueError) as exc:
        return [], f"{type(exc).__name__}: {exc}"
    found: list[tuple[str, int]] = []
    in_scripts = filename.startswith("scripts/")

    def register(symbol: str, line: int) -> None:
        if symbol in inventory:
            found.append((symbol, line))

    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                parts = alias.name.split(".")
                if len(parts) >= 2 and parts[0] == "scripts":
                    register(parts[1], node.lineno)
        elif isinstance(node, ast.ImportFrom):
            module = node.module or ""
            if node.level == 0 and module == "scripts":
                for alias in node.names:
                    register(alias.name, node.lineno)
            elif node.level == 0 and module.startswith("scripts."):
                register(module.split(".")[1], node.lineno)
            elif in_scripts and node.level == 1 and not module:
                for alias in node.names:
                    register(alias.name, node.lineno)
            elif in_scripts and node.level == 1 and module:
                register(module.split(".")[0], node.lineno)
        elif isinstance(node, ast.Call) and node.args:
            first = node.args[0]
            if not isinstance(first, ast.Constant) or not isinstance(first.value, str):
                continue
            name = first.value
            if not name.startswith("scripts."):
                continue
            func = node.func
            if (isinstance(func, ast.Name) and func.id in {"__import__", "import_module", "run_module"}) or (
                isinstance(func, ast.Attribute) and func.attr in {"import_module", "run_module"}
            ):
                register(name.split(".")[1], node.lineno)
    return found, None


def examine(root: Path, *, max_bytes: int = 8 * 1024 * 1024) -> dict:
    scripts_path = root / "scripts"
    if not scripts_path.is_dir():
        raise ValueError("No scripts/ directory at " + str(root))
    script_files = sorted(p for p in scripts_path.glob("*.py") if p.name != "__init__.py")
    if not script_files:
        raise ValueError("No script modules found")
    inventory = {p.stem for p in script_files}
    # A text mention is deliberately more conservative than a Python import.
    choices = "|".join(re.escape(name) for name in sorted(inventory, key=lambda n: (-len(n), n)))
    pattern = re.compile(r"(?<![A-Za-z0-9_])(?:scripts[./\\])?(?P<name>" + choices + r")(?:\.py)?(?![A-Za-z0-9_])")
    references: list[dict] = []
    seen: set[tuple[str, int, str]] = set()
    source_errors: list[dict] = []
    parse_errors: list[dict] = []
    skipped_large: list[dict] = []
    scanned = 0
    source_count = 0

    input_paths: list[Path] = []
    for folder, dirs, files in os.walk(root):
        parent = Path(folder)
        dirs[:] = sorted(d for d in dirs if d not in SKIP_DIRS
                         and not rel(parent / d, root).startswith("output/cleanup_audit"))
        input_paths.extend(parent / name for name in files)

    for path in sorted(input_paths):
        if not path.is_file():
            continue
        relative = rel(path, root)
        if relative.startswith("output/cleanup_audit/"):
            continue  # Prevent scanner results from becoming references on reruns.
        if relative == "tools/AuditScriptDependencies.py":
            continue  # Its declared PRESERVE_ENTRYPOINTS list is not evidence.
        if path.suffix.lower() not in SCANNED_EXTENSIONS:
            continue
        try:
            size = path.stat().st_size
            if size > max_bytes:
                skipped_large.append({"path": relative, "bytes": size})
                continue
            raw = path.read_bytes()
            source = raw.decode("utf-8-sig", errors="replace")
        except OSError as exc:
            source_errors.append({"path": relative, "error": f"{type(exc).__name__}: {exc}"})
            continue
        scanned += 1
        lines = source.splitlines()
        if path.suffix.lower() == ".py":
            source_count += 1
            imported, error = ast_imports(source, relative, inventory)
            if error:
                parse_errors.append({"path": relative, "error": error})
            for name, number in imported:
                if relative == f"scripts/{name}.py":
                    continue
                key = (relative, number, name)
                if key not in seen:
                    references.append({"script": name, "source": relative,
                                       "line": number, "kind": "STATIC_IMPORT",
                                       "excerpt": lines[number - 1][:240] if number <= len(lines) else ""})
                    seen.add(key)

        default_kind = classify_source(relative)
        for number, line in enumerate(lines, start=1):
            for match in pattern.finditer(line):
                name = match.group("name")
                if relative == f"scripts/{name}.py":
                    continue
                key = (relative, number, name)
                if key not in seen:
                    references.append({"script": name, "source": relative,
                                       "line": number, "kind": default_kind,
                                       "excerpt": line.strip()[:240]})
                    seen.add(key)

    counts: dict[str, Counter] = defaultdict(Counter)
    for reference in references:
        counts[reference["script"]][reference["kind"]] += 1
    rows = []
    for path in script_files:
        name = path.stem
        count = counts[name]
        hard = count["STATIC_IMPORT"] + count["TEST_REFERENCE"] + count["CODE_REFERENCE"]
        note = count["DOC_REFERENCE"] + count["SHELL_REFERENCE"] + count["CONFIG_REFERENCE"] + count["TEXT_REFERENCE"]
        saved = count["SAVED_RESULT_REFERENCE"]
        if name in PRESERVE_ENTRYPOINTS:
            classification = "RETAIN_DECLARED_WORKFLOW"
        elif hard:
            classification = "RETAIN_STATIC_CODE_OR_TEST_REFERENCE"
        elif note:
            classification = "REVIEW_DOCUMENTED_ENTRYPOINT"
        elif saved:
            classification = "REVIEW_FROZEN_PROVENANCE"
        else:
            classification = "REVIEW_NO_STATIC_REFERENCE"
        rows.append({"script": name, "project_path": rel(path, root),
                     "classification": classification,
                     "sha256": sha256(path),
                     **{kind.lower(): count[kind] for kind in REF_CLASSES},
                     "references_total": sum(count.values())})

    references.sort(key=lambda r: (r["script"], r["kind"], r["source"], r["line"]))
    return {
        "format": "t3-script-dependency-audit-v1",
        "head": head(root), "scripts_total": len(rows),
        "scanned_text_files": scanned, "scanned_python_files": source_count,
        "source_read_errors": source_errors, "python_parse_errors": parse_errors,
        "large_text_files_skipped": skipped_large,
        "classifications": dict(Counter(row["classification"] for row in rows)),
        "rows": rows, "references": references,
        "warning": "Read-only static analysis; no classification authorises deleting or archiving a script.",
    }


def render_summary(result: dict) -> str:
    lines = [
        "# T3 script dependency inventory (read-only)", "",
        f"- HEAD: `{result['head']}`",
        f"- Script modules inventoried: **{result['scripts_total']}**",
        f"- Text files scanned: **{result['scanned_text_files']}**",
        f"- Python files parsed or attempted: **{result['scanned_python_files']}**",
        f"- Python parse errors: **{len(result['python_parse_errors'])}**",
        f"- Text read errors: **{len(result['source_read_errors'])}**",
        f"- Oversize text files skipped: **{len(result['large_text_files_skipped'])}**", "",
        "**No scripts have been archived. No 'unreferenced' script is automatically safe to archive.**", "",
        "| Classification | Modules |", "|---|---:|",
    ]
    for key, value in sorted(result["classifications"].items()):
        lines.append(f"| `{key}` | {value} |")
    lines += ["", "## Inventory", "",
              "| Script | Classification | Static imports | Tests | Other source | Docs | Saved results |",
              "|---|---|---:|---:|---:|---:|---:|"]
    for row in result["rows"]:
        lines.append("| `{script}` | `{classification}` | {static_import} | {test_reference} | {code_reference} | {doc_reference} | {saved_result_reference} |".format(**row))
    lines += ["", "## First-order dependencies of declared workflows", ""]
    protected = {r["script"] for r in result["rows"] if r["classification"] == "RETAIN_DECLARED_WORKFLOW"}
    for row in result["rows"]:
        if row["script"] not in protected:
            continue
        found = sorted({r["script"] for r in result["references"]
                        if r["kind"] == "STATIC_IMPORT" and r["source"] == row["project_path"]})
        lines.append(f"- `{row['script']}` imports from scripts/: " + (", ".join(f"`{n}`" for n in found) if found else "none detected"))
    lines += ["", "## Interpretation", "",
              "- `STATIC_IMPORT` is a resolvable Python import or literal dynamic import; retain both ends while dependencies exist.",
              "- Tests, code, shell commands, and documented entrypoints can still require scripts without formal imports.",
              "- `SAVED_RESULT_REFERENCE` denotes immutable result provenance. Do not rewrite saved fit JSONs to make a file seem unused.",
              "- Absence of references cannot detect external invocations or dynamically constructed module names.",
              "- A separate review of runtime invocation and reproducibility is required before any retirement.",
              "- `references.csv` lists exact file/line evidence. `audit.json` stores the full inventory and SHA-256 hashes."]
    if result["source_read_errors"] or result["python_parse_errors"] or result["large_text_files_skipped"]:
        lines += ["", "## Coverage exceptions", "",
                  "Inspect `audit.json` for unreadable files, parse errors and skipped large sources before making changes."]
    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", type=Path, default=ROOT,
                    help="Repository checkout (defaults to source script parents[1])")
    ap.add_argument("--output", type=Path, default=DEFAULT_OUT,
                    help="Output directory (relative to --root by default)")
    ap.add_argument("--max-text-mb", type=int, default=8,
                    help="Maximum input text size per file in MiB (default: 8)")
    args = ap.parse_args(argv)
    root = args.root.resolve()
    if args.max_text_mb < 1:
        ap.error("--max-text-mb must be >= 1")
    output = (root / args.output).resolve() if not args.output.is_absolute() else args.output.resolve()
    # Avoid accidental replacement of the checkout through a bad --output.
    if output == root or output == root / "scripts":
        ap.error("--output must be a dedicated report directory")
    result = examine(root, max_bytes=args.max_text_mb * 1024 * 1024)
    output.mkdir(parents=True, exist_ok=True)
    (output / "audit.json").write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    with (output / "inventory.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(result["rows"][0]))
        writer.writeheader()
        writer.writerows(result["rows"])
    with (output / "references.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=["script", "source", "line", "kind", "excerpt"])
        writer.writeheader()
        writer.writerows(result["references"])
    (output / "summary.md").write_text(render_summary(result), encoding="utf-8")
    print(f"SCRIPT MODULES: {result['scripts_total']}")
    print(f"TEXT FILES SCANNED: {result['scanned_text_files']}")
    for key, value in sorted(result["classifications"].items()):
        print(f"{key}: {value}")
    print(f"COVERAGE WARNINGS: read={len(result['source_read_errors'])}, parse={len(result['python_parse_errors'])}, large={len(result['large_text_files_skipped'])}")
    print("READ-ONLY: no scripts or source data moved or deleted.")
    print("REPORT:", output / "summary.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
