"""Read-only, second-pass classification of scripts/*.py references.

Separates AST-resolved imports and literal execution calls from ordinary text
mentions. Uses the installed first-pass scanner to scan the *current checkout*;
only writes reports below output/cleanup_audit/script_dependencies_refined_v2/.

No import or execution of scripts/ modules, no archival, no modifications to
source, tests, or saved numerical results.

Usage: python -X utf8 -m tools.RefineScriptDependencyAudit
"""
from __future__ import annotations

import argparse
import ast
import csv
from collections import Counter, defaultdict, deque
import json
import os
from pathlib import Path
import re

from tools.AuditScriptDependencies import ROOT, SKIP_DIRS, examine

OUT_REL = Path("output/cleanup_audit/script_dependencies_refined_v2")
SELF_SOURCES = {
    "tools/RefineScriptDependencyAudit.py",
    "docs/SCRIPT_DEPENDENCY_REFINEMENT.md",
    "tests/python/test_t3_refined_script_dependencies.py",
}
SOURCE_EDGE_KINDS = frozenset({"PYTHON_IMPORT", "PYTHON_LITERAL_MODULE_LOAD", "PYTHON_LITERAL_EXECUTION"})
TEST_EDGE_KINDS = frozenset({"TEST_IMPORT", "TEST_LITERAL_MODULE_LOAD", "TEST_LITERAL_EXECUTION"})


def _qualified_func(node: ast.expr) -> str:
    if isinstance(node, ast.Name):
        return node.id
    if isinstance(node, ast.Attribute):
        left = _qualified_func(node.value)
        return f"{left}.{node.attr}" if left else node.attr
    return ""


def _literal_strings(node: ast.AST) -> list[str]:
    """String fragments of a *call argument*, not arbitrary source lines."""
    return [n.value for n in ast.walk(node)
            if isinstance(n, ast.Constant) and isinstance(n.value, str)]


def _target_from_module(text: str, inventory: set[str]) -> str | None:
    m = re.fullmatch(r"scripts\.([A-Za-z][A-Za-z0-9_]*)", text.strip())
    return m.group(1) if m and m.group(1) in inventory else None


def _targets_from_executable_arg(text: str, inventory: set[str]) -> set[str]:
    # These are explicit module or script-path tokens in a known execution API.
    ans = set()
    for m in re.finditer(r"(?<![A-Za-z0-9_])scripts\.([A-Za-z][A-Za-z0-9_]*)(?![A-Za-z0-9_])", text):
        if m.group(1) in inventory:
            ans.add(m.group(1))
    for m in re.finditer(r"(?<![A-Za-z0-9_])(?:scripts[/\\])?([A-Za-z][A-Za-z0-9_]*)\.py(?![A-Za-z0-9_])", text):
        if m.group(1) in inventory:
            ans.add(m.group(1))
    return ans


def classify_python_evidence(source: str, filename: str, inventory: set[str]) -> tuple[list[dict], str | None]:
    """Return genuine syntactic dependencies; calls with nonliteral args are not guessed.

    Only known import APIs and execution APIs constitute strong evidence.
    Documentation strings, comments, ordinary literals and function *names*
    do not become dependency edges merely because they mention a script.
    """
    try:
        tree = ast.parse(source, filename=filename)
    except (SyntaxError, ValueError) as exc:
        return [], f"{type(exc).__name__}: {exc}"

    is_test = filename.startswith("tests/") or "/tests/" in filename
    prefix = "TEST" if is_test else "PYTHON"
    found: list[dict] = []
    seen: set[tuple[str, int, str]] = set()

    def add(name: str, line: int, kind: str) -> None:
        if name not in inventory or filename == f"scripts/{name}.py":
            return
        key = (name, line, kind)
        if key not in seen:
            found.append({"script": name, "source": filename, "line": line, "kind": kind})
            seen.add(key)

    inside_scripts = filename.startswith("scripts/")
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                parts = alias.name.split(".")
                if len(parts) >= 2 and parts[0] == "scripts":
                    add(parts[1], node.lineno, f"{prefix}_IMPORT")
        elif isinstance(node, ast.ImportFrom):
            mod = node.module or ""
            if node.level == 0 and mod == "scripts":
                for alias in node.names:
                    add(alias.name, node.lineno, f"{prefix}_IMPORT")
            elif node.level == 0 and mod.startswith("scripts."):
                add(mod.split(".")[1], node.lineno, f"{prefix}_IMPORT")
            elif inside_scripts and node.level == 1 and not mod:
                for alias in node.names:
                    add(alias.name, node.lineno, f"{prefix}_IMPORT")
            elif inside_scripts and node.level == 1 and mod:
                add(mod.split(".")[0], node.lineno, f"{prefix}_IMPORT")
        elif isinstance(node, ast.Call) and node.args:
            qual = _qualified_func(node.func)
            # A literal module load is a dependency, even if not a normal import.
            if qual in {"importlib.import_module", "import_module", "__import__", "runpy.run_module", "run_module"}:
                if isinstance(node.args[0], ast.Constant) and isinstance(node.args[0].value, str):
                    target = _target_from_module(node.args[0].value, inventory)
                    if target:
                        add(target, node.lineno, f"{prefix}_LITERAL_MODULE_LOAD")
            # A literal execution call references scripts as a runtime dependency.
            elif qual in {
                "runpy.run_path", "run_path", "subprocess.run", "subprocess.call",
                "subprocess.check_call", "subprocess.check_output", "subprocess.Popen",
                "Popen", "os.system", "os.popen", "os.execv", "os.execvp",
            }:
                for literal in _literal_strings(node.args[0]):
                    for target in _targets_from_executable_arg(literal, inventory):
                        add(target, node.lineno, f"{prefix}_LITERAL_EXECUTION")
    return sorted(found, key=lambda d: (d["source"], d["line"], d["script"], d["kind"])), None


def _source_text(root: Path, source: str) -> str:
    return (root / source).read_text(encoding="utf-8-sig", errors="replace")


def refined_report(root: Path, base: dict) -> dict:
    if base["python_parse_errors"] or base["source_read_errors"] or base["large_text_files_skipped"]:
        raise RuntimeError("First-pass audit coverage is incomplete. Inspect warnings before interpreting refinement.")

    inventory = {r["script"] for r in base["rows"]}
    # 1st pass is conservative: its text references are *not* imported as strong edges.
    first_refs = [r for r in base["references"] if r["source"] not in SELF_SOURCES]
    path_to_refs: dict[str, list[dict]] = defaultdict(list)
    for r in first_refs:
        path_to_refs[r["source"]].append(r)

    strong_edges: list[dict] = []
    parse_errors: list[dict] = []
    # Rescan *every* Python file, including tests, even if the first scanner
    # found no textual name; this catches AST modules built from literals.
    py_paths: set[str] = {r["source"] for r in first_refs if r["source"].endswith(".py")}
    for directory, dirs, files in os.walk(root):
        parent = Path(directory)
        dirs[:] = sorted(d for d in dirs if d not in SKIP_DIRS
                         and not (parent / d).relative_to(root).as_posix().startswith("output/cleanup_audit"))
        for filename in files:
            if filename.endswith(".py"):
                path = (parent / filename).relative_to(root).as_posix()
                if not path.startswith("output/cleanup_audit/"):
                    py_paths.add(path)
    pyfiles = sorted(py_paths)
    for source in pyfiles:
        if source in SELF_SOURCES or not (root / source).is_file():
            continue
        text = _source_text(root, source)
        edges, error = classify_python_evidence(text, source, inventory)
        if error:
            parse_errors.append({"source": source, "error": error})
        else:
            lines = text.splitlines()
            for e in edges:
                e["excerpt"] = lines[e["line"]-1].strip()[:240] if e["line"] <= len(lines) else ""
                strong_edges.append(e)
    if parse_errors:
        raise RuntimeError("Python parse failures in second-pass audit: " + str(parse_errors[:5]))

    by_script: dict[str, list[dict]] = defaultdict(list)
    for e in strong_edges:
        by_script[e["script"]].append(e)
    weak_by_script: dict[str, list[dict]] = defaultdict(list)
    for r in first_refs:
        weak_by_script[r["script"]].append(r)

    # Active workflow closure follows only syntactically identified dependencies
    # from current declared active scripts, never textual mentions.
    declared = {r["script"] for r in base["rows"] if r["classification"] == "RETAIN_DECLARED_WORKFLOW"}
    imported_by_script: dict[str, set[str]] = defaultdict(set)
    for edge in strong_edges:
        src = edge["source"]
        if src.startswith("scripts/") and src.endswith(".py") and edge["kind"] in SOURCE_EDGE_KINDS:
            imported_by_script[Path(src).stem].add(edge["script"])
    closure = set(declared)
    queue = deque(sorted(declared))
    while queue:
        for child in sorted(imported_by_script[queue.popleft()] - closure):
            closure.add(child)
            queue.append(child)

    rows = []
    for old in base["rows"]:
        name = old["script"]
        edges = by_script[name]
        refs = weak_by_script[name]
        non_test_import = sum(e["kind"] == "PYTHON_IMPORT" for e in edges)
        non_test_dynamic = sum(e["kind"] in {"PYTHON_LITERAL_MODULE_LOAD", "PYTHON_LITERAL_EXECUTION"} for e in edges)
        test_import = sum(e["kind"] == "TEST_IMPORT" for e in edges)
        test_dynamic = sum(e["kind"] in {"TEST_LITERAL_MODULE_LOAD", "TEST_LITERAL_EXECUTION"} for e in edges)
        code_mentions = sum(r["kind"] == "CODE_REFERENCE" for r in refs)
        test_mentions = sum(r["kind"] == "TEST_REFERENCE" for r in refs)
        shell_mentions = sum(r["kind"] == "SHELL_REFERENCE" for r in refs)
        documentation_mentions = sum(r["kind"] == "DOC_REFERENCE" for r in refs)
        saved_mentions = sum(r["kind"] == "SAVED_RESULT_REFERENCE" for r in refs)
        config_mentions = sum(r["kind"] in {"CONFIG_REFERENCE", "TEXT_REFERENCE"} for r in refs)
        if name in declared:
            status = "RETAIN_DECLARED_WORKFLOW"
        elif name in closure:
            status = "RETAIN_TRANSITIVE_ACTIVE_DEPENDENCY"
        elif non_test_import or non_test_dynamic:
            status = "RETAIN_OTHER_PYTHON_DEPENDENCY"
        elif test_import or test_dynamic:
            status = "RETAIN_TEST_EXECUTION_OR_IMPORT"
        elif code_mentions:
            status = "REVIEW_PYTHON_TEXT_MENTION_ONLY"
        elif test_mentions:
            status = "REVIEW_TEST_TEXT_MENTION_ONLY"
        elif shell_mentions:
            status = "REVIEW_SHELL_ENTRYPOINT"
        elif documentation_mentions or config_mentions:
            status = "REVIEW_DOC_OR_CONFIG_ONLY"
        elif saved_mentions:
            status = "REVIEW_FROZEN_PROVENANCE_ONLY"
        else:
            status = "REVIEW_NO_DETECTED_REFERENCE"
        rows.append({
            "script": name,
            "project_path": old["project_path"],
            "v1_classification": old["classification"],
            "refined_classification": status,
            "in_declared_active_closure": name in closure,
            "python_imports": non_test_import,
            "python_literal_import_or_exec": non_test_dynamic,
            "test_imports": test_import,
            "test_literal_import_or_exec": test_dynamic,
            "code_text_mentions": code_mentions,
            "test_text_mentions": test_mentions,
            "shell_mentions": shell_mentions,
            "doc_mentions": documentation_mentions,
            "config_other_mentions": config_mentions,
            "saved_result_mentions": saved_mentions,
            "sha256": old["sha256"],
        })
    old_group = Counter(r["v1_classification"] for r in rows)
    new_group = Counter(r["refined_classification"] for r in rows)
    old_77 = [r for r in rows if r["v1_classification"] == "RETAIN_STATIC_CODE_OR_TEST_REFERENCE"]
    review_77 = [r for r in old_77 if r["refined_classification"].startswith("REVIEW_")]
    return {
        "format": "t3-script-dependency-refinement-v2",
        "head": base["head"],
        "scanned_text_files": base["scanned_text_files"],
        "scanned_python_files": len(pyfiles),
        "v1_categories": dict(old_group),
        "refined_categories": dict(new_group),
        "original_77_count": len(old_77),
        "original_77_review_count": len(review_77),
        "active_workflow_closure": sorted(closure),
        "rows": rows,
        "strong_evidence": sorted(strong_edges, key=lambda e: (e["script"], e["source"], e["line"], e["kind"])),
        "weak_evidence": first_refs,
        "warnings": [
            "Only static AST imports and literal import/exec calls are treated as strong executable dependencies.",
            "String-built dynamic imports, external shell usage and generated entrypoints may not be detected.",
            "REVIEW does not mean safe to archive or delete; provenance and historical reproducibility still require inspection.",
            "No physics code or stored results are modified; scanner writes reports only.",
        ],
    }


def render_summary(report: dict) -> str:
    lines = [
        "# Refined T3 script dependency audit (read-only)", "",
        f"- Checkout HEAD: `{report['head']}`",
        f"- Script modules: **{len(report['rows'])}**",
        f"- Python files parsed: **{report['scanned_python_files']}**",
        f"- Text files scanned by first pass: **{report['scanned_text_files']}**",
        f"- Originally in broad `RETAIN_STATIC_CODE_OR_TEST_REFERENCE`: **{report['original_77_count']}**",
        f"- From that group now needing *manual* review rather than presumed code/test retention: **{report['original_77_review_count']}**", "",
        "## All scripts by refined status", "", "| Status | Count |", "|---|---:|",
    ]
    for key, value in sorted(report["refined_categories"].items()):
        lines.append(f"| `{key}` | {value} |")
    lines += ["", "## Breakdown of the original 77", "", "| Refined status | Count |", "|---|---:|"]
    counts = Counter(r["refined_classification"] for r in report["rows"]
                     if r["v1_classification"] == "RETAIN_STATIC_CODE_OR_TEST_REFERENCE")
    for key, value in sorted(counts.items()):
        lines.append(f"| `{key}` | {value} |")
    lines += ["", "## Review candidates from the original 77", "",
              "These are **not** approved for archival. See `inventory.csv` and `strong_evidence.csv`.", "",
              "| Script | New status | Code text | Test text | Docs | Saved result |",
              "|---|---|---:|---:|---:|---:|"]
    subset = [r for r in report["rows"] if r["v1_classification"] == "RETAIN_STATIC_CODE_OR_TEST_REFERENCE"
              and r["refined_classification"].startswith("REVIEW_")]
    for r in subset:
        lines.append("| `{script}` | `{refined_classification}` | {code_text_mentions} | {test_text_mentions} | {doc_mentions} | {saved_result_mentions} |".format(**r))
    if not subset:
        lines.append("| None detected | — | — | — | — | — |")
    lines += ["", "## Declared-workflow dependency closure", "",
              "Only AST imports and literal execution calls contribute to the closure:", "",
              ", ".join(f"`{x}`" for x in report["active_workflow_closure"]), "",
              "## Evidence interpretation", "",
              "- `PYTHON_IMPORT` and `TEST_IMPORT`: syntax-level imports (including relative imports in scripts/).",
              "- `*_LITERAL_MODULE_LOAD`: literal importlib/runpy/__import__ module loads.",
              "- `*_LITERAL_EXECUTION`: literal script/module arguments of known execution APIs.",
              "- `CODE_REFERENCE` and `TEST_REFERENCE` from the first pass are **text mentions**, not evidence of an actual import.",
              "- `SHELL_REFERENCE` remains a manual executable entrypoint review, not proof of inactivity.",
              "- Saved-result references are preserved as provenance; do not rewrite historical JSON to remove them.",
              "- This tool does not move, delete, rename or execute any script. No `--apply` mode exists.", ""]
    return "\n".join(lines)


def run(root: Path) -> dict:
    root = root.resolve()
    base = examine(root)
    report = refined_report(root, base)
    out = root / OUT_REL
    out.mkdir(parents=True, exist_ok=True)
    (out / "audit.json").write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    with (out / "inventory.csv").open("w", encoding="utf-8", newline="") as h:
        writer = csv.DictWriter(h, fieldnames=list(report["rows"][0].keys()))
        writer.writeheader()
        writer.writerows(report["rows"])
    with (out / "strong_evidence.csv").open("w", encoding="utf-8", newline="") as h:
        writer = csv.DictWriter(h, fieldnames=["script", "source", "line", "kind", "excerpt"])
        writer.writeheader()
        writer.writerows(report["strong_evidence"])
    (out / "summary.md").write_text(render_summary(report), encoding="utf-8")
    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    args = parser.parse_args(argv)
    report = run(args.root)
    print(f"SCRIPT MODULES: {len(report['rows'])}")
    print(f"ORIGINAL BROAD GROUP: {report['original_77_count']}")
    print(f"NOW NEED MANUAL REVIEW: {report['original_77_review_count']}")
    for status, n in sorted(report["refined_categories"].items()):
        print(f"{status}: {n}")
    print("READ-ONLY: no scripts or scientific results moved, changed or executed.")
    print("REPORT:", args.root.resolve() / OUT_REL / "summary.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
