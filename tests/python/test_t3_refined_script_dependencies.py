"""Offline synthetic tests of the local read-only script dependency refinement."""
from __future__ import annotations

from pathlib import Path
import tempfile

from tools.RefineScriptDependencyAudit import classify_python_evidence, refined_report, run


def test_ast_import_versus_plain_mentions():
    inventory = {"Alpha", "Beta", "Gamma", "Delta"}
    src = '''# from scripts import Delta\nDOC = "scripts.Gamma"\nfrom scripts import Alpha\nimport scripts.Beta\nimport importlib\nmod = importlib.import_module("scripts.Gamma")\n'''
    e, error = classify_python_evidence(src, "Numerical/demo.py", inventory)
    assert error is None
    assert {(x["script"], x["kind"]) for x in e} == {
        ("Alpha", "PYTHON_IMPORT"),
        ("Beta", "PYTHON_IMPORT"),
        ("Gamma", "PYTHON_LITERAL_MODULE_LOAD"),
    }


def test_relative_import_and_explicit_runtime_launch():
    src = '''from . import Alpha\nimport subprocess\nsubprocess.run(["python", "scripts/Beta.py"])\nsubprocess.run(["python", variable_name])\n'''
    e, error = classify_python_evidence(src, "scripts/Gamma.py", {"Alpha", "Beta", "Gamma"})
    assert error is None
    assert {(x["script"], x["kind"]) for x in e} == {
        ("Alpha", "PYTHON_IMPORT"), ("Beta", "PYTHON_LITERAL_EXECUTION")
    }


def test_test_import_classification():
    e, err = classify_python_evidence('from scripts import Gamma\n', "tests/test_example.py", {"Gamma"})
    assert err is None
    assert [(x["script"], x["kind"]) for x in e] == [("Gamma", "TEST_IMPORT")]


def test_synthetic_full_scan_read_only_and_transitive_closure():
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        (root / "scripts").mkdir()
        (root / "tests").mkdir()
        (root / "docs").mkdir()
        (root / "tools").mkdir()
        # The first pass declares 11 real workflow names; the synthetic checkout
        # exercises a detected source->source closure without relying on them.
        (root / "scripts" / "Alpha.py").write_text("from scripts import Beta\n", encoding="utf-8")
        (root / "scripts" / "Beta.py").write_text("print('working')\n", encoding="utf-8")
        (root / "scripts" / "Gamma.py").write_text("print('working')\n", encoding="utf-8")
        (root / "tests" / "test_one.py").write_text("from scripts import Gamma\n", encoding="utf-8")
        (root / "docs" / "old.md").write_text("Run scripts/Alpha.py manually\n", encoding="utf-8")
        before = {p.relative_to(root).as_posix(): p.read_bytes() for p in root.rglob("*") if p.is_file()}
        result = run(root)
        row = {r["script"]: r for r in result["rows"]}
        assert row["Beta"]["refined_classification"] == "RETAIN_OTHER_PYTHON_DEPENDENCY"
        assert row["Gamma"]["refined_classification"] == "RETAIN_TEST_EXECUTION_OR_IMPORT"
        assert row["Alpha"]["refined_classification"] == "REVIEW_DOC_OR_CONFIG_ONLY"
        assert result["original_77_count"] == 2
        assert (root / "output/cleanup_audit/script_dependencies_refined_v2/summary.md").is_file()
        for path, data in before.items():
            assert (root / path).read_bytes() == data


def test_active_workflow_transitive_closure():
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        (root / "scripts").mkdir()
        (root / "tests").mkdir()
        (root / "scripts" / "RunT3ReducedScalarSobolComparison.py").write_text(
            "from scripts import Alpha\n", encoding="utf-8")
        (root / "scripts" / "Alpha.py").write_text(
            "from scripts import Beta\n", encoding="utf-8")
        (root / "scripts" / "Beta.py").write_text(
            "DATA = 1\n", encoding="utf-8")
        report = run(root)
        rows = {r["script"]: r for r in report["rows"]}
        assert rows["RunT3ReducedScalarSobolComparison"]["refined_classification"] == "RETAIN_DECLARED_WORKFLOW"
        assert rows["Alpha"]["refined_classification"] == "RETAIN_TRANSITIVE_ACTIVE_DEPENDENCY"
        assert rows["Beta"]["refined_classification"] == "RETAIN_TRANSITIVE_ACTIVE_DEPENDENCY"
        assert set(report["active_workflow_closure"]) == {
            "RunT3ReducedScalarSobolComparison", "Alpha", "Beta"}


def test_plain_test_strings_do_not_become_executable_dependencies():
    inventory = {"Alpha"}
    e, err = classify_python_evidence(
        '"""Run scripts/Alpha.py manually"""\nPATH = "scripts/Alpha.py"\n',
        "tests/test_doc.py", inventory)
    assert err is None
    assert e == []
