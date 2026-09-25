"""Documentation math should use Markdown-rendered dollar delimiters."""

from __future__ import annotations

from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def test_all_repository_markdown_uses_dollar_math_delimiters() -> None:
    markdown_files = sorted(PROJECT_ROOT.glob("*.md"))
    assert {path.name for path in markdown_files} == {
        "INFO_LAGRANGIAN.md",
        "INFO_PIPELINE.md",
        "INFO_RGE.md",
    }

    for path in markdown_files:
        source = path.read_text(encoding="utf-8")
        assert r"\[" not in source, path
        assert r"\]" not in source, path
        assert r"\(" not in source, path
        assert r"\)" not in source, path
