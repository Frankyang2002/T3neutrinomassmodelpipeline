"""Canonical consolidated guides use consistent paired Markdown/LaTeX delimiters.

Both dollar-delimited maths and paired \( \) / \[ \] are accepted, because the
consolidated physics documents use the latter as well. A mismatched pair is a
formatting regression.
"""
from __future__ import annotations

from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
CANONICAL_GUIDES = ("README.md", "docs/PHYSICS.md", "docs/NUMERICAL.md", "docs/RESULTS.md")


def test_all_repository_markdown_uses_dollar_math_delimiters() -> None:
    for name in CANONICAL_GUIDES:
        path = PROJECT_ROOT / name
        assert path.is_file(), f"Missing canonical guide: {name}"
        source = path.read_text(encoding="utf-8-sig")
        for start, end in ((r"\[", r"\]"), (r"\(", r"\)")):
            assert source.count(start) == source.count(end), (name, start, end)
