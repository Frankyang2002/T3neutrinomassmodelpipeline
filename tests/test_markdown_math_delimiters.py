"""Canonical project guides use Markdown dollar math delimiters.

Supplementary technical notes may legitimately use other LaTeX delimiters;
this test intentionally covers only the three canonical guides.
"""

from __future__ import annotations

from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
CANONICAL_GUIDES = (
    "INFO_LAGRANGIAN.md",
    "INFO_PIPELINE.md",
    "INFO_RGE.md",
)


def test_all_repository_markdown_uses_dollar_math_delimiters() -> None:
    for name in CANONICAL_GUIDES:
        path = PROJECT_ROOT / name
        assert path.is_file(), f"Missing canonical guide: {name}"
        source = path.read_text(encoding="utf-8-sig")
        for obsolete in (r"\[", r"\]", r"\(", r"\)"):
            assert obsolete not in source, f"{name} contains {obsolete!r}"
