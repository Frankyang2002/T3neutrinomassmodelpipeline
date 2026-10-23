from __future__ import annotations

from Reports.C12GroupFactorPDF import _extract_c12_sections, transform_group_factor_tex


SAMPLE_TEX = r'''
\documentclass{article}
\usepackage{booktabs}
\begin{document}
\section*{Existing C5 tables}
Some C5 content here.
\subsection*{T3-A-p0}
Representation: $(d_{S_1}, d_{S_2}, d_F, \alpha) = (1,3,2,0)$. Matched independent $C_{12}$ components: 20.
\begin{tabular}{lll}
\toprule
One-loop source & Group factor & Alignment \\
\midrule
Yukawa wavefunction & 0 & yes \\
Scalar-pair gauge/quartic & -2\lambda_{12,\mathrm{Inv1}} & yes \\
Mixed Yukawa and gauge & not proportional to C_{12} & no \\
Fermion anomalous dimension & -\frac{3 g_1^2}{2} - \frac{9 g_2^2}{2} + y_e \overline{y_e} & yes \\
\bottomrule
\end{tabular}
\subsection*{T3-B-m3}
Representation: $(d_{S_1}, d_{S_2}, d_F, \alpha) = (2,2,1,-3)$. Matched independent $C_{12}$ components: 48.
\begin{tabular}{lll}
\toprule
One-loop source & Group factor & Alignment \\
\midrule
Yukawa wavefunction & 0 & yes \\
Scalar anomalous dimension & -3 g_2^2 & yes \\
\bottomrule
\end{tabular}
\subsection*{T3-B-m1}
Representation: $(d_{S_1}, d_{S_2}, d_F, \alpha) = (2,2,1,-1)$. Matched independent $C_{12}$ components: 48.
\begin{tabular}{lll}
\toprule
One-loop source & Group factor & Alignment \\
\midrule
Yukawa wavefunction & 0 & yes \\
Scalar anomalous dimension & -3 g_2^2 - 4\sqrt{3}\lambda_{12Inv2}/3 & yes \\
\bottomrule
\end{tabular}
\subsection*{T3-B-p1}
Representation: $(d_{S_1}, d_{S_2}, d_F, \alpha) = (2,2,1,1)$. Matched independent $C_{12}$ components: 48.
\begin{tabular}{lll}
\toprule
One-loop source & Group factor & Alignment \\
\midrule
Yukawa wavefunction & 0 & yes \\
Scalar anomalous dimension & -3 g_2^2 & yes \\
\bottomrule
\end{tabular}
\end{document}
'''


def test_extracts_existing_c12_sections():
    sections = _extract_c12_sections(SAMPLE_TEX)
    assert [s.model for s in sections] == ["T3-A-p0", "T3-B-m3", "T3-B-m1", "T3-B-p1"]
    assert sections[0].matched_components == 20
    assert any(label == "Scalar-pair gauge/quartic" for label, _, _ in sections[0].rows)


def test_transform_replaces_individual_sections_with_run_block_tables():
    out = transform_group_factor_tex(SAMPLE_TEX)
    assert "{12}$ representation-specific group factors" in out
    assert "Run block 1: T3-A-p0--T3-A-p0" in out
    assert "Run block 2: T3-B-m3--T3-B-p1" in out
    assert "T3 -- B -- m3" in out
    assert "T3 -- B -- m1" in out
    assert "T3 -- B -- p1" in out
    assert "\\textit{not proportional}" in out
    # Original per-model section headings should be gone from the rewritten body.
    assert out.count("\\subsection*{T3-A-p0}") == 0
    assert out.count("\\subsection*{T3-B-m3}") == 0


def test_existing_non_c12_content_is_preserved():
    out = transform_group_factor_tex(SAMPLE_TEX)
    assert "Existing C5 tables" in out
    assert "Some C5 content here." in out
    assert out.strip().endswith(r"\end{document}")
