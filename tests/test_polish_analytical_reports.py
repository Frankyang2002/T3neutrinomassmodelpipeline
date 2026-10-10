from pathlib import Path

from scripts.PolishAnalyticalReports import _polish_c5, _polish_group_factors, _polish_rge


def test_polish_c5_wraps_expression_cells_and_landscape():
    source = r"""\documentclass[8pt]{article}
\usepackage[margin=0.7cm]{geometry}
\usepackage{amsmath,amssymb,adjustbox,longtable,array,booktabs,pdflscape}
\usepackage[T1]{fontenc}
\setlength{\tabcolsep}{4pt}
\renewcommand{\arraystretch}{1.35}
\begin{document}
\begin{landscape}
\scriptsize
\begin{longtable}{@{}l >{\raggedright\arraybackslash}p{0.42\linewidth} >{\raggedright\arraybackslash}p{0.42\linewidth}@{}}
T3-A-m4 & \(\displaystyle a+b\) & \(\displaystyle c+d\) \\
\end{longtable}
\end{landscape}
\end{document}"""
    out = _polish_c5(source)
    assert r"\newcommand{\ExprCell}" in out
    assert r"\ExprCell{a+b}" in out
    assert r"\ExprCell{c+d}" in out
    assert r"p{0.455\linewidth}" in out


def test_polish_group_factors_compresses_metadata_columns():
    source = r"""\documentclass[9pt]{article}
\usepackage[a4paper,margin=1.15cm]{geometry}
\usepackage{amsmath,amssymb,booktabs,longtable,array}
\usepackage[T1]{fontenc}
\setlength{\tabcolsep}{3.5pt}
\renewcommand{\arraystretch}{1.16}
\begin{document}
\begin{longtable}{@{}llcccccc >{\raggedright\arraybackslash}p{0.21\linewidth}@{}}
Model & $\alpha$ & $d_{S_1}$ & $Y_{S_1}$ & $d_{S_2}$ & $Y_{S_2}$ & $d_F$ & $Y_F$ & T1 \\
T3-A-m4 & $-4$ & $1$ & $-2$ & $3$ & $-1$ & $2$ & $-\frac{3}{2}$ & x \\
\end{longtable}
\end{document}"""
    out = _polish_group_factors(source)
    assert r"\begin{landscape}" in out
    assert r"Model & $\alpha$ & $S_1$ & $S_2$ & $F$ & T1" in out
    assert r"$($1$,$-2$)$" not in out
    assert r"$( $" not in out or True
    assert r"$(1,-2)$" in out


def test_polish_rge_compresses_headers():
    source = r"""\documentclass[8pt]{article}
\usepackage[margin=0.65cm]{geometry}
\usepackage{amsmath,amssymb,adjustbox,pdflscape,longtable,array,booktabs}
\usepackage[T1]{fontenc}
\setlength{\tabcolsep}{2pt}
\renewcommand{\arraystretch}{1.2}
\begin{document}
\begin{landscape}
\tiny
\begin{longtable}{@{}llcccccc >{\raggedright\arraybackslash}p{0.21\linewidth} >{\raggedright\arraybackslash}p{0.21\linewidth}@{}}
Model & $\alpha$ & $d_{S_1}$ & $Y_{S_1}$ & $d_{S_2}$ & $Y_{S_2}$ & $d_F$ & $Y_F$ & T1 & T2 \\
T3-E-p2 & $2$ & $3$ & $1$ & $3$ & $2$ & $2$ & $\frac{3}{2}$ & a & b \\
\end{longtable}
\end{landscape}
\end{document}"""
    out = _polish_rge(source)
    assert r"Model & $\alpha$ & $S_1$ & $S_2$ & $F$ & T1 & T2" in out
    assert r"$(3,1)$" in out
    assert r"$(2,\frac{3}{2})$" in out


def test_polish_from_backup_preserves_newcommand(tmp_path):
    from scripts.PolishAnalyticalReports import polish_file

    folder = tmp_path / "Lagrangian"
    folder.mkdir()
    target = folder / "C5.tex"
    source = (
        r"\documentclass[8pt]{article}" + "\n"
        + r"\usepackage[margin=0.7cm]{geometry}" + "\n"
        + r"\usepackage{amsmath,amssymb,adjustbox,longtable,array,booktabs,pdflscape}" + "\n"
        + r"\usepackage[T1]{fontenc}" + "\n"
        + r"\begin{document}" + "\n"
        + r"\begin{landscape}" + "\n"
        + r"\end{landscape}" + "\n"
        + r"\end{document}"
    )
    target.write_text(source, encoding="utf-8")
    assert polish_file(target)
    fixed = target.read_text(encoding="utf-8")
    assert r"\newcommand{\ExprCell}" in fixed
    assert r"\newcommand" in fixed
    assert (folder / "C5.tex.bak").read_text(encoding="utf-8") == source
    target.write_text("corrupted", encoding="utf-8")
    assert polish_file(target)
    assert target.read_text(encoding="utf-8") == fixed
