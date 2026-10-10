"""Regenerate a legible, C5-only comparison from the saved Matchete report.

Never modifies physics calculations. Refuses to combine nonidentical kernels.
"""
from __future__ import annotations

import argparse
import re
import shutil
import subprocess
from pathlib import Path

FLAVOUR = r"\overline{y_{1,i,r}} \overline{y_{2,j,r}}"
MODEL_RE = re.compile(r"^T3-[A-E]-(?:m|p)\d+$")
PREFIX = {"A": ("+", r"\sqrt{3}"), "B": ("+", r"\sqrt{3}"),
          "C": ("-", ""), "D": ("+", r"\sqrt{3}"),
          "E": ("-", r"\sqrt{6}")}


def parse_ordered(source: str) -> tuple[list[tuple[str, str]], str]:
    rows: list[tuple[str, str]] = []
    kernels: set[str] = set()
    for line in source.splitlines():
        if not line.startswith("T3-") or " & " not in line:
            continue
        cells = line.rsplit(r" \\", 1)[0].split(" & ")
        if len(cells) != 3:
            raise ValueError("Expected original three-column Matchete C5 table")
        model, ordered, _symmetric = cells
        if not MODEL_RE.fullmatch(model):
            raise ValueError(f"Unexpected model name: {model}")
        cls = model.split("-")[1]
        sign, radical = PREFIX[cls]
        if not ordered.startswith(r"\(\displaystyle A_5^{ij}=") or not ordered.endswith(r"\)"):
            raise ValueError(f"Missing original ordered Matchete expression for {model}")
        body = ordered[len(r"\(\displaystyle A_5^{ij}="):-2].strip()
        leading = ("-" if sign == "-" else "") + r"\frac{2}{3}"
        if not body.startswith(leading):
            raise ValueError(f"Unexpected overall class factor for {model}")
        body = body[len(leading):].strip()
        if body.count(FLAVOUR) != 18:
            raise ValueError(f"Unexpected flavour structure for {model}")
        if radical:
            if body.count(radical) != 18:
                raise ValueError(f"Unexpected square-root factor for {model}")
            body = body.replace(radical + " ", "")
        if r"\sqrt{3}" in body or r"\sqrt{6}" in body:
            raise ValueError(f"Unfactored radical for {model}")
        body = body.replace(FLAVOUR, "").strip()
        kernels.add(body)
        rows.append((model, cls))
    if len(rows) != 16 or len({m for m, _ in rows}) != 16:
        raise ValueError(f"Expected exactly 16 models, found {len(rows)}")
    if len(kernels) != 1:
        raise ValueError("Mass-dependent Matchete kernels are not identical; refusing to combine")
    return rows, kernels.pop()


def build_tex(rows: list[tuple[str, str]]) -> str:
    lines = [
        r"\documentclass[10pt,a4paper]{article}",
        r"\usepackage[margin=2.0cm]{geometry}",
        r"\usepackage{amsmath,amssymb,booktabs,longtable,array}",
        r"\usepackage[T1]{fontenc}",
        r"\setlength{\parindent}{0pt}",
        r"\renewcommand{\arraystretch}{1.18}",
        r"\begin{document}",
        r"\section*{T3: physical symmetric Weinberg coefficient $C_5$}",
        r"The physical coefficient is displayed without the redundant ordered Matchete column.",
        r"For all 16 models the saved ordered coefficients have the same mass-dependent",
        r"structure after removing the flavour tensor and exact class prefactor:",
        r"\[",
        r" C_{5,m}^{ij}=\kappa_m\,\mathcal F_{\rm match}(M_F,m_1,m_2,\bar\mu,\epsilon)\,",
        r" \left(\overline{y_{1,i,r}}\,\overline{y_{2,j,r}}+",
        r" \overline{y_{1,j,r}}\,\overline{y_{2,i,r}}\right).",
        r"\]",
        r"The repeated fermion index $r$ is summed as in the saved output.",
        r"The common factor $\mathcal F_{\rm match}$ includes $\lambda_5$, $\hbar$, the",
        r"Matchete integral $I_3$, and its original mass and regulator dependence.",
        r"It is defined exactly by the common kernel in",
        r"\texttt{C5\_common\_kernel.txt}; no loop-integral simplification or",
        r"renormalisation has been assumed here.",
        r"\medskip",
        r"\begin{longtable}{@{}lllc@{}}",
        r"\toprule",
        r"Model & Class & $(d_{S_1},d_{S_2},d_F)$ & $\kappa_m$ \\",
        r"\midrule\endfirsthead",
        r"\toprule",
        r"Model & Class & $(d_{S_1},d_{S_2},d_F)$ & $\kappa_m$ \\",
        r"\midrule\endhead",
    ]
    dims = {"A": "(1,3,2)", "B": "(2,2,1)", "C": "(2,2,3)",
            "D": "(3,1,2)", "E": "(3,3,2)"}
    for model, cls in rows:
        sign, radical = PREFIX[cls]
        coefficient = ("-" if sign == "-" else "") + rf"\frac{{2{radical}}}{{3}}" if radical else ("-" if sign == "-" else "") + r"\frac{2}{3}"
        lines.append(f"{model} & {cls} & ${dims[cls]}$ & ${coefficient}$ " + r"\\")
    lines += [
        r"\bottomrule\end{longtable}",
        r"\paragraph{Scope.} The prefactors and common kernel are obtained",
        r"by literal factor extraction from the previously generated Matchete",
        r"expressions, not from an independent matching calculation.",
        r"The full original expression remains preserved separately.",
        r"\end{document}", "",
    ]
    return "\n".join(lines)


def run(root: Path, compile_pdf: bool) -> None:
    root.mkdir(parents=True, exist_ok=True)
    tex = root / "C5.tex"
    candidates = [root / "C5_original_matchete.tex", root / "C5.tex.bak", tex]
    original = next((p for p in candidates if p.is_file() and "Ordered coefficient extracted from Matchete" in p.read_text(encoding="utf-8")), None)
    if original is None:
        raise FileNotFoundError("Cannot find original three-column C5.tex; restore the generated source first")
    source = original.read_text(encoding="utf-8")
    rows, kernel = parse_ordered(source)
    archived = root / "C5_original_matchete.tex"
    if not archived.exists():
        archived.write_text(source, encoding="utf-8")
    (root / "C5_common_kernel.txt").write_text(
        "Exact common Matchete kernel, after removing the ordered flavour product, "
        "class overall prefactor, and per-term class radical.\n"
        "The C5 report multiplies this kernel by kappa_m and the symmetric flavour combination.\n\n"
        + kernel + "\n", encoding="utf-8",
    )
    tex.write_text(build_tex(rows), encoding="utf-8")
    print(f"Wrote {tex} ({len(rows)} model rows; identical kernels verified)")
    if compile_pdf:
        engine = shutil.which("pdflatex")
        if not engine:
            raise RuntimeError("pdflatex was not found on PATH")
        result = subprocess.run([engine,"-interaction=nonstopmode","-halt-on-error",tex.name],
                                cwd=root, capture_output=True, text=True)
        if result.returncode:
            print(result.stdout[-6000:])
            raise RuntimeError("C5 compilation failed")
        print(f"Compiled {tex.with_suffix('.pdf')}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("Reports/output/analytical/Lagrangian"))
    parser.add_argument("--compile", action="store_true")
    args = parser.parse_args()
    run(args.root, args.compile)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
