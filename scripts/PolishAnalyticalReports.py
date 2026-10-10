from __future__ import annotations

import argparse
import re
import shutil
import subprocess
from pathlib import Path


def _replace_once(text: str, old: str, new: str) -> str:
    if old not in text:
        return text
    return text.replace(old, new, 1)


def _wrap_landscape_document(text: str) -> str:
    if r"\begin{landscape}" in text:
        return text
    text = _replace_once(text, r"\begin{document}", "\\begin{document}\n\\begin{landscape}")
    text = text.replace(r"\end{document}", "\\end{landscape}\n\\end{document}")
    return text


def _modernise_general_preamble(
    text: str,
    *,
    landscape: bool = False,
    margin: str = "0.8cm",
    font_size: str = "8pt",
) -> str:
    text = re.sub(
        r"\\documentclass\[[^\]]*\]\{article\}",
        lambda _m: f"\\documentclass[{font_size}]{{article}}",
        text,
        count=1,
    )
    text = re.sub(
        r"\\usepackage\[[^\]]*\]\{geometry\}",
        lambda _m: f"\\usepackage[a4paper,margin={margin}]{{geometry}}",
        text,
        count=1,
    )
    package_line = (
        r"\usepackage{amsmath,amssymb,mathtools,booktabs,longtable,array,"
        r"adjustbox,pdflscape,microtype,ragged2e}"
    )
    if "amsmath" in text and "usepackage" in text:
        text = re.sub(
            r"\\usepackage\{[^\n]*amsmath[^\n]*\}",
            lambda _m: package_line,
            text,
            count=1,
        )
    else:
        text = _replace_once(text, r"\usepackage[T1]{fontenc}", package_line + "\n\\usepackage[T1]{fontenc}")

    insert_lines = [
        r"\setlength{\parindent}{0pt}",
        r"\setlength{\LTleft}{0pt}",
        r"\setlength{\LTright}{0pt}",
        r"\allowdisplaybreaks",
    ]
    for line in insert_lines:
        if line not in text:
            text = _replace_once(text, r"\begin{document}", line + "\n\\begin{document}")

    if landscape:
        text = _wrap_landscape_document(text)
    return text


def _split_data_row(line: str) -> list[str] | None:
    if " & " not in line or not line.rstrip().endswith(r"\\"):
        return None
    body = line.rstrip()[:-2].strip()
    return [part.strip() for part in body.split(" & ")]


def _strip_math(token: str) -> str:
    token = token.strip()
    if token.startswith("$") and token.endswith("$"):
        return token[1:-1].strip()
    return token


def _compress_model_metadata_tables(text: str) -> str:
    lines = text.splitlines()
    out: list[str] = []
    for line in lines:
        if r"\begin{longtable}{@{}llcccccc" in line:
            line = (
                r"\begin{longtable}{@{}lcccc "
                r">{\RaggedRight\arraybackslash}p{0.36\linewidth} "
                r">{\RaggedRight\arraybackslash}p{0.36\linewidth} "
                r">{\RaggedRight\arraybackslash}p{0.36\linewidth} "
                r">{\RaggedRight\arraybackslash}p{0.36\linewidth} "
                r">{\RaggedRight\arraybackslash}p{0.36\linewidth}@{}}"
            )
        line = line.replace(
            r"Model & $\alpha$ & $d_{S_1}$ & $Y_{S_1}$ & $d_{S_2}$ & $Y_{S_2}$ & $d_F$ & $Y_F$ & ",
            r"Model & $\alpha$ & $S_1$ & $S_2$ & $F$ & ",
        )
        parts = _split_data_row(line)
        if parts and len(parts) >= 9 and parts[0].startswith("T3-"):
            fixed = [
                parts[0],
                parts[1],
                f"$({_strip_math(parts[2])},{_strip_math(parts[3])})$",
                f"$({_strip_math(parts[4])},{_strip_math(parts[5])})$",
                f"$({_strip_math(parts[6])},{_strip_math(parts[7])})$",
            ] + parts[8:]
            line = " & ".join(fixed) + r" \\" 
        out.append(line)
    return "\n".join(out)


def _polish_c5(text: str) -> str:
    text = _modernise_general_preamble(text, landscape=True, margin="0.75cm", font_size="7pt")
    text = text.replace(r"\setlength{\tabcolsep}{4pt}", r"\setlength{\tabcolsep}{2.5pt}")
    text = text.replace(r"\renewcommand{\arraystretch}{1.35}", r"\renewcommand{\arraystretch}{1.18}")
    macro = (
        r"\newcommand{\ExprCell}[1]{"
        r"\parbox[t]{\linewidth}{\vspace{0pt}\centering"
        r"\adjustbox{max width=\linewidth,max totalheight=0.78\textheight,"
        r"keepaspectratio}{$\displaystyle #1$}}}"
    )
    if macro not in text:
        text = _replace_once(text, r"\begin{document}", macro + "\n\\begin{document}")
    text = text.replace(
        r"\begin{longtable}{@{}l >{\raggedright\arraybackslash}p{0.42\linewidth} >{\raggedright\arraybackslash}p{0.42\linewidth}@{}}",
        r"\begin{longtable}{@{}l >{\RaggedRight\arraybackslash}p{0.455\linewidth} >{\RaggedRight\arraybackslash}p{0.455\linewidth}@{}}",
    )
    lines: list[str] = []
    for line in text.splitlines():
        parts = _split_data_row(line)
        if parts and len(parts) == 3 and parts[0].startswith("T3-"):
            new_parts = [parts[0]]
            for expr in parts[1:]:
                expr = expr.strip()
                if expr.startswith(r"\(\displaystyle ") and expr.endswith(r"\)"):
                    expr = expr[len(r"\(\displaystyle ") : -len(r"\)")]
                elif expr.startswith(r"\(") and expr.endswith(r"\)"):
                    expr = expr[2:-2]
                new_parts.append(rf"\ExprCell{{{expr}}}")
            line = " & ".join(new_parts) + r" \\" 
        lines.append(line)
    return "\n".join(lines)


def _polish_group_factors(text: str) -> str:
    text = _modernise_general_preamble(text, landscape=True, margin="0.9cm", font_size="8pt")
    text = text.replace(r"\setlength{\tabcolsep}{3.5pt}", r"\setlength{\tabcolsep}{2.5pt}")
    text = text.replace(r"\renewcommand{\arraystretch}{1.16}", r"\renewcommand{\arraystretch}{1.10}")
    text = _compress_model_metadata_tables(text)
    text = text.replace(r"\tiny", r"\scriptsize")
    return text


def _polish_rge(text: str) -> str:
    text = _modernise_general_preamble(text, landscape=True, margin="0.7cm", font_size="7pt")
    text = text.replace(r"\setlength{\tabcolsep}{2pt}", r"\setlength{\tabcolsep}{1.6pt}")
    text = text.replace(r"\renewcommand{\arraystretch}{1.2}", r"\renewcommand{\arraystretch}{1.10}")
    text = _compress_model_metadata_tables(text)
    return text


def _normalise_latex_control_sequences(text: str) -> str:
    return re.sub(r"\\\\(?=[A-Za-z])", r"\\", text)


def polish_file(path: Path) -> bool:
    backup = path.with_suffix(path.suffix + ".bak")
    original = backup.read_text(encoding="utf-8") if backup.exists() else path.read_text(encoding="utf-8")
    if path.name == "C5.tex":
        updated = _polish_c5(original)
    elif "GroupFactors" in path.parts:
        updated = _polish_group_factors(original)
    elif "RGE" in path.parts:
        updated = _polish_rge(original)
    else:
        return False
    updated = _normalise_latex_control_sequences(updated)
    if updated == original:
        return False
    if not backup.exists():
        backup.write_text(original, encoding="utf-8")
    path.write_text(updated, encoding="utf-8")
    return True


def compile_tex(path: Path) -> None:
    engine = shutil.which("pdflatex")
    if not engine:
        print(f"[skip] pdflatex not found: {path}")
        return
    for _ in range(2):
        completed = subprocess.run(
            [engine, "-interaction=nonstopmode", "-halt-on-error", path.name],
            cwd=path.parent,
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        if completed.returncode:
            diagnostic = "\n".join(completed.stdout.splitlines()[-35:])
            raise RuntimeError(
                f"LaTeX compilation failed for {path}.\n"
                f"Log: {path.with_suffix('.log')}\n{diagnostic}"
            )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Polish combined analytical T3 LaTeX reports in place.")
    parser.add_argument("--root", type=Path, default=Path("Reports/output/analytical"))
    parser.add_argument("--compile", action="store_true")
    args = parser.parse_args(argv)

    tex_files = sorted(args.root.rglob("*.tex"))
    if not tex_files:
        raise SystemExit(f"No .tex files found under {args.root}")

    touched: list[Path] = []
    for path in tex_files:
        if polish_file(path):
            touched.append(path)
            print(f"[updated] {path}")
            if args.compile:
                compile_tex(path)
    print(f"Updated {len(touched)} files.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
