"""Assemble the two verified SM+Weinberg running figures into one vector PDF.

Usage from repository root:
    python -m Numerical.analysis.MakeT3ThesisMassMixingPanel

No RGE recalculation. No modification of physics outputs.
"""
from pathlib import Path
import fitz

ROOT = Path(__file__).resolve().parents[2]
FIGURES = ROOT / 'Reports' / 'output' / 'thesis_assets' / 'figures'

def main() -> None:
    first = FIGURES / 'fig08_t3b_weinberg_mass_splitting_running.pdf'
    second = FIGURES / 'fig09_t3b_weinberg_angle_running.pdf'
    for path in (first, second):
        if not path.is_file():
            raise FileNotFoundError(path)
    with fitz.open(first) as a, fitz.open(second) as b:
        out = fitz.open()
        page = out.new_page(width=934, height=370)
        page.show_pdf_page(fitz.Rect(16, 14, 466, 332), a, 0,
                           keep_proportion=True)
        page.show_pdf_page(fitz.Rect(475, 14, 925, 342), b, 0,
                           keep_proportion=True)
        destination = FIGURES / 'fig10_mass_vs_mixing_running_combined.pdf'
        out.save(destination, garbage=3, deflate=True)
        out[0].get_pixmap(matrix=fitz.Matrix(2, 2), alpha=False).save(
            FIGURES / 'fig10_mass_vs_mixing_running_combined.png')
        out.close()
    print(f'Written: {destination}')

if __name__ == '__main__':
    main()
