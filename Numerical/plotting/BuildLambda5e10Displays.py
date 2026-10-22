"""Rebuild 5e-10 saved-result overview and interactive RGE-running dashboards.

Run from the project root:
    python -m Numerical.plotting.BuildLambda5e10Displays

Does not integrate RGEs or optimise parameters. Requires saved successful
running_diagnostics.json from RunLambda5e10Comparisons.
"""
from __future__ import annotations

from pathlib import Path

from Numerical.plotting import BuildLambda5e9Displays as overview
from Numerical.plotting.InteractiveModelComparison import write_dashboard

ROOT = Path(__file__).resolve().parents[2]
STUDY = ROOT / 'output' / 'studies' / 'lambda_5e-10_refit' / 'comparison'
OUTPUT = ROOT / 'Reports' / 'output' / 'display_lambda_5e10'


def main() -> None:
    overview.STUDY = STUDY
    overview.DEFAULT_DEST = OUTPUT
    original_page = overview.page

    def corrected_page(title: str, contents: str) -> str:
        for old, new in (
            ('5×10⁻⁹', '5×10⁻¹⁰'),
            ('5×10<sup>−9</sup>', '5×10<sup>−10</sup>'),
            ('2.5×10<sup>−9</sup>', '2.5×10<sup>−10</sup>'),
        ):
            title = title.replace(old, new)
            contents = contents.replace(old, new)
        return original_page(title, contents)

    overview.page = corrected_page
    # Build overview through existing CLI with its default --output.
    # Avoid consuming arbitrary user arguments in the delegated parser.
    import sys
    if len(sys.argv) > 1:
        raise SystemExit('Usage: python -m Numerical.plotting.BuildLambda5e10Displays')
    overview.main()

    # Separately render trajectory dashboards, including model selectors.
    for mode, label in (
        ('optimised', 'Individually optimised'),
        ('shared_a_m2', 'Common A(-2) Yukawas'),
    ):
        source = STUDY / mode / 'models'
        if not source.is_dir():
            print(f'[SKIP] No models directory: {source}')
            continue
        destination = OUTPUT / (mode + '_rge.html')
        try:
            write_dashboard(source, destination, f'T3 5e-10: {label} RGE running')
        except FileNotFoundError as exc:
            print(f'[SKIP] {exc}')
        else:
            print('Interactive RGE dashboard:', destination)


if __name__ == '__main__':
    main()
