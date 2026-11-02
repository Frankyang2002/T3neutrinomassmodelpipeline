"""Restore the established coloured T3 RGE report layout from saved diagnostics.

Presentation-only wrapper around the repository's existing
Numerical.plotting.InteractiveModelComparison.write_dashboard.
Does not integrate RGEs, match Wilson coefficients, or manufacture LLSS curves.
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

from Numerical.plotting.InteractiveModelComparison import write_dashboard

ROOT = Path(__file__).resolve().parents[1]
SCENARIOS = ("smallY_smallL", "smallY_largeL", "largeY_smallL", "largeY_largeL")


def regroup_dashboard(html_text: str) -> str:
    """Change only the variable selector; preserve the existing coloured chart."""
    groups = re.findall(r'<optgroup label="[^"]*">.*?</optgroup>', html_text, re.S)
    if len(groups) != 3:
        raise ValueError("Unexpected legacy RGE selector. Refusing to edit report HTML.")
    prefix, llss_group, others = groups
    if 'Weinberg operator' not in prefix or 'LLSS' not in llss_group or 'Other running observables' not in others:
        raise ValueError("Unexpected RGE dashboard sections; source may have changed.")
    opt = re.compile(r'<option\s+value="([^"]+)"[^>]*>.*?</option>', re.S)
    uv, obs = [], []
    for item in opt.finditer(others):
        if item.group(1) in {'y1norm', 'y2norm', 'lambdaT3'}:
            uv.append(item.group(0))
        else:
            obs.append(item.group(0))
    if len(uv) != 3 or len(obs) < 5:
        raise ValueError("Expected UV Yukawa/scalar and low-energy observable selectors were not found.")
    new_groups = (
        prefix.replace('Weinberg operator — complete matching history',
                       'C5 running — EFT1 continuation and SM')
        + llss_group.replace('Weinberg-like LLSS — analytic diagnostic',
                             'CLLSS in EFT1 — analytic kernel only (not Wilson trajectory)')
        + '<optgroup label="UV — Yukawa and scalar running">' + ''.join(uv) + '</optgroup>'
        + '<optgroup label="Observables">' + ''.join(obs) + '</optgroup>'
    )
    original = ''.join(groups)
    if original not in html_text:
        raise ValueError("Dashboard optgroups are not contiguous; source layout changed.")
    result = html_text.replace(original, new_groups, 1)
    if result == html_text:
        raise RuntimeError('Regrouping did not change the original HTML')
    return result


def build_one(study_dir: Path, out_file: Path, title: str) -> Path:
    # Existing upstream writer provides identical established plots, colours,
    # region shading, matched-continuation warnings, and saved source data.
    out_file.parent.mkdir(parents=True, exist_ok=True)
    write_dashboard(study_dir, out_file, title)
    original = out_file.read_text(encoding='utf-8')
    restored = regroup_dashboard(original)
    out_file.write_text(restored, encoding='utf-8')
    return out_file


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--scenario', choices=(*SCENARIOS, 'all'), default='all')
    ap.add_argument('--input-root', type=Path, default=ROOT / 'output/full/comparison')
    ap.add_argument('--output-root', type=Path, default=ROOT / 'Reports/output/display/comparison')
    args = ap.parse_args()
    for scenario in SCENARIOS if args.scenario == 'all' else (args.scenario,):
        source = args.input_root / scenario
        output = args.output_root / scenario / 'interactive_comparison.html'
        print(build_one(source, output, f'T3 normalized comparison: {scenario}'))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
