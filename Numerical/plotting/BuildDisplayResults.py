"""Regenerate a compact comparison display from saved numerical diagnostics.

Presentation only: does not run RGEs, matching, fits, or symbolic tools.
The --report-input argument is retained for CLI/API compatibility but is not read.
The source of truth for this display is output/full/comparison/*/running_diagnostics.json.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
import shutil
from pathlib import Path
from typing import Any

from Numerical.plotting.InteractiveModelComparison import write_dashboard
from Numerical.plotting.RunningResultFigures import generate_running_result_figures

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DATA_INPUT = PROJECT_ROOT / 'output' / 'full'
DEFAULT_REPORT_INPUT = PROJECT_ROOT / 'Reports' / 'output' / 'full'  # compatibility only
DEFAULT_OUTPUT = PROJECT_ROOT / 'Reports' / 'output' / 'display'

KEEP_FIGURES = (
    'uv_coupling_running.png', 'intermediate_direct_weinberg_running.png',
    'c5_threshold_contributions.png', 'c5_running.png',
    'neutrino_mass_splitting_running.png', 'neutrino_mixing_running.png',
)
SCENARIOS = {
    'smallY_smallL': (0.005, 0.1), 'smallY_largeL': (0.005, 1.0),
    'largeY_smallL': (0.5, 0.1), 'largeY_largeL': (0.5, 1.0),
}
SCENARIO_ORDER = tuple(SCENARIOS)
MODEL_CLASS_BY_DIMS = {
    (1, 3, 2): 'A', (2, 2, 1): 'B', (2, 2, 3): 'C',
    (3, 1, 2): 'D', (3, 3, 2): 'E',
}
LAMBDA_T3_NORMALISATION = {
    'A': 1.0, 'B': 1.0, 'C': math.sqrt(3.0),
    'D': 1.0, 'E': 1.0 / math.sqrt(2.0),
}


def _load_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding='utf-8-sig'))
    if not isinstance(payload, dict):
        raise ValueError(f'{path}: expected JSON object')
    return payload


def _model_name(path: Path, payload: dict[str, Any]) -> str:
    search = payload.get('benchmark_search', {})
    model = search.get('model') if isinstance(search, dict) else None
    if isinstance(model, dict) and model.get('model_key'):
        return str(model['model_key'])
    for part in reversed(path.parts):
        if part.startswith('T3_dS1_'):
            return part
    return path.parent.parent.name


def _model_class(path: Path, payload: dict[str, Any]) -> str:
    search = payload.get('benchmark_search', {})
    model = search.get('model') if isinstance(search, dict) else None
    if isinstance(model, dict):
        try:
            dimensions = (int(model['d_s1']), int(model['d_s2']), int(model['d_f']))
        except (KeyError, TypeError, ValueError):
            pass
        else:
            if dimensions in MODEL_CLASS_BY_DIMS:
                return MODEL_CLASS_BY_DIMS[dimensions]
    match = re.search(r'T3_dS1_(\d+)_dS2_(\d+)_dF_(\d+)_', _model_name(path, payload))
    if match:
        dimensions = tuple(map(int, match.groups()))
        if dimensions in MODEL_CLASS_BY_DIMS:
            return MODEL_CLASS_BY_DIMS[dimensions]
    raise ValueError(f'Cannot determine representation for {path}')


def _last(value: Any) -> Any:
    return value[-1] if isinstance(value, list) and value else ''


def _summary_row(scenario: str, path: Path, payload: dict[str, Any]) -> dict[str, Any]:
    low = payload.get('low_energy', {})
    final = payload.get('final_running', {})
    masses = low.get('masses_ev', [])
    search = payload.get('benchmark_search', {})
    model = search.get('model', {}) if isinstance(search, dict) else {}
    if not isinstance(model, dict):
        model = {}
    y, reference = SCENARIOS[scenario]
    model_class = _model_class(path, payload)
    factor = LAMBDA_T3_NORMALISATION[model_class]
    return {
        'scenario': scenario, 'model': _model_name(path, payload),
        'model_class': model_class, 'dS1': model.get('d_s1', ''),
        'dS2': model.get('d_s2', ''), 'dF': model.get('d_f', ''),
        'alpha': model.get('alpha', ''), 'yukawa_scale': y,
        'lambdaT3_reference': reference, 'lambdaT3_effective': reference * factor,
        'lambdaT3_normalisation_factor': factor,
        'm1_eV': masses[0] if len(masses) > 0 else '',
        'm2_eV': masses[1] if len(masses) > 1 else '',
        'm3_eV': masses[2] if len(masses) > 2 else '',
        'delta_m21_sq_eV2': low.get('delta_m21_sq_ev2', _last(final.get('delta_m21_sq_ev2'))),
        'delta_m3l_sq_eV2': low.get('delta_m31_sq_ev2', _last(final.get('delta_m3l_sq_ev2'))),
        'sin2_theta12': _last(final.get('sin2_theta12')),
        'sin2_theta13': _last(final.get('sin2_theta13')),
        'sin2_theta23': _last(final.get('sin2_theta23')),
        'takagi_residual': low.get('takagi_residual', ''),
    }


def _diagnostics(data_input: Path, scenario: str) -> list[tuple[Path, dict[str, Any]]]:
    root = data_input / 'comparison' / scenario
    successful = []
    for path in sorted(root.rglob('running_diagnostics.json')) if root.is_dir() else []:
        try:
            payload = _load_json(path)
        except (OSError, ValueError):
            continue
        if payload.get('status') == 'Success':
            successful.append((path, payload))
    return successful


def _collect_rows(data_root: Path) -> list[dict[str, Any]]:
    rows = [
        _summary_row(scenario, path, payload)
        for scenario in SCENARIO_ORDER
        for path, payload in _diagnostics(data_root, scenario)
    ]
    rows.sort(key=lambda row: (row['scenario'], row['model']))
    return rows


def _write_csv(rows: list[dict[str, Any]], path: Path) -> None:
    if not rows:
        return
    with path.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _write_readme(path: Path, *, rows: list[dict[str, Any]], report_model_count: int) -> None:
    lines = [
        '# T3 numerical study — regenerated comparison display', '',
        'Generated solely from saved `output/full/comparison` diagnostics.',
        'No RGEs, matching, fits or numerical trajectories were recomputed.',
        'The legacy Reports/output/full directory is not an input.', '',
        'The six per-model plots are produced by RunningResultFigures.py;',
        'cross-model HTML dashboards by InteractiveModelComparison.py.',
        'A missing curve is not interpolated or reconstructed by this builder.', '',
        '## Operator-running convention', '',
        'Intermediate plots show the implemented direct one-loop LLSS -> Weinberg',
        'contribution. No LLSS self-running Wilson trajectory is inferred.', '',
        '## Comparison benchmark normalisation', '',
        'lambda_T3^eff = f_R lambda_T3^ref, with',
        'f_R(A,B,C,D,E) = (1, 1, sqrt(3), 1, 1/sqrt(2)).',
        'This acts on benchmark inputs, not physical matching group factors.',
        'The fixed comparison textures are not oscillation-data fits.', '',
    ]
    for scenario, (y, lam) in SCENARIOS.items():
        lines.append(f'- `{scenario}`: Y = {y:g}, lambda_T3^ref = {lam:g}')
    lines += ['', f'Model directories generated: {report_model_count}.',
              f'Successful numerical rows: {len(rows)}.', '']
    path.write_text('\n'.join(lines), encoding='utf-8')


def _ensure_safe_output(data_input: Path, report_input: Path, output: Path) -> None:
    for protected in (data_input, report_input):
        if output == protected or output in protected.parents or protected in output.parents:
            raise ValueError(f'Output must not contain or be inside protected source: {protected}')
    if output == PROJECT_ROOT or PROJECT_ROOT in output.parents and output.name == 'output':
        raise ValueError('Refusing ambiguous output directory')


def build_display(data_input: Path, report_input: Path, output: Path) -> Path:
    """Regenerate selected report figures directly from saved comparison JSON.

    ``report_input`` remains in the signature for existing callers, but is unused.
    The existing output is replaced only after prevalidating available inputs.
    """
    data_input = Path(data_input).resolve()
    report_input = Path(report_input).resolve()
    output = Path(output).resolve()
    _ensure_safe_output(data_input, report_input, output)
    scenarios = [(name, _diagnostics(data_input, name)) for name in SCENARIO_ORDER]
    scenarios = [(name, items) for name, items in scenarios if items]
    if not scenarios:
        raise FileNotFoundError(f'No successful comparison diagnostics under {data_input / "comparison"}')

    # The old display is replaced only after confirming that source data exists.
    # Regeneration is presentation-only; the source diagnostics are never modified.
    if output.exists():
        shutil.rmtree(output)
    output.mkdir(parents=True)
    report_model_count = 0
    for scenario, diagnostics in scenarios:
        source = data_input / 'comparison' / scenario
        dest = output / 'comparison' / scenario
        dest.mkdir(parents=True, exist_ok=True)
        write_dashboard(source, dest / 'interactive_comparison.html',
                        f'T3 normalized comparison: {scenario}')
        for path, payload in diagnostics:
            model_name = _model_name(path, payload)
            figures = dest / 'models' / model_name
            # Generate all supported plots using the original plotting functions;
            # discard non-selected files to keep the display compact.
            generate_running_result_figures(path, figures)
            for item in figures.iterdir():
                if item.is_file() and item.name not in KEEP_FIGURES:
                    item.unlink()
            report_model_count += 1
    rows = _collect_rows(data_input)
    _write_csv(rows, output / 'model_summary.csv')
    _write_readme(output / 'README.md', rows=rows, report_model_count=report_model_count)
    print(f'Display regenerated from diagnostics: {output}')
    print(f'Model figures: {report_model_count}; successful rows: {len(rows)}')
    return output


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-input', type=Path, default=DEFAULT_DATA_INPUT)
    parser.add_argument('--report-input', type=Path, default=DEFAULT_REPORT_INPUT,
                        help='Compatibility argument; no longer read')
    parser.add_argument('--output', type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    build_display(args.data_input, args.report_input, args.output)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
