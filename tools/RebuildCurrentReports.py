"""Rebuild current T3 figures and displays from SAVED numerical data only.

Project path: tools/RebuildCurrentReports.py
Run from project root: python -m tools.RebuildCurrentReports

No RGE integration, fitting, Matchete, RGBeta, or symbolic reconstruction.
The superseded four-scenario workflow, Study A/B HTML, and standalone
interactive_correlations.html are excluded. Historical scan directories are
retained; figures are presentation-only, with canonical T3A-m2 style labels.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import sys
from collections.abc import Callable

from tools.RedirectStudyFigures import redirect_matplotlib_saves
from tools.LabelSensitivityDirections import label_html
from Numerical.plotting.DisplayModelLabels import (
    presentation_labels_on_save, presentation_envelope_labels,
)

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / 'Reports' / 'output'
FULL = REPORT / 'full'
DISPLAY = REPORT / 'display'
SCENARIO_DEFAULT = 'smallY_smallL'


def read_json(path: Path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def scan_rows(path: Path) -> list[dict]:
    value = read_json(path)
    if not isinstance(value, list):
        raise ValueError(f'Expected saved JSON list: {path}')
    return value


def successful_diagnostics(folder: Path) -> list[Path]:
    found = []
    for path in sorted(folder.glob('*/running_diagnostics.json')):
        if read_json(path).get('status') == 'Success':
            found.append(path)
    return found


def rebuild_running_figures(scenario: str) -> None:
    from Numerical.plotting.RunningResultFigures import generate_running_result_figures

    total = 0
    for kind in ('fitted_rge', 'common_input_rge'):
        folder = FULL / kind / scenario
        for diagnostic in successful_diagnostics(folder):
            generate_running_result_figures(diagnostic, DISPLAY / 'full_studies' / kind / scenario / diagnostic.parent.name / 'figures')
            total += 1
    if not total:
        raise FileNotFoundError('No successful saved fitted/common-input RGE diagnostics')
    print(f'Regenerated running plots for {total} model trajectories')


def rebuild_fitted_rge_display(scenario: str) -> None:
    from Numerical.orchestration.studies.fitted_rge.RunFittedRGEDisplay import (
        add_parameters_link, build_display_from_completed, write_parameter_report,
    )
    from Numerical.orchestration.studies.fitted_rge.RunFittedSensitivityComparison import best_candidate

    folder = FULL / 'fitted_rge' / scenario
    records = []
    for path in successful_diagnostics(folder):
        model = path.parent.name
        row = {'model': model, 'status': 'AlreadyComplete', 'diagnostic_path': str(path),
               'fit_config': '', 'fit_chi2': None, 'fit_stage': ''}
        selected = best_candidate(scenario, model)
        if selected is not None:
            chi2, stage, config, _ = selected
            row.update(fit_config=str(config), fit_chi2=chi2, fit_stage=stage)
        records.append(row)
    if not records:
        raise FileNotFoundError(f'No successful saved fitted-RGE data under {folder}')
    dashboard = DISPLAY / f'optimised_rge_{scenario}.html'
    dashboard.parent.mkdir(parents=True, exist_ok=True)
    build_display_from_completed(records, folder, dashboard, f'Optimised T3 RGE: {scenario}')
    if any(row['fit_config'] for row in records):
        write_parameter_report(records, folder, DISPLAY, scenario, dashboard.name)
        add_parameters_link(dashboard, scenario)
    else:
        print('[NOTE] Parameter table skipped: fitted configs not found locally')


def rebuild_common_rge_display(scenario: str) -> None:
    from Numerical.orchestration.studies.common_inputs.RunCommonInputRGE import make_dashboard

    folder = FULL / 'common_input_rge' / scenario
    rows = [{'model': path.parent.name, 'status': 'AlreadyComplete',
             'diagnostics': str(path)} for path in successful_diagnostics(folder)]
    if not rows:
        raise FileNotFoundError(f'No saved common-input RGE trajectories: {folder}')
    make_dashboard(rows, DISPLAY / f'common_input_rge_{scenario}.html')


def fitted_sensitivity_plots(scenario: str) -> None:
    from Numerical.orchestration.studies.fitted_rge.RunFittedSensitivityComparison import plot_comparison

    folder = FULL / 'fitted_sensitivity' / scenario
    rows = scan_rows(folder / 'comparison_summary.json')
    # Saved summaries may contain absolute paths from a different checkout.
    for row in rows:
        if 'model' in row:
            row['scan_json'] = str(folder / row['model'] / 'sensitivity_diagnostics.json')
    plot_comparison(rows, folder)


def fitted_sensitivity_thesis_plots(scenario: str) -> None:
    folder = FULL / 'fitted_sensitivity' / scenario
    if not any(folder.glob('*/sensitivity_diagnostics.json')):
        raise FileNotFoundError(f'No saved sensitivity scans: {folder}')
    # A separate Python process would otherwise bypass the figure-label context.
    subprocess.run([sys.executable, '-m', 'tools.PlotDisplayWrapper',
                    'Numerical.plotting.PlotFittedSensitivityFigures',
                    '--input', str(folder), '--output', str(DISPLAY / 'full_studies' / 'fitted_sensitivity' / scenario / 'thesis_figures')],
                   cwd=ROOT, check=True)


def fitted_sensitivity_envelopes(scenario: str) -> None:
    from Numerical.plotting.BuildFittedModelEnvelopes import build

    output = DISPLAY / 'fitted_model_envelopes_all.html'
    with presentation_envelope_labels():
        build(FULL / 'fitted_sensitivity' / scenario, output)
    label_html(output)


def common_input_envelopes(scenario: str) -> None:
    from Numerical.orchestration.studies.common_inputs.RunCommonInputEnvelopes import make_common_envelopes

    output = DISPLAY / 'common_input_model_envelopes_all.html'
    with presentation_envelope_labels():
        make_common_envelopes(FULL / 'common_input_sensitivity' / scenario, output)
    label_html(output)


def normalized_sensitivity(scenario: str) -> None:
    from Numerical.plotting.NormalizedFittedSensitivity import analyze

    analyze(FULL / 'fitted_sensitivity' / scenario,
            FULL / 'normalized_fitted_sensitivity' / scenario)


def generic_saved_scan(folder: Path, filename: str, fn: Callable) -> None:
    fn(scan_rows(folder / filename), folder)


def fermion_threshold_variation_plots(scenario: str) -> None:
    from Numerical.plotting.FermionThresholdVariationFigures import rebuild

    rebuild('fermion-threshold', scenario)


def weinberg_running_plots(scenario: str) -> None:
    from Numerical.plotting.FermionThresholdVariationFigures import rebuild

    rebuild('weinberg', scenario)


def fitted_fermion_plots(scenario: str) -> None:
    from Numerical.plotting.FermionThresholdVariationFigures import rebuild

    rebuild('fitted-mass', scenario)


def common_coupling_plots(scenario: str) -> None:
    from Numerical.orchestration.studies.common_inputs.RunCommonInputCouplingScan import save_plots

    folder = FULL / 'common_input_coupling_scan' / scenario
    save_plots(scan_rows(folder / 'scan.json'), folder / 'plots')


def common_fermion_plots(scenario: str) -> None:
    from Numerical.orchestration.studies.common_inputs.RunCommonInputFermionScale import save_plots

    folder = FULL / 'common_input_fermion_scale' / scenario
    save_plots(scan_rows(folder / 'scan.json'), folder / 'plots')


def common_scalar_plots(scenario: str) -> None:
    from Numerical.fitting.OscillationFit import OscillationFitTarget
    from Numerical.orchestration.studies.scalar_rge.RunCommonScalarComparison import make_plots

    folder = FULL / 'common_scalar_comparison' / scenario
    rows = scan_rows(folder / 'comparison_summary.json')
    for row in rows:
        if 'model' in row:
            row['scan_json'] = str(folder / row['model'] / 'common_scalar_scan.json')
    target = OscillationFitTarget.from_json(ROOT / 'configs' / 'nufit_6_1_ic24_no.json')
    make_plots(rows, folder, target)


def common_scalar_control_plots(scenario: str) -> None:
    from Numerical.orchestration.studies.scalar_rge.RunCommonScalarRGControls import create_plots

    folder = FULL / 'common_scalar_rg_controls' / scenario
    create_plots(scan_rows(folder / 'running_vs_fixed.json'), folder / 'plots')


def common_scalar_beta_plots(scenario: str) -> None:
    from Numerical.orchestration.studies.scalar_rge.RunCommonScalarBetaDominance import plot_results

    folder = FULL / 'common_scalar_beta_dominance' / scenario
    payload = read_json(folder / 'beta_decomposition.json')
    plot_results(payload['rows'], folder / 'plots')


def heavy_fermion_control_plots(scenario: str) -> None:
    from Numerical.orchestration.studies.fermion_scales.RunHeavyFermionScaleComparison import save_plots

    folder = FULL / 'heavy_fermion_scale_comparison' / scenario
    payload = read_json(folder / 'running_vs_fixed_mass_scan.json')
    save_plots(payload['points'], folder / 'plots')


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--scenario', default=SCENARIO_DEFAULT,
                    help='Retained saved-data folder name, not the obsolete scenario sweep')
    ap.add_argument('--only', choices=('all', 'displays', 'plots'), default='all',
                    help='Select HTML displays, static figures, or both')
    ap.add_argument('--strict', action='store_true',
                    help='Fail on optional missing saved input as well as errors')
    args = ap.parse_args(argv)
    if not (ROOT / 'pipeline.py').is_file():
        ap.error(f'Unexpected project root: {ROOT}')

    tasks: list[tuple[str, str, Callable[[], None]]] = [
        ('plots', 'Fitted/common-input per-model running figures', lambda: rebuild_running_figures(args.scenario)),
        ('displays', 'Fitted RGE interactive dashboard and UV parameter table', lambda: rebuild_fitted_rge_display(args.scenario)),
        ('displays', 'Common A-m2 RGE interactive dashboard', lambda: rebuild_common_rge_display(args.scenario)),
        ('displays', 'Individually fitted gamma-envelope HTML', lambda: fitted_sensitivity_envelopes(args.scenario)),
        ('displays', 'Common A-m2 gamma-envelope HTML', lambda: common_input_envelopes(args.scenario)),
        ('plots', 'Fitted gamma observable-pair plots', lambda: fitted_sensitivity_plots(args.scenario)),
        ('plots', 'Fitted gamma thesis figures', lambda: fitted_sensitivity_thesis_plots(args.scenario)),
        ('plots', 'Fitted sigma-normalized sensitivity', lambda: normalized_sensitivity(args.scenario)),
        ('plots', 'Fitted fermion-scale figures', lambda: fitted_fermion_plots(args.scenario)),
        ('plots', 'Fermion-threshold variation: stage-on vs segment-frozen plots', lambda: fermion_threshold_variation_plots(args.scenario)),
        ('plots', 'Fermion-threshold variation: isolated Weinberg-running plots', lambda: weinberg_running_plots(args.scenario)),
        ('plots', 'Common-input coupling scan figures', lambda: common_coupling_plots(args.scenario)),
        ('plots', 'Common-input fermion-scale figures', lambda: common_fermion_plots(args.scenario)),
        ('plots', 'Common-scalar coupling figures', lambda: common_scalar_plots(args.scenario)),
        ('plots', 'Common-scalar running/frozen control figures', lambda: common_scalar_control_plots(args.scenario)),
        ('plots', 'Common-scalar beta decomposition figures', lambda: common_scalar_beta_plots(args.scenario)),
        ('plots', 'Heavy-fermion running/frozen control figures', lambda: heavy_fermion_control_plots(args.scenario)),
    ]
    ok, missing, failed = 0, [], []
    for group, name, action in tasks:
        if args.only != 'all' and group != args.only:
            continue
        print(f'\n=== {name} ===', flush=True)
        try:
            # Applies to any static plot produced by the current rebuild workflow.
            # A no-op for HTML builders; it changes only figure text on save.
            with presentation_labels_on_save():
                with redirect_matplotlib_saves():
                    action()
            print('[OK]', name, flush=True)
            ok += 1
        except FileNotFoundError as exc:
            print('[SKIP] Missing saved input:', exc, flush=True)
            missing.append(name)
        except Exception as exc:
            print(f'[FAIL] {type(exc).__name__}: {exc}', flush=True)
            failed.append(name)
    print(f'\nRebuild summary: {ok} succeeded; {len(missing)} missing-input skips; {len(failed)} failed.')
    if missing:
        print('Missing input:', '; '.join(missing))
    if failed:
        print('Failures:', '; '.join(failed))
    print('Numerical source: Reports/output/full/; figures: Reports/output/display/full_studies/')
    print('No physics calculations or symbolic generators were executed.')
    return 1 if failed or args.strict and missing else 0


if __name__ == '__main__':
    raise SystemExit(main())
