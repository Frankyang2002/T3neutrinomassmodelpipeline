import csv
import importlib.util
import json
from pathlib import Path

MODULE = Path(__file__).resolve().parents[1] / 'Numerical/plotting/BuildFrozenGeneralScalarLambdaReports.py'
spec = importlib.util.spec_from_file_location('frozen_reports', MODULE)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def setup_runs(tmp_path: Path):
    base = tmp_path / 'output/studies/general_scalar_extended/bulk_48'
    for run_id, case, model, chi2 in [
        ('frozen_massnorm_a0_v1', 'small_uv_frozen_5e-9', 'A_+0', 1e-18),
        ('frozen_massnorm_remaining_v1', 'small_uv_frozen_5e-9', 'E_-2', 2e-18),
    ]:
        path = base / run_id / case / mod.safe(model)
        path.mkdir(parents=True)
        state = {
            'status': 'Accepted',
            'best': {
                'chi2': chi2,
                'prediction': {
                    'mass_hi_ev': 0.05,
                    'mass_mid_ev': 0.009,
                    'mass_lo_ev': 0.001,
                    'delta_m21_sq_ev2': 7.5e-5,
                    'delta_m3l_sq_ev2': 2.5e-3,
                    'sin2_theta12': 0.30,
                    'sin2_theta13': 0.022,
                    'sin2_theta23': 0.57,
                },
            },
        }
        (path / 'state.json').write_text(json.dumps(state), encoding='utf-8')
    plotting = tmp_path / 'Numerical/plotting'
    plotting.mkdir(parents=True)
    (plotting / 'InteractiveModelComparison.py').write_text('<html>__TITLE__<script>const DATA=__DATA__;</script>__NOTE__</html>')
    (plotting / 'BuildFittedModelEnvelopes.py').write_text('<html>__TITLE__<script>const DATA=__DATA__;</script>__NOTE__</html>')


def test_summary_and_report_build(tmp_path, monkeypatch):
    setup_runs(tmp_path)
    monkeypatch.setattr(mod, 'ROOT', tmp_path)
    out = tmp_path / 'Reports/output/display_frozen_general_scalar/test'
    manifest = mod.build(['frozen_massnorm_a0_v1', 'frozen_massnorm_remaining_v1'], out)
    assert manifest['accepted_row_count'] == 2
    assert (out / 'index.html').is_file()
    assert (out / 'interactive_rge.html').is_file()
    assert (out / 'sensitivity_report.html').is_file()
    text = (out / 'sensitivity_report.html').read_text(encoding='utf-8').lower()
    assert 'no real frozen t1/t2/t3 gamma scan csv was supplied' in text


def test_sensitivity_csv_parsing(tmp_path, monkeypatch):
    setup_runs(tmp_path)
    monkeypatch.setattr(mod, 'ROOT', tmp_path)
    csv_path = tmp_path / 'gamma.csv'
    with csv_path.open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=['model', 'case', 'direction', 'gamma', 'sin2_theta12', 'delta_m21_sq_ev2'])
        writer.writeheader()
        for direction in ('T1', 'T2', 'T3'):
            for gamma in (0.8, 1.0, 1.2):
                writer.writerow({
                    'model': 'A_+0',
                    'case': 'small_uv_frozen_5e-9',
                    'direction': direction,
                    'gamma': gamma,
                    'sin2_theta12': 0.3,
                    'delta_m21_sq_ev2': 7.5e-5,
                })
    models, skipped, state = mod.load_sensitivity(csv_path, {})
    assert state == 'SAVED_FROZEN_GAMMA_CSV'
    assert len(models) == 1
    assert skipped == []


def test_missing_absolute_masses_are_not_invented(tmp_path):
    cfg = tmp_path / 'fitted_config.json'
    cfg.write_text(json.dumps({'conditional_general_scalar_fit': {'chi2': 1e-18}}))
    assert mod._extract_masses({'delta_m21_sq_ev2': 7.5e-5, 'delta_m3l_sq_ev2': 2.5e-3}, cfg) == (None, None, None)
    assert mod._extract_masses({'mass_hi_ev': 0.06, 'mass_mid_ev': 0.01, 'mass_lo_ev': 0.001}, cfg) == (0.06, 0.01, 0.001)
    assert mod._extract_masses({'mass_hi_ev': 0.06, 'mass_mid_ev': None, 'mass_lo_ev': 0.001}, cfg) == (None, None, None)


def test_report_with_realistic_prediction_missing_masses(tmp_path, monkeypatch):
    setup_runs(tmp_path)
    monkeypatch.setattr(mod, 'ROOT', tmp_path)
    path = tmp_path / 'output/studies/general_scalar_extended/bulk_48/frozen_massnorm_a0_v1/small_uv_frozen_5e-9/A_p0/state.json'
    state = json.loads(path.read_text())
    for key in ('mass_hi_ev', 'mass_mid_ev', 'mass_lo_ev'):
        del state['best']['prediction'][key]
    path.write_text(json.dumps(state))
    dest = tmp_path / 'report'
    manifest = mod.build(['frozen_massnorm_a0_v1', 'frozen_massnorm_remaining_v1'], dest)
    assert manifest['accepted_row_count'] == 2
    with (dest / 'all_results.csv').open(newline='', encoding='utf-8') as stream:
        rows = list(csv.DictReader(stream))
    a0 = next(row for row in rows if row['model'] == 'A_+0')
    assert a0['mass_hi_ev'] == a0['mass_mid_ev'] == a0['mass_lo_ev'] == ''
    assert a0['delta_m21_sq_ev2'] == '7.5e-05'
