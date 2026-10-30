import csv
import importlib.util
import json
from pathlib import Path

mod_path = Path(__file__).resolve().parents[1] / 'Numerical/plotting/BuildFrozenGeneralScalarLambdaReports.py'
spec = importlib.util.spec_from_file_location('frozen_build', mod_path)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def setup_repo(tmp_path: Path):
    base = tmp_path / 'output/studies/general_scalar_extended/bulk_48/frozen_massnorm_a0_v1/small_uv_frozen_5e-9/A_p0'
    base.mkdir(parents=True)
    state = {
        'status': 'Accepted',
        'best': {
            'chi2': 1e-18,
            'prediction': {'delta_m21_sq_ev2': 7.5e-5, 'delta_m3l_sq_ev2': 2.5e-3, 'sin2_theta12': 0.3, 'sin2_theta13': 0.022, 'sin2_theta23': 0.57},
            'parameters': {f'p{i}': i for i in range(18)},
        },
    }
    (base / 'state.json').write_text(json.dumps(state), encoding='utf-8')
    sens_dir = tmp_path / 'output/studies/general_scalar_extended/frozen_display/sensitivity_v1'
    sens_dir.mkdir(parents=True)
    with (sens_dir / 'frozen_gamma_points.csv').open('w', newline='', encoding='utf-8') as f:
        wr = csv.DictWriter(f, fieldnames=['model', 'case', 'direction', 'gamma', 'status', 'chi2', 'sin2_theta13', 'sin2_theta23'])
        wr.writeheader()
        for d in ('T1', 'T2', 'T3'):
            for g in (0.8, 1.0, 1.2):
                wr.writerow({'model': 'A_+0', 'case': 'small_uv_frozen_5e-9', 'direction': d, 'gamma': g, 'status': 'Accepted', 'chi2': 1.0, 'sin2_theta13': 0.022, 'sin2_theta23': 0.57})
    traj_dir = tmp_path / 'output/studies/general_scalar_extended/frozen_display/rge_v1/trajectories_v1'
    traj_dir.mkdir(parents=True)
    traj = {'model': 'A_+0', 'case': 'small_uv_frozen_5e-9', 'note': 'x', 'uv_trajectory': [{'scale_gev': 1e7, 'y1_fro': 1.0}], 'eft1_trajectory': [{'scale_gev': 1e3, 'delta_m3l_sq_ev2': 2.5e-3}]}
    (traj_dir / 'A_p0__small_uv_frozen_5e-9.json').write_text(json.dumps(traj), encoding='utf-8')
    return sens_dir / 'frozen_gamma_points.csv', traj_dir


def test_missing_masses_ok(tmp_path, monkeypatch):
    sens_csv, traj_dir = setup_repo(tmp_path)
    monkeypatch.setattr(mod, 'ROOT', tmp_path)
    monkeypatch.setattr(mod, 'OUT_BASE', tmp_path / 'output/studies/general_scalar_extended/bulk_48')
    monkeypatch.setattr(mod, 'REPORT_BASE', tmp_path / 'Reports/output/display_frozen_general_scalar')
    rows, _ = mod.summary_rows(['frozen_massnorm_a0_v1'])
    assert len(rows) == 1
    assert rows[0]['mass_hi_ev'] is None
    man = mod.build(['frozen_massnorm_a0_v1'], 'out', str(traj_dir), str(sens_csv))
    assert man['accepted_row_count'] == 1
    assert man['trajectory_model_count'] == 1
    assert man['sensitivity_model_count'] == 1
    assert (tmp_path / 'Reports/output/display_frozen_general_scalar/out/interactive_rge.html').is_file()
    assert (tmp_path / 'Reports/output/display_frozen_general_scalar/out/sensitivity_report.html').is_file()
