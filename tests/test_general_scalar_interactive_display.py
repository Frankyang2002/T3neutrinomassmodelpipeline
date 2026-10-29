"""Read-only offline general-scalar reporting tests (synthetic data only)."""
from __future__ import annotations

import csv
import json
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.BuildGeneralScalarInteractiveDisplay import (
    MODELS, build, choose_preferred, convert_row, load_parameter_scans,
    read_continuous, _script_json,
)


def csv_file(path: Path, rows: list[dict]):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def make_fixture(root: Path):
    root.mkdir(parents=True, exist_ok=True)
    for model in MODELS:
        tag = model.replace('+', 'p').replace('-', 'm')
        file = root / tag / 'continuous_trajectory.json'
        file.parent.mkdir(parents=True, exist_ok=True)
        def row(mu):
            return {'mu_gev':mu,'mS1Sq_gev2':1_000_000+mu/100,
                    'mS2Sq_gev2':1_200_000+mu/100,
                    'mu2_gev2':{'real':40,'imag':0},
                    'extra_quartics':{'lambdaExample':{'real':.01+mu/1e11,'imag':0}}}
        file.write_text(json.dumps({
            'format':'t3-continuous-general-scalar-candidate-v1',
            'status':'CandidateContinuousUVFermionThresholdEFT1Integrated',
            'model':model, 'quadratic_kind':'mixed',
            'uv_trajectory':[row(1e7),row(1e6),row(1e5)],
            'eft1_trajectory':[row(1e5),row(1e4),row(1e3)],
            'scalar_endpoint_spectrum':{'minimum_mass_squared_gev2':1000000},
            'scalar_endpoint_physically_admissible':True,
        }),encoding='utf-8')
    csv_file(root / 'neutrino_v1/summary.csv',[
        {'model':m,'status':'Baseline','largest_neutrino_mass_ev':'.05',
         'middle_neutrino_mass_ev':'.01','smallest_neutrino_mass_ev':'.005',
         'delta_m21_sq_ev2':'0.000075','sin2_theta12':'.3'}
        for m in MODELS
    ])
    csv_file(root / 'self_anomalous_massdep_v1/combined_16_summary.csv',[
        {'model':m,'matching_category': 'self-bilinear-one-scale-massdependent-conditional'
         if m in ('A_+0','A_-2','D_+0','D_-2','E_+0','E_-2') else
         'existing-candidate-matching-different-scheme',
         'mass1_ev_descending':'.051','mass2_ev_descending':'.01','mass3_ev_descending':'.005'}
        for m in MODELS
    ])
    csv_file(root / 'direct_wilson_sequential_v1/summary.csv',[
        {'model':m,'direct_wilson_source_evolved':'False' if m in ('A_-4','B_-3','C_-3','E_-4') else 'True',
         'max_mnu_ev':'.052' if m not in ('A_-4','B_-3','C_-3','E_-4') else '',
         'middle_mnu_ev':'.01' if m not in ('A_-4','B_-3','C_-3','E_-4') else '',
         'min_mnu_ev':'.005' if m not in ('A_-4','B_-3','C_-3','E_-4') else ''}
        for m in MODELS
    ])
    csv_file(root / 'bc_extra_wilson_sequential_v1/summary.csv',[
        {'model':m,'updated_largest_mnu_ev':'.053','updated_middle_mnu_ev':'.02',
         'updated_smallest_mnu_ev':'.004'} for m in ('B_-1','C_-1')
    ])
    csv_file(root / 'bc_coupled_wilson_v1/summary.csv',[
        {'model':m,'new_largest_mnu_ev':'.054','new_middle_mnu_ev':'.03',
         'new_smallest_mnu_ev':'.004'} for m in ('B_-1','C_-1')
    ])


def test_generate_16_model_report(tmp_path):
    src = tmp_path / 'input'
    make_fixture(src)
    dest = tmp_path / 'reports'
    result = build(src,dest,run_id='continuous_v1')
    assert result['model_count'] == 16
    assert result['selected_matching_sources'] == {
        'baseline':4,'self_massdep':6,'direct':4,'extra':0,'coupled':2,
    }
    assert result['matching_variant_count'] == 16+6+12+2+2
    assert result['general_scalar_parameter_sensitivity_scan_loaded'] is False
    assert 'ovals are intentionally absent' in (dest/'sensitivity_report.html').read_text()
    assert all((dest/f).is_file() for f in ('index.html','interactive_rge.html','all_results.html',
                                               'all_results.csv','matching_variants.csv','provenance.json'))
    with (dest/'all_results.csv').open(newline='',encoding='utf-8') as f:
        records = {row['model']:row for row in csv.DictReader(f)}
    assert records['B_-1']['source'] == 'coupled'
    assert records['A_+0']['source'] == 'self_massdep'
    assert records['B_+1']['source'] == 'direct'
    assert records['A_-4']['source'] == 'baseline'
    assert records['B_-1']['sin2_theta12'] == ''  # never borrow mixing angles
    assert float(records['B_-1']['gap_small_mid_sq_ev2']) == pytest.approx(.03**2-.004**2)
    assert 'lambdaExample' in (dest/'interactive_rge.html').read_text()


def test_preflight_missing_baseline_does_not_create_report(tmp_path):
    src=tmp_path/'input'
    make_fixture(src)
    (src/'neutrino_v1/summary.csv').unlink()
    out=tmp_path/'reports'
    with pytest.raises(ValueError,match='benchmark mass table is required'):
        build(src,out,run_id='continuous_v1')
    assert not out.exists()


def test_preflight_missing_trajectory_does_not_create_report(tmp_path):
    src=tmp_path/'input'
    make_fixture(src)
    (src/'A_p0/continuous_trajectory.json').unlink()
    out=tmp_path/'reports'
    with pytest.raises(ValueError,match='missing 16-model scalar trajectory'):
        build(src,out,run_id='continuous_v1')
    assert not out.exists()


def test_parameter_points_are_not_invented(tmp_path):
    scan=tmp_path/'scan.csv'
    csv_file(scan,[{'model':'B_-1','direction':'T1','gamma':g,
                    'sin2_theta13':.02*g,'delta_m3l_sq_ev2':.002*g} for g in (.5,1,2)])
    points,note=load_parameter_scans(scan,set(MODELS))
    assert len(points)==3 and note is None
    with pytest.raises(ValueError,match='missing gamma=1'):
        csv_file(scan,[{'model':'B_-1','direction':'T1','gamma':g,
                        'sin2_theta13':.02*g,'delta_m3l_sq_ev2':.002*g} for g in (.5,.7,2)])
        load_parameter_scans(scan,set(MODELS))


def test_fake_self_massdep_pass_through_rejected():
    r={'matching_category':'existing-candidate-matching-different-scheme',
       'mass1_ev_descending':.05,'mass2_ev_descending':.02,'mass3_ev_descending':.01}
    assert convert_row('C_-1','self_massdep',r) is None


def test_json_script_escaping():
    assert '<script>' not in _script_json({'x':'<script>'})


def test_unavailable_missing_optional_sources_reported(tmp_path):
    src=tmp_path/'input'
    make_fixture(src)
    (src/'bc_extra_wilson_sequential_v1/summary.csv').unlink()
    (src/'bc_coupled_wilson_v1/summary.csv').unlink()
    manifest=build(src,tmp_path/'reports',run_id='continuous_v1')
    assert manifest['selected_matching_sources']['coupled']==0
    assert len(manifest['unavailable_optional_sources'])==2


def test_actual_saved_gamma_points_are_rendered_only_when_supplied(tmp_path):
    src=tmp_path/'input'
    make_fixture(src)
    scan=tmp_path/'gamma.csv'
    csv_file(scan,[{'model':'B_-1','direction':'T1','gamma':g,
                    'sin2_theta13':.02*g,'delta_m3l_sq_ev2':.002*g}
                   for g in (.5,1.,2.)])
    out=tmp_path/'report'
    manifest=build(src,out,run_id='continuous_v1',sensitivity_csv=scan)
    assert manifest['general_scalar_parameter_sensitivity_scan_loaded'] is True
    assert manifest['parameter_scan_record_count']==3
    assert (out/'parameter_sensitivity_points.csv').exists()
    display=(out/'sensitivity_report.html').read_text()
    assert 'Show descriptive sample covariance oval' in display
    assert '"gamma":0.5' in display


def test_output_root_cannot_contain_inputs(tmp_path):
    src=tmp_path/'run'/'input'
    make_fixture(src)
    with pytest.raises(ValueError,match='must not contain'):
        build(src,src.parent,run_id='continuous_v1')
