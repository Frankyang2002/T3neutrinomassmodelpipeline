"""Numerical extra-only Majorana mass and relative-phase diagnostics."""
from __future__ import annotations
import copy
import json
from pathlib import Path

import numpy as np
import pytest

from Numerical.experimental.GeneralT3NeutrinoImpact import analyze_extra_mass_response, load_and_analyze

A = "lambdaT3ExtraQ020020I1"
B = "lambdaT3ExtraQ020200I1"


def block(m):
    m=np.asarray(m, dtype=complex)
    return {"real": m.real.tolist(), "imag": m.imag.tolist(),
            "frobenius_norm": float(np.linalg.norm(m, 'fro'))}


def sample(model="B_m1", conjugated_a=True):
    a=np.diag([1e-17,2e-17,3e-17]).astype(complex)
    b=np.diag([1j*2e-17,1e-17,-1e-17]).astype(complex)
    return {"model":model,"status":"ExperimentalHardBoundaryC5",
            "units":"GeV^-1", "finite_mass_ratio_logs_preserved":True,
            "delta_c5":block(a+b), "contributions":{A:block(a),B:block(b)},
            "kernels":[{"coupling":A,"conjugated":conjugated_a,"yukawa":"y2"},
                       {"coupling":B,"conjugated":False,"yukawa":"y1"}]}


def test_contribution_matrix_norm_and_takagi():
    r=analyze_extra_mass_response(sample())
    assert r["status"]=="ExtraOnlyNeutrinoMassImpact"
    original=r["scenarios"]["both_as_input"]
    mass=np.array(original["extra_neutrino_mass_matrix_ev"]["real"])+1j*np.array(original["extra_neutrino_mass_matrix_ev"]["imag"])
    assert np.allclose(sorted(np.linalg.svd(mass, compute_uv=False)), original["extra_only_takagi_singular_values_ev"])
    assert np.isclose(original["c5"]["frobenius_norm"], np.linalg.norm(np.diag([1e-17+2j*1e-17, 3e-17, 2e-17])))
    assert original["extra_neutrino_mass_matrix_ev"]["frobenius_norm"]>0


def test_single_contribution_and_zero():
    r=analyze_extra_mass_response(sample())
    assert r["scenarios"]["extra_couplings_zero"]["extra_neutrino_mass_matrix_ev"]["frobenius_norm"]==0.0
    assert r["scenarios"]["extra_couplings_zero"]["extra_only_takagi_singular_values_ev"] == [0.0,0.0,0.0]
    assert r["scenarios"]["only_first"]["c5"]["frobenius_norm"]>0
    assert r["scenarios"]["only_second"]["c5"]["frobenius_norm"]>0


def test_interference_decomposition_exact():
    r=analyze_extra_mass_response(sample())
    d=r["c5_norm_squared_decomposition_gev_minus2"]
    assert np.isclose(d["combined_norm_squared"],d["sum_individual_norm_squared"]+d["interference_2ReTrAstarB"])


def test_phase_grid_and_conjugation_orientation():
    r=analyze_extra_mass_response(sample())
    assert len(r["phase_scan"])==16
    assert r["conjugated_coupling_atoms"][A] is True
    one=next(x for x in r["phase_scan"] if x["physical_coupling_phase_shifts_deg"]=={A:90,B:0})
    # Physical +90 degree rotation on A acts as -i on Bar[lambda_A].
    mm=np.diag([-1j*1e-17+2j*1e-17,-2j*1e-17+1e-17,-3j*1e-17-1e-17])
    ref=np.linalg.norm(mm)*(246.22**2/2)*1e9
    assert np.isclose(one["mnu_norm_ev"],ref)


def test_phase_scan_min_max_within_bounds():
    r=analyze_extra_mass_response(sample())
    all_vals=[t["mnu_norm_ev"] for t in r["phase_scan"]]
    assert r["minimum_mass_norm_phase"]["mnu_norm_ev"]==min(all_vals)
    assert r["maximum_mass_norm_phase"]["mnu_norm_ev"]==max(all_vals)
    assert min(all_vals)<=r["scenarios"]["both_as_input"]["extra_neutrino_mass_matrix_ev"]["frobenius_norm"]<=max(all_vals)


def test_contribution_sum_corruption_rejected():
    data=sample()
    data["delta_c5"]=block(np.zeros((3,3)))
    with pytest.raises(ValueError,match="do not sum"):
        analyze_extra_mass_response(data)


def test_missing_or_repeated_metadata_rejected():
    data=sample()
    data["kernels"][1]["coupling"]=A
    with pytest.raises(ValueError,match="orientation"):
        analyze_extra_mass_response(data)


def test_bad_matrix_units_status_or_vev_rejected():
    for key,value in (("units","eV^-1"),("status","Failed"),("model","A_m2")):
        data=sample();data[key]=value
        with pytest.raises(ValueError): analyze_extra_mass_response(data)
    with pytest.raises(ValueError,match="vev_gev"): analyze_extra_mass_response(sample(),vev_gev=-1)


def test_nonsymmetric_matrix_rejected():
    data=sample()
    m=np.zeros((3,3),complex);m[0,1]=1e-17
    data["contributions"][A]=block(m)
    with pytest.raises(ValueError,match="not symmetric"):
        analyze_extra_mass_response(data)


def test_json_load_and_cli(tmp_path):
    p=tmp_path/'pilot.json';p.write_text(json.dumps(sample('C_m1')))
    assert load_and_analyze(p)["model"]=="C_m1"
    import subprocess,sys
    root=Path(__file__).resolve().parents[2]
    inp=tmp_path/'inputs';inp.mkdir()
    for model in ('B_m1','C_m1'):
        (inp/f'{model}_matchete_hard_c5.json').write_text(json.dumps(sample(model)))
    out=tmp_path/'outputs'
    task=subprocess.run([sys.executable,str(root/'scripts/RunGeneralT3ExtraC5PhaseResponse.py'),'--input-root',str(inp),'--output-root',str(out)],cwd=root,capture_output=True,text=True)
    assert task.returncode==0,task.stderr
    assert (out/'extra_mass_impact_summary.csv').is_file()
    assert json.loads((out/'extra_mass_impact_summary.json').read_text())['status']=='Success'


def test_cli_fails_closed_on_missing_input(tmp_path):
    import subprocess,sys
    root=Path(__file__).resolve().parents[2]
    inp=tmp_path/'inputs';inp.mkdir()
    out=tmp_path/'outputs'
    proc=subprocess.run([sys.executable,str(root/'scripts/RunGeneralT3ExtraC5PhaseResponse.py'),'--input-root',str(inp),'--output-root',str(out)],cwd=root,capture_output=True,text=True)
    assert proc.returncode==2
    assert json.loads((out/'extra_mass_impact_summary.json').read_text())['status']=='Incomplete'
