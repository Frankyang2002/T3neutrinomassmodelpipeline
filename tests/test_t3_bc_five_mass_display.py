"""Unit tests for the read-only Workstream B UV/EFT1 display tool."""
from __future__ import annotations

from copy import deepcopy
import csv
import hashlib
import json
import math
from pathlib import Path

import pytest

from scripts import BuildT3BCFiveMassDisplay as display


def make_row(model: str, stage: str, scale: float, step: int, *, shift: float = 0.0):
    # Purely synthetic fixtures; no RGE computation.
    m1 = 1000000.0 + 12000*step + shift
    m2 = 1200000.0 + 4000*step
    mixed = complex(10000+400*step, 3000+250*step)
    centre = (m1 + m2)/2
    width = math.hypot((m1-m2)/2, abs(mixed))
    low = centre-width
    high = centre+width
    return {
        "mu_gev": float(scale),
        "mHSq_gev2": -7800.0 - 100*step - shift,
        "mS1Sq_gev2": m1,
        "mS2Sq_gev2": m2,
        "muT3GeneralS1S2_gev2": {"real": mixed.real, "imag": mixed.imag},
        "scalar_spectrum": {
            "m_minus_squared_gev2": low,
            "m_plus_squared_gev2": high,
            "real_degeneracy_per_eigenvalue": 4,
            "mass_squared_eigenvalues_gev2": [low]*4+[high]*4,
            "tachyonic": False,
        },
    }


def payload_for(model: str) -> dict:
    shift = 0 if model == "B_-1" else 1000
    uv = [make_row(model, "UV", 1e7, 0),
          make_row(model, "UV", 1e6, 1, shift=shift),
          make_row(model, "UV", 1e5, 2, shift=shift)]
    eft = [deepcopy(uv[-1]),
           make_row(model, "EFT1", 1e4, 3, shift=shift),
           make_row(model, "EFT1", 1e3, 4, shift=shift)]
    # Values deliberately maintain perfect UV/EFT1 boundary continuity.
    return {
        "format": display.SOURCE_FORMAT,
        "status": display.SOURCE_STATUS,
        "model": model,
        "solver_success": True,
        "one_loop_convention": display.CONVENTION,
        "loop_factor_applied_exactly_once": True,
        "previous_mass_betas_replaced_in_both_stages": True,
        "physical_matchete_rgbeta_mapping_certified": False,
        "scalar_mass_eigenstate_threshold_matching_performed": False,
        "physical_neutrino_mass_prediction": False,
        "five_real_mass_coordinates": list(display.MASS_COORDS),
        "uv_range_gev": [1e7, 1e5],
        "eft1_range_gev": [1e5, 1e3],
        "initial_higgs_mass_squared_gev2_at_uv": -7800,
        "initial_mu2_fraction_at_uv": 0.01,
        "initial_mu2_phase_at_uv": 0.6,
        "rtol": 1e-8,
        "atol": 1e-11,
        "threshold": {
            "scale_gev": 1e5,
            "tree_level_renormalisable_only": True,
            "finite_mass_matching_performed": False,
            "wilson_coefficient_matching_performed": False,
            "heavy_fermion_singular_masses_gev": [1e5, 1.1e5, 1.2e5],
        },
        "uv_beta_keys": list(display.MASS_SYMBOLS) + ["gY", "MF", "y1", "y2"],
        "eft1_beta_keys": list(display.MASS_SYMBOLS) + ["gY"],
        "uv_solver_nfev": 212,
        "eft1_solver_nfev": 47,
        "uv_trajectory": uv,
        "eft1_trajectory": eft,
        "uv_initial": deepcopy(uv[0]),
        "uv_threshold": deepcopy(uv[-1]),
        "eft1_threshold": deepcopy(eft[0]),
        "eft1_final": deepcopy(eft[-1]),
    }


def fixture_files(root: Path) -> dict[str, Path]:
    root.mkdir(parents=True, exist_ok=True)
    paths = {}
    for model in display.MODELS:
        tag = model.replace("-", "m")
        path = root / f"{tag}_uv_to_eft1_coupled_five_mass_v1.json"
        path.write_text(json.dumps(payload_for(model)), encoding="utf-8")
        paths[model] = path
    return paths


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_seven_components_and_boundary(tmp_path):
    path = fixture_files(tmp_path)["B_-1"]
    report = display.read_model(path, "B_-1")
    assert len(display.KEYS) == 7
    assert report["uv"][-1]["mS12"] == report["eft1"][0]["mS12"]
    assert report["uv"][-1]["im_mu"] == report["eft1"][0]["im_mu"]
    assert report["uv"][0]["mH2"] == -7800


def test_reject_threshold_mass_jump(tmp_path):
    path = fixture_files(tmp_path)["C_-1"]
    payload = json.loads(path.read_text())
    payload["eft1_trajectory"][0]["mHSq_gev2"] += 50
    payload["eft1_threshold"]["mHSq_gev2"] += 50
    path.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="mass state discontinuity"):
        display.read_model(path, "C_-1")


def test_reject_incorrect_spectrum_degeneracy(tmp_path):
    path = fixture_files(tmp_path)["B_-1"]
    payload = json.loads(path.read_text())
    payload["uv_trajectory"][1]["scalar_spectrum"]["mass_squared_eigenvalues_gev2"][0] += 500
    path.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="eight-fold spectrum"):
        display.read_model(path, "B_-1")


def test_reject_claimed_physical_certificate(tmp_path):
    path = fixture_files(tmp_path)["B_-1"]
    payload = json.loads(path.read_text())
    payload["physical_matchete_rgbeta_mapping_certified"] = True
    path.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="source metadata mismatch"):
        display.read_model(path, "B_-1")


def test_read_only_export_and_csv_counts(tmp_path):
    sources = fixture_files(tmp_path / "inputs")
    before = {model: digest(path) for model,path in sources.items()}
    outdir = tmp_path / "display"
    summary = display.build(input_dir=tmp_path / "inputs", output_dir=outdir,
                            plots=False)
    assert {model: digest(path) for model,path in sources.items()} == before
    assert summary["status"] == "PresentationOnlyNoPhysicsRecomputed"
    assert len(summary["outputs"]) == 4
    assert all((outdir / name).is_file() for name in summary["outputs"])
    with (outdir / "bc_uv_eft1_mass_trajectories.csv").open(newline="") as fh:
        rows = list(csv.DictReader(fh))
    assert len(rows) == 12  # 2 models × (3 UV + 3 EFT1)
    with (outdir / "bc_uv_eft1_mass_endpoints.csv").open(newline="") as fh:
        endpoints = list(csv.DictReader(fh))
    assert len(endpoints) == 6
    assert "No RGEs were recomputed" in (outdir / "bc_uv_eft1_benchmark_summary.md").read_text()
    assert summary["model_summary"]["B_-1"]["uv_solver_nfev"] == 212
    assert summary["eft1_final_C_minus_B"]["mH2"] == pytest.approx(-1000)


def test_plots_from_saved_states(tmp_path):
    pytest.importorskip("matplotlib")
    fixture_files(tmp_path / "input")
    summary = display.build(input_dir=tmp_path / "input",
                            output_dir=tmp_path / "display", plots=True)
    figures = [p for p in summary["outputs"] if p.endswith(".png")]
    assert len(figures) == 7
    assert all((tmp_path / "display" / f).stat().st_size > 1000 for f in figures)
