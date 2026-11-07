"""Preflight/resume regression tests; these do NOT run Wolfram/Matchete."""
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts/RunT3GeneralOrdinaryScalarHard.py"
# Allow stand-alone ZIP tests outside the installed repository.
if not SCRIPT.is_file():
    SCRIPT = Path(__file__).resolve().parents[2] / "scripts/RunT3GeneralOrdinaryScalarHard.py"
spec = importlib.util.spec_from_file_location("t3_ordinary_hard_test", SCRIPT)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def fixture(tmp_path, *, model="D_+2", wrong_chirality=False):
    root = tmp_path / "source"
    base = root / mod.slug(model)
    data = base / "data"
    data.mkdir(parents=True)
    (data / "eft1_after_F_tree_inputform.txt").write_text(
        "Field[NewScalar1, Scalar, {}, {}]*Field[NewScalar2, Scalar, {}, {}] + "
        "Coupling[MS1, {}, 1]^2 + Coupling[MS2, {}, 1]^2 + "
        "Coupling[y1, {Index[p, Flavor], Index[r, NFlavor]}, 0]*"
        "Coupling[y2, {Index[q, Flavor], Index[r, NFlavor]}, 0]/"
        "Coupling[MF, {Index[r, NFlavor]}, 0] + "
        "Bar[Coupling[y1, {Index[p, Flavor], Index[r, NFlavor]}, 0]]*"
        "Bar[Coupling[y2, {Index[q, Flavor], Index[r, NFlavor]}, 0]]/"
        "Coupling[MF, {Index[r, NFlavor]}, 0] " + "X" * 1200,
        encoding="utf-8")
    (data / "eft1_after_F_wilson_seed.json").write_text(json.dumps({
        "IntegratedFields": ["F"], "ActiveHeavyFields": ["S1", "S2"],
        "TreeWilsonTermCount": 2,
        "TreeWilsonTerms": [{"Chirality": "PR"},
                            {"Chirality": "PR" if wrong_chirality else "PL"}]
    }), encoding="utf-8")
    (data / "eft1_after_F_scalar_quartic_seed.json").write_text(json.dumps({
        "IntegratedFields": ["F"], "ActiveHeavyFields": ["S1", "S2"],
        "ScalarQuarticTermCount": 8}), encoding="utf-8")
    (data / "eft1_after_F_rge_export.json").write_text(json.dumps({
        "IntegratedFields": ["F"], "ActiveHeavyFields": ["S1", "S2"]
    }), encoding="utf-8")
    (base / "native_status.json").write_text(json.dumps({
        "model": model, "status": "FThresholdTreeWilsonExported", "loop_order": 0,
        "operator_basis_is_full_general_requested": True,
        "integrated_fields": ["F"], "light_bsm_fields": ["S1", "S2"],
        "heavy_flavour_count": 3,
        "full_general_extra_yukawa_count": 0,
        "full_general_quartic_count": 0,
        "full_general_quadratic_count": 0
    }), encoding="utf-8")
    (base / "status.json").write_text(json.dumps({
        "model": model, "status": "FermionThresholdTreeWilsonExported"
    }), encoding="utf-8")
    wl = tmp_path / "test.wl"
    wl.write_text("(* no wolfram in unit tests *)", encoding="utf-8")
    return root, wl


def test_two_ordinary_models_have_unique_slugs():
    assert mod.slug("D_+2") == "D_p2"
    assert mod.slug("A_-4") == "A_m4"
    assert mod.select_models("D_+2,A_-4") == ("D_+2", "A_-4")


def test_reject_other_mass_classes():
    import pytest
    with pytest.raises(mod.PreflightError):
        mod.select_models("D_-2,B_-1")
    with pytest.raises(mod.PreflightError):
        mod.select_models("D_+2,D_+2")


def test_source_preflight(tmp_path):
    root, wl = fixture(tmp_path)
    row = mod.preflight("D_+2", root, wl)
    assert row["status"] == "Ready"
    assert row["tree_wilson_terms"] == 2
    assert len(row["source_fingerprints"]) == len(mod.SOURCES) + 1


def test_source_changes_detected(tmp_path):
    root, wl = fixture(tmp_path)
    before = mod.preflight("D_+2", root, wl)
    with (root / "D_p2" / "data/eft1_after_F_tree_inputform.txt").open("a") as f:
        f.write(" changed")
    after = mod.preflight("D_+2", root, wl)
    assert before["source_fingerprints"] != after["source_fingerprints"]


def test_missing_pl_is_blocked(tmp_path):
    import pytest
    root, wl = fixture(tmp_path, wrong_chirality=True)
    with pytest.raises(mod.PreflightError, match="chiral"):
        mod.preflight("D_+2", root, wl)


def test_missing_source_is_blocked(tmp_path):
    import pytest
    root, wl = fixture(tmp_path)
    (root / "D_p2" / "data/eft1_after_F_wilson_seed.json").unlink()
    with pytest.raises(mod.PreflightError, match="Missing"):
        mod.preflight("D_+2", root, wl)


def test_preflight_only_does_not_run_wolfram(tmp_path, capsys):
    root, wl = fixture(tmp_path)
    a = SimpleNamespace(models="D_+2", source_root=root,
                        wolfram_script=wl, preflight_only=True,
                        run_id="ignored", resume=False, timeout_seconds=1200,
                        wolframscript="wolframscript")
    assert mod.run(a) == 0
    assert "no fermion rematching" in capsys.readouterr().out


def test_wolfram_safe_alpha_argument():
    assert mod.alpha_argument("A_-4") == "m4"
    assert mod.alpha_argument("D_+2") == "p2"
    assert mod.alpha_argument("D_+0") == "p0"


def test_reject_unexpected_mass_contractions(tmp_path):
    import pytest
    root, wl = fixture(tmp_path)
    src = root / "D_p2" / "data" / "eft1_after_F_tree_inputform.txt"
    src.write_text(src.read_text().replace("Coupling[MF, {Index[r, NFlavor]}, 0]",
                                              "Coupling[MF, {}, 0]", 1))
    with pytest.raises(mod.PreflightError, match="indexed diagonal MF"):
        mod.preflight("D_+2", root, wl)


def test_saved_failure_diagnosis_addressed_by_wolfram_template():
    wl = (ROOT / "Lagrangian/RunT3GeneralOrdinaryScalarHard.wl").read_text()
    assert "RelabelIndices[perEigenvalueTree]" in wl
    assert "Matchete`Coupling[MF, {}, 0]" in wl
    assert "PerEigenvalueHeavyMassSumRequired" in wl
    assert "m\" ~~ DigitCharacter" in wl
