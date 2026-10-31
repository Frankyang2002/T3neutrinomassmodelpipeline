"""Static contract tests for the v2 Wolfram runner; not native-Wolfram evidence."""
from pathlib import Path
import re

RUNNER = (Path(__file__).resolve().parents[1]
          / "RGE/running/rgbeta/RunT3BCFullScalarMassRGBeta.wl")


def code_without_comments(source: str) -> str:
    """Drop Wolfram (* ... *) blocks for source-level assertions."""
    return re.sub(r"\(\*.*?\*\)", "", source, flags=re.S)


def test_frontend_context_and_non_message_based_beta_evaluation():
    text = code_without_comments(RUNNER.read_text(encoding="utf-8"))
    assert 'Catch[UsingFrontEnd[T3BHMRun[]],"T3BHM"]' in text
    assert 'Quiet[BetaTerm[sym,1]]' in text
    assert 'Check[BetaTerm[sym,1]' not in text
    assert 'Check[Catch[T3BHMRun[]' not in text


def test_mass_sector_and_quartic_registration_unchanged():
    text = code_without_comments(RUNNER.read_text(encoding="utf-8"))
    assert 'T3RGBetaBuild[2,2,dF,-1]' in text
    assert 'T3RGBetaBuildEFT1[2,2,dF,-1]' in text
    assert 'AddScalarMass[mHSq,{Bar@H,H}' in text
    assert 'AddScalarMass[muT3GeneralS1S2,{Bar@S1,Bar@S2}' in text
    assert 'Get[extension]' in text
    assert 'T3GeneralQuarticReusedCount =!= 6' in text
    assert 'lambdaT3GeneralQ110101I1' in text
    assert 'massSymbols = {mHSq,mS1Sq,mS2Sq,muT3GeneralS1S2}' in text


def test_blocked_diagnostics_written_without_false_success():
    text = code_without_comments(RUNNER.read_text(encoding="utf-8"))
    assert '"beta_diagnostics"' in text
    assert '"unresolved" -> MemberQ[failedBetas,name]' in text
    assert '"status"->"Blocked"' in text
    assert 'If[Lookup[result,"status",""] === "Blocked"' in text
    assert 'Exit[5]' in text
    assert '"physical_mapping_certified"->False' in text
    assert '"format"->"workstream-b-bc-full-scalar-mass-live-rgbeta-v1"' in text


def test_runner_has_no_overwrite_and_no_production_model_writes():
    text = code_without_comments(RUNNER.read_text(encoding="utf-8"))
    assert 'If[FileExistsQ[out]' in text
    assert 'T3RGBetaModel.wl' in text
    assert 'T3RGBetaWriteJSON' not in text  # export remains local opt-in
    assert 'Export[out,result,"RawJSON"]' in text
