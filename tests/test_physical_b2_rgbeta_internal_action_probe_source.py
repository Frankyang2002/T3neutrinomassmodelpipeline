"""Static regression guards for isolated B2 Wolfram native probe.

These checks do not execute Mathematica and do not certify physical action
normalisations.  Native evaluation still requires wolframscript on Windows.
"""
from pathlib import Path

SOURCE = (Path(__file__).resolve().parents[1] / "scripts" /
          "RunPhysicalB2RGBetaInternalActionProbe.wl")


def source() -> str:
    return SOURCE.read_text(encoding="utf-8")


def test_downvalues_gets_resolved_symbol_not_unevaluated_toexpression():
    script = source()
    assert 'B4GroupInvDownValuesString[names_List]' in script
    assert 'With[{symbol=ToExpression[First[names]]}, B4String[DownValues[symbol]]]' in script
    assert 'DownValues[ToExpression[' not in script.split('B4GroupInvDownValuesString[names_List] :=', 1)[1]


def test_missing_or_ambiguous_internal_function_does_not_fabricate_definition():
    script = source()
    assert 'If[Length[names]===1,' in script
    assert '"unresolved"];' in script
    assert '"GroupInvBar_symbols"->$B4GroupInv' in script


def test_all_expected_registry_operators_still_probed():
    script = source()
    for symbol in ('y1', 'y2', 'MF', 'lambdaT3', 'lambdaH1', 'lambdaH2'):
        assert 'B4CouplingRow[$' in script
        assert f',{symbol},' in script
    assert 'AllSixRegistryRatiosMeasured' in script


def test_no_physical_certificate_from_internal_factors():
    script = source()
    for token in (
        '"internal_prefactor_is_physical_action_factor"->False',
        '"numeric_component_axes_verified"->False',
        '"physical_source_certified"->False',
        '"physical_Lagrangian_coefficient_verified"->False',
    ):
        assert token in script
