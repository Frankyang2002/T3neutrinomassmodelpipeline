"""Static regression of the RGBeta live probe's *verified source interface*.

These checks cannot replace Wolfram/RGBeta execution. They reject the
previous held-symbol, package-key and Wolfram pattern-index regressions; they
verify that blocked native evaluations retain per-model diagnostic evidence.
"""
from __future__ import annotations

from pathlib import Path

PROBE = (Path(__file__).resolve().parents[1] / "scripts" /
         "ProbeT3BCMixedBilinearRGBeta.wl")


def text() -> str:
    return PROBE.read_text(encoding="utf-8")


def test_downvalues_resolves_a_literal_package_symbol():
    source = text()
    assert 'lamDefinitionCount=Length[DownValues[RGBeta`PackageScope`Lam]]' in source
    assert 'DownValues[Symbol[lamName]]' not in source.replace(
        '(* DownValues has a holding attribute: DownValues[Symbol[name]] is invalid.', '')
    assert 'RGBeta`PackageScope`Lam[t3BMa,t3BMb,t3BMc,t3BMd,True]' in source


def test_native_symbol_keys_use_unique_runtime_key_resolution():
    source = text()
    for name in ('Fields', 'UniqueArrangements', 'SelfConjugate', 'Invariant', 'Coupling'):
        assert f'T3BRNativeGet[label,info,"{name}"]' in source
    for unsafe in ('info[Fields]', 'info[UniqueArrangements]', 'info[Invariant]',
                   'info[Coupling]', 'info[SelfConjugate]'):
        assert unsafe not in source
    assert 'Length[matches]=!=1' in source
    assert 'Lookup[info,First[matches]' in source


def test_fail_closed_and_separate_native_evidence():
    source = text()
    assert 'FileExistsQ[out]' in source
    assert 'KeyExistsQ[$scalarMasses,muT3GeneralS1S2]' in source
    assert '"native_lam_mass_sparse_array_rules_count"->Length[relevant]' in source
    assert '"native_rgbeta_epsilon_component_phase_verified"->False' in source
    assert '"physical_mapping_certified"->False' in source
    assert '"B_-1"->T3BRSafeOne["B_-1",1]' in source
    assert '"C_-1"->T3BRSafeOne["C_-1",3]' in source
    assert 'report=CheckAbort[<|' in source
    assert 'result=CheckAbort[' in source
    assert 'Check[Catch[T3BROne[label,dF],"T3BRMixed"],$Failed]' in source
    assert 'model diagnostic JSON saved' in source


def test_live_lam_arguments_are_plain_symbols_not_wolfram_patterns():
    """Underscores create patterns in Wolfram, which are not valid RHS indices."""
    import re
    source = text()
    call = re.search(r"RGBeta`PackageScope`Lam\[([^]]+)\]", source)
    assert call is not None
    args = [a.strip() for a in call.group(1).split(",")]
    assert args == ["t3BMa", "t3BMb", "t3BMc", "t3BMd", "True"]
    assert all(re.fullmatch(r"[A-Za-z][A-Za-z0-9]*", a) for a in args)
    # Verify the negative control: the exact old spelling must fail the check.
    assert re.fullmatch(r"[A-Za-z][A-Za-z0-9]*", "t3BM_a") is None


def test_native_lam_messages_are_contained_per_model():
    source = text()
    assert 'T3BRSafeOne[label_String,dF_Integer]' in source
    assert '"model"->label,"status"->"Blocked"' in source
    assert '"physical_mapping_certified"->False' in source
    assert '"probe_revision"->"inspect-native-sparse-rules-with-block-diagnostics-v4"' in source
    assert 'report=CheckAbort[Check[<|' not in source


def test_sparse_array_native_rules_are_the_only_positive_lam_witness():
    source = text()
    assert 'arrays=Cases[live,_SparseArray,{0,Infinity}]' in source
    assert 'rules=Flatten[ArrayRules/@arrays,1]' in source
    assert 'relevant=Select[rules,' in source
    assert 'StringContainsQ[ToString[#,InputForm],"muT3GeneralS1S2"]' in source
    assert 'If[Length[relevant]<1,' in source
    # The previous check incorrectly rejected opaque native TStructure.
    assert 'FreeQ[live,muT3GeneralS1S2,{0,Infinity}] ||' not in source
    assert 'If[live===$Failed || live===$Aborted ||' not in source


def test_blocked_model_preserves_runtime_sparse_diagnostics():
    source = text()
    assert 'T3BRFail[model_,why_,details_:<||>]' in source
    assert '"native_lam_sparse_array_count"->Length[arrays]' in source
    assert '"native_lam_total_sparse_array_rules_count"->Length[rules]' in source
    assert '"native_lam_mass_sparse_array_rules_count"->Length[relevant]' in source
    assert '"native_lam_excerpt_inputform"->' in source
    assert '"native_lam_runtime_diagnostics"->diag' in source
    assert '"native_mass_registry_keys_inputform"->' in source
    assert 'T3BRFail[label,\n    "Native massive Lam has no mixed-mass coupling in inspected sparse rules",\n    diag]' in source
