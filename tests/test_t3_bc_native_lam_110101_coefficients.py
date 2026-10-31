"""Schema regressions, NOT a substitute for a live wolframscript run."""
from __future__ import annotations

from copy import deepcopy
from pathlib import Path
import pytest

from scripts.VerifyT3BCNativeLam110101Coefficients import (
    FORMAT, MODELS, STATUS, SYMBOLS, verify,
)


def fixture():
    return {
        'format': FORMAT,
        'status':'LiveLam110101SourceCoefficientComponentAudit',
        'source':'installed RGBeta`PackageScope`Lam',
        'action_level_matchete_independently_evaluated':False,
        'physical_coupling_mapping_certified':False,
        'physical_matching_certified':False,
        'models':{
            model:{
                'model':model,'status':STATUS,
                'live_native_lam_called':True,
                'symbol':'RGBeta`PackageScope`Lam',
                'source_symbols':list(SYMBOLS),
                'source_scalar_field_positions':[[1,2,4,6],[2,1,3,5]],
                'registered_complex_quartic_normalization':'6*Matrix[...]',
                'unique_arrangements_per_coupling':24,
                'component_tests':64,'expected_component_tests':64,
                'expected_lam_over_registered_invariant':'1/4',
                'coefficient_extraction':'symbolic Coefficient on evaluated native Lam',
                'native_evaluated_coupling_forms':[
                    {'orientation':orientation, 'coupling':symbol,
                     'evaluated_coupling_form':(
                        'Bar['+symbol+']' if orientation=='conjugate' else symbol)}
                    for orientation in ('direct','conjugate') for symbol in SYMBOLS
                ],
                'failed_component_count':0,
                'first_mismatches':[],
                'maximum_absolute_component_residual':0.,
                'matchete_runtime_evaluated':False,
                'physical_coupling_mapping_certified':False,
                'physical_matching_certified':False,
            } for model in MODELS
        }
    }


def test_valid_fixture():
    result=verify(fixture())
    assert result['total_components_checked']==128
    assert result['physical_coupling_mapping_certified'] is False


@pytest.mark.parametrize('change',[
    lambda x: x['models']['B_-1'].update(status='LiveLam110101SourceComponentCoefficientsMismatch'),
    lambda x: x['models']['B_-1'].update(component_tests=32),
    lambda x: x['models']['B_-1'].update(failed_component_count=1),
    lambda x: x['models']['B_-1'].update(maximum_absolute_component_residual=1e-7),
    lambda x: x['models']['B_-1'].update(registered_complex_quartic_normalization='24*Matrix[...]'),
    lambda x: x['models']['C_-1'].update(unique_arrangements_per_coupling=12),
    lambda x: x['models']['C_-1'].update(expected_lam_over_registered_invariant='6'),
    lambda x: x['models']['C_-1'].update(physical_coupling_mapping_certified=True),
    lambda x: x.update(action_level_matchete_independently_evaluated=True),
    lambda x: x['models']['C_-1'].update(live_native_lam_called=False),
    lambda x: x['models']['C_-1'].update(coefficient_extraction='unverified text search'),
    lambda x: x['models']['B_-1'].update(native_evaluated_coupling_forms=[]),
])
def test_reject_bad_or_overclaimed_evidence(change):
    d=deepcopy(fixture());change(d)
    with pytest.raises(ValueError):verify(d)


def test_script_uses_live_lambda_and_component_checks():
    root=Path(__file__).resolve().parents[1]
    path=root/'scripts'/'ProbeT3BCNativeLam110101Coefficients.wl'
    script=path.read_text(encoding='utf-8')
    for text in ('RGBeta`PackageScope`Lam','T3BCSelectCoupling',
                 'Extract[array,positions[[slot]]]','T3BCComponent',
                 't3a,t3b,t3c,t3d','"physical_matching_certified"->False',
                 '"failed_component_count"->Length[failures]'):
        assert text in script
    assert 'If[FileExistsQ[outFile]' in script


def test_source_runtime_extraction_v2_is_fail_closed_and_saves_diagnostics():
    script=(Path(__file__).resolve().parents[1]/'scripts'/
            'ProbeT3BCNativeLam110101Coefficients.wl').read_text(encoding='utf-8')
    assert 'T3BCLookupDiagnostic' in script
    assert '"requested_entry_inputform"->T3BCTruncate[entry,4000]' in script
    assert '"coupling_bearing_sparse_rules_excerpt"' in script
    assert '"matrix_symbol_names"->Names["*`Matrix"]' in script
    assert 'Coefficient[Expand[expression],#,1]' in script
    assert 'T3BCSelectCoupling[entry,symbols[[k]],conjFlag]' in script
    assert 'Expand[entry-Total[selectedTerms]] == 0' in script
    assert '"status"->"BlockedNativeLamComponentExtraction"' in script
    assert 'Catch[T3BCProbeModel[spec[[1]],spec[[2]]],"T3BCCoeff"]' in script
    assert 'Export[outFile,result,"RawJSON"]' in script
    assert 'Exit[7]' in script
    assert '"physical_coupling_mapping_certified"->False' in script
    assert '"physical_matching_certified"->False' in script


def test_source_extracts_live_bare_coefficients_and_enforces_hc():
    script=(Path(__file__).resolve().parents[1]/'scripts'/
            'ProbeT3BCNativeLam110101Coefficients.wl').read_text(encoding='utf-8')
    assert 'Coefficient[Expand[expression],#,1]' in script
    assert '{Bar[target],Conjugate[target]}' in script
    assert 'Length[which] =!= 1' in script
    assert 'native_evaluated_coupling_forms' in script
    assert 'Native Lam selected entry has extra/nonlinear coupling terms' in script
    assert '"physical_matching_certified"->False' in script


def test_wolfram_outer_do_iterator_is_an_argument_not_separate_statement():
    """Catch the actual v3 failure: `; {slot,1,2}` leaves slot$ unbound."""
    import re
    script=(Path(__file__).resolve().parents[1]/'scripts'/
            'ProbeT3BCNativeLam110101Coefficients.wl').read_text(encoding='utf-8')
    assert 'Do[\n    conjFlag=(slot==2);' in script
    assert re.search(r'\]\],\s*\{slot,1,2\}\]',script), (
        'The outer Do body MUST end in a comma before its iterator. '
        'A semicolon leaves slot symbolic in positions[[slot]] and Extract.'
    )
    assert not re.search(r';\s*\{slot,1,2\}\]',script)
    assert 'entry=Extract[array,positions[[slot]]]' in script
    assert 'ok=(count==64 && failures=={})' in script
