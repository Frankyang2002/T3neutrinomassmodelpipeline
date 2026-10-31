"""Schema/claim tests. Actual installed RGBeta runtime must be run locally."""
import copy
import pytest
from scripts.VerifyT3BCNativeLamProbe import verify, EXPECTED, FORMAT


def fixture():
    doc = {'format':FORMAT,'source':'installed RGBeta runtime Lam',
           'physical_coupling_mapping_certified':False,
           'physical_matching_certified':False,'models':{}}
    for model in ('B_-1','C_-1'):
        doc['models'][model] = {
            'model':model,'actual_installed_rgbeta_lam_called':True,
            'native_registered_symbols':list(EXPECTED),
            'lam_action_potential_reconstructed':False,
            'physical_coupling_mapping_certified':False,
            'matchete_runtime_evaluated':False,
            'lam_expression_sha256':'0'*64,
            'native_coupling_factors':[f'6*Matrix[{name}][##1] & ' for name in EXPECTED],
            'native_self_conjugate':[False,False],
            'status':'NativeLamExpressionAndSparseTensorCaptured',
            'sparse_array_count':1,
            'lam_contains_couplings':{name:True for name in EXPECTED},
        }
    return doc


def test_fixture_schema():
    out = verify(fixture())
    assert len(out['results']) == 2
    assert out['physical_coupling_mapping_certified'] is False


@pytest.mark.parametrize('edit', [
    lambda d: d['models']['B_-1'].update({'actual_installed_rgbeta_lam_called':False}),
    lambda d: d['models']['C_-1'].update({'sparse_array_count':0}),
    lambda d: d['models']['B_-1'].update({'lam_contains_couplings':{EXPECTED[0]:True}}),
    lambda d: d['models']['C_-1'].update({'native_self_conjugate':[False,True]}),
    lambda d: d['models']['B_-1'].update({'native_coupling_factors':['24*Matrix[x]', '6*Matrix[y]']}),
    lambda d: d['models']['B_-1'].update({'physical_coupling_mapping_certified':True}),
    lambda d: d['models']['B_-1'].update({'status':'NativeLamEvaluationRequiresInspection'}),
    lambda d: d.update({'physical_matching_certified':True}),
])
def test_reject_changed_or_incomplete_evidence(edit):
    d=copy.deepcopy(fixture()); edit(d)
    with pytest.raises(ValueError): verify(d)


def test_wolfram_probe_uses_declared_package_scope():
    """Prevent regression to the broken main/private-context lookup."""
    from pathlib import Path
    source = (Path(__file__).resolve().parents[1] / 'scripts' /
              'ProbeT3BCNativeRGBetaLam.wl').read_text(encoding='utf-8')
    assert 'lamName = "RGBeta`PackageScope`Lam";' in source
    assert 'Names[lamName] =!= {lamName}' in source
    assert 'Length[DownValues[f]]' in source
    assert 'expr=CheckAbort[Check[lamFunction[' in source
    assert '"lam_symbol"->lamName' in source
    assert 'Names["RGBeta`Private`*"]' not in source
