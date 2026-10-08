import json
from pathlib import Path
from Numerical.diagnostics.VerifyT3CGCovariance import audit, generators, tensor_residual

DATA = Path(__file__).resolve().parents[1] / 'Reports/output/full/multiplicity'

def test_generator_algebra():
    gs = generators(json.loads((DATA/'matchete_generator_probe.json').read_text()))
    assert set(gs) == {1,2,3}

def test_all_exported_cg_tensors_are_covariant():
    component = json.loads((DATA/'component_support_audit.json').read_text())
    probe = json.loads((DATA/'matchete_generator_probe.json').read_text())
    result = audit(component, probe)
    assert result['count'] == 16
    assert all(x['all_covariant'] for x in result['models'])
    assert all(x['verified_loop_multiplicity'] is None for x in result['models'])

def test_corrupt_tensor_is_not_covariant():
    component = json.loads((DATA/'component_support_audit.json').read_text())
    probe = json.loads((DATA/'matchete_generator_probe.json').read_text())
    t = dict(component['cases'][0]['vertices']['T3MixCG'])
    t['components'] = [dict(x) for x in t['components']]
    t['components'][0]['cg'] = '99'
    assert tensor_residual(t, (False,False,False,True), generators(probe))
