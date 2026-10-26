"""Stage 28 tests. Synthetic native metadata never count as live evidence."""
from copy import deepcopy
import json
from pathlib import Path

import pytest

from scripts.Stage27TripletIntertwiner import audit as stage27_audit
from scripts.Stage28TripletNativeStructureAudit import ALL, EXPECTED, report_for, run

ROOT = Path(__file__).resolve().parents[2]
FOLDER = ROOT/'output/general_t3_stage22'


def witness():
    a = FOLDER/'matchete_tensors.json'
    b = FOLDER/'rgbeta_native_invariants.json'
    if not a.is_file() or not b.is_file():
        pytest.skip('Native Stage22 Matchete and RGBeta tensor fixture absent')
    return stage27_audit(json.loads(a.read_text()), json.loads(b.read_text()))


def fake_native():
    # A strict schema fixture: this does not assert any live Wolfram outcome.
    rows=[]
    for model in ALL:
        dims, names = EXPECTED[model]
        rows.append(dict(model=model, scalar_su2_dimensions=list(dims),
                         triplet_fields=list(names),native_S2_representation='SU2L[S2]',
                         dimension=3, quadratic_casimir=2, trace_normalisation=2,
                         native_generator_inputform='tGen[SU2L[S2], t28A, t28x, t28y]',
                         native_generator_is_symbolic=True,
                         native_conjugate_generator_minus_transpose=True,
                         native_delS2_orthogonality=True,
                         native_complex_conjugate_metric_two=True,
                         native_numerical_generator_components_exported=False))
    return dict(format='stage28-live-rgbeta-s2-structure-v1', models_checked=6,
                models=rows, all_group_representation_checks_passed=True,
                all_native_barred_generator_rules_passed=True,
                all_native_delS2_contractions_passed=True,
                native_numerical_triplet_generators_available=False,
                physical_cross_package_component_basis_certified=False,
                physical_coupling_mapping_certified=False)


def test_metric_and_generator_invariance_exact():
    result=report_for(witness(),fake_native())
    assert result['matchete_triplet_epsilon_invariant_metric']==[
        ['-1','0','0'],['0','-1','0'],['0','0','-1']]
    assert result['matchete_triplet_invariant_metric_and_generators_exact'] is True


def test_six_model_and_all_triplet_field_counts():
    result=report_for(witness(),fake_native())
    assert result['models_checked']==6
    assert result['native_triplet_fields_in_six_models']==8


def test_never_promotes_numerical_component_certificate():
    result=report_for(witness(),fake_native())
    assert result['physical_coupling_mapping_certified'] is False
    assert result['native_rgbeta_numerical_triplet_generators_exported'] is False
    assert result['barred_triplet_component_intertwiner_physically_certified'] is False


@pytest.mark.parametrize('field,value',[
    ('dimension',4),('quadratic_casimir',3),('trace_normalisation',1),
    ('native_generator_is_symbolic',False),
    ('native_conjugate_generator_minus_transpose',False),
    ('native_delS2_orthogonality',False),
    ('native_complex_conjugate_metric_two',False),
    ('native_numerical_generator_components_exported',True),
])
def test_reject_tampered_native_source_claim(field,value):
    data=fake_native()
    data['models'][0][field]=value
    with pytest.raises(ValueError,match='Stage28 fail-closed'):
        report_for(witness(),data)


def test_reject_wrong_triplet_species():
    data=fake_native(); data['models'][0]['triplet_fields']=['S1']
    with pytest.raises(ValueError,match='assignment'):
        report_for(witness(),data)


def test_reject_false_certificate():
    data=fake_native();data['physical_coupling_mapping_certified']=True
    with pytest.raises(ValueError,match='falsely certified'):
        report_for(witness(),data)


def test_reject_duplicate_model():
    data=fake_native(); data['models'][1]['model']=data['models'][0]['model']
    with pytest.raises(ValueError,match='repeated'):
        report_for(witness(),data)


def test_fail_closed_when_live_export_missing(tmp_path):
    with pytest.raises(ValueError,match='native Stage28 JSON'):
        run(tmp_path,None,check_stage27=False)


def test_real_live_stage28_and_stage27_if_present():
    if not (FOLDER/'stage28_live_rgbeta_triplet_structure.json').is_file() or not \
       (FOLDER/'stage27_triplet_generator_witness.json').is_file() or not \
       (FOLDER/'stage26_conditional_quartic_action.json').is_file():
        pytest.skip('User Wolfram output and complete Stage25-27 evidence needed')
    result=run(FOLDER,None)
    assert result['stage27_fresh_evidence_checked'] is True
    assert result['physical_coupling_mapping_certified'] is False
