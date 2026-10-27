"""Offline contracts for the experimental direct one-loop probe (no live Matchete)."""
from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from scripts import RunGeneralT3MassEigenOneLoop as probe


def passing_report(model: str, stem: str) -> dict:
    return {
        'status': probe.STATUS, 'model': model,
        'mass_basis_ordinary_check_passed': True,
        'heavy_mass_basis_detailed_check_passed': True,
        'one_loop_match_stage_completed': True,
        'weinberg_extraction_status': 'Success',
        'matching_EFT_order': 5, 'matching_loop_order': 1,
        'matching_approach': 'SimultaneousAllHeavyOneScaleExploratory',
        'weinberg_holomorphic_term_count': 2,
        'numerical_c5_computed': False,
        'production_sequential_matching_performed': False,
        'physical_matching_performed': False,
        'ordered_coefficient_is_canonical_flavour_symmetric_c5': False,
        'matchete_rgbeta_quartic_normalization_certified': False,
        'matched_eft_inputform_file': stem + '_matched_eft.inputform.txt',
        'weinberg_ordered_coefficient_inputform_file': stem + '_weinberg_ordered_coefficient.inputform.txt',
    }


@pytest.mark.parametrize('model',['B_-1','C_-1'])
def test_stem(model):
    assert probe._stem(model,'final') == model.replace('-','m') + '_final'


@pytest.mark.parametrize('which',['bad','',None])
def test_bad_endpoint(which):
    with pytest.raises(ValueError,match='Unsupported'):
        probe._stem('B_-1',which)


def test_probe_report_requires_both_exported_artifacts(tmp_path):
    stem='B_m1_final_one_loop_probe'
    report=passing_report('B_-1',stem)
    for suffix in ('_matched_eft.inputform.txt',
                   '_weinberg_ordered_coefficient.inputform.txt'):
        (tmp_path/(stem+suffix)).write_text('Nonempty symbolic Matchete result',encoding='utf-8')
    a,b=probe.validate_probe_report(report,'B_-1',tmp_path,stem)
    assert a.name.endswith('_matched_eft.inputform.txt')
    assert b.name.endswith('_weinberg_ordered_coefficient.inputform.txt')


@pytest.mark.parametrize('field,value',[
    ('status','MatcheteMassEigenHeavyBasisPassedPhysicalMatchingNotPerformed'),
    ('model','C_-1'),
    ('one_loop_match_stage_completed',False),
    ('weinberg_extraction_status','NotFound'),
    ('weinberg_holomorphic_term_count',0),
    ('numerical_c5_computed',True),
    ('production_sequential_matching_performed',True),
    ('physical_matching_performed',True),
    ('ordered_coefficient_is_canonical_flavour_symmetric_c5',True),
    ('matchete_rgbeta_quartic_normalization_certified',True),
    ('matching_EFT_order',6),
])
def test_probe_report_refuses_bad_certifications(tmp_path,field,value):
    stem='B_m1_final_one_loop_probe'
    report=passing_report('B_-1',stem)
    report[field]=value
    with pytest.raises(ValueError):
        probe.validate_probe_report(report,'B_-1',tmp_path,stem)


def test_probe_report_refuses_missing_artifact(tmp_path):
    stem='B_m1_final_one_loop_probe'
    report=passing_report('B_-1',stem)
    with pytest.raises(ValueError,match='missing/empty'):
        probe.validate_probe_report(report,'B_-1',tmp_path,stem)


def test_wolfram_probe_contains_matching_and_safe_status():
    content=Path(__file__).resolve().parents[2]/'Lagrangian/RunGeneralT3MassEigenOneLoop.wl'
    text=content.read_text(encoding='utf-8')
    assert 'RunT3Matching[lag,5,1]' in text
    assert 'ExtractWeinbergCoefficient[matchedEFT]' in text
    assert 'CanonicalCoefficient' in text  # Explicitly not used for masses MSX,MSY.
    assert '"SimultaneousAllHeavyOneScaleExploratory"' in text
    assert '"numerical_c5_computed"->False' in text
    assert '"production_sequential_matching_performed"->False' in text
    assert '"status"->"OneLoopDirectSymbolicWeinbergProbePassedNumericalC5Pending"' in text
    assert 'If[!TrueQ[check] || !heavy' in text
    # The proven native operator reconstructions must not be dropped.
    for symbol in ('T3MEBuildNativeYukawa1','T3MEBuildNativeYukawa2',
                   'T3MEBuildNativeExtraQuartics','T3MEBuildNativeS2ScalarQuartics',
                   'T3MEBuildNativeH2Portals','T3MEBuildNativeTopology'):
        assert symbol in text


def test_run_model_refuses_uncertified_previous_check(tmp_path,monkeypatch):
    monkeypatch.setattr(probe,'prepare_model',lambda *a,**k: {
        'payload':{'unitary_U':{}}, 'parameters': tmp_path/'parameters.json'})
    monkeypatch.setattr(probe,'_prior_certificate',lambda *a,**k: (_ for _ in ()).throw(ValueError('missing certificate')))
    with pytest.raises(ValueError, match='missing certificate'):
        probe.run_model('B_-1',outdir=tmp_path/'out')


def test_run_model_mock_success_respects_output_and_provenance(tmp_path,monkeypatch):
    root=tmp_path/'repo'
    (root/'Lagrangian').mkdir(parents=True)
    (root/'output').mkdir()
    p=root/'Lagrangian/RunGeneralT3MassEigenOneLoop.wl';p.write_text('WL content')
    m=root/'Lagrangian/RunMatching.wl';m.write_text('production matching content')
    inputpath=root/'output/parameters.json';inputpath.write_text('{}')
    monkeypatch.setattr(probe,'ROOT',root)
    monkeypatch.setattr(probe,'WOLFRAM',p)
    monkeypatch.setattr(probe,'MATCHING_SOURCE',m)
    monkeypatch.setattr(probe,'prepare_model',lambda *a,**k:{'parameters':inputpath,'payload':{'unitary_U':{}},'masses':[1.,2.]})
    certificate=root/'output/prior_matchete.json';certificate.write_text('{}')
    monkeypatch.setattr(probe,'_prior_certificate',lambda *a,**k:certificate)
    monkeypatch.setattr(probe.shutil,'which',lambda *a,**k:'wolframscript')
    def fake_run(command,**kwargs):
        target=Path(command[-1]);stem=target.stem
        report=passing_report('B_-1',stem)
        target.write_text(json.dumps(report))
        for suffix in ('_matched_eft.inputform.txt',
                       '_weinberg_ordered_coefficient.inputform.txt'):
            (target.parent/(stem+suffix)).write_text('symbolic coefficient')
        return SimpleNamespace(returncode=0,stdout='Matched',stderr='')
    monkeypatch.setattr(probe.subprocess,'run',fake_run)
    result=probe.run_model('B_-1',outdir=root/'output/probe')
    assert result['numerical_c5_computed'] is False
    assert result['production_sequential_matching_performed'] is False
    assert result['status']==probe.STATUS
    assert (root/'output/probe/B_m1_final_one_loop_probe.log').is_file()
