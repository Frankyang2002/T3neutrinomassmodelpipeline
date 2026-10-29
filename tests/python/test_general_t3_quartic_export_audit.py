"""Stage21 complete catalogue export auditing using controlled fixtures."""
from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import sys

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from scripts import RunGeneralT3QuarticExportAudit as stage21


@pytest.fixture
def registry():
    return json.loads((ROOT/'RGE/general/generated/general_rgbeta_quartic_registry.json').read_text())


def _make_export(model,stage,entry,*,extra_betas=None):
    added = [v['name'] for v in entry['new_rgbeta_symbols']]
    reused = [v['existing_name'] for v in entry['preexisting_rgbeta_symbols']]
    basic = ['gY','g2','g3','lambdaH','lambdaS1','lambdaS2','lambdaH1','lambdaH2','lambda12','lambdaT3']
    allnames = basic + added + reused + (['y1','y2','MF'] if stage == 'UV' else [])
    allnames = list(dict.fromkeys(allnames))
    betas = {k:'0' for k in allnames}
    betas.update(extra_betas or {})
    meta = {
        'GeneralScalarModel':model, 'GeneralScalarStage':stage,
        'GeneralScalarBasisStatus':'CandidateNativeRGBetaUnmappedToMatchete',
        'GeneralScalarQuarticBetaBasis':'RGBetaFullQuarticProjection',
        'ReportBetaConvention':stage21._CONVENTION,
        'GeneralScalarQuadraticMassMixing':'NotIncluded',
        'NativeRGBetaCheckInvarianceUsed':False,
        'GeneralScalarAdditionalQuartics':entry['total_additional_quartics'],
        'GeneralScalarExistingStructuralCount':len(reused),
        'GeneralScalarNewBetaKeys':added,
        'GeneralScalarExistingQuarticBetasProjected':basic[3:] + reused,
        'dS1':entry['dimensions'][0],'dS2':entry['dimensions'][1],
        'dF':entry['dimensions'][2],'alpha':entry['alpha'],
    }
    return {'status':'Success','metadata':meta,'betas':betas,'report_betas':dict(betas)}


@pytest.fixture
def exports(tmp_path,registry):
    p=tmp_path/'exports'
    p.mkdir()
    for name,model in registry['models'].items():
        for stage in stage21.STAGES:
            (p/f'{name}_{stage.lower()}_rgbeta.json').write_text(json.dumps(_make_export(name,stage,model)))
    return p


def _write(p,obj):
    p.write_text(json.dumps(obj))


def test_registry_has_all_50_model_specific_quartic_slots(registry):
    assert len(registry['models'])==16
    assert sum(v['total_additional_quartics'] for v in registry['models'].values())==50
    assert sum(len(v['new_rgbeta_symbols']) for v in registry['models'].values())==38
    assert sum(len(v['preexisting_rgbeta_symbols']) for v in registry['models'].values())==12


def test_synthetic_32_exports_structurally_pass(exports,tmp_path):
    result=stage21.audit_catalogue(export_dir=exports,output_dir=tmp_path/'report')
    assert result['export_count']==32
    assert result['stage_export_counts']=={'UV':16,'EFT1':16}
    assert result['inventory']['matchete_extra_quartics']==50
    assert result['matchete_rgbeta_normalisation_certified'] is False
    assert result['quadratic_mass_mixing_included'] is False
    assert result['numeric_probe_status']=='NotRequested'
    assert (tmp_path/'report/stage21_summary.json').is_file()


def test_missing_export_is_fatal(exports,registry):
    (exports/'B_-1_eft1_rgbeta.json').unlink()
    with pytest.raises(ValueError,match='Missing exported'):
        stage21.audit_catalogue(export_dir=exports,output_dir=None)


def test_missing_reused_quartic_is_fatal(exports,registry):
    p=exports/'B_-1_uv_rgbeta.json'
    payload=json.loads(p.read_text())
    key=registry['models']['B_-1']['preexisting_rgbeta_symbols'][0]['existing_name']
    del payload['report_betas'][key]
    del payload['betas'][key]
    _write(p,payload)
    with pytest.raises(ValueError,match='missing required general quartic betas'):
        stage21.audit_catalogue(export_dir=exports,output_dir=None)


def test_unresolved_rgbeta_beta_is_fatal(exports):
    p=exports/'A_+0_uv_rgbeta.json'
    payload=json.loads(p.read_text())
    key=next(iter(payload['betas']))
    payload['betas'][key]=payload['report_betas'][key]='BetaTerm[lambdaX, 1]'
    _write(p,payload)
    with pytest.raises(ValueError,match='unresolved RGBeta expressions'):
        stage21.audit_catalogue(export_dir=exports,output_dir=None)


def test_projection_flag_is_required(exports):
    p=exports/'A_+0_eft1_rgbeta.json'
    payload=json.loads(p.read_text())
    payload['metadata']['GeneralScalarQuarticBetaBasis']='IndependentBetaTerm'
    _write(p,payload)
    with pytest.raises(ValueError,match='GeneralScalarQuarticBetaBasis mismatch'):
        stage21.audit_catalogue(export_dir=exports,output_dir=None)


def test_beta_basis_must_have_same_uv_eft1_scalar_names(exports):
    p=exports/'A_+0_eft1_rgbeta.json'
    payload=json.loads(p.read_text())
    payload['betas']['extraVar']='0'
    payload['report_betas']['extraVar']='0'
    _write(p,payload)
    with pytest.raises(ValueError,match='unexpected EFT1-only'):
        stage21.audit_catalogue(export_dir=exports,output_dir=None)


def test_wrong_model_and_representation_are_rejected(exports):
    p=exports/'B_-1_uv_rgbeta.json'
    payload=json.loads(p.read_text())
    payload['metadata']['alpha']=1
    _write(p,payload)
    with pytest.raises(ValueError,match='alpha mismatch'):
        stage21.audit_catalogue(export_dir=exports,output_dir=None)


def test_synthetic_parser_probe_passes_clean_zeros(exports,tmp_path):
    result=stage21.audit_catalogue(export_dir=exports,output_dir=tmp_path,
        numeric_probe=True,evaluator=lambda expression, environment: 0.)
    assert result['numeric_parser_failures']==0
    assert all(row['numerical_parser_probe']['synthetic_inputs_only'] for row in result['rows'])


def test_parser_probe_reports_errors_without_misreporting_physical_matching(exports):
    result=stage21.audit_catalogue(export_dir=exports,output_dir=None,numeric_probe=True,
        evaluator=lambda expression, environment: (_ for _ in ()).throw(ValueError('synthetic error')))
    assert result['numeric_parser_failures']>0
    assert result['numeric_probe_status']=='ParserDiagnosticsRequireReview'
    assert result['matchete_rgbeta_normalisation_certified'] is False


def test_expressions_need_same_raw_and_report_keys(exports):
    p=exports/'B_-1_uv_rgbeta.json'
    payload=json.loads(p.read_text())
    del payload['betas']['lambdaH']
    _write(p,payload)
    with pytest.raises(ValueError,match='raw/report beta keys differ'):
        stage21.audit_catalogue(export_dir=exports,output_dir=None)
