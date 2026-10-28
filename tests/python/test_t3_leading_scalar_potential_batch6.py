"""Batch 6: self-contained tests of cross-audit synthesis and fail-closed inputs."""
from __future__ import annotations

from copy import deepcopy
import hashlib
import json

import pytest

from scripts import FinalizeT3LeadingScalarPotential as a

EXTRA = {'A_-4':0,'A_-2':5,'A_+0':5,'B_-3':1,'B_-1':7,'B_+1':1,
         'C_-3':1,'C_-1':7,'C_+1':1,'D_-2':5,'D_+0':5,'D_+2':0,
         'E_-4':0,'E_-2':6,'E_+0':6,'E_+2':0}
MASS = {'A_+0':'self_s1','A_-2':'self_s2','B_-1':'mixed','C_-1':'mixed',
        'D_+0':'self_s1','D_-2':'self_s2','E_+0':'self_s1','E_-2':'self_s2'}


def fixture_reports():
    ops=[]
    per={}
    rows3=[]
    rows4={}
    rows5={}
    for m in a.MODELS:
        n=EXTRA[m]
        for i in range(1):
            ops.append({'model':m,'kind':'baseline_topology','coupling':'lambdaT3'})
        for i in range(2):
            ops.append({'model':m,'kind':'baseline_quadratic','coupling':'MS'+str(i)})
        for i in range(5):
            ops.append({'model':m,'kind':'baseline_quartic_family','coupling':'base'+str(i)})
        for i in range(n):
            direct=m in a.DIRECT and i<2
            ops.append({'model':m,'kind':'additional_quartic','coupling':'additionalQ'+str(i),
                        'higgs_legs':2 if direct else 0,'observed_in_saved_C5':direct})
        kind=MASS.get(m,'none')
        if kind!='none':
            ops.append({'model':m,'kind':'additional_quadratic','coupling':'mu2_extra'})
        per[m]={'new_two_H_vertices':a.EXPECTED_2H[m],
                'wick_eligible':2 if m in a.DIRECT else 0,
                'observed_in_saved_C5':2 if m in a.DIRECT else 0}
        rows3.append({'model':m,'quadratic_kind':kind,'proven_nonzero_C5_contribution':False})
        rows4[m]={'extra_quadratic_exists':kind!='none'}
        direct=['q1','q2'] if m in a.DIRECT else []
        uv={'direct_C5_quartics':direct,'direct_C5_quartic_orders':{q:1 for q in direct},
            'immediately_generated_quartics':['q1','q2','q3'] if m in a.DIRECT else [],
            'quadratic_first_generation_order':None}
        eft={'immediately_generated_quartics':[], 'quadratic_first_generation_order':None}
        rows5[m]={'extra_quartics':n,'extra_quadratic':kind!='none',
                  'stages':{'uv':uv,'eft1':eft}}
    return {
        'batch1':{'summary':{'strict_complete':True,'models_checked':16,
                            'extra_quartic_slots':50,'extra_quadratic_slots':8,
                            'observed_additional_quartic_C5_slots':4},'operators':ops},
        'batch2':{'status':'WickNecessaryConditionsReproduceSavedC5Inventory',
                  'source_provenance_verified_on_local_checkout':True,
                  'new_matching_calculation_performed':False,
                  'counts':{'models_checked':16,'additional_quartics_total':50,'new_two_H_vertices':24,
                            'new_direct_wick_eligible':4,'saved_C5_observed_new':4,
                            'excluded_H_Hdag_by_higgs_orientation':10,
                            'excluded_Hdag2_by_scalar_orientation':10,'new_no_H_vertex':26,
                            'extra_quadratics_not_screened':8},'per_model':per},
        'batch3':{'status':'CandidateQuadraticMassSpectraValidatedNotC5','models_checked':16,
                  'additional_quadratic_slots':8,'self_bilinears':6,'mixed_bilinears':2,
                  'no_tachyons_in_reported_candidate_spectra':True,
                  'general_neutrino_C5_matching_completed':False,
                  'physical_matchete_rgbeta_mapping_certified':False,'models':rows3},
        'batch4':{'model_count':16,'quadratic_slots':8,
                  'uv_and_eft1_beta_exports_checked':True,
                  'source_vanishes_if_all_additional_quartics_zero_at_one_loop':True,
                  'not_a_zero_mu2_trajectory':True,'not_a_neutrino_mass_error_estimate':True,
                  'models':rows4},
        'batch5':{'status':'CandidateReducedScalarOneLoopRGGenerationOrderAudit',
                  'models_checked':16,'quartic_slots_checked':50,'quadratic_slots_checked':8,
                  'uv_immediately_sourced_quartic_slots':6,
                  'uv_quadratic_first_generated_at_order_2_count':0,
                  'model_results':rows5},
    }


def test_synthesis_counts_and_physics_flags():
    out=a.build(fixture_reports(),{'test':'hash'})
    assert out['status']=='ConditionalLeadingOneLoopT3ScalarSectorAudit'
    assert out['quartic_result']['not_direct_at_selected_order']==46
    assert out['quartic_result']['direct_extra_first_generated_in_candidate_UV_running']==4
    assert out['quadratic_result']['allowed_complex_bilinear_slots']==8
    assert out['conditions']['general_C5_neutrino_mass_error_computed'] is False


def test_bc_direct_and_mass_order():
    x={r['model']:r for r in a.build(fixture_reports(),{})['model_results']}
    assert x['B_-1']['extra_direct_C5_quartics_if_nonzero']==2
    assert x['C_-1']['uv_quartics_sourced_at_order_one']==3
    assert x['A_-2']['extra_bilinear_kind']=='self_s2'
    assert x['B_+1']['extra_bilinear_kind']=='none'
    assert x['A_-4']['extra_quartic_slots']==0
    assert all(r['candidate_uv_bilinear_rg_order'] is None for r in x.values())


def test_summary_prose_does_not_claim_physical_error():
    doc=a.render_markdown(a.build(fixture_reports(),{}))
    assert 'not a numerical neutrino-mass error bound' in doc
    assert 'zero at the specified UV boundary' in doc
    assert 'mass-basis' not in doc.split('## Verified operator classification')[0]


def test_missing_source_rejected(tmp_path):
    with pytest.raises(a.FinalAuditError,match='missing report'):
        a.load_reports(tmp_path)


def test_fingerprints_and_io(tmp_path):
    f=fixture_reports()
    for name,relative in a.INPUTS.items():
        p=tmp_path/relative
        p.parent.mkdir(parents=True,exist_ok=True)
        p.write_text(json.dumps(f[name]),encoding='utf-8')
    reports,hashes=a.load_reports(tmp_path)
    assert len(hashes)==5
    assert set(hashes) == {relative.as_posix() for relative in a.INPUTS.values()}
    assert all(len(digest) == 64 and all(ch in '0123456789abcdef' for ch in digest)
               for digest in hashes.values())
    final=a.build(reports,hashes)
    js,md=a.write(final,tmp_path/'out')
    assert json.loads(js.read_text())['status']==final['status']
    assert 'Batch 6' in md.read_text()


def test_matching_2H_mutation_rejected():
    f=fixture_reports()
    f['batch2']['per_model']['B_-1']['wick_eligible']=1
    with pytest.raises(a.FinalAuditError,match='Wick'):
        a.build(f,{})


def test_extra_mass_inventory_mutation_rejected():
    f=fixture_reports()
    f['batch3']['models'][0]['quadratic_kind']='mixed'
    with pytest.raises(a.FinalAuditError,match='quadratic kind'):
        a.build(f,{})


def test_nonzero_C5_mass_claim_rejected():
    f=fixture_reports()
    f['batch3']['models'][0]['proven_nonzero_C5_contribution']=True
    with pytest.raises(a.FinalAuditError,match='falsely claims C5'):
        a.build(f,{})


def test_mass_rge_generation_rejected():
    f=fixture_reports()
    f['batch5']['model_results']['B_-1']['stages']['uv']['quadratic_first_generation_order']=2
    with pytest.raises(a.FinalAuditError,match='quadratic can be generated'):
        a.build(f,{})


def test_missing_direct_quartic_generation_rejected():
    f=fixture_reports()
    f['batch5']['model_results']['C_-1']['stages']['uv']['direct_C5_quartic_orders']['q1']=None
    with pytest.raises(a.FinalAuditError,match='direct quartic generation'):
        a.build(f,{})


def test_unverified_source_provenance_rejected():
    f=fixture_reports()
    f['batch2']['source_provenance_verified_on_local_checkout']=False
    with pytest.raises(a.FinalAuditError,match='Batch 2'):
        a.build(f,{})


def test_physical_mapping_claim_rejected():
    f=fixture_reports()
    f['batch3']['physical_matchete_rgbeta_mapping_certified']=True
    with pytest.raises(a.FinalAuditError,match='Batch 3'):
        a.build(f,{})


def test_inconsistent_UV_quartic_sources_rejected():
    f=fixture_reports()
    f['batch5']['model_results']['B_-1']['stages']['uv']['immediately_generated_quartics'].pop()
    with pytest.raises(a.FinalAuditError,match='source quartic generations'):
        a.build(f,{})


def test_missing_model_rejected():
    f=fixture_reports()
    del f['batch4']['models']['A_-4']
    with pytest.raises(a.FinalAuditError,match='same 16 models'):
        a.build(f,{})
