"""Offline fixture and analytic checks, not native Mathematica/RGBeta validation."""
import copy
import json
import math
from pathlib import Path

import pytest

from scripts.AuditT3QuadraticImpact import (
    MODELS, REPS, AuditError, analytic_spectrum, audit_model, kind_from_manifest,
    run, stage_diagnostics,
)


def _kind(model):
    if model in {'A_+0','D_+0','E_+0'}:
        return 'self_s1'
    if model in {'A_-2','D_-2','E_-2'}:
        return 'self_s2'
    if model in {'B_-1','C_-1'}:
        return 'mixed'
    return 'none'


def _record_manifest(kind):
    counts = {'self_s1':[0,0,0,2,0,0], 'self_s2':[0,0,0,0,0,2],
              'mixed':[0,0,0,1,0,1]}
    structures=[] if kind=='none' else [{'counts':counts[kind], 'invariants':1}]
    return {'quadratic_complex_couplings_not_enabled':len(structures),
            'quadratic_structures_not_enabled':structures}


def _spect(model,kind,m1,m2,mu):
    dims=REPS[model[0]]
    eig=analytic_spectrum(kind,m1,m2,mu,*dims)
    return {'model':model,'quadratic_kind':kind,
            'analytic_spectrum_verified':True,
            'included_in_threshold_matching':False,
            'matchete_mass_basis_verified':False,
            'rgbeta_mass_normalisation_verified':False,
            'extra_mass_coefficient_gev2':{'real':mu.real,'imag':mu.imag},
            'dimensions':list(dims),
            'mass_squared_eigenvalues_gev2':eig,
            'analytic_mass_squared_eigenvalues_gev2':eig,
            'tachyonic':min(eig)<0,
            'minimum_mass_squared_gev2':min(eig)}


def _record(model, frac=0.01):
    kind=_kind(model)
    frac=frac if kind!='none' else 0.0
    m1,m2=1.3e6,1e6
    a=frac*min(m1,m2)*complex(math.cos(0.3), math.sin(0.3))
    b=frac*1.2e6*complex(math.cos(0.3), math.sin(0.3))
    symbols={'self_s1':'muT3GeneralS1Sq','self_s2':'muT3GeneralS2Sq',
             'mixed':'muT3GeneralS1S2','none':None}
    record={'format':'general-t3-candidate-quadratic-eft1-v1',
            'status':'CandidateQuadraticAndFullQuarticEFT1Integrated',
            'model':model, 'quadratic_kind':kind,
            'mu2_rgbeta_symbol':symbols[kind],
            'mu2_beta_included':kind!='none',
            'complete_one_loop_quartic_class_projected':True,
            'physical_matchete_rgbeta_mapping_certified':False,
            'scalar_mass_threshold_matching_performed':False,
            'wilson_coefficient_matching_performed':False,
            'physical_neutrino_mass_prediction':False,
            'quadratic_operators_evolved':int(kind!='none'),
            'mu2_fraction_initial':frac}
    for s,masses,mu in [('initial',(m1,m2),a),('final',(1.5e6,1.2e6),b)]:
        record[s+'_baseline']={'mS1Sq':masses[0],'mS2Sq':masses[1]}
        record['mu2_'+s+'_gev2']={'real':mu.real,'imag':mu.imag}
        record[s+'_scalar_mass_spectrum']=_spect(model,kind,*masses,mu)
    return record


def _make_files(root: Path):
    manifest={'model_count':16,'models':{m:_record_manifest(_kind(m)) for m in MODELS}}
    batch2={'status':'WickNecessaryConditionsReproduceSavedC5Inventory',
            'counts':{'models_checked':16,'additional_quartics_total':50,
                      'new_direct_wick_eligible':4,'extra_quadratics_not_screened':8}}
    data_dir=root/'output/general_t3_quadratic_eft1'
    data_dir.mkdir(parents=True)
    (root/'docs').mkdir()
    (root/'docs/general_scalar_trial_operator_manifest.json').write_text(json.dumps(manifest))
    (root/'output/scalar_relevance').mkdir(parents=True)
    (root/'output/scalar_relevance/leading_c5_vertex_audit.json').write_text(json.dumps(batch2))
    summary={'status':'Candidate16ModelQuadraticEFT1Integrated',
             'models_checked':16,'new_quartics_evolved':38,
             'reused_quartics_evolved':12,'quadratic_model_slots_evolved':8,
             'physical_mapping_certified':False,'matching_extended':False,
             'mu2_fraction':0.01,'models':{}}
    for m in MODELS:
        r=_record(m)
        name=m.replace('+','p').replace('-','m')
        (data_dir/f'{name}_eft1_quadratic_running.json').write_text(json.dumps(r))
        summary['models'][m]={'quadratic_kind':_kind(m),'tachyonic':False}
    (data_dir/'all_models_quadratic_eft1_summary.json').write_text(json.dumps(summary))
    return manifest


def test_spectrum_self_split_and_multiplicity():
    eigen=analytic_spectrum('self_s1',100.,400.,3.+4.j,1,3)
    assert eigen==[90.,110.]+[400.]*6
    eigen=analytic_spectrum('self_s2',100.,400.,3.+4.j,3,1)
    assert eigen==[100.]*6+[390.,410.]


def test_mixed_analytic_spectrum_and_mixing_zero_limit():
    values=analytic_spectrum('mixed',120.,100.,3+4j,2,2)
    lower,upper=110-math.sqrt(125),110+math.sqrt(125)
    assert values==[lower]*4+[upper]*4
    assert analytic_spectrum('mixed',120.,100.,0j,2,2)==[100.]*4+[120.]*4


def test_nonzero_forbidden_coefficient_rejected():
    with pytest.raises(AuditError):
        analytic_spectrum('none',100.,200.,1.+0j,1,3)


def test_manifest_identifies_eight_bilinears():
    kinds=[kind_from_manifest(_record_manifest(_kind(m)),m) for m in MODELS]
    assert kinds.count('none')==8
    assert kinds.count('self_s1')==3
    assert kinds.count('self_s2')==3
    assert kinds.count('mixed')==2


def test_candidate_mixed_effect_at_one_loop_not_declared_proven():
    model='B_-1'; row=audit_model(model,_record(model),{'models':{model:_record_manifest('mixed')}},0.01)
    assert row['potential_leading_one_loop_C5_effect_if_nonzero']
    assert not row['proven_nonzero_C5_contribution']
    assert row['final']['sin_2theta_mixed_doublets']>0
    assert row['final']['max_relative_mass_shift_at_fixed_diagonals']>0


def test_zero_allowed_quadratic_is_a_boundary_choice_not_symmetry():
    model='A_+0'
    r=_record(model,frac=0.0)
    row=audit_model(model,r,{'models':{model:_record_manifest('self_s1')}},0.0)
    assert row['direct_effect_on_scalar_propagators_at_tree_level']
    assert row['initial']['max_relative_mass_squared_shift_at_fixed_diagonals']==0.0


def test_none_uses_zero_fraction_even_if_other_models_nonzero():
    model='A_-4'
    row=audit_model(model,_record(model),{'models':{model:_record_manifest('none')}},0.01)
    assert row['quadratic_kind']=='none'
    assert row['final']['mu2_abs_gev2']==0.0
    assert not row['potential_leading_one_loop_C5_effect_if_nonzero']


def test_mass_mismatch_fails_closed():
    model='B_-1';r=_record(model)
    r['final_scalar_mass_spectrum']['mass_squared_eigenvalues_gev2'][0]+=300
    with pytest.raises(AuditError,match='mass-spectrum mismatch'):
        audit_model(model,r,{'models':{model:_record_manifest('mixed')}},0.01)


def test_bad_matchete_flag_fails_closed():
    model='B_-1';r=_record(model)
    r['final_scalar_mass_spectrum']['matchete_mass_basis_verified']=True
    with pytest.raises(AuditError,match='provenance'):
        audit_model(model,r,{'models':{model:_record_manifest('mixed')}},0.01)


def test_bad_initial_fraction_fails_closed():
    model='B_-1';r=_record(model)
    r['mu2_initial_gev2']={'real':0,'imag':0}
    with pytest.raises(AuditError):
        audit_model(model,r,{'models':{model:_record_manifest('mixed')}},0.01)


def test_full_offline_fixture_writes_three_reports(tmp_path):
    _make_files(tmp_path)
    result=run(tmp_path)
    assert result['models_checked']==16
    assert result['additional_quadratic_slots']==8
    assert not result['general_neutrino_C5_matching_completed']
    assert result['no_tachyons_in_reported_candidate_spectra']
    for name in ('quadratic_impact_audit.json','quadratic_impact_audit.csv',
                 'quadratic_impact_summary.md'):
        assert (tmp_path/'output/scalar_relevance'/name).is_file()


def test_missing_model_fails_closed(tmp_path):
    _make_files(tmp_path)
    (tmp_path/'output/general_t3_quadratic_eft1/B_m1_eft1_quadratic_running.json').unlink()
    with pytest.raises(AuditError,match='Missing input'):
        run(tmp_path)


def test_bad_batch2_counts_fails_closed(tmp_path):
    _make_files(tmp_path)
    path=tmp_path/'output/scalar_relevance/leading_c5_vertex_audit.json'
    x=json.loads(path.read_text());x['counts']['new_direct_wick_eligible']=5
    path.write_text(json.dumps(x))
    with pytest.raises(AuditError,match='Batch 2'):
        run(tmp_path)


def test_zero_mode_can_be_tachyon_free_but_not_neutrino_bound():
    m='B_-1'
    data=_record(m)
    row=audit_model(m,data,{'models':{m:_record_manifest('mixed')}},0.01)
    assert row['final']['all_scalar_mass_squared_positive']
    assert row['final']['max_relative_mass_shift_at_fixed_diagonals']<0.01
    assert row['final']['matching_neutrino_mass_shift_evaluated'] is False
