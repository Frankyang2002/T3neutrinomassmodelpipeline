"""Offline Batch8 checks: source contracts, replay, C5 norm, scope/markdown."""
import hashlib
import json
import math
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts import RunT3ZeroBoundaryLowEnergyC5 as a


def z(x):
    return {'real': x, 'imag': 0.0}


def fixture_snap():
    return {'mu_gev': 1e7, 'mS1Sq_gev2': 1e6, 'mS2Sq_gev2': 1.2e6,
            'lambdaT3': z(.1), 'quadratic_mu2_gev2': z(0),
            'new_quartics': {'n': z(0)},
            'reused_quartics': {f'r{i}': z(0) for i in range(6)}}


def fixture_saved(root, model='B_-1'):
    source = ['RGE/general/generated/general_rgbeta_quartic_registry.json',
              'configs/source.json', 'beta/UV.json', 'beta/EFT1.json',
              'stage31/report.json']
    for n in source:
        p = root / n
        p.parent.mkdir(parents=True, exist_ok=True)
        if n.endswith('registry.json'):
            p.write_text(json.dumps({'models': {model: {'new_rgbeta_symbols': [{'name':'n'}],
                'preexisting_rgbeta_symbols':[{'existing_name':f'r{i}'} for i in range(6)]}}}))
        else:
            p.write_text(n)
    hashes = {n: hashlib.sha256((root/n).read_bytes()).hexdigest() for n in source}
    snap = fixture_snap()
    return {'format':'t3-zero-extra-boundary-bc-renormalisable-running-v1',
            'status':'CandidateBCZeroBoundaryTwoSegmentRunning','model':model,
            'boundary_scale_gev':1e7,'fermion_threshold_gev':1e5,
            'scalar_threshold_gev':1e3,
            'boundary_all_additional_couplings_zero':True,
            'uv_extra_quartic_slots':7,'uv_quadratic_slots':1,
            'physical_matchete_rgbeta_map_certified':False,
            'full_sequential_C5_computed':False,
            'physical_neutrino_mass_error_bound':False,
            'source_sha256': hashes,
            **{k:snap for k in ['uv_initial','uv_final','eft1_initial','eft1_final']}}


def test_model_tags():
    assert a.safe_tag('B_-1') == 'B_m1'
    assert a.safe_tag('C_-1') == 'C_m1'
    with pytest.raises(ValueError): a.safe_tag('A_+0')


def test_saved_valid_and_hashes(tmp_path):
    saved = fixture_saved(tmp_path)
    a.validate_saved_batch7(saved, 'B_-1', root=tmp_path)


def test_saved_rejects_mutated_hash(tmp_path):
    saved = fixture_saved(tmp_path)
    (tmp_path/'beta/UV.json').write_text('changed')
    with pytest.raises(ValueError, match='changed'):
        a.validate_saved_batch7(saved, 'B_-1', root=tmp_path)


def test_saved_rejects_wrong_model(tmp_path):
    saved = fixture_saved(tmp_path)
    with pytest.raises(ValueError): a.validate_saved_batch7(saved, 'C_-1', root=tmp_path)


def test_saved_rejects_nonzero_boundary(tmp_path):
    saved = fixture_saved(tmp_path)
    saved['uv_initial']['reused_quartics']['r4']['real'] = 0.01
    with pytest.raises(ValueError, match='not exactly zero'):
        a.validate_saved_batch7(saved, 'B_-1', root=tmp_path)


def test_saved_rejects_false_certification(tmp_path):
    saved = fixture_saved(tmp_path)
    saved['physical_matchete_rgbeta_map_certified'] = True
    with pytest.raises(ValueError): a.validate_saved_batch7(saved, 'B_-1', root=tmp_path)


def test_saved_rejects_unsafe_path(tmp_path):
    saved = fixture_saved(tmp_path)
    saved['source_sha256']['../escape'] = '0'*64
    with pytest.raises(ValueError, match='unsafe'):
        a.validate_saved_batch7(saved, 'B_-1', root=tmp_path)


def test_replay_accepts_exact_and_perturbation():
    snap = fixture_snap()
    saved = {name:snap for name in ('uv_initial','uv_final','eft1_initial','eft1_final')}
    actual = json.loads(json.dumps(saved))
    actual['eft1_final']['lambdaT3']['real'] += 1e-11
    a.validate_replay(saved, actual)


def test_replay_rejects_wrong_endpoint():
    snap = fixture_snap()
    saved = {name:snap for name in ('uv_initial','uv_final','eft1_initial','eft1_final')}
    actual = json.loads(json.dumps(saved))
    actual['eft1_final']['new_quartics']['n']['real'] = 0.005
    with pytest.raises(ValueError, match='replay differs'):
        a.validate_replay(saved, actual)


def test_fractional_c5():
    ref = np.eye(3,dtype=complex) * 1e-16
    assert a.fractional_matrix_change(ref, 1.1*ref) == pytest.approx(.1)
    assert a.fractional_matrix_change(ref, ref) == 0


def test_fractional_c5_zero_rejected():
    with pytest.raises(ValueError, match='zero reference'):
        a.fractional_matrix_change(np.zeros((3,3)), np.eye(3))


def test_complex_matrix_check():
    zmat = {'real':np.eye(3).tolist(), 'imag':np.zeros((3,3)).tolist()}
    assert np.array_equal(a._matrix_from_json(zmat,'x'),np.eye(3))
    zmat['real'][0][1] = .1
    with pytest.raises(ValueError, match='nonsymmetric'):
        a._matrix_from_json(zmat,'x')


def test_reject_wrong_prior(monkeypatch,tmp_path):
    import types
    fake = types.ModuleType('scripts.Stage36FinalWeinbergRunning')
    fake.FORMAT = 'stage36-t3-existing-c5-sm-weinberg-v1'
    fake.STATUS = 'SuccessExistingC5WeinbergRunning'
    monkeypatch.setitem(sys.modules, 'scripts.Stage36FinalWeinbergRunning',fake)
    artifact = tmp_path/'source.json'
    artifact.write_text('matched')
    prior = {'format':fake.FORMAT,'status':fake.STATUS,'model':'B_-1',
             'source_existing_c5_used':True,'matching_source_is_hierarchical':False,
             'new_general_quartics_included_in_matching':False,
             'sm_weinberg_running_performed':True,
             'physical_coupling_mapping_certified':False,
             'mu_scalar_threshold_gev':1000,'mu_low_gev':100,
             'source_matching_artifact_sha256':a.b7.digest(artifact),
             'source_matching_artifact':'source.json'}
    with pytest.raises(ValueError, match='unsupported provenance'):
        a.verify_prior_stage36(prior,'B_-1',artifact,root=tmp_path)
    prior['matching_source_is_hierarchical'] = True
    a.verify_prior_stage36(prior,'B_-1',artifact,root=tmp_path)
    prior['source_matching_artifact_sha256'] = '0'*64
    with pytest.raises(ValueError, match='changed'):
        a.verify_prior_stage36(prior,'B_-1',artifact,root=tmp_path)


def test_summary_preserves_caveat():
    row = {'model':'B_-1','neutrino_observables':{
            'masses_ev':[.001,.008,.05],'delta_m21_sq_ev2':6.3e-5,
            'delta_m31_sq_ev2':.0025},
           'comparison_to_saved_stage36_nonzero_stage29_boundary':{
             'fractional_C5_fro_change_at_100gev':.03}}
    summary = a.make_summary_md([row])
    assert '0.05' in summary
    assert 'distinct UV boundary scenarios' in summary
    assert 'not a physical' in summary


def test_frozen_root_not_required_for_unit_tests():
    assert a.OUT.relative_to(a.ROOT).as_posix() == 'output/scalar_relevance/zero_boundary_low_energy'
    assert a.MODELS == ('B_-1','C_-1')
