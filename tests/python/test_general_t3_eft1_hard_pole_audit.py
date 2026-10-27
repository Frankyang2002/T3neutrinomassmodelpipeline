"""Offline lexical/operator factorisation regression tests; not EFT1 pole-RGE matching."""
from __future__ import annotations

from pathlib import Path
import json

import pytest

from scripts import RunGeneralT3EFT1HardPoleAudit as task

MF = 'Coupling[MF, {}, 0]'
LOG = task.MATCH_LOG


def basic_term(chirality: int, *, quartic: str = 'lambdaT3') -> str:
    return (
        f'(hbar*Coupling[{quartic}, {{}}, 0]'
        '*Field[H, Scalar, {Index[a, SU2L]}, {}]'
        '*Field[H, Scalar, {Index[b, SU2L]}, {}]'
        '*Field[l, Fermion, {Index[a, Flavor]}, {}]'
        '*Field[l, Fermion, {Index[b, Flavor]}, {}]'
        f'*NCM[DiracProduct[GammaCC, Proj[{chirality}]]])/{MF}'
    )


def fixture_terms() -> list[tuple[int, str]]:
    all_terms = []
    for sign, chirality, quartic in (
        (1, -1, 'lambdaT3'),
        (-1, 1, 'lambdaT3ExtraQ020020I1'),
    ):
        finite = basic_term(chirality, quartic=quartic)
        log = finite.replace(f'])/{MF}', f']*{LOG})/{MF}')
        pole = finite.replace(f'/{MF}', f'/(ϵ*{MF})')
        all_terms.extend(((sign, finite), (sign, log), (sign, pole)))
    return all_terms


def sample() -> str:
    return task.render_terms(fixture_terms())


def report(model: str = 'B_-1') -> dict:
    stem = task.model_stem(model)
    return {
        'format': task.INPUT_FORMAT,
        'status': task.INPUT_STATUS,
        'model': model,
        'matching_approach': 'FermionFirstTreeAndOneLoopScalarsRemainLight',
        'matching_EFT_order': 5, 'matching_loop_order': 1,
        'source_quartic_count': 17, 'source_yukawa_count': 2,
        'fermion_heavy_mass_basis_detailed_check_passed': True,
        'tree_matching_completed': True, 'one_loop_matching_completed': True,
        'raw_hard_llhh_no_fermion_field': True,
        'loop_difference_fermion_absent': True,
        'raw_hard_llhh_scheme_renormalized': False,
        'raw_hard_llhh_pole_rge_certified': False,
        'eft1_rge_performed': False,
        'scalar_threshold_matching_performed': False,
        'one_loop_hard_fermion_matching_performed': False,
        'physical_C5_computed': False,
        'sequential_EFT1_completed': False,
        'physical_matching_performed': False,
        'raw_hard_llhh_candidate_term_count': 204,
        'loop_difference_term_count': 6,
        'raw_hard_llhh_file': stem + '_raw_hard_llhh_sector.inputform.txt',
        'one_loop_difference_file': stem + '_eft1_one_loop_correction.inputform.txt',
        'llss_inputform_file': stem + '_llss_tree.inputform.txt',
    }


def write_sources(folder: Path, model: str = 'B_-1') -> None:
    folder.mkdir(parents=True, exist_ok=True)
    stem = task.model_stem(model)
    (folder / (stem + '.json')).write_text(json.dumps(report(model)))
    (folder / (stem + '_raw_hard_llhh_sector.inputform.txt')).write_text(sample())
    (folder / (stem + '_eft1_one_loop_correction.inputform.txt')).write_text(sample())
    (folder / (stem + '_llss_tree.inputform.txt')).write_text(
        'GammaCC Field[NewScalarX] Field[NewScalarY]')


@pytest.mark.parametrize('model',task.MODELS)
def test_model_stem(model):
    assert task.model_stem(model).startswith(model.replace('-', 'm'))


@pytest.mark.parametrize('model,which',[('A_+0','final'),('B_-1','other')])
def test_unsupported_stems(model,which):
    with pytest.raises(ValueError,match='Unsupported'):
        task.model_stem(model,which)


def test_valid_triplets_and_chirality_split():
    result = task.factorize_hard_sector(sample(), expected_count=6)
    assert result['term_count'] == 6
    assert result['independent_signed_operator_structures'] == 2
    assert result['sector_counts'] == {'finite': 2,'log': 2,'pole': 2}
    assert result['chirality_counts']['holomorphic_PL.finite'] == 1
    assert result['chirality_counts']['conjugate_PR.log'] == 1
    assert 'Proj[-1]' in result['holomorphic_finite_operator']
    assert 'Proj[1]' not in result['holomorphic_finite_operator']
    assert 'ϵ' not in result['coefficient_sector_strings']['pole']
    assert 'Log[' not in result['coefficient_sector_strings']['log']
    assert 'GammaCC' in result['coefficient_sector_strings']['finite']
    assert 'CG[' not in result['coefficient_sector_strings']['finite']  # synthetic only


def test_exact_source_roundtrip():
    raw=sample()
    assert task.render_terms(task._balanced_terms(raw))==raw
    assert task._balanced_terms(raw)==fixture_terms()


def test_log_and_pole_coefficients_are_same_multiset():
    out=task.factorize_hard_sector(sample(),expected_count=6)
    from collections import Counter
    for key in ('pole','log'):
        assert Counter(task._balanced_terms(out['coefficient_sector_strings'][key])) == Counter(
            task._balanced_terms(out['coefficient_sector_strings']['finite'])
        )


@pytest.mark.parametrize('replacement', [
    'Log[other]', 'Log[μbar2/MSX^2]', 'Log[μbar2/MF^2]^2',
])
def test_unknown_log_fails(replacement):
    corrupted=sample().replace(LOG, replacement)
    with pytest.raises(ValueError):
        task.factorize_hard_sector(corrupted,expected_count=6)


@pytest.mark.parametrize('replacement', ['ε*','epsilon*','ϵ^2*','ϵ*ϵ*'])
def test_unknown_pole_symbol_or_power_fails(replacement):
    broken=sample().replace('ϵ*',replacement)
    with pytest.raises(ValueError):
        task.factorize_hard_sector(broken, expected_count=6)


def test_mismatched_pole_coupling_fails():
    raw=sample().replace('ϵ*Coupling[MF, {}, 0]',
                          'ϵ*Coupling[MF2, {}, 0]', 1)
    with pytest.raises(ValueError):
        task.factorize_hard_sector(raw, expected_count=6)


def test_missing_term_or_sign_mismatch_fails():
    parts=fixture_terms()
    with pytest.raises(ValueError,match='Expected 6'):
        task.factorize_hard_sector(task.render_terms(parts[:-1]), expected_count=6)
    index=next(i for i,(_,t) in enumerate(parts) if 'ϵ*' in t)
    parts[index]=(-parts[index][0],parts[index][1])
    with pytest.raises(ValueError):
        task.factorize_hard_sector(task.render_terms(parts),expected_count=6)


@pytest.mark.parametrize('targ,new', [
    ('GammaCC','Gamma'),('hbar','hbar^2'),
    ('Field[H, Scalar','Field[NewScalarX, Scalar'),
    ('Proj[-1]','Proj[7]'),
])
def test_wrong_operator_is_rejected(targ,new):
    with pytest.raises(ValueError):
        task.factorize_hard_sector(sample().replace(targ,new),expected_count=6)


def test_corrupted_top_level_parser_fails():
    with pytest.raises(ValueError):
        task._balanced_terms(sample() + ']')


@pytest.mark.parametrize('key,value', [
    ('status','Success'), ('model','C_-1'),
    ('raw_hard_llhh_candidate_term_count', 0),
    ('raw_hard_llhh_pole_rge_certified',True),
    ('physical_C5_computed',True),
    ('scalar_threshold_matching_performed',True),
])
def test_certificate_not_trusted_if_inconsistent(tmp_path,key,value):
    write_sources(tmp_path)
    cert=report();cert[key]=value
    (tmp_path/(task.model_stem('B_-1')+'.json')).write_text(json.dumps(cert))
    with pytest.raises(ValueError):
        task.audit_model('B_-1',input_dir=tmp_path,output_dir=tmp_path/'out')


def test_raw_sector_must_be_exact_subset_of_full_difference(tmp_path, monkeypatch):
    write_sources(tmp_path)
    rawpath=tmp_path/(task.model_stem('B_-1')+'_raw_hard_llhh_sector.inputform.txt')
    rawpath.write_text(sample().replace('lambdaT3','lambdaT3ExtraQ020200I1'))
    # Valid quartic in hard sector, but not present in one-loop difference.
    monkeypatch.setattr(task,'factorize_hard_sector',lambda *_a,**_k:{} )
    with pytest.raises(ValueError, match='not an exact signed subset'):
        task.audit_model('B_-1',input_dir=tmp_path,output_dir=tmp_path/'out')


def test_certificate_source_read_and_validate(tmp_path):
    write_sources(tmp_path)
    task._read_report(tmp_path/(task.model_stem('B_-1')+'.json'),
                      'B_-1',task.model_stem('B_-1'))


def test_residue_normalization_for_B_and_C():
    mf=MF
    assert task._operator_coefficient_term(f'(X)/(ϵ*{mf})','pole') == f'(X)/{mf}'
    assert task._operator_coefficient_term(f'(X)/(Sqrt[3]*ϵ*{mf})','pole') == f'(X)/(Sqrt[3]*{mf})'


def test_code_does_not_claim_physical_c5():
    s=(Path(__file__).resolve().parents[2]/'scripts/RunGeneralT3EFT1HardPoleAudit.py').read_text()
    assert "'physical_majorana_c5_computed': False" in s
    assert "'independent_eft1_llss_to_c5_rge_check_performed': False" in s
    assert "'flavour_su2_c5_operator_projection_performed': False" in s
