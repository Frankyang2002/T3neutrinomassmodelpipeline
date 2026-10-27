"""Regression checks: exact raw two-threshold pole cancellation, not C5."""
from __future__ import annotations

from collections import Counter
from pathlib import Path
import json

import pytest

from scripts import RunGeneralT3TwoThresholdPoleAssembly as m


def test_balanced_split_exact():
    raw = 'f[x + y]*a + (f[z] - 3)*b - (f[t] + (q - v))'
    values = m.split_terms(raw)
    assert values == [(1, 'f[x + y]*a'), (1, '(f[z] - 3)*b'),
                      (-1, '(f[t] + (q - v))')]
    assert m.render_terms(values) == raw


@pytest.mark.parametrize('raw', [
    '', '0', 'x + ', '(x + y', 'x + ]y', 'x - ',
])
def test_reject_malformed_source(raw):
    with pytest.raises(ValueError, match='fail-closed'):
        m.split_terms(raw)


def test_unary_minus_extraction_preserves_sign():
    assert m.strip_outer_unary_minus(1, '-((x+y)/(ϵ*M))') == (-1, '(x+y)/(ϵ*M)')
    assert m.strip_outer_unary_minus(-1, '-((x+y)/(ϵ*M))') == (1, '(x+y)/(ϵ*M)')
    assert m.strip_outer_unary_minus(1, '-x + y') == (1, '-x + y')


def test_exact_operator_subset_normalizes_leading_unary_minus():
    subset = m.normalized_subset([(1, '-((hbar*x)/(ϵ*MF))')])
    source = m.normalized_subset([(-1, '(hbar*x)/(ϵ*MF)'), (1, 'background')])
    assert not subset - source


@pytest.mark.parametrize('model', m.MODELS)
def test_model_stem(model):
    assert m.stem(model, 'final').endswith('_final')


@pytest.mark.parametrize('model', ['A_0','invalid','C_+1'])
def test_other_models_rejected(model):
    with pytest.raises(ValueError, match='fail-closed'):
        m.stem(model, 'final')


def test_source_kind_uses_only_known_logs():
    assert m.kind('a*' + m.LOG_X, hard=False) == 'log'
    assert m.kind('a*' + m.LOG_Y, hard=False) == 'log'
    assert m.kind('a*' + m.LOG_F, hard=True) == 'log'
    assert m.kind('a/(ϵ*MF)', hard=True) == 'pole'
    assert m.kind('a', hard=False) == 'finite'
    with pytest.raises(ValueError, match='fail-closed'):
        m.kind('a*' + m.LOG_X, hard=True)
    with pytest.raises(ValueError, match='fail-closed'):
        m.kind('a/(ϵ*MF)*' + m.LOG_X, hard=False)
    with pytest.raises(ValueError, match='fail-closed'):
        m.kind('a*' + m.LOG_X + '*' + m.LOG_Y, hard=False)


def _epsilon(ia, ib, barred=False):
    name = 'Bar[eps[SU2L]]' if barred else 'eps[SU2L]'
    def idx(i):
        s = f'Index[d$${i}, SU2L[fund]]'
        return f'Bar[{s}]' if barred else s
    return f'CG[{name}, {{{idx(ia)}, {idx(ib)}}}]'


def _closed():
    return (_epsilon(1, 3) + '*' + _epsilon(2, 4) + '*'
            'Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*'
            'Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*'
            'NCM[Field[l, Fermion, {Index[d$$3, SU2L[fund]]}, {}], '
            'Field[l, Fermion, {Index[d$$4, SU2L[fund]]}, {}]]')


def test_canonicalize_epsilon_swaps_sign_once():
    src = _closed()
    assert m.canonicalize_su2_epsilon(1, src) == (1, src)
    swapped = src.replace(_epsilon(1, 3), _epsilon(3, 1))
    sign, normalized = m.canonicalize_su2_epsilon(1, swapped)
    assert sign == -1 and normalized == src


def test_epsilon_twice_restores_sign():
    src = _closed()
    swapped = src.replace(_epsilon(1, 3), _epsilon(3, 1)).replace(
        _epsilon(2, 4), _epsilon(4, 2))
    assert m.canonicalize_su2_epsilon(1, swapped) == (1, src)


def test_closed_dummy_relabel_preserves_double_contractions():
    original = _closed()
    permuted = m.relabel_su2_dummies(original)
    assert permuted != original
    assert m.relabel_su2_dummies(permuted) == original
    assert Counter(m.SU2_INDEX.findall(permuted)) == Counter(
        {'1':2,'2':2,'3':2,'4':2})


def test_fail_on_open_su2_leg():
    with pytest.raises(ValueError, match='fail-closed'):
        m.canonicalize_su2_epsilon(1, _epsilon(1, 3) + '*' + _epsilon(2, 4))


def test_fail_on_missing_epsilon():
    with pytest.raises(ValueError, match='fail-closed'):
        m.canonicalize_su2_epsilon(1, _closed().replace('CG[eps[SU2L]', 'CG[other]'))


def test_mu_f_substitution_and_hard_log_zero():
    assert m._at_mu_f(1, 'foo*'+m.LOG_F) is None
    result = m._at_mu_f(-1, 'a*'+m.LOG_X+'*'+m.LOG_Y)
    assert result is not None and result[0] == -1
    assert 'Coupling[MF, {}, 0]^2/Coupling[MSX, {}, 0]^2' in result[1]
    assert 'Coupling[MF, {}, 0]^2/Coupling[MSY, {}, 0]^2' in result[1]
    assert 'μbar2' not in result[1]


def test_source_certificates_and_real_pole_matching(tmp_path):
    # Runs automatically on a real project with locally saved Matchete data;
    # remains a unit-test suite on a clean CI checkout.
    if not (m.HARD_DIR / 'B_m1_final_eft1_fermion_first_oneloop.json').exists():
        pytest.skip('Locally generated Matchete matching outputs are absent')
    for model in m.MODELS:
        report = m.assemble(model, output_dir=tmp_path / model)
        assert report['matched_canonical_pole_term_count'] == 68
        assert report['raw_scalar_sector_counts'] == {'pole':68,'log':100,'finite':64}
        assert report['scalar_logarithm_terms_by_mass'] == {'MSX':50,'MSY':50}
        assert not report['physical_sequential_C5_computed']
        a, b = report['outputs']['hard_pole_residue_canonical'], report['outputs']['scalar_negative_pole_residue_canonical']
        assert (tmp_path / model / a).read_bytes() == (tmp_path / model / b).read_bytes()
        for file in report['outputs'].values():
            data = (tmp_path / model / file).read_text(encoding='utf-8')
            assert m.POLE not in data
        assert 'μbar2' not in (tmp_path / model / report['outputs']['fixed_muF_unprojected']).read_text(encoding='utf-8')


def test_tampering_source_pole_breaks_matching(tmp_path):
    if not (m.HARD_DIR / 'B_m1_final_eft1_fermion_first_oneloop.json').exists():
        pytest.skip('Locally generated Matchete matching outputs are absent')
    soft_dir = tmp_path / 'modified_soft'
    soft_dir.mkdir()
    source = m.SOFT_DIR / 'B_m1_final_eft1_scalar_threshold_probe_c5_relevant'
    for suffix in ('.json', '_scalar_one_loop_correction.inputform.txt', '_raw_scalar_llhh.inputform.txt'):
        target = soft_dir / (source.name + suffix)
        target.write_bytes((source.parent / (source.name + suffix)).read_bytes())
    raw = soft_dir / (source.name + '_raw_scalar_llhh.inputform.txt')
    text = raw.read_text(encoding='utf-8')
    text = text.replace('Coupling[lambdaT3, {}, 0]', 'Coupling[lambdaT3Tampered, {}, 0]', 1)
    raw.write_text(text, encoding='utf-8')
    with pytest.raises(ValueError, match='fail-closed'):
        m.assemble('B_-1', soft_dir=soft_dir, output_dir=tmp_path / 'out')


def test_permuted_epsilon_without_sign_flip_is_detected():
    original = _closed()
    a = m.canonicalize_su2_epsilon(1, original)
    b = m.canonicalize_su2_epsilon(1, original.replace(_epsilon(1,3),_epsilon(3,1)))
    assert a != b


def test_real_c_requires_dummy_relabelling():
    if not (m.HARD_DIR / 'C_m1_final_eft1_fermion_first_oneloop.json').exists():
        pytest.skip('Locally generated Matchete matching outputs are absent')
    hstem = 'C_m1_final_eft1_fermion_first_oneloop'
    sstem = 'C_m1_final_eft1_scalar_threshold_probe_c5_relevant'
    hsec = m.as_sectors(m.load_text(m.HARD_DIR / (hstem + '_raw_hard_llhh_sector.inputform.txt')), hard=True)
    ssec = m.as_sectors(m.load_text(m.SOFT_DIR / (sstem + '_raw_scalar_llhh.inputform.txt')), hard=False)
    hc = Counter(m.canonical_pole(s,t,model='C_-1',soft=False) for s,t in hsec['pole'])
    # Without the closed-dummy map, exactly half the C_-1 residues are unmatched.
    unrelabeled = Counter((-s, op) for s,t in ssec['pole'] for s,op in [
        m.canonical_pole(s,t,model='B_-1',soft=True)])
    relabeled = Counter((-s, op) for s,t in ssec['pole'] for s,op in [
        m.canonical_pole(s,t,model='C_-1',soft=True)])
    assert sum((hc & unrelabeled).values()) == 34
    assert hc == relabeled


def test_modified_matching_scope_is_rejected(tmp_path):
    if not (m.HARD_DIR / 'B_m1_final_eft1_fermion_first_oneloop.json').exists():
        pytest.skip('Locally generated Matchete matching outputs are absent')
    soft_dir = tmp_path / 'modified_soft'
    soft_dir.mkdir()
    base = 'B_m1_final_eft1_scalar_threshold_probe_c5_relevant'
    cert = m.load_json(m.SOFT_DIR / (base + '.json'))
    cert['matching_scope'] = 'CompleteGeneralT3EFT'
    (soft_dir / (base + '.json')).write_text(json.dumps(cert), encoding='utf-8')
    with pytest.raises(ValueError, match='fail-closed'):
        m.assemble('B_-1', soft_dir=soft_dir, output_dir=tmp_path / 'out')
