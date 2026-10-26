"""Stage17 polynomial group tensors and operator-mixing sanity tests."""
from pathlib import Path
import sys
import sympy as sp
import pytest
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from Numerical.experimental.GeneralT3ExtraEFT1Mixing import (
    ACTIVE_RGBETA, _real_fields, extra_potential_polynomial, real_quartic_entries,
    extra_scalar_pair_component, _assert_linear_in_quartics, _classify_flavor,
    calculate_extra_eft1_mixing,
)


def test_all_quartic_components_real_and_symmetric():
    d=real_quartic_entries()
    assert d
    assert all(tuple(sorted(k))==k and len(k)==4 for k in d)
    for v in d.values():
        assert sp.expand(sp.conjugate(v)-v)==0


def test_all_terms_mix_two_higgs_and_two_s1_or_s2():
    for k in real_quartic_entries():
        assert sum(1<=i<=4 for i in k)==2
        assert sum(5<=i<=8 for i in k)==2 or sum(9<=i<=12 for i in k)==2
        assert not (any(5<=i<=8 for i in k) and any(9<=i<=12 for i in k))


def test_quartic_separation_and_additivity():
    both=real_quartic_entries()
    s2=real_quartic_entries(include_s1=False)
    s1=real_quartic_entries(include_s2=False)
    assert set(s1).isdisjoint(s2)
    assert both==dict(s1,**{}) | s2
    assert real_quartic_entries(include_s1=False,include_s2=False)=={}


def test_lexicographic_derivative_tensor_matches_potential():
    r=_real_fields()
    poly=extra_potential_polynomial()
    tensor=real_quartic_entries()
    for k in [(1,1,9,9),(1,1,7,7),(2,2,5,5),(1,2,9,10)]:
        expected=sp.diff(poly,*[r[i-1] for i in k])
        actual=tensor.get(tuple(sorted(k)),sp.S.Zero)
        assert sp.expand(actual-expected)==0


def test_identity_s2_and_epsilon_s1_signs():
    s2=real_quartic_entries(include_s1=False)
    s1=real_quartic_entries(include_s2=False)
    l2,l1=map(sp.Symbol,ACTIVE_RGBETA)
    # 4 real fields factors of (1/sqrt(2))^4, plus HC.
    assert sp.simplify(s2[(1,1,9,9)] - (l2+sp.conjugate(l2)))==0
    assert sp.simplify(s1[(1,1,7,7)] - (l1+sp.conjugate(l1)))==0
    assert sp.simplify(s1.get((1,1,5,5),sp.S.Zero))==0  # antisymmetric SU2 contraction


def test_scalar_pair_delta_on_synthetic_wilson():
    tensor=real_quartic_entries(include_s1=False)
    y2=sp.Symbol('y2')
    class Wilson:
        def __call__(self,i,j,c,d):
            return sp.conjugate(y2)**2 if (i,j,c,d)==(1,2,9,9) else sp.S.Zero
    result=extra_scalar_pair_component(tensor,Wilson(),1,2,1,1)
    l2=sp.Symbol(ACTIVE_RGBETA[0])
    assert sp.simplify(result-(l2+sp.conjugate(l2))*sp.conjugate(y2)**2)==0
    assert _classify_flavor(result)==['y2']
    assert _assert_linear_in_quartics(result)


def test_quartic_linearity_checker_rejects_square():
    lam=sp.Symbol(ACTIVE_RGBETA[0]); y=sp.Symbol('y1')
    assert _assert_linear_in_quartics(lam*sp.conjugate(y)**2)
    assert not _assert_linear_in_quartics(lam**2)
    assert not _assert_linear_in_quartics(1+lam)


@pytest.mark.parametrize('model',['B_m1','C_m1'])
def test_missing_real_certification_rejected(model):
    with pytest.raises(ValueError,match='certified Stage15'):
        calculate_extra_eft1_mixing(model,{}, {}, {})


def test_project_seed_analysis_if_available():
    for model, d in [('B_m1','B'),('C_m1','C')]:
        folder=ROOT/f'output/T3_{d}_alpha_m1/data'
        stage15=ROOT/'output/general_t3_stage15/source_lagrangian_convention.json'
        if not (folder/'eft1_after_F_wilson_seed.json').exists() or not stage15.exists():
            pytest.skip('Matchete Wilson seeds and Stage15 certificate are in the user repository')
        import json
        result=calculate_extra_eft1_mixing(model,
          json.loads((folder/'eft1_after_F_wilson_seed.json').read_text()),
          json.loads((folder/'eft1_rgbeta_rge.json').read_text()),
          json.loads(stage15.read_text()))
        assert result['nonzero_weinberg_real_components']>0
        assert result['extra_low_energy_neutrino_mass'] is None
