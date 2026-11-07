"""Offline checks for the Batch 4 source audit and analytic fixed-order lift."""
from __future__ import annotations

import math
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from Numerical.running.T3OrdinaryHardFlavor import (
    HardFlavorError, LOOP, mass_kernel, source_pole_matrix,
)
from Numerical.running.T3OrdinaryInheritedFixedOrder import (
    audit_inherited_source, matching_kernel_at_heavy_mass, top_level_terms,
)


def mock_source() -> str:
    """Six structurally faithful archived one-loop source terms."""
    common = 'Coupling[MF, {}, 0]'
    terms=[]
    for group in ('holo','anti'):
        if group=='holo':
            yuk='Bar[Coupling[y1, {Index[r, Flavor], Index[x, NFlavor]}, 0]]*Bar[Coupling[y2, {Index[s, Flavor], Index[x, NFlavor]}, 0]]*Coupling[lambdaT3, {}, 0]'
        else:
            yuk='Coupling[y1, {Index[r, Flavor], Index[x, NFlavor]}, 0]*Coupling[y2, {Index[s, Flavor], Index[x, NFlavor]}, 0]*Bar[Coupling[lambdaT3, {}, 0]]'
        pre=f'2*hbar*{yuk}*Field[H, Scalar, {{}}, {{}}]*Field[H, Scalar, {{}}, {{}}]'
        tail='NCM[Field[l, Fermion, {}, {}], DiracProduct[GammaCC], Field[l, Fermion, {}, {}]]'
        for kind in ('finite','pole','log'):
            log=f'Log[μbar2/{common}^2]*' if kind=='log' else ''
            eps='ϵ*' if kind=='pole' else ''
            terms.append(f'+ ({pre}*{log}{tail})/(Sqrt[3]*{eps}{common})')
    return '\n'.join(terms)


def test_parser_splits_only_outer_plus():
    expr=mock_source()
    terms=top_level_terms(expr)
    assert len(terms)==6
    assert all('DiracProduct[GammaCC]' in t for t in terms)


@pytest.mark.parametrize('model', ['D_+2','A_-4'])
def test_archived_triplets_pass(tmp_path,model):
    p=tmp_path/'F_one_loop.txt'
    p.write_text(mock_source(),encoding='utf-8')
    a=audit_inherited_source(p,model)
    assert a['pure_sm_weinberg_term_count']==6
    assert a['pole_cancellation_at_common_couplings']
    assert not a['complete_counterterm_and_scheme_certified']


def test_missing_pole_fails(tmp_path):
    p=tmp_path/'F_one_loop.txt'
    p.write_text(mock_source().replace('ϵ*',''),encoding='utf-8')
    with pytest.raises(HardFlavorError,match='duplicate|incomplete'):
        audit_inherited_source(p,'D_+2')


def test_bad_relative_sign_fails(tmp_path):
    p=tmp_path/'F_one_loop.txt'
    p.write_text(mock_source().replace('+ (','- (',1),encoding='utf-8')
    with pytest.raises(HardFlavorError,match='negative'):
        audit_inherited_source(p,'A_-4')


def test_mismatched_log_coefficient_fails(tmp_path):
    p=tmp_path/'F_one_loop.txt'
    p.write_text(mock_source().replace('2*hbar*','3*hbar*',1),encoding='utf-8')
    with pytest.raises(HardFlavorError,match='differ|normalization'):
        audit_inherited_source(p,'A_-4')


def test_tripled_index_mass_fails(tmp_path):
    p=tmp_path/'F_one_loop.txt'
    p.write_text(mock_source().replace('Coupling[MF, {}, 0]',
      'Coupling[MF, {Index[x, NFlavor]}, 0]'),encoding='utf-8')
    with pytest.raises(HardFlavorError,match='heavy-mass structure|denominator'):
        audit_inherited_source(p,'D_+2')


def test_two_matrix_pole_cancellation_with_complex_yukawa():
    rng=np.random.default_rng(4112)
    y1=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
    y2=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
    mf=np.array([8.5e4,1.7e5,2.3e5])
    lam=complex(.015,-.001)
    hard,components=source_pole_matrix(y1,y2,mf,lam)
    inherited=sum((-c for c in components),np.zeros((3,3),complex))
    assert np.linalg.norm(hard+inherited) < 1e-25
    assert np.linalg.norm(hard-hard.T) < 1e-25


@pytest.mark.parametrize('a,b,mf',[(1e6,1.3e6,1e5),(2e6,2e6,1e5),
                                    (1e6,1.00000001e6,2e5)])
def test_inherited_log_identity_and_scale_invariance(a,b,mf):
    muF=1e5
    for muS in (700.,1000.,3000.):
        explicit=(mass_kernel(a,b,muS)-2*math.log(muS/muF)
                  -(1+2*math.log(muF/mf)))
        compact=matching_kernel_at_heavy_mass(a,b,mf)
        assert math.isclose(explicit,compact,rel_tol=1e-12,abs_tol=2e-12)


def test_evolved_coupling_pole_does_not_cancel_by_fiat():
    y=np.eye(3,dtype=complex)*1e-3
    mf=np.array([1e5,1.2e5,1.5e5])
    pf,_=source_pole_matrix(y,y,mf,.01)
    ps,_=source_pole_matrix(y,y,mf,.012)
    assert np.linalg.norm(ps-pf) > 0
    assert math.isclose(np.linalg.norm(ps-pf)/np.linalg.norm(pf),.2,rel_tol=1e-12)


def test_source_is_non_empty_and_balanced():
    with pytest.raises(HardFlavorError):
        top_level_terms('')
    with pytest.raises(HardFlavorError):
        top_level_terms('('+'x'*150)


@pytest.mark.parametrize('model,tag', [('D_+2','T3_D_alpha_p2'),
                                        ('A_-4','T3_A_alpha_m4')])
def test_repository_archived_f_loop_if_present(model,tag):
    """Also test the real legacy Matchete archive on a user's local checkout."""
    path=ROOT/'output/full16_symbolic'/tag/'data/eft1_after_F_one_loop_inputform.txt'
    if not path.is_file():
        pytest.skip('Real archived Matchete first threshold is not shipped in this small code ZIP')
    result=audit_inherited_source(path,model)
    assert result['pure_sm_weinberg_term_count']==6
    assert result['holomorphic_triplet']==['finite','pole','log']
