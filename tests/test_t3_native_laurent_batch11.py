"""Batch 11 native-source algebra tests. Saved-source tests run when present."""
from __future__ import annotations

from pathlib import Path
import math

import numpy as np
import pytest
import sympy as sp

from Numerical.parallel_t3.complete.two_model_native_laurent_batch11 import (
    ONELOOP, SOFT_PATHS, _canonical_soft, parse_laurent,
    evaluate_soft, combine_native, hard_from_pole_at_scale,
)

C12='Coupling[T3EFT1W12, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]'
CEX='Coupling[T3EFT1WExtra, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]'
M1='Coupling[MS1, {}, 0]'
M2='Coupling[MS2, {}, 0]'
LAM='Coupling[lambdaT3, {}, 0]'
LE='Coupling[lambdaT3ExtraQ020101I1, {}, 0]'


def a_native_soft_text():
    """Independent rewriting of the 10-term saved A soft native formula."""
    delta=f'({M1}^2-{M2}^2)'
    return (f'(-2*Sqrt[3]*hbar*Bar[{C12}]*{LAM}/3)*('
            f'(-2*{M1}^2-{M2}^2)/{delta}'
            f'+({M1}^2-{M2}^2)/(ϵ*{delta})'
            f'+({M1}^2*Log[μbar2/{M1}^2]-'
            f'{M2}^2*Log[μbar2/{M2}^2])/{delta})')


def test_exact_A_pole_coefficient_and_scale_dependence():
    regular,pole=parse_laurent(a_native_soft_text(),'A_-4')
    W=(np.array([[1+2j,2j,3], [4,5,6], [7,8j,9]],complex))*1e-9
    lam=0.13+0.017j
    params=dict(w12=W,wextra=np.zeros_like(W),lam=lam,extra=0,
                m1=1290.,m2=1100.,mminus=1100.,mplus=1100.)
    reg,P=evaluate_soft(regular,pole,**params,mu_bar=1000.)
    expected=(-2/math.sqrt(3))*ONELOOP*lam*W.conj()
    np.testing.assert_allclose(P,expected,rtol=2e-12,atol=1e-24)
    reg2,_=evaluate_soft(regular,pole,**params,mu_bar=777.)
    np.testing.assert_allclose(reg2-reg,2*P*math.log(.777),rtol=2e-12,atol=1e-23)


def test_pole_closed_source_split_is_matching_scale_invariant():
    P=np.array([[1+1j,2,3],[2,4+4j,5],[3,5,6]],complex)*1e-15
    h=P.copy()
    s=-P.copy()
    hreg=np.array([[3,1j,0],[1j,4,2j],[0,2j,5]],complex)*1e-15
    sreg=np.array([[1,2j,0],[2j,0,1j],[0,1j,2]],complex)*1e-15
    a=combine_native(hreg,h,sreg,s,1e5,1e3)
    assert a['source_pole_closure_relative']==0.
    f1=1.2e5;s1=810.
    hreg1=hard_from_pole_at_scale(hreg,h,1e5,f1)
    sreg1=sreg+2*s*math.log(s1/1e3)
    b=combine_native(hreg1,h,sreg1,s,f1,s1)
    np.testing.assert_allclose(a['native_finite_source_sector'],
                               b['native_finite_source_sector'],rtol=1e-12,atol=1e-27)
    assert np.linalg.norm(a['source_derived_eft1_scale_log'])>0


def test_D_extra_complex_phase_is_not_silently_discarded():
    mm='Coupling[T32MinusMass, {}, 0]'
    mp='Coupling[T32PlusMass, {}, 0]'
    # Small source-shaped test of the D scalar coupling and mass lexer.
    source=(f'Sqrt[3]*hbar*Bar[{LE}]*Bar[{CEX}]*'
            f'({M1}^2+{mm}^2+{mp}^2)/(ϵ*{M1}^2)'
            f'+Sqrt[3]*hbar*Bar[{LE}]*Bar[{CEX}]*'
            f'Log[μbar2/{M1}^2]+0*Bar[{C12}]*{LAM}')
    reg,pole=parse_laurent(source,'D_-2')
    x=np.eye(3,dtype=complex)*(1+2j)*1e-10
    kw=dict(w12=np.zeros((3,3),complex),wextra=x,lam=0.1,
            extra=0.03-0.02j,m1=1400.,m2=1000.,
            mminus=980.,mplus=1020.,mu_bar=1000.)
    r,p=evaluate_soft(reg,pole,**kw)
    assert np.linalg.norm(r)>0 and np.linalg.norm(p)>0
    assert abs(p[0,0].imag)>0
    np.testing.assert_allclose(p[0,1],0.)


def test_unknown_native_coupling_is_rejected():
    malformed=a_native_soft_text()+ '+Coupling[SomeOtherPortal,{},0]'
    with pytest.raises(ValueError,match='unrecognised'):
        _canonical_soft(malformed,'A_-4')


@pytest.mark.parametrize('model',['A_-4','D_-2'])
def test_saved_native_source_parses_if_present(model):
    # With the full Windows checkout these are the actual archived 10/68-term
    # Matchete soft sources, not synthetic loop functions.
    project=Path(__file__).resolve().parents[1]
    file=project/SOFT_PATHS[model]
    if not file.is_file():
        pytest.skip('Native source not in standalone test environment')
    r,p=parse_laurent(file.read_text(encoding='utf-8-sig'),model)
    assert p != 0 and r != 0
    assert not p.has(sp.Symbol('EPS'))
