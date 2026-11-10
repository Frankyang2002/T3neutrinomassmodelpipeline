"""Batch 5 arithmetic and native-subtraction scenario regressions.

The repo's original hard_matching.py is exercised when installed; the tests
never require or mutate any saved 16-model numerical output.
"""
from __future__ import annotations
import numpy as np
import pytest

from Numerical.parallel_t3.wilson.matching_soft_subtraction_batch5 import (
    _rel,native_laurent_from_evaluated,comparison_with_native,
    triangle_soft_ledger)


def _mock_native(mass=(90.,100.,110.),mu=100.):
    m=np.asarray(mass,float)
    coeff=np.zeros((2,3,3,3),complex)
    coeff[0,0,0,0]=.002+.0002j
    coeff[0,1,1,1]=-.0004+.0001j
    coeff[1,2,2,2]=.0012-.0007j
    # The source's *exact* raw pole=constant=log Laurent structure.
    factor=1+2*np.log(mu/m)
    return coeff*factor[None,:,None,None],coeff


def test_native_source_Laurent_constant_log_and_pole_close():
    per,raw=_mock_native()
    ans=native_laurent_from_evaluated(per,np.asarray([90.,100.,110.]),100.)
    np.testing.assert_allclose(ans['per_native_kernel_and_fermion_constant'],raw,rtol=1e-13,atol=1e-15)
    np.testing.assert_allclose(ans['pole_residue'],raw.sum(axis=(0,1)),rtol=1e-13,atol=1e-15)
    np.testing.assert_allclose(ans['finite_constant']+ans['finite_scale_log'],per.sum(axis=(0,1)),rtol=1e-13,atol=1e-15)
    assert ans['native_finite_is_MSbar_hard'] is False


def test_native_laurent_rejects_near_zero_denominator_and_nonfinite():
    m=np.array([100*np.exp(.5),100.,110.])
    with pytest.raises(ValueError,match='vanishing'):
        native_laurent_from_evaluated(np.zeros((1,3,3,3),complex),m,100.)
    bad=np.zeros((1,3,3,3),complex);bad[0,0,0,0]=np.nan
    with pytest.raises(ValueError,match='geometry'):
        native_laurent_from_evaluated(bad,[90.,100.,110.],100.)


def test_residual_pole_not_silently_set_to_zero():
    z=np.eye(3,dtype=complex)*1e-9
    analytic={'scales':{1.:{'hard_local_dim5':z,'EFT1_soft_finite':z/4}},
              'eft_soft_pole_at_MF':-z/2}
    native={'formal_pole_dropped':z*2,'pole_residue':z/3}
    result=comparison_with_native(analytic,native)
    np.testing.assert_allclose(result['hypothetical_hard_MF'],z*1.75)
    np.testing.assert_allclose(result['unaccounted_pole_or_counterterm_residual'],-z/6)
    assert result['pole_cancellation_relative']>0
    assert result['hypothesis_is_MSbar_hard'] is False
    assert result['counterterms_applied'] is False


def test_native_source_does_not_infer_wrong_tensor_geometry():
    with pytest.raises(ValueError,match='geometry'):
        native_laurent_from_evaluated(np.ones((4,2,3,3)),[90,100,110],100)


def test_real_source_kernel_renormalization_scale_cancellation():
    pytest.importorskip('Numerical.lo_certification.completeness.hard_matching')
    pytest.importorskip('Numerical.parallel_t3.threshold.matching')
    # Synthetic flavour-symmetric LLSS and real HHSS tensor. These are not
    # production gauge-invariant multiplets, so tolerate their SU2 projection;
    # scale cancellation is independent of this model geometry.
    C=np.zeros((3,3,3,2,2,2,2),complex)
    Q=np.zeros((4,4,2,2),float)
    for k in range(3):
        C[k,k,k,0,0,0,0]=(0.04+0.01j)/(100+k*10)
        C[k,k,k,1,1,1,1]=(0.04+0.01j)/(100+k*10)
    Q[2,2,0,0]=0.2
    Q[0,0,1,1]=0.2
    masses=np.array([90.,100.,110.])
    mass2=np.diag([40.**2,50.**2])
    ans=triangle_soft_ledger(C,masses,Q,mass2,100.,projection_tolerance=2.)
    assert ans['fermion_count']==3
    assert ans['scalar_pairs_evaluated']==3*3  # per fermion, unordered 2x2 scalar pairs
    assert ans['matching_scale_cancellation_relative']<1e-10
    assert ans['full_UV_scale_invariance_relative']<1e-10
    assert ans['dimension7plus_relative_to_full_triangle']>=0
    for scale in (0.5,1.,2.):
        row=ans['scales'][scale]
        np.testing.assert_allclose(
            row['hard_local_dim5']+row['EFT1_soft_finite']+
            row['power_suppressed_dim7plus_remainder'],
            row['full_UV_triangle'],rtol=1e-10,atol=1e-15)


def test_native_formal_kernel_need_not_equal_analytic_UV_hard():
    per,_=_mock_native()
    native=native_laurent_from_evaluated(per,[90.,100.,110.],100.)
    assert _rel(native['formal_pole_dropped'],native['pole_residue'])>0
    assert not native['native_finite_is_MSbar_hard']
