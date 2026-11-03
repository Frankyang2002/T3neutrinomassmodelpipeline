"""Batch 8: exact derivative of finite-mixing kernel, relative tensor residuals,
complex fit derivative, fail-closed all-16 integration and source protection."""
from __future__ import annotations
import copy
import json
import math
import os
from pathlib import Path
import numpy as np
import pytest
from Numerical.lo_patterns_16 import splitting_convergence as s8
from Numerical.lo_patterns_16 import factorization_theorem as t7
from Numerical.lo_patterns_16 import kernel_factorization as k6
from Numerical.lo_patterns_16 import texture_synthesis as t5
from Numerical.lo_patterns_16.analytic_explanation import scalar_kernels


@pytest.mark.parametrize('m2,a,b,delta',[
    (2.0e4,1.2e6,8.7e5,1.1e4),
    (1.5e7,1.3e6,9.1e5,-3200+1000j),
    (8.1e6,4.1e5,9.1e5,2200-1700j),
    (6.4e9,8.1e5,1.4e6,-5300j),
])
def test_exact_kernel_partials_match_symmetric_finite_differences(m2,a,b,delta):
    k,da,db=s8.exact_kernel_derivatives(m2,a,b,delta)
    assert k==pytest.approx(scalar_kernels(m2,a,b,delta)[0],rel=1e-12)
    h=3e-5
    dda=(scalar_kernels(m2,a*(1+h),b,delta)[0]-
         scalar_kernels(m2,a*(1-h),b,delta)[0])/(2*a*h)
    ddb=(scalar_kernels(m2,a,b*(1+h),delta)[0]-
         scalar_kernels(m2,a,b*(1-h),delta)[0])/(2*b*h)
    assert abs(da-dda)<1e-5*max(abs(da),abs(dda),1e-30)
    assert abs(db-ddb)<1e-5*max(abs(db),abs(ddb),1e-30)


def test_exact_derivative_contains_mixing_correction_absent_from_batch7():
    m2,a,b,d=7e6,1.5e6,9e5,110000.
    _,fa,fb=s8.exact_kernel_derivatives(m2,a,b,d)
    _,da,db=t7.derivative_D(m2,a,b)
    assert abs(fa-d*da)>1e-6*abs(fa)
    assert abs(fb-d*db)>1e-6*abs(fb)


def test_mixing_phase_and_sign_covariance():
    a=s8.exact_kernel_derivatives(7e6,1e6,8e5,2500+800j)
    b=s8.exact_kernel_derivatives(7e6,1e6,8e5,2500-800j)
    c=s8.exact_kernel_derivatives(7e6,1e6,8e5,-2500-800j)
    for x,y,z in zip(a,b,c):
        assert y==pytest.approx(x.conjugate(),rel=1e-13)
        assert z==pytest.approx(-x,rel=1e-13)


@pytest.mark.parametrize('bad',(0.,-10.,float('nan')))
def test_kernel_invalid_masses_rejected(bad):
    with pytest.raises(s8.SplittingConvergenceError):
        s8.exact_kernel_derivatives(1e6,bad,1e5,2000)


def test_tachyonic_scalar_rejected():
    with pytest.raises(s8.SplittingConvergenceError,match='Tachyonic'):
        s8.exact_kernel_derivatives(1e6,1e4,1e4+50,2e4)


def test_exact_tensor_linearization_quadratic_remainder():
    rng=np.random.default_rng(67)
    a=rng.normal(size=(3,3,3))+1j*rng.normal(size=(3,3,3))
    b=rng.normal(size=(3,3,3))+1j*rng.normal(size=(3,3,3))
    c=rng.normal(size=(3,3,3))+1j*rng.normal(size=(3,3,3))
    errs=[]
    for t in (.5,.25,.125):
        exact=a+t*b+t*t*c
        errs.append(s8.norm(exact-a-t*b))
    assert errs[0]/errs[1]==pytest.approx(4.,rel=1e-12)
    assert errs[1]/errs[2]==pytest.approx(4.,rel=1e-12)


def test_complex_projection_derivative_matches_finite_difference():
    rng=np.random.default_rng(59)
    ref=rng.normal(size=(3,3,3))+1j*rng.normal(size=(3,3,3))
    ref1=.1*(rng.normal(size=(3,3,3))+1j*rng.normal(size=(3,3,3)))
    other=(1.7+.2j)*ref+0.2*(rng.normal(size=(3,3,3))+1j*rng.normal(size=(3,3,3)))
    other1=.2*(rng.normal(size=(3,3,3))+1j*rng.normal(size=(3,3,3)))
    analysis=s8.fit_with_derivative(ref,other,ref1,other1)
    h=1e-5
    plus=s8.point_fit(ref+h*ref1,other+h*other1)
    minus=s8.point_fit(ref-h*ref1,other-h*other1)
    dk=(plus['kappa']-minus['kappa'])/(2*h)
    dd=(plus['transverse']-minus['transverse'])/(2*h)
    assert dk==pytest.approx(analysis['kappa_derivative'],rel=2e-8)
    assert s8.norm(dd-analysis['transverse_derivative'])<2e-8*s8.norm(analysis['transverse_derivative'])
    assert abs(np.vdot(s8.flatten(ref),analysis['baseline_transverse']))<1e-11*s8.norm(ref)*s8.norm(other)


def test_proportional_tensors_have_no_transverse_signal():
    rng=np.random.default_rng(88)
    a=rng.normal(size=(3,3,3))+1j*rng.normal(size=(3,3,3))
    a1=.1*rng.normal(size=(3,3,3))
    d=s8.fit_with_derivative(a,(3+1j)*a,a1,(3+1j)*a1)
    assert s8.norm(d['baseline_transverse'])<1e-13*s8.norm(a)
    assert s8.norm(d['transverse_derivative'])<1e-13*s8.norm(a1)


def test_slope_guard_and_noise_floor():
    assert s8._slope({.125:1e-6,.25:4e-6},.125,.25,1e-12)==pytest.approx(2.)
    assert s8._slope({.125:1e-16,.25:4e-16},.125,.25,1e-12) is None


def test_split_input_respects_endpoints_and_never_mutates_source():
    from test_lo_kernel_factorization_16 import base_toy
    inp,_,d1=base_toy('A_-2')
    before=json.dumps(inp,sort_keys=True)
    for t in [0,.125,1.]:
        p=s8.split_input(inp,d1,t)
        A=t5.decode(inp['scalar_A_gev2'],'A')
        S=t5.decode(k6.deform_scalars(inp,d1,'charge_split_diagnostic')['scalar_A_gev2'],'S')
        assert np.allclose(t5.decode(p['scalar_A_gev2'],'split'),A+t*(S-A))
    assert json.dumps(inp,sort_keys=True)==before
    with pytest.raises(s8.SplittingConvergenceError,match='outside'):
        s8.split_input(inp,d1,2)


def test_output_dir_refuses_overwrite_or_outside(tmp_path):
    d=tmp_path/s8.TARGET;d.mkdir(parents=True)
    (d/'keep').write_text('safe')
    with pytest.raises(s8.SplittingConvergenceError,match='already populated'):
        s8.analyse(tmp_path)
    assert (d/'keep').read_text()=='safe'
    with pytest.raises(s8.SplittingConvergenceError,match='inside Chat A'):
        s8.analyse(tmp_path,output_dir=tmp_path/'outside')


def make_synthetic_results(tmp_path):
    from test_lo_factorization_theorem_16 import _synthetic_result
    root=_synthetic_result(tmp_path)
    t7.analyse(root)
    return root


def test_synthetic_all16_full_workflow_and_readonly(tmp_path):
    root=make_synthetic_results(tmp_path)
    sha_before={str(p):k6.digest(p) for p in root.rglob('*') if p.is_file()}
    p=s8.analyse(root)
    assert p['model_count']==16
    assert len(p['points'])==16*len(s8.TS)
    assert len(p['derivative_checks'])==16
    assert max(z['baseline_saved_reconstruction_relative_error'] for z in p['derivative_checks'])<5e-8
    assert max(z['endpoint_saved_reconstruction_relative_error'] for z in p['derivative_checks'])<5e-8
    assert all(z['exact_remainder_order'] is None or 1.5<z['exact_remainder_order']<2.5
               for z in p['model_summary'])
    assert all((r['transverse_increment_fractional_error'] is None or
                r['transverse_increment_fractional_error']>=0) for r in p['points'])
    for path,old_sha in sha_before.items():
        assert k6.digest(Path(path))==old_sha
    out=root/s8.TARGET
    for name in ('BATCH08_SPLITTING_CONVERGENCE.md','batch08_convergence.json',
                 'batch08_scaling_points.csv','batch08_model_summary.csv',
                 'batch08_derivative_validation.csv','SOURCE_MANIFEST.json','HANDOFF.json'):
        assert (out/name).is_file(),name
    assert '0/16' in (out/'BATCH08_SPLITTING_CONVERGENCE.md').read_text()


def test_upstream_batch7_tamper_blocks_and_does_not_write(tmp_path):
    root=make_synthetic_results(tmp_path)
    p=root/s8.BATCH07/'SOURCE_MANIFEST.json'
    p.write_text(p.read_text()+'\n')
    with pytest.raises(s8.SplittingConvergenceError,match='Handoff manifest'):
        s8.analyse(root)
    assert not (root/s8.TARGET).exists()


def test_source_input_tamper_blocks_report(tmp_path):
    root=make_synthetic_results(tmp_path)
    p=root/'output/general_observable/lo_direct_one_loop/provisional_inputs/C_+1.json'
    p.write_text(p.read_text()+'\n')
    with pytest.raises(s8.SplittingConvergenceError,match='Source SHA'):
        s8.analyse(root)
    assert not (root/s8.TARGET).exists()


def test_fake_physical_certification_blocks(tmp_path):
    root=make_synthetic_results(tmp_path)
    p=root/s8.BATCH07/'batch07_theorem.json'
    obj=json.loads(p.read_text());obj['physical_normalisation_certified']=True
    p.write_text(json.dumps(obj))
    with pytest.raises(s8.SplittingConvergenceError,match='Physical certification'):
        s8.analyse(root)
    assert not (root/s8.TARGET).exists()


@pytest.mark.skipif(os.environ.get('RUN_LO_PATTERN_INTEGRATION')!='1',
                    reason='Opt-in real local checkout; synthetic 16-model integration runs by default')
def test_real_checkout_integration():
    root=Path(__file__).resolve().parents[1]
    p=s8.analyse(root,output_dir=root/s8.BASE/'analysis_batch08_pytest_manual')
    assert p['model_count']==16
