"""Minimal fixed-order and source-composition regressions for Batch 2."""
import numpy as np
import pytest

from Numerical.parallel_t3.wilson.strict_one_loop_batch2 import first_order_endpoint,LOOP


def test_first_order_retains_hard_exactly_once():
    t=np.zeros((2,2,3,3),complex)
    t[0,1,1,2]=t[1,0,2,1]=complex(2,-1)
    h=np.zeros_like(t)
    h[0,0,0,0]=1.4e-4j
    b1=np.zeros_like(t);b1[0,1,1,2]=3+2j
    b2=np.zeros_like(t);b2[1,0,2,1]=-1.6+0.6j
    log=-np.log(100.)
    endpoint,parts=first_order_endpoint(t,h,{'channel_a':b1,'channel_b':b2},log)
    np.testing.assert_allclose(endpoint,t+h+(b1+b2)*log/LOOP,rtol=1e-15,atol=1e-15)
    np.testing.assert_array_equal(endpoint[0,0,0,0],h[0,0,0,0])
    assert set(parts)=={'channel_a','channel_b'}


def test_first_order_rejects_bad_threshold_and_beta():
    x=np.zeros((2,2,2,2),complex)
    with pytest.raises(ValueError,match='descending'):
        first_order_endpoint(x,x,{'term':x},0.1)
    with pytest.raises(ValueError,match='ADM channel'):
        first_order_endpoint(x,x,{'term':np.zeros((3,3,2,2))},-1)
    with pytest.raises(ValueError,match='Missing'):
        first_order_endpoint(x,x,{},-1)


def test_quartic_channels_have_correct_shape_and_finite_values():
    from pathlib import Path
    from Numerical.parallel_t3.wilson.rough_stream2 import load_canonical
    from Numerical.parallel_t3.wilson.quartic_channels_batch2 import named_quartics
    project=Path(__file__).resolve().parents[1]
    manifest=(project/'output/parallel_t3/rough_inputs/fhard_triangle_all16_batch1_v1/a/A_m4/boundary.json')
    if not manifest.is_file():
        pytest.skip('Batch1 Windows output not installed on this machine')
    case=load_canonical(project,manifest.parent,'A_-4',c5_boundary='require')
    names,channels,report=named_quartics(project,'A_-4','MF',case.background.quartic[-1])
    assert names and channels.shape[0]==len(names)
    assert np.isfinite(channels).all()
    assert report['reassembly_relative']<1e-10


def test_hard_mf_quartic_channels_reconstruct_batch1():
    from pathlib import Path
    from Numerical.parallel_t3.wilson.quartic_channels_batch2 import (
        named_quartics,split_existing_f_hard)
    from Numerical.parallel_t3.wilson.rough_stream2 import load_canonical
    project=Path(__file__).resolve().parents[1]
    batch1=project/'output/parallel_t3/fhard_triangle_all16_batch1_v1'
    if not (batch1/'summary.json').is_file():
        pytest.skip('Batch1 Windows output not installed')
    case=load_canonical(project,
        project/'output/parallel_t3/rough_inputs/fhard_triangle_all16_batch1_v1/a/A_m4',
        'A_-4',c5_boundary='require')
    names,channels,_=named_quartics(project,'A_-4','MF',case.background.quartic[-1])
    pieces,audit=split_existing_f_hard(project,'A_-4',batch1,names,channels)
    assert len(pieces)==len(names)
    assert audit['relative_reproduction_batch1']<1e-9
