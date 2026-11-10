"""Fast source-format and strict-perturbative-order tests for Batch13."""
from __future__ import annotations

import importlib.util
from pathlib import Path
import numpy as np
import pytest

P=Path(__file__).resolve().parents[1]/'Numerical/parallel_t3/complete/two_model_master_wilson_batch13.py'
spec=importlib.util.spec_from_file_location('candidate_batch13',P)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def test_source_to_full_master_pair_tensor():
    labels=[f'L[{p},{i}]' for p in range(1,4) for i in (1,2)]
    raw=np.zeros((3,3,2,2,8,8),complex)
    raw[0,1,0,0,0,1]=raw[0,1,0,0,1,0]=2+3j
    raw[1,0,0,0,0,1]=raw[1,0,0,0,1,0]=2+3j
    c=m.source_to_master(raw,tuple(labels),np.eye(8),'A_-4')
    assert c.shape==(6,6,12,12)
    assert c[0,2,4,5]==2+3j
    assert np.allclose(c,c.swapaxes(0,1))
    assert np.allclose(c,c.swapaxes(2,3))
    assert not np.any(c[:,:,:4,:4])


def test_rotation_is_orthogonal_and_does_not_change_scalar_norm():
    labels=[f'L[{p},{i}]' for p in range(1,4) for i in (1,2)]
    raw=np.zeros((3,3,2,2,8,8),complex)
    raw[0,0,0,0,1,1]=1j
    rot=np.eye(8)
    rot[[0,1],:]=rot[[1,0],:]
    out=m.source_to_master(raw,tuple(labels),rot,'D_-2')
    assert out[0,0,4,4]==1j
    assert np.linalg.norm(out)==pytest.approx(np.linalg.norm(raw))
    with pytest.raises(ValueError,match='orthogonal'):
        m.source_to_master(raw,tuple(labels),2*np.eye(8),'D_-2')


def test_log_replacement_never_double_counts_native_leading_log():
    a=np.eye(3,dtype=complex)*(10+2j)
    log=np.eye(3,dtype=complex)*(4-1j)
    master=np.eye(3,dtype=complex)*(3+2j)
    got=m.replace_eft1_log(a,log,master,2.)
    assert np.allclose(got,a-log+2*master)
    assert not np.allclose(got,a+log+2*master)
    assert np.allclose(m.replace_eft1_log(a,log,log,1),a)


def test_invalid_candidate_normalization_blocks():
    a=np.eye(3,dtype=complex)
    for x in (0.,float('nan'),float('inf')):
        with pytest.raises(ValueError,match='normalization'):
            m.replace_eft1_log(a,a,a,x)


def test_candidate_aliases_cover_native_archived_a_and_d():
    assert m.CANDIDATE_SOURCE_ALIASES['λ']=='lambdaH'
    assert m.CANDIDATE_SOURCE_ALIASES['lambdaS1Inv2']=='lambdaS1Adj'
    assert m.CANDIDATE_SOURCE_ALIASES['lambdaT3']=='lambdaT3'
    assert m.BASIS_SIZES=={'A_-4':(1,3),'D_-2':(3,1)}


def test_hhss_projection_uses_neutral_higgs_real_index_2():
    labels=tuple(f'L[{p},{i}]' for p in range(1,4) for i in (1,2))
    C=np.zeros((6,6,12,12),complex)
    C[0,2,2,2]=C[2,0,2,2]=2j
    result=m.projection(C,labels)
    assert result[0,1]==8j
    assert result[1,0]==8j
    assert result[1,1]==0


def test_optional_repo_complete_integration_contract():
    # The real source files exist only in the user's Windows project checkout.
    root=Path(__file__).resolve().parents[1]
    if not (root/m.B9/'A_-4/result.json').is_file():
        pytest.skip('Full project Batch9/11 source outputs only in local checkout')
    assert (root/m.B11/'A_-4/result.json').is_file()


def test_d2_three_additional_scalar_only_quartics_retained():
    q=np.zeros((12,)*4,float)
    coefficients={k:{'real':0.01*(i+1),'imag':0.002*i}
                  for i,k in enumerate((
                      'lambdaT3GeneralQ000004I1',
                      'lambdaT3GeneralQ000013I1',
                      'lambdaT3GeneralQ001102I1'))}
    info=m.add_d2_extra_scalar_only(q,coefficients)
    assert info['D2_extra_scalar_only_count']==3
    assert not info['D2_scalar_only_extra_RGBeta_to_Matchete_map_certified']
    assert np.linalg.norm(q)>0
    assert np.allclose(q,q.swapaxes(0,1),atol=1e-14)
    assert np.allclose(q,q.swapaxes(0,2),atol=1e-14)
    assert np.allclose(q,q.swapaxes(0,3),atol=1e-14)
    assert np.all(q[:4]==0)


def test_native_batch1_trajectory_numpy_coupling_matrix():
    """Real Batch1.at_fermion_threshold gives complex NumPy matrices."""
    matrix=np.array([[1+2j,0,0],[0,3-4j,0],[0,0,5j]],dtype=np.complex128)
    actual=m.decode_complex(matrix)
    assert actual.shape==(3,3)
    assert np.array_equal(actual,matrix)


def test_native_batch1_trajectory_numpy_scalar_and_json_forms():
    """Both NPZ 0d scalars and the two native JSON layouts are accepted."""
    scalar=np.asarray(0.125-0.02j,dtype=np.complex128)
    assert m.decode_complex(scalar).shape==()
    assert m.decode_complex(scalar).item()==scalar.item()
    assert m.decode_complex(np.complex128(0.375+0.05j)).item()==0.375+0.05j
    obj={'real':[[1,2],[3,4]],'imag':[[0,1],[-1,0]]}
    assert np.array_equal(m.decode_complex(obj),np.array([[1,2+1j],[3-1j,4]]))
    nested=[[{'real':1,'imag':2}, {'real':3,'imag':0}]]
    # general_batch1_api is available in the user's checkout, not necessarily
    # in this standalone release staging.
    if (P.parents[3]/'Numerical/parallel_t3/background/general_batch1_api.py').is_file():
        assert np.array_equal(m.decode_complex(nested),np.array([[1+2j,3]]))


def test_sm_yukawa_tensor_accepts_numpy_threshold_state():
    """Regression for the exact line that blocked both A(-4) and D(-2)."""
    labels=[]
    for species in ('Q','uc','dc','L','ec'):
        for gen in range(1,4):
            if species=='Q':
                labels.extend(f'Q[{gen},{c},{w}]' for c in range(1,4) for w in range(1,3))
            elif species in ('uc','dc'):
                labels.extend(f'{species}[{gen},{c}]' for c in range(1,4))
            elif species=='L':
                labels.extend(f'L[{gen},{w}]' for w in range(1,3))
            else:
                labels.append(f'ec[{gen}]')
    assert len(labels)==45
    state={'yu':np.eye(3,dtype=np.complex128)*0.9,
           'yd':np.eye(3,dtype=np.complex128)*0.03,
           'ye':np.eye(3,dtype=np.complex128)*0.01}
    Y=m.sm_yukawa_tensor(state,tuple(labels))
    assert Y.shape==(45,45,12)
    assert np.isfinite(Y).all()
    assert np.array_equal(Y,Y.swapaxes(0,1))
    assert np.linalg.norm(Y)>0


def test_malformed_native_complex_input_fails_closed():
    with pytest.raises(ValueError,match='Mismatched'):
        m.decode_complex({'real':[1,2],'imag':[1]})
    with pytest.raises(ValueError,match='Nonfinite'):
        m.decode_complex(np.asarray([float('nan')+1j]))
