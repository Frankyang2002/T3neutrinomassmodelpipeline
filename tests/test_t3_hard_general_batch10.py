"""Numerical and source-inventory checks independent of old Windows outputs."""

import numpy as np
import pytest

from Numerical.parallel_t3.complete.two_model_hard_general_batch10 import (
    EXPECTED_D, FORMAT, extra_hhss_from_source, matrix, matrix_json,
    matrix_metrics, reference_hhss_tensor, symmetric, _reference_number,
)


def test_exact_reference_bosonic_orbits_and_phase():
    rows=[{'H1_H2_S1_S2_0based':[0,1,0,7], 'value_exact':'sqrt(2)/2'}]
    re=reference_hhss_tensor(rows, 3+4j, kind='Re')
    im=reference_hhss_tensor(rows, 3+4j, kind='Im')
    for v in (re, im):
        assert v.shape == (4,4,8,8)
        np.testing.assert_allclose(v, v.swapaxes(0,1))
        np.testing.assert_allclose(v, v.swapaxes(2,3))
    assert np.isclose(re[0,1,0,7],3/np.sqrt(2))
    assert np.isclose(im[1,0,7,0],4/np.sqrt(2))
    assert np.count_nonzero(re)==4
    assert _reference_number('-sqrt(3)/2') == pytest.approx(-np.sqrt(3)/2)
    with pytest.raises(ValueError):
        _reference_number('garbage+')


def test_exact_d2_hhss_both_extra_channels_not_silently_zero():
    model = {'models': {'D_-2': {'new_rgbeta_symbols': [
        {'counts':list(k), 'name':v} for k,v in EXPECTED_D.items()]}}}
    sectors=[]
    for counts,name in EXPECTED_D.items():
        ref={'id':'reference_invariant_1',
             'Re_c_unit_tensor_components':[
                {'H1_H2_S1_S2_0based':[0,0,1,6], 'value_exact':'2'}],
             'Im_c_unit_tensor_components':[]}
        sectors.append({'counts_H_BarH_S1_BarS1_S2_BarS2':list(counts),
                        'source_registry_symbols':[{'rgbeta_symbol':name,'source':'new'}],
                        'exact_reference_real_component_tensor_basis':[ref]})
    doc={'format':'t3-ws2-extra-hhss-exact-su2-reference-space-v1',
         'model':'D_-2','monomial_sectors':sectors}
    couplings={k:0.1 for k in EXPECTED_D.values()}
    out,meta=extra_hhss_from_source(doc,model,couplings)
    assert set(out)==set(couplings)
    for name,q in out.items():
        assert np.isclose(q[0,0,1,6],0.2)
        assert np.isclose(q[0,0,6,1],0.2)
        assert meta[name]['reference_coordinate_basis_is_canonical_SU2_not_certified_RGBeta']
    with pytest.raises(ValueError):
        extra_hhss_from_source(doc,model,{'lambdaT3GeneralQ020101I1':0.1})


def test_formal_poles_are_not_counterterms():
    w=np.array([[1j,0,0],[0,2j,0],[0,0,-1j]],complex)
    np.testing.assert_allclose(symmetric(w),2*w)
    d=matrix_metrics(2*w,w)
    assert d['relative_to_B']==pytest.approx(1)
    assert d['fit_is_not_applied'] is True
    assert FORMAT.endswith('batch10-v1')


def test_real_batch9_and_native_hard_json_matrix_encodings_agree():
    """Regression: old tests missed the actual hard-source native encoding."""
    # Byte-for-byte schema produced by the repository's encode_native(array).
    def encode_native(a):
        return [[{'real': float(z.real), 'imag': float(z.imag)} for z in row]
                for row in a]

    # An ordered, NOT Hermitian or symmetric, matrix: preserve chirality/phase.
    ordered = np.array([
        [1+2j, 3-4j, 5j],
        [-2+0.5j, 0j, 9-1j],
        [7+3j, 8j, -6-2j],
    ], dtype=complex)
    hard_source = encode_native(ordered)
    batch9_source = matrix_json(ordered)
    assert isinstance(hard_source, list) and isinstance(hard_source[0][0], dict)
    assert isinstance(batch9_source, dict)
    np.testing.assert_array_equal(matrix(hard_source, 'hard', symmetric=False), ordered)
    np.testing.assert_array_equal(matrix(batch9_source, 'Batch9', symmetric=False), ordered)
    np.testing.assert_array_equal(symmetric(matrix(hard_source, symmetric=False)),
                                  ordered + ordered.T)
    with pytest.raises(ValueError, match='not complex-symmetric'):
        matrix(hard_source)


def test_complex_json_matrix_decoding_fails_closed():
    # Byte-for-byte schema produced by the repository's encode_native(array).
    def encode_native(a):
        return [[{'real': float(z.real), 'imag': float(z.imag)} for z in row]
                for row in a]

    correct = np.eye(3, dtype=complex) * (2+3j)
    np.testing.assert_array_equal(matrix(encode_native(correct)), correct)
    np.testing.assert_array_equal(matrix(matrix_json(correct)), correct)

    for malformed in (
        [[], [], []],
        [[{'real': 1, 'imag': 0}] * 3] * 2,
        [[{'real': 1, 'imag': 0}] * 3] * 2 + [[{'real': 1}]*3],
        np.eye(3),
        {'real': np.eye(3).tolist(), 'imag': np.eye(2).tolist()},
        {'real': np.eye(3).tolist(), 'imag': (np.eye(3)*np.nan).tolist()},
    ):
        with pytest.raises((ValueError, TypeError)):
            matrix(malformed, symmetric=False)
