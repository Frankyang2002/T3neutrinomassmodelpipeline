"""Pure numerical regressions and fail-closed contract for the threshold preparer."""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

from Numerical.experimental.GeneralT3MixedScalarMassBasis import MixedDoubletMassBasis
from scripts.RunT3BCSeparatedScalarThresholds import (
    PHASE_Y2, STATUS, FORMAT, exact_scale_crossing,
    yukawa_weights, boundary, validate_inputs, check_native_mixing,
)


def constant_basis(model: str = 'B_-1', *, mX: float = 1200., mY: float = 1500.,
                   angle: float = 0.28) -> MixedDoubletMassBasis:
    """Synthesize a nontrivially mixed, positive Hermitian scalar mass matrix."""
    c, s = math.cos(angle), math.sin(angle)
    U = np.array([[c, -s], [s, c]], dtype=complex)
    mass = U @ np.diag([mX*mX, mY*mY]) @ U.conj().T
    return MixedDoubletMassBasis.build(model, mass[0,0].real, mass[1,1].real, mass[0,1])


@pytest.mark.parametrize('which,expected', [(0,1200.), (1,1500.)])
def test_crossing_exact_constant_masses(which, expected):
    fn = lambda t: constant_basis()
    t = exact_scale_crossing(fn, math.log(1000), math.log(1e5), which)
    assert math.exp(t) == pytest.approx(expected, abs=1e-9)


def test_mass_crossing_fails_when_unbracketed():
    with pytest.raises(ValueError, match='no sign-changing'):
        exact_scale_crossing(lambda t: constant_basis(),
                             math.log(2000), math.log(1e5), 0)


def test_crossing_detects_tachyon():
    tachyon = MixedDoubletMassBasis.build('B_-1', -1., 2., 0.1)
    with pytest.raises(ValueError, match='tachyonic'):
        exact_scale_crossing(lambda t: tachyon,
                             math.log(1000), math.log(1e5), 1)


@pytest.mark.parametrize('model,alpha', [('B_-1',-1), ('C_-1',1)])
def test_yukawa_weights_native_phase(model, alpha):
    basis = constant_basis(model)
    info = yukawa_weights(model, basis.U, 1, alpha)
    w1 = np.array([complex(z['real'],z['imag']) for z in info['y1_CConj_fermion_S1_source']])
    w2 = np.array([complex(z['real'],z['imag']) for z in info['y2_unconjugated_fermion_BarS2_source']])
    assert np.allclose(w1, basis.U[0])
    assert np.allclose(w2, -alpha*basis.U[1])
    assert abs(np.vdot(w1,w2)) < 1e-12


@pytest.mark.parametrize('model,wrong', [('B_-1',1), ('C_-1',-1)])
def test_wrong_y2_native_sign_rejected(model, wrong):
    with pytest.raises(ValueError, match='wrong native'):
        yukawa_weights(model, constant_basis(model).U, 1, wrong)


def test_yukawa_weight_non_unitary_rejected():
    with pytest.raises(ValueError, match='mixing matrix invalid'):
        yukawa_weights('B_-1', 2 * np.eye(2), 1, -1)


def test_threshold_projector_is_rank_one():
    basis = constant_basis()
    entry = boundary('B_-1',math.log(1500),basis,1,-1,'Y')
    p = np.array(entry['projector_in_phi_multiplet']['real']) + \
        1j*np.array(entry['projector_in_phi_multiplet']['imag'])
    assert entry['mu_gev'] == pytest.approx(1500)
    assert np.allclose(p@p,p)
    assert np.linalg.matrix_rank(p) == 1


def test_threshold_off_shell_rejected():
    with pytest.raises(ValueError, match='off-shell'):
        boundary('B_-1',math.log(1800),constant_basis(),1,-1,'Y')


def test_unsupported_model_rejected():
    with pytest.raises(ValueError, match='wrong native'):
        yukawa_weights('A_+0',np.eye(2),1,-1)


def valid_payloads():
    from scripts import RunGeneralT3QuadraticEFT1 as qeft
    q = {'format':qeft.FORMAT,'status':'CandidateQuadraticAndFullQuarticEFT1Integrated',
         'model':'B_-1','quadratic_kind':'mixed', 'mu2_beta_included':True,
         'new_quartics_evolved':1,'reused_quartics_evolved':6,
         'physical_matchete_rgbeta_mapping_certified':False,
         'scalar_mass_threshold_matching_performed':False,
         'wilson_coefficient_matching_performed':False}
    ref={'status':'ConditionalEFT1ScalarEndpointOneScaleC5NotPhysicalSequential',
         'model':'B_-1','endpoint':'final','running_sha256':'abc',
         'scalar_threshold_matching_performed':False}
    return q,ref


def test_valid_inputs():
    q,ref=valid_payloads()
    validate_inputs('B_-1',q,ref)


@pytest.mark.parametrize('flag', [
    'physical_matchete_rgbeta_mapping_certified',
    'scalar_mass_threshold_matching_performed',
    'wilson_coefficient_matching_performed'])
def test_falsely_certified_input_rejected(flag):
    q,ref=valid_payloads();q[flag]=True
    with pytest.raises(ValueError,match='incompatible'):
        validate_inputs('B_-1',q,ref)


def test_sequential_claims_deliberately_false_in_source():
    src = Path(__file__).resolve().parents[1] / 'scripts/RunT3BCSeparatedScalarThresholds.py'
    text = src.read_text(encoding='utf-8')
    for key in ('Y_threshold_Wilson_matching_performed',
                'X_only_EFT2_beta_functions_available_here',
                'X_threshold_Wilson_matching_performed',
                'finite_scalar_threshold_matching_performed',
                'physical_sequential_C5_certified'):
        assert f"'{key}': False" in text
    assert FORMAT.endswith('-v1') and 'NotPerformed' in STATUS


def test_complex_mixing_phases_are_kept():
    matrix = np.array([[1200.**2, 4e4 * np.exp(0.77j)],
                       [4e4 * np.exp(-0.77j), 1500.**2]], dtype=complex)
    basis = MixedDoubletMassBasis.build('C_-1', matrix[0,0].real,
                                       matrix[1,1].real, matrix[0,1])
    result = yukawa_weights('C_-1', basis.U, 1, 1)
    w2 = np.array([complex(v['real'],v['imag'])
                   for v in result['y2_unconjugated_fermion_BarS2_source']])
    assert np.allclose(w2, -basis.U[1])
    assert np.max(np.abs(w2.imag)) > 1e-3


def test_beta_source_registry_hash_checks_remain_fail_closed():
    source = (Path(__file__).resolve().parents[1] /
              'scripts/RunT3BCSeparatedScalarThresholds.py').read_text(encoding='utf-8')
    assert 'saved qeft beta source differs' in source
    assert 'quadratic EFT1 registry/manifest source changed' in source
    assert 'one-scale replay used a different quadratic EFT1 source' in source
    assert 'X_only_EFT2_beta_functions_available_here' in source
