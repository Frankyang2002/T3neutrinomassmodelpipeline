"""Offline interface and algebra checks; actual Matchete test is on Windows."""
from __future__ import annotations

from pathlib import Path
import json
import math
import re
from types import SimpleNamespace

import numpy as np
import pytest

from scripts import RunT3GeneralMatcheteMatching as matching

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / 'Lagrangian' / 'RunT3GeneralSelfMassBasis.wl'
SELF = {'A_-2', 'A_+0', 'D_-2', 'D_+0', 'E_-2', 'E_+0'}


def test_dispatch_six_self_and_eight_other_plus_two_mixed():
    assert len(matching.ALL) == 16
    assert matching.SELF == SELF
    assert set(matching.ALL) - matching.SELF - matching.MIXED == {
        'A_-4', 'B_-3', 'B_+1', 'C_-3', 'C_+1', 'D_+2', 'E_-4', 'E_+2'
    }
    for model in SELF:
        runner, name = matching.source_spec(model)
        assert runner.name == 'RunT3GeneralSelfMassBasis.wl'
        assert 'SelfBilinear' in name
    for model in matching.MIXED:
        with pytest.raises(ValueError, match='mixed-scalar'):
            matching.source_spec(model)


def test_real_bilinear_diagonalization_and_nonzero_splitting():
    # Source Lagrangian convention: Lmass = -m^2 |S|^2 + mu*S*^2 + mu*S^2.
    # With S=(R+iI)/sqrt(2), the Hessian is diag(m^2-2mu, m^2+2mu).
    msq, mu = 2.25e6, 10000.
    hessian = np.diag([msq-2*mu, msq+2*mu])
    eig = np.linalg.eigvalsh(hessian)
    np.testing.assert_allclose(eig, [2230000., 2270000.])
    assert eig[1] - eig[0] == pytest.approx(4*mu)
    # Complex-to-real field map is unitary: S=(R+iI)/sqrt(2).
    U = np.array([[1, 1j], [1, -1j]]) / math.sqrt(2)
    np.testing.assert_allclose(U.conj().T@U, np.eye(2), atol=1e-14)


def test_wolfram_source_retains_yukawas_and_quartics_and_checks_matching():
    s = SCRIPT.read_text(encoding='utf-8')
    assert 'BuildT3Lagrangian[model]' in s
    assert 'AdditionalYukawasIncluded' in s
    assert 'AdditionalQuarticsIncluded' in s
    assert 'RemoveField[NewScalar2]' in s
    assert 'RemoveField[NewScalar1]' in s
    assert 'SelfConjugate -> True' in s
    assert 'Mass -> {Heavy, mass}' in s
    assert 'FreeLag[NewFermion, rest, a, b]' in s
    assert 'CheckLagrangian[newUV]' in s
    assert 'RunT3Matching[newUV,5,1]' in s
    assert 'ExtractWeinbergCoefficient[matched]' in s
    assert 'self_bilinear_nonzero_by_mass_splitting' in s
    assert 'complex_bilinear_phase_supported' in s
    assert 'c5_coefficient.txt' in s
    assert s.count('[') == s.count(']')


def test_matching_process_uses_new_wolfram_runner_with_all_flags(tmp_path, monkeypatch):
    monkeypatch.setattr(matching, 'ROOT', tmp_path)
    dest = tmp_path/'matching'
    def fake_run(command, **kwargs):
        assert 'RunT3GeneralSelfMassBasis.wl' in command[2]
        assert command[3:] == [str(dest/'A_m2'), 'A', 'm2']
        assert all(kwargs['env'][name] == '1' for name in (
            'T3_EXPERIMENTAL_GENERAL_FERMIONS',
            'T3_EXPERIMENTAL_GENERAL_QUARTICS',
            'T3_EXPERIMENTAL_GENERAL_QUADRATICS'))
        folder = dest/'A_m2'
        (folder/'c5_coefficient.txt').write_text('Coupling[y1] Coupling[y2] Coupling[y2Conj] T32PlusMass T32MinusMass')
        (folder/'self_mass_basis_metadata.json').write_text(json.dumps({
            'model':'A_-2', 'self_bilinear_nonzero_by_mass_splitting': True,
            'complex_bilinear_phase_supported': False}))
        return SimpleNamespace(returncode=0,
            stdout='GENERAL SELF-BILINEAR MASS-BASIS C5: SUCCESS', stderr='')
    monkeypatch.setattr(matching.subprocess, 'run', fake_run)
    status = matching.run_model('A_-2', dest, 'wolframscript', reuse=False)
    assert status['status'] == 'Matched'
    assert status['self_bilinear_in_mass_eigenvalues'] is True
    assert status['matching_scheme'] == matching.SCHEME_SELF


def test_mixed_is_explicitly_blocked(tmp_path):
    result = matching.run_model('B_-1', tmp_path, 'wolframscript', reuse=False)
    assert result['status'] == 'BlockedMixedScalarMassBasis'
    assert not (tmp_path/'B_m1'/'c5_coefficient.txt').exists()
