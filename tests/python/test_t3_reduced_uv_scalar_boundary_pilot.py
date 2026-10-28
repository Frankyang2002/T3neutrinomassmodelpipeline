"""Offline numerical tests; do not require Wolfram or saved T3 running files."""
import importlib.util
from pathlib import Path
import numpy as np
import pytest

SRC = Path(__file__).resolve().parents[2] / 'scripts/ProbeT3ReducedUVScalarBoundaries.py'
spec = importlib.util.spec_from_file_location('t3_uv_mass_pilot', SRC)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def sample_map():
    a = np.array([[1.08, 0.02], [0.01, 0.95]])
    b = np.array([1600., 1400.])
    return a, b


def test_affine_reconstruction_and_checks():
    a, b = sample_map()
    x0 = np.array([1., 1.2])
    step = 4.
    f = lambda x: a@x+b
    aa, bb = m.build_affine_map(x0, f(x0), f(x0 + [step, 0]),
                                 f(x0 + [0, step]), step)
    assert np.allclose(aa, a, atol=1e-12)
    assert np.allclose(bb, b, atol=1e-12)
    assert m.check_affinity(aa, bb, [(np.zeros(2), f(np.zeros(2))),
                                   (x0+[step, step], f(x0+[step, step]))]) < 1e-9


def test_fails_closed_non_affine():
    a, b = sample_map()
    with pytest.raises(ValueError, match='not numerically validated'):
        m.check_affinity(a,b,[(np.zeros(2),b+[.2,0])])


def test_signed_solution_for_heavy_positive_offset():
    a, b = sample_map()
    ans = m.signed_target_solution(a,b,target_tev=2.)
    assert ans['status'] == 'SolvedAlgebraically'
    assert not ans['positive_uv_mass_squared']
    x = np.array(ans['uv_mass_sq_tev2'])
    assert np.allclose(a@x+b,[4.,4.],atol=1e-9)


def test_positive_uv_infeasible_for_heavy_offset():
    a,b = sample_map()
    assert m.find_mass_window(a,b)['feasible'] is False
    assert m.find_mass_window(a,b,.5,5.)['feasible'] is False


def test_positive_uv_feasible_when_offset_small():
    a = np.eye(2)
    b = np.array([1., 3.])
    candidate = m.find_mass_window(a,b)
    assert candidate['feasible']
    x = np.array(candidate['candidate_uv_mass_sq_tev2'])
    y = a@x+b
    assert np.all(y >= .25) and np.all(y <= 25.)


def test_light_uv_candidate_with_constraints():
    candidate = m.find_mass_window(np.eye(2), np.zeros(2), .5, 5.)
    assert candidate['feasible']
    x = np.array(candidate['candidate_uv_mass_sq_tev2'])
    assert np.all(x >= .25-1e-9) and np.all(x <= 25.+1e-9)


def test_coupled_feasibility_not_independent_bounds():
    # Both scalar targets cannot be simultaneously achieved by positive UV masses.
    a = np.array([[1.,1.],[-1.,-1.]])
    b = np.array([100.,100.])
    assert not m.find_mass_window(a,b)['feasible']


def test_reject_bad_bounds():
    with pytest.raises(ValueError):
        m.find_mass_window(np.eye(2), np.zeros(2), 5., .5)
    with pytest.raises(ValueError):
        m.build_affine_map([0,0],[0,0],[1,0],[0,1],0.)
    with pytest.raises(ValueError):
        m.signed_target_solution(np.eye(2),np.zeros(2),10.)


def test_rank_deficient_map_does_not_claim_exact_root():
    a = np.ones((2,2))
    ans = m.signed_target_solution(a, [0,0], 2.)
    assert ans['status'] == 'IllConditionedMap'
