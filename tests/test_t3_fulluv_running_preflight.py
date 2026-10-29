from scripts.ProbeT3FullUVRunning import candidate_vectors, finite_or_none, BLOCKED_LEGACY, REUSED
import numpy as np

def test_numeric_sanitation():
    assert finite_or_none(float('inf')) is None
    assert finite_or_none('abc') is None
    assert finite_or_none(2) == 2.0

def test_special_categories():
    assert 'A_-4' in BLOCKED_LEGACY
    assert 'C_-1' in REUSED
    assert not BLOCKED_LEGACY.intersection(REUSED)

def test_seed_always_includes_historical():
    hist={'seed':np.arange(18)}
    seeds=list(candidate_vectors('A_+0','',hist))
    assert len(seeds)==1
    assert seeds[0][0]=='historical_fulluv'
    assert np.array_equal(seeds[0][1],hist['seed'])
