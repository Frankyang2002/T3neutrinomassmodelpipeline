from scripts.RunT3FrozenMassNormalizedRemaining import accepted_chi2, plan, PARITY_BLOCKED

def test_acceptance_is_numerical():
    assert accepted_chi2({'chi2': 1e-18}, 5.86) == 1e-18
    assert accepted_chi2({'chi2': 6},5.86) is None
    assert accepted_chi2({'chi2': float('nan')},5.86) is None

def test_skip_accepted_and_blocked():
    rows=[{'model':'A_+0','study':'small_uv_frozen_5e-9','status':'CandidateConstrainedPoorFit', 'chi2':600},
          {'model':'A_-4','study':'small_uv_frozen_5e-9','status':'FailedOrUnavailable','chi2':None},
          {'model':'C_+1','study':'small_uv_frozen_5e-9','status':'CandidateConstrainedAccepted','chi2':1}]
    pending,skipped=plan(rows,['small_uv_frozen_5e-9'],['A_+0','A_-4','C_+1'],5.86,{('A_+0','small_uv_frozen_5e-9'):{'source':'previous','chi2':1e-18}})
    assert pending == []
    assert len(skipped) == 3
    assert 'A_-4' in PARITY_BLOCKED

def test_unresolved_selected():
    rows=[{'model':'D_+0','study':'large_uv_frozen_0p5','status':'CandidateConstrainedPoorFit','chi2':200}]
    pending,_=plan(rows,['large_uv_frozen_0p5'],['D_+0'],5.86,{})
    assert pending == [('D_+0','large_uv_frozen_0p5')]
