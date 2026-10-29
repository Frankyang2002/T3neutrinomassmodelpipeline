"""Focused checks for frozen scenario selection and candidate retention."""
from scripts import ContinueT3GeneralScalarSobol as mod


def test_only_requested_frozen_cases():
    rows = {(m,c): {'model':m,'study':c,'status':'CandidateConstrainedPoorFit','chi2':100.}
            for m in ('A_-2',) for c in mod.CASES}
    result=mod.selected_cells(rows,('A_-2',),8,{},('small_uv_frozen_5e-9',))
    assert result==[('A_-2','small_uv_frozen_5e-9')]


def test_existing_accepted_not_searched():
    rows={('A_-2','large_uv_frozen_0p5'):{'status':'CandidateConstrainedAccepted','chi2':0.01}}
    assert not mod.selected_cells(rows,('A_-2',),8,{},('large_uv_frozen_0p5',))


def test_retains_best_poor_fit():
    old={'status':'CandidateConstrainedPoorFit','chi2':100.0}
    new={'status':'CandidateConstrainedPoorFit','chi2':200.0}
    assert mod.choose(old,new) is old
    assert mod.choose(new,old) is old


def test_acceptance_prevents_overwrite():
    good={'status':'CandidateConstrainedAccepted','chi2':1.0}
    alt={'status':'CandidateConstrainedAccepted','chi2':0.1}
    assert mod.choose(good,alt) is good
