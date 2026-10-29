from scripts.ContinueT3GeneralScalarSobol import acceptable, choose, selected_cells


def row(status,chi2):
    return {'status':status,'chi2':chi2}


def test_accepted_never_overwritten():
    accepted=row('CandidateConstrainedAccepted',1.0)
    assert choose(accepted,row('CandidateConstrainedAccepted',0.1)) is accepted


def test_new_acceptable_beats_poor():
    good=row('CandidateConstrainedAccepted',4.0)
    assert choose(row('CandidateConstrainedPoorFit',100),good) is good


def test_reject_no_mass_safe_and_nan():
    assert not acceptable(row('NoMassSafeNeutrinoFit',None))
    assert not acceptable(row('CandidateConstrainedAccepted',float('nan')))


def test_better_poor_fit_preserved():
    poor=row('CandidateConstrainedPoorFit',500)
    better=row('CandidateConstrainedPoorFit',200)
    assert choose(poor,better) is better


def test_respects_round_caps():
    rows={('A_-2','low'):row('CandidateConstrainedPoorFit',30),('A_-2','high'):row('CandidateConstrainedAccepted',1)}
    import scripts.ContinueT3GeneralScalarSobol as s
    original=s.CASES
    s.CASES=('low','high')
    try:
        assert selected_cells(rows,['A_-2'],3,{('A_-2','low'):2})==[('A_-2','low')]
        assert selected_cells(rows,['A_-2'],3,{('A_-2','low'):3})==[]
    finally:
        s.CASES=original
