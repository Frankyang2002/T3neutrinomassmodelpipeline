from Numerical.diagnostics.AuditT3ChargeChannels import audit, charge_channels


def case(f, s1, s2):
    return {"model": {"dS1": 2, "dS2": 2, "dF": 1, "alpha": -1},
            "fields": {k: {"weights_and_charges": [{"Q": str(q)} for q in vals]}
                       for k, vals in {"F": f, "S1": s1, "S2": s2}.items()}}


def test_single_intersection():
    row = charge_channels(case([0], [0, -1], [0, 1]))
    assert row["internal_charge_channels"] == ["0"]
    assert row["charge_allowed_channel_count"] == 1
    assert row["verified_loop_multiplicity"] is None


def test_two_intersections_preserve_no_compensation():
    row = charge_channels(case([0, -1], [0, -1, -2], [1, 0, -1]))
    assert row["internal_charge_channels"] == ["-1", "0"]
    assert row["charge_allowed_channel_count"] == 2
    assert row["compensation_factor"] is None


def test_empty_intersection():
    assert audit({"models": [case([1], [2], [3])]})["models"][0]["charge_allowed_channel_count"] == 0
