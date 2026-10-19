"""Isolated tests for Yukawa grouping, without running the physics pipeline."""
import numpy as np

from Numerical.fitting.YukawaSimilarity import group_collection, make_candidate, pair_distance


def _case(name, x):
    matrix = np.eye(3, dtype=complex) * x
    return {"model": name, "sha256": name, "pair": (matrix, matrix),
            "payload": {"representation": {"alpha": len(name)},
                        "base_state": {"ordinary": {
                            "y1": {"real": (matrix.real).tolist(), "imag": (matrix.imag).tolist()},
                            "y2": {"real": (matrix.real).tolist(), "imag": (matrix.imag).tolist()},
                            "lambdaT3": {"real": 0.1, "imag": 0.0}}}}}


def test_distance_and_clustering():
    a, b, c = [_case(*p) for p in (("a", 1.0), ("b", 1.01), ("c", 4.0))]
    assert pair_distance(a["pair"], a["pair"]) == 0
    assert np.isclose(pair_distance(a["pair"], b["pair"]), pair_distance(b["pair"], a["pair"]))
    result = group_collection([c, a, b], threshold=0.1)
    assert [g["size"] for g in result["groups"]] == [2, 1]
    assert result["groups"][0]["reference_model"] == "a"


def test_candidate_only_changes_yukawas():
    a, b = _case("a", 1.0), _case("b", 2.0)
    candidate = make_candidate(a, b)
    assert candidate["base_state"]["ordinary"]["y1"] == b["payload"]["base_state"]["ordinary"]["y1"]
    assert candidate["base_state"]["ordinary"]["lambdaT3"] == a["payload"]["base_state"]["ordinary"]["lambdaT3"]
    assert candidate["representation"] == a["payload"]["representation"]
    assert candidate["yukawa_grouping"]["status"] == "UNVALIDATED_REFERENCE_REPLACEMENT"
