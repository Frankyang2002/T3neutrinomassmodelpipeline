from __future__ import annotations

import pytest

from common.T3Model import (
    T3_CLASSES,
    neutral_t3_class_points,
    t3_has_neutral_bsm_component,
    t3_neutral_component_fields,
    valid_t3_dimensions,
    valid_t3_topology_dimensions,
)


EXPECTED_NEUTRAL_POINTS = (
    ("A", -4),
    ("A", -2),
    ("A", 0),
    ("B", -3),
    ("B", -1),
    ("B", 1),
    ("C", -3),
    ("C", -1),
    ("C", 1),
    ("D", -2),
    ("D", 0),
    ("D", 2),
    ("E", -4),
    ("E", -2),
    ("E", 0),
    ("E", 2),
)


def test_alpha_minus4_to_2_scan_has_expected_16_neutral_points() -> None:
    points = neutral_t3_class_points(range(-4, 3))
    assert points == EXPECTED_NEUTRAL_POINTS
    assert len(points) == 16


@pytest.mark.parametrize("model_class,alpha", EXPECTED_NEUTRAL_POINTS)
def test_selected_scan_points_really_contain_a_neutral_bsm_state(
    model_class: str,
    alpha: int,
) -> None:
    dims = T3_CLASSES[model_class]
    assert t3_has_neutral_bsm_component(*dims, alpha)
    assert t3_neutral_component_fields(*dims, alpha)


def test_non_neutral_benchmark_point_is_rejected_by_neutrality_predicate() -> None:
    assert not t3_has_neutral_bsm_component(*T3_CLASSES["B"], 0)


def test_force_topology_validator_extends_dimension_scope_without_relaxing_t3() -> None:
    # (4, 4, 3) lies outside current d<=3 production support but satisfies
    # dS=dF+1 for both scalars and 4⊗4 contains the triplet.
    assert not valid_t3_dimensions(4, 4, 3)
    assert valid_t3_topology_dimensions(4, 4, 3)

    # Positive dimensions alone are not enough: this violates the T3 Yukawa
    # representation condition and must remain forbidden even in forced mode.
    assert not valid_t3_topology_dimensions(5, 5, 1)
