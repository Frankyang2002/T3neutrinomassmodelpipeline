from __future__ import annotations

import numpy as np
import pytest

from Numerical.FinalC5ContributionDiagnostics import (
    FinalC5ContributionBreakdown,
)


def test_contribution_breakdown_requires_hard_plus_direct_equal_combined() -> None:
    hard = np.eye(3, dtype=complex)
    direct = 0.1 * np.eye(3, dtype=complex)

    result = FinalC5ContributionBreakdown(
        hard=hard,
        direct_running=direct,
        combined=hard + direct,
    ).validated()

    payload = result.as_json_dict()

    assert np.isclose(
        payload["direct_to_hard_norm_ratio"],
        0.1,
    )


def test_contribution_breakdown_rejects_inconsistent_combined() -> None:
    hard = np.eye(3, dtype=complex)
    direct = 0.1 * np.eye(3, dtype=complex)

    with pytest.raises(ValueError, match=r"hard \+ direct-running"):
        FinalC5ContributionBreakdown(
            hard=hard,
            direct_running=direct,
            combined=hard,
        ).validated()
