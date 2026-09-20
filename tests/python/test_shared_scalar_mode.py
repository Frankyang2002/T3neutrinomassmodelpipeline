from __future__ import annotations

from pathlib import Path
import sys
import tempfile

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from common.T3Model import (
    shared_scalar_formal_dimensions,
    valid_shared_scalar_dimensions,
)
from common.Thresholds import (
    build_eft_stage_records,
    threshold_plan_for_wolfram,
    validate_threshold_plan,
)


def test_shared_scalar_dimensions():
    assert valid_shared_scalar_dimensions(2, 1)
    assert valid_shared_scalar_dimensions(2, 3)
    assert not valid_shared_scalar_dimensions(2, 2)
    assert shared_scalar_formal_dimensions(2, 1) == (2, 2, 1)


def test_shared_threshold_plan_is_physical_but_expands_for_matchete(tmp_path: Path):
    plan = validate_threshold_plan([["F"], ["S"]], shared_scalar=True)
    assert plan == (("F",), ("S",))
    assert threshold_plan_for_wolfram(plan, shared_scalar=True) == (
        ("F",),
        ("S1", "S2"),
    )

    stages = build_eft_stage_records(
        plan,
        tmp_path,
        shared_scalar=True,
    )
    assert stages[0].active_heavy_fields == ("F", "S")
    assert stages[1].label == "EFT_1_after_F"
    assert stages[1].active_heavy_fields == ("S",)
    assert stages[2].label == "EFT_2_after_S"
    assert stages[2].active_heavy_fields == ()


def test_shared_rge_model_counts_one_scalar_block():
    import sympy as sp
    from RGE.general.ScalarBasis import ScalarBasis

    model = ScalarBasis.t3_shared(2, sp.Rational(1, 2), include_higgs=True)
    assert tuple(model.blocks) == ("H", "S")
    assert model.total_real_scalar_dimension == 8
    assert tuple(model.block("S").indices) == (5, 6, 7, 8)


def main() -> None:
    print("=" * 72)
    print("SHARED-SCALAR MODE REGRESSION")
    print("=" * 72)

    test_shared_scalar_dimensions()
    print("PASS: shared-scalar dimensions")

    with tempfile.TemporaryDirectory() as directory:
        test_shared_threshold_plan_is_physical_but_expands_for_matchete(
            Path(directory)
        )
    print("PASS: shared-scalar threshold plan")

    test_shared_rge_model_counts_one_scalar_block()
    print("PASS: one physical scalar block")

    print("PASS: shared-scalar mode")


if __name__ == "__main__":
    main()
