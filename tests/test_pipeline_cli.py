"""Regression tests for the pipeline command-line configuration boundary."""

from __future__ import annotations

import pytest

from common.PipelineCLI import (
    build_argument_parser,
    resolve_model_mode,
    resolve_pipeline_plan,
    study_name,
)


def test_three_field_mode_keeps_alpha_zero_default() -> None:
    parser = build_argument_parser()
    args = parser.parse_args(["--dims", "2", "2", "1"])

    shared = resolve_model_mode(parser, args)

    assert shared is False
    assert args.alpha == 0
    assert study_name(args) == "single"


def test_shared_scalar_mode_fixes_alpha_minus_one() -> None:
    parser = build_argument_parser()
    args = parser.parse_args(["--dims", "2", "1"])

    shared = resolve_model_mode(parser, args)

    assert shared is True
    assert args.alpha == -1


def test_shared_scalar_mode_rejects_other_alpha() -> None:
    parser = build_argument_parser()
    args = parser.parse_args(["--dims", "2", "1", "--alpha", "0"])

    with pytest.raises(SystemExit) as exc_info:
        resolve_model_mode(parser, args)

    assert exc_info.value.code == 2


def test_scalar_first_threshold_plan_requires_explicit_truncation_opt_in() -> None:
    parser = build_argument_parser()
    args = parser.parse_args(
        [
            "--dims", "2", "2", "1",
            "--threshold", "S1",
            "--threshold", "F",
            "--threshold", "S2",
        ]
    )
    shared = resolve_model_mode(parser, args)

    with pytest.raises(SystemExit) as exc_info:
        resolve_pipeline_plan(parser, args, shared_scalar_mode=shared)

    assert exc_info.value.code == 2


def test_scalar_first_threshold_plan_is_preserved_after_explicit_opt_in() -> None:
    parser = build_argument_parser()
    args = parser.parse_args(
        [
            "--dims", "2", "2", "1",
            "--threshold", "S1",
            "--threshold", "F",
            "--threshold", "S2",
            "--allow-truncated-scalar-first",
        ]
    )
    shared = resolve_model_mode(parser, args)
    plan = resolve_pipeline_plan(parser, args, shared_scalar_mode=shared)

    assert plan.threshold_plan == (("S1",), ("F",), ("S2",))
    assert plan.threshold_scales == ("MS1", "MF", "MS2")
    assert plan.scalar_first_truncates_leading_path is True
    assert [interval.active_heavy_fields for interval in plan.running_intervals] == [
        frozenset({"F", "S2"}),
        frozenset({"S2"}),
    ]


def test_study_name_preserves_nested_full_study_path() -> None:
    parser = build_argument_parser()
    args = parser.parse_args(["--smoke", "--study", "full/my study"])

    assert study_name(args) == "full/my_study"


def test_existing_primary_cli_modes_remain_mutually_exclusive() -> None:
    parser = build_argument_parser()

    with pytest.raises(SystemExit) as exc_info:
        parser.parse_args(["--smoke", "--dimension-comparison"])

    assert exc_info.value.code == 2
