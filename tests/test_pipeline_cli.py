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


def test_verified_fermion_first_plan_is_accepted() -> None:
    parser = build_argument_parser()
    args = parser.parse_args(
        [
            "--dims", "2", "2", "1",
            "--threshold", "F",
            "--threshold", "S1", "S2",
        ]
    )
    shared = resolve_model_mode(parser, args)
    plan = resolve_pipeline_plan(parser, args, shared_scalar_mode=shared)

    assert plan.threshold_plan == (("F",), ("S1", "S2"))
    assert plan.threshold_scales == ("MF", "MS")
    assert plan.is_verified_fermion_first_hierarchy is True
    assert plan.is_supported_production_order is True


def test_scalar_first_threshold_plan_is_rejected_for_production() -> None:
    parser = build_argument_parser()
    args = parser.parse_args(
        [
            "--dims", "2", "2", "1",
            "--threshold", "S1",
            "--threshold", "F", "S2",
        ]
    )
    shared = resolve_model_mode(parser, args)

    with pytest.raises(SystemExit) as exc_info:
        resolve_pipeline_plan(parser, args, shared_scalar_mode=shared)

    assert exc_info.value.code == 2


def test_partially_split_scalar_hierarchy_is_rejected_for_production() -> None:
    parser = build_argument_parser()
    args = parser.parse_args(
        [
            "--dims", "2", "2", "1",
            "--threshold", "F",
            "--threshold", "S1",
            "--threshold", "S2",
        ]
    )
    shared = resolve_model_mode(parser, args)

    with pytest.raises(SystemExit) as exc_info:
        resolve_pipeline_plan(parser, args, shared_scalar_mode=shared)

    assert exc_info.value.code == 2


def test_retired_scalar_first_cli_options_are_absent() -> None:
    parser = build_argument_parser()
    assert "--eft-max-dimension" not in parser._option_string_actions
    assert "--allow-truncated-scalar-first" not in parser._option_string_actions


def test_study_name_preserves_nested_full_study_path() -> None:
    parser = build_argument_parser()
    args = parser.parse_args(["--smoke", "--study", "full/my study"])

    assert study_name(args) == "full/my_study"


def test_existing_primary_cli_modes_remain_mutually_exclusive() -> None:
    parser = build_argument_parser()

    with pytest.raises(SystemExit) as exc_info:
        parser.parse_args(["--smoke", "--dimension-comparison"])

    assert exc_info.value.code == 2
