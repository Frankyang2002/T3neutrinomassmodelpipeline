"""Architecture and boundary regressions for the scalar-first d=6 foundation."""

from __future__ import annotations

import argparse
from pathlib import Path

from common.EFT import EFTTruncation
from common.PipelineCLI import (
    build_argument_parser,
    resolve_model_mode,
    resolve_pipeline_plan,
)
from common.PipelinePlan import PipelinePlan
from common.RunRecords import RunRecord
from RGE.running.intermediate.ScalarFirstDimensionSixSeed import (
    extract_scalar_first_dimension_six_seed,
)
from studies.FullT3Study import _forward_common_cli_args


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def test_scalar_first_d6_cli_does_not_need_truncation_override() -> None:
    parser = build_argument_parser()
    args = parser.parse_args(
        [
            "--dims", "2", "2", "1",
            "--threshold", "S1",
            "--threshold", "F",
            "--threshold", "S2",
            "--eft-max-dimension", "6",
        ]
    )
    shared = resolve_model_mode(parser, args)
    plan = resolve_pipeline_plan(parser, args, shared_scalar_mode=shared)

    assert plan.truncation == EFTTruncation(max_operator_dimension=6)
    assert plan.has_scalar_first_threshold is True
    assert plan.scalar_first_truncates_leading_path is False


def test_full_study_forwards_dimension_six_matching_order() -> None:
    args = argparse.Namespace(
        debug_reports=False,
        force=False,
        numerical=None,
        threshold=None,
        threshold_scale=None,
        allow_truncated_scalar_first=False,
        eft_max_dimension=6,
    )
    assert _forward_common_cli_args(args) == ["--eft-max-dimension", "6"]


def test_matching_boundary_threads_pipeline_truncation_to_runner() -> None:
    source = (
        PROJECT_ROOT / "Lagrangian" / "T3ModelMatching.py"
    ).read_text(encoding="utf-8")
    assert "eft_order = pipeline_plan.truncation.max_operator_dimension" in source
    assert source.count("eft_order=eft_order") == 3

    runner = (PROJECT_ROOT / "Lagrangian" / "Runner.py").read_text(
        encoding="utf-8"
    )
    assert "str(eft_order)" in runner
    assert "Current T3 matching supports EFT order 5 or 6" in runner


def test_exact_matchete_scalar_first_candidate_is_exported_verbatim(
    tmp_path: Path,
) -> None:
    plan = PipelinePlan.from_threshold_configuration(
        [["S1"], ["F"], ["S2"]],
        ["MS1", "MF", "MS2"],
        truncation=EFTTruncation(max_operator_dimension=6),
    )
    interval = plan.running_intervals[0]

    record = RunRecord(
        name="T3-B",
        alpha=-1,
        d_s1=2,
        d_s2=2,
        d_f=1,
        return_code=0,
        summary={},
        output_dir=tmp_path,
        shared_scalar=False,
    )
    record.eft_stages = plan.build_stage_records(tmp_path)
    stage = record.stage_with_active_fields("F", "S2")

    data_dir = tmp_path / "data"
    data_dir.mkdir()
    tree_path = data_dir / "eft1_after_S1_tree_inputform.txt"
    candidate = (
        "kappa*Field[l, Fermion, {Index[i, SU2L[fund]]}, {}]*"
        "Field[NewFermion, Fermion, {Index[r, NFlavor]}, {}]*"
        "Field[H, Scalar, {Index[a, SU2L[fund]]}, {}]*"
        "Field[H, Scalar, {Index[b, SU2L[fund]]}, {}]*"
        "Field[NewScalar2, Scalar, {Index[c, SU2L[fund]]}, {}]"
    )
    spectator = (
        "lambdaH*Field[H, Scalar, {Index[a, SU2L[fund]]}, {}]*"
        "Field[H, Scalar, {Index[b, SU2L[fund]]}, {}]"
    )
    tree_path.write_text(candidate + " + " + spectator, encoding="utf-8")

    payload = extract_scalar_first_dimension_six_seed(record, stage, interval)

    assert payload["status"] == "Success"
    assert payload["operator_dimension"] == 6
    assert payload["candidate_term_count"] == 1
    assert payload["candidates"][0]["TermInputForm"] == candidate
    assert payload["normalization_status"] == "RawMatcheteTermsOnly"
    assert payload["tensor_adapter_status"] == "NotImplemented"
    assert (
        data_dir / "eft1_after_S1_psi2phi3_seed.json"
    ).is_file()
