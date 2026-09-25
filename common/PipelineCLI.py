"""Command-line configuration for the T3 neutrino-mass pipeline.

This module owns only user-facing configuration: CLI options, shared-scalar
mode selection, threshold-plan construction, and study-directory naming.
It deliberately contains no matching, RGE, neutrino-mass, report, or model-
building calculations.  ``pipeline.py`` imports these helpers so the central
backbone can remain a short, readable description of the calculation order.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from common.EFT import EFTTruncation
from common.PipelinePlan import PipelinePlan


def build_argument_parser() -> argparse.ArgumentParser:
    """Build the command-line interface for the T3 pipeline."""

    parser = argparse.ArgumentParser(
        description=(
            "Run T3 matching for known benchmark models "
            "or arbitrary valid SU(2) irreps."
        )
    )
    mode = parser.add_mutually_exclusive_group()

    mode.add_argument(
        "--smoke",
        action="store_true",
        help="five T3 regression models",
    )
    mode.add_argument(
        "--hypercharge-comparison",
        action="store_true",
        help=(
            "scan alpha=-4,-3,...,2 for every T3-A...E SU(2) assignment; "
            "normal mode runs only points with at least one neutral BSM state"
        ),
    )
    mode.add_argument(
        "--dimension-comparison",
        action="store_true",
        help=(
            "compare T3-A...E at fixed alpha=0 to isolate SU(2) "
            "representation dependence; normal mode keeps only neutral-compatible points"
        ),
    )
    mode.add_argument(
        "--full",
        action="store_true",
        help=(
            "run smoke, hypercharge comparison, and dimension comparison "
            "sequentially under output/full and Reports/output/full"
        ),
    )
    mode.add_argument(
        "--dims",
        nargs="+",
        type=int,
        metavar="D",
        help=(
            "three numbers DS1 DS2 DF give the ordinary T3 model; "
            "two numbers DS DF give one physical shared scalar, e.g. "
            "--dims 2 1 for the scotogenic singlet-fermion model"
        ),
    )

    parser.add_argument(
        "--study",
        type=str,
        default=None,
        help=(
            "output/report study folder. Defaults to smoke, extended, "
            "interesting, or single depending on the selected run mode; "
            "examples: hypercharge, dimensions"
        ),
    )
    parser.add_argument(
        "--numerical",
        type=Path,
        default=None,
        help=(
            "optional JSON parameter point for numerical matched-EFT running; "
            "e.g. python pipeline.py --numerical examples/t3_numerical_example.json"
        ),
    )
    parser.add_argument(
        "--alpha",
        type=int,
        default=None,
        help=(
            "hypercharge parameter for three-number --dims mode (default 0). "
            "Two-number shared-scalar mode fixes alpha=-1."
        ),
    )
    parser.add_argument(
        "--threshold",
        action="append",
        nargs="+",
        metavar="FIELD",
        default=None,
        help=(
            "ordered threshold group; ordinary mode uses F,S1,S2 and "
            "shared-scalar mode uses F,S. Repeat the option for successive "
            "thresholds. Fields in one group are integrated out together. "
            "If omitted, F S1 S2 are integrated out together."
        ),
    )
    parser.add_argument(
        "--threshold-scale",
        action="append",
        metavar="SCALE",
        default=None,
        help=(
            "matching scale for each --threshold occurrence, in the same "
            "order. Values may be symbolic (MF, MS1, ...) or numeric. "
            "For F -> (S1,S2), defaults are MF and MS."
        ),
    )
    parser.add_argument(
        "--eft-max-dimension",
        type=int,
        choices=(5, 6),
        default=5,
        help=(
            "maximum EFT operator dimension retained by Matchete matching. "
            "The default is 5. Use 6 for scalar-first studies that must retain "
            "the leading psi^2 phi^3 intermediate operator."
        ),
    )
    parser.add_argument(
        "--allow-truncated-scalar-first",
        action="store_true",
        help=(
            "allow a scalar-first threshold plan with the current d<=5 "
            "truncation. This is diagnostic/experimental only: scalar-first "
            "matching can require an intermediate dimension-six operator for "
            "the leading Weinberg contribution, and unsupported intermediate "
            "EFT running will still block an authoritative low-energy result."
        ),
    )
    parser.add_argument(
        "--debug-reports",
        action="store_true",
        help=(
            "also keep raw Wolfram logs and generate the full "
            "UV/EFT expression reports"
        ),
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help=(
            "bypass the neutral-component requirement and the current "
            "SU(2)-dimension production-support limit. Genuine T3 topology "
            "conditions are still enforced."
        ),
    )
    return parser


def resolve_model_mode(
    parser: argparse.ArgumentParser,
    args: argparse.Namespace,
) -> bool:
    """Validate dimension/alpha options and return shared-scalar mode."""

    shared_scalar_mode = False
    if args.dims is not None:
        if len(args.dims) not in (2, 3):
            parser.error("--dims requires either DS DF or DS1 DS2 DF.")
        shared_scalar_mode = len(args.dims) == 2

    if shared_scalar_mode:
        if args.alpha is not None and args.alpha != -1:
            parser.error("Two-number shared-scalar mode requires alpha=-1.")
        args.alpha = -1
    elif args.alpha is None:
        args.alpha = 0

    return shared_scalar_mode


def resolve_pipeline_plan(
    parser: argparse.ArgumentParser,
    args: argparse.Namespace,
    *,
    shared_scalar_mode: bool,
) -> PipelinePlan:
    """Build the physical EFT execution plan requested on the command line."""

    try:
        plan = PipelinePlan.from_threshold_configuration(
            args.threshold,
            args.threshold_scale,
            shared_scalar=shared_scalar_mode,
            truncation=EFTTruncation(
                max_operator_dimension=args.eft_max_dimension,
            ),
        )
    except (TypeError, ValueError) as exc:
        parser.error(str(exc))
        raise AssertionError("argparse.error() should not return") from exc

    if (
        plan.scalar_first_truncates_leading_path
        and not args.allow_truncated_scalar_first
    ):
        parser.error(
            "Scalar-first threshold ordering with the current d<=5 truncation "
            "can discard an intermediate dimension-six operator that feeds the "
            "leading T3 Weinberg coefficient. Use "
            "--allow-truncated-scalar-first only for an explicit experimental "
            "truncated run; the verified fermion-first/common-threshold route "
            "remains authoritative."
        )

    return plan


def study_name(args: argparse.Namespace) -> str:
    """Return the output/report study directory name for one pipeline run."""

    if args.study:
        # Keep '/' so --full can intentionally create nested study roots such
        # as full/hypercharge. Spaces are normalised for CLI convenience.
        return args.study.strip().replace(" ", "_")
    if args.smoke:
        return "smoke"
    if args.hypercharge_comparison:
        return "hypercharge"
    if args.dimension_comparison:
        return "dimensions"
    if args.dims:
        return "single"
    return
