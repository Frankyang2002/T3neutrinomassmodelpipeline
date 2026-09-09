# Main script for the T3 model pipeline.
#
# This file:
# 1. Reads the requested T3 representation assignment / benchmark model.
# 2. Builds the UV model and performs EFT matching.
# 3. Organises the matched Weinberg coefficient C5.
# 4. Runs the UV and EFT RGE stages.
# 5. Constructs the neutrino-mass matrix.
# 6. Optionally performs numerical running and calculates neutrino observables.
# 7. Generates the final reports and run summaries.
#
from __future__ import annotations

import argparse
import json
from pathlib import Path

from RGE.matching.MatchedEFTRGE import run_matched_eft_rge
from RGE.stages.FlavorMatchedRGEStage import run_flavor_matched_rge
from RGE.stages.NeutrinoMassStage import run_neutrino_mass_stage
from RGE.stages.NumericalPipelineStage import run_numerical_pipeline_stage
from RGE.phenomenology.NeutrinoObservables import run_neutrino_observables_stage
from RGE.running.RGBetaT3Running import run_rgbeta_t3

from common.Paths import (
    OUTPUT_DIR,
    REPORT_OUTPUT_DIR,
)

from common.Records import RunRecord
from common.T3Model import EXTENDED, INTERESTING, SMOKE
from Lagrangian.Runner import validate_dimensions, obtain_class_dimensions
from Reports.ReportGeneration import (
    report_output_dir_for,
    write_reports,
)
from Reports.RGEReport import write_and_compile_rge_report
from Reports.RGEComparison import write_and_compile_rge_comparison
from Reports.EFTRGEComparison import write_and_compile_weinberg_rge_comparison
from Reports.EFTRenormalisableRGE import write_and_compile_eft_renormalisable_rge


def run_uv_rgbeta_stage(record: RunRecord) -> bool:
    """ Generate the one-loop beta functions for the renormalisable UV T3 model
    using RGBeta.

    This stage describes running above the heavy-particle matching threshold,
    before the T3 fields are integrated out."""
    summary = record.summary

    # Fails if we dont have an actual build for our model
    if summary.get("BuildStatus") != "Success":
        summary["UVRGEStatus"] = "NotRun"
        return False

    data_dir = record.output_dir / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    output_path = data_dir / "uv_rgbeta_rge.json"

    print(f"  {record.name}: starting RGBeta UV-RGE stage...", flush=True)

    # Use run_rgbeta_t3
    try:
        result = run_rgbeta_t3(
            record.d_s1,
            record.d_s2,
            record.d_f,
            record.alpha,
        )
    except Exception as exc:
        summary["UVRGEStatus"] = "Failed"
        summary["UVRGEError"] = str(exc)
        print(f"  {record.name}: RGBeta UV-RGE failed: {exc}")
        return False

    # Pass the T3 SU(2) representation dimensions and hypercharge parameter
    # to the RGBeta Wolfram runner, which constructs the UV model beta functions.
    output_path.write_text(
        json.dumps(result.raw, indent=2),
        encoding="utf-8",
    )

    summary["UVRGEStatus"] = result.status
    summary["UVRGEFile"] = output_path.relative_to(record.output_dir).as_posix()
    summary["UVRGEBetaCount"] = len(result.betas)

    print(
        f"  {record.name}: RGBeta UV-RGE={result.status} "
        f"({len(result.betas)} beta functions) -> {output_path}"
    )

    return result.status == "Success"

def print_summary(records: list[RunRecord]) -> int:
    """Print the results of all completed T3 runs."""

    print(
        "\n"
        + "=" * 72
        + "\nT3 MODEL SUMMARY\n"
        + "=" * 72
    )

    successful = 0

    for record in records:
        summary = record.summary

        build_ok = summary.get("BuildStatus") == "Success"
        match_ok = summary.get("MatchingStatus") == "Success"
        uv_rge_ok = summary.get("UVRGEStatus") == "Success"

        successful += int(build_ok and match_ok and uv_rge_ok)

        print(
            f"{record.name} alpha={record.alpha}: "
            f"build={summary.get('BuildStatus')}, "
            f"UV-RGE={summary.get('UVRGEStatus')}, "
            f"match={summary.get('MatchingStatus')}, "
            f"T3={summary.get('T3IngredientsPresent')}, "
            f"Weinberg={summary.get('WeinbergOperatorPresent')}"
        )

        # If matching did not generate the Weinberg operator, there is
        # no C5 coefficient to extract.
        if not summary.get("WeinbergOperatorPresent"):
            continue

        extraction = summary.get(
            "WeinbergExtractionStatus",
            "Unknown",
        )

        if extraction == "Success":
            n_holo = summary.get(
                "WeinbergHolomorphicTermCount",
                0,
            )
            n_hc = summary.get(
                "WeinbergConjugateTermCount",
                0,
            )

            print(
                f"  C5 extraction: Success "
                f"({n_holo} holomorphic + {n_hc} HC terms)"
            )

            if coefficient_file := summary.get("WeinbergCoefficientFile"):
                print(
                    f"  C5: "
                    f"{record.output_dir / coefficient_file}"
                )

        else:
            print(
                f"  Weinberg terms: "
                f"{summary.get('WeinbergTermCount', 0)}; "
                f"C5: {extraction}"
            )

    # Save all model summaries together so we can easily compare runs
    # or use them for regression testing.
    aggregate = OUTPUT_DIR / "t3_model_comparison.json"

    aggregate.write_text(
        json.dumps(
            [record.summary for record in records],
            indent=2,
        ),
        encoding="utf-8",
    )

    print(
        f"\n{successful}/{len(records)} completed build+matching."
        f"\nAggregate: {aggregate}"
    )

    # Standard command-line convention:
    #   0 = everything succeeded
    #   1 = at least one model failed
    return 0 if successful == len(records) else 1


def organise_c5_input(record: RunRecord) -> Path | None:
    """
    Put the matched Weinberg coefficient in its standard machine-readable
    location before any EFT RGE stage starts.
    From this point onward, c5_coefficient.txt represents C5 at the
    matching scale M and is the main input to the EFT-side calculations.
    """
    
    summary = record.summary
    coefficient_file = summary.get("WeinbergCoefficientFile")

    if not coefficient_file:
        return None

    c5_path = record.output_dir / coefficient_file

    if not c5_path.is_file():
        return None

    data_dir = record.output_dir / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    organised_path = data_dir / c5_path.name

    if c5_path != organised_path:
        if organised_path.exists():
            raise FileExistsError(
                f"Refusing to overwrite existing RGE input: {organised_path}"
            )

        c5_path.replace(organised_path)
        c5_path = organised_path

    summary["WeinbergCoefficientFile"] = (
        c5_path.relative_to(record.output_dir).as_posix()
    )
    return c5_path


def matched_c5_path(record: RunRecord) -> Path | None:
    """Return the organised matched C5 coefficient for the RGE stages."""

    coefficient_file = record.summary.get("WeinbergCoefficientFile")

    if not coefficient_file:
        return None

    c5_path = record.output_dir / coefficient_file

    if not c5_path.is_file():
        return None

    return c5_path


def run_matched_eft_rge_stage(
    record: RunRecord,
    debug_reports: bool = False,
) -> bool:
    """
    Evaluate the matched Weinberg coefficient using the general psi^2 phi^2
    one-loop RGE machinery and verify that it reduces to the known
    one-generation SMEFT Weinberg RGE. 
    Just for testing as one generation is simple to test
    """

    summary = record.summary
    c5_path = matched_c5_path(record)

    if c5_path is None:
        summary["RGEStatus"] = "NotRun"
        return False

    print(f"  {record.name}: starting matched-EFT RGE stage...", flush=True)

    try:
        rge_summary = run_matched_eft_rge(
            c5_path=c5_path,
            output_dir=record.output_dir,
            debug_outputs=debug_reports,
        )
    except Exception as exc:
        summary["RGEStatus"] = "Failed"
        summary["RGEError"] = str(exc)
        print(
            f"  {record.name}: matched-EFT RGE failed: {exc}"
        )
        return False

    summary.update(rge_summary)

    print(
        f"  {record.name}: matched-EFT RGE=Success"
        f" -> {record.output_dir / rge_summary['C5BetaFile']}"
    )

    return summary.get("RGEStatus") == "Success"


def run_flavor_rge_stage(
    record: RunRecord,
    debug_reports: bool = False,
) -> bool:
    """Run the symbolic full-flavor Weinberg RGE for one model."""

    summary = record.summary
    c5_path = matched_c5_path(record)

    if c5_path is None:
        summary["FlavorRGEStatus"] = "NotRun"
        return False

    print(
        f"  {record.name}: starting symbolic full-flavor RGE stage...",
        flush=True,
    )

    try:
        flavor_summary = run_flavor_matched_rge(
            c5_path=c5_path,
            output_dir=record.output_dir,
            debug_outputs=debug_reports,
        )
    except Exception as exc:
        summary["FlavorRGEStatus"] = "Failed"
        summary["FlavorRGEError"] = str(exc)
        print(
            f"  {record.name}: full-flavor RGE failed: {exc}"
        )
        return False

    summary.update(flavor_summary)

    print(
        f"  {record.name}: full-flavor RGE=Success"
        f" -> "
        f"{record.output_dir / flavor_summary['C5FlavorBetaMatrixFile']}"
    )

    return summary.get("FlavorRGEStatus") == "Success"


def run_symbolic_neutrino_mass_stage(record: RunRecord) -> bool:
    """Construct the symbolic neutrino-mass matrix from the matched C5."""

    summary = record.summary
    c5_path = matched_c5_path(record)

    if c5_path is None:
        summary["NeutrinoMassStatus"] = "NotRun"
        return False

    print(
        f"  {record.name}: starting symbolic neutrino mass stage...",
        flush=True,
    )

    try:
        mass_summary = run_neutrino_mass_stage(
            c5_path=c5_path,
            output_dir=record.output_dir,
        )
    except Exception as exc:
        summary["NeutrinoMassStatus"] = "Failed"
        summary["NeutrinoMassError"] = str(exc)
        print(
            f"  {record.name}: neutrino mass stage failed: {exc}"
        )
        return False

    summary.update(mass_summary)

    print(
        f"  {record.name}: neutrino mass=Success"
        f" -> "
        f"{record.output_dir / mass_summary['NeutrinoMassMatrixFile']}"
    )

    return summary.get("NeutrinoMassStatus") == "Success"


def run_numerical_rge_stage(
    record: RunRecord,
    numerical_config: Path,
) -> bool:
    """Run numerical EFT evolution and obtain the resulting neutrino observables."""

    summary = record.summary
    c5_path = matched_c5_path(record)

    if c5_path is None:
        summary["NumericalRGEStatus"] = "NotRun"
        return False

    print(f"  {record.name}: starting numerical RGE stage...", flush=True)

    try:
        numerical_summary = run_numerical_pipeline_stage(
            c5_path=c5_path,
            output_dir=record.output_dir,
            config_path=numerical_config,
        )
    except Exception as exc:
        summary["NumericalRGEStatus"] = "Failed"
        summary["NumericalRGEError"] = str(exc)
        print(
            f"  {record.name}: numerical RGE failed: {exc}"
        )
        return False

    summary.update(numerical_summary)

    print(
        f"  {record.name}: numerical RGE calculation finished.",
        flush=True,
    )

    mass_matrix_path = (
        record.output_dir
        / numerical_summary["NeutrinoMassMatrixLowScaleFile"]
    )

    numerical_payload = json.loads(
        numerical_config.read_text(encoding="utf-8")
    )
    ordering = numerical_payload.get("ordering", "NO")

    print(f"  {record.name}: starting neutrino observables...", flush=True)

    try:
        observable_summary = run_neutrino_observables_stage(
            mass_matrix_path=mass_matrix_path,
            output_dir=record.output_dir,
            ordering=ordering,
        )
    except Exception as exc:
        summary["NeutrinoObservableStatus"] = "Failed"
        summary["NeutrinoObservableError"] = str(exc)
        print(
            f"  {record.name}: neutrino observables failed: {exc}"
        )
        return False

    summary.update(observable_summary)

    print(
        f"  {record.name}: neutrino observables=Success"
        f" -> "
        f"{record.output_dir / observable_summary['NeutrinoObservablesFile']}"
    )

    print(
        f"  {record.name}: numerical RGE=Success"
        f" -> "
        f"{record.output_dir / numerical_summary['NeutrinoMassMatrixLowScaleFile']}"
    )

    return (
        summary.get("NumericalRGEStatus") == "Success"
        and summary.get("NeutrinoObservableStatus") == "Success"
    )


def finish_runs(
    records: list[RunRecord],
    debug_reports: bool = False,
    physics_failed: bool = False,
) -> int:
    """Print final summaries and generate reports after the physics pipeline."""

    status = print_summary(records)

    if physics_failed:
        status = 1

    # Generate the Lagrangian/matching reports only after all physics stages have
    # updated the run summaries.
    write_reports(records, debug_reports)

    # Compare all successful UV one-loop beta functions coupling-by-coupling.
    write_and_compile_rge_comparison(records)

    # Keep the EFT reporting split into the dimension-five Weinberg sector and
    # the renormalisable SM sector that survives below the heavy threshold.
    write_and_compile_weinberg_rge_comparison(records)
    write_and_compile_eft_renormalisable_rge(records)

    for record in records:
        coefficient_pdf = report_output_dir_for(record) / "c5_coefficient.pdf"
        record.summary["WeinbergCoefficientPDFFile"] = (
            coefficient_pdf.name if coefficient_pdf.exists() else ""
        )

        # The RGE report is also a final reporting step.  By this point all
        # symbolic and optional numerical RGE stages have already run.
        write_and_compile_rge_report(record)

    # Save all model summaries together so we can easily compare runs
    # or use them for regression testing.
    aggregate = OUTPUT_DIR / "t3_model_comparison.json"

    aggregate.write_text(
        json.dumps(
            [record.summary for record in records],
            indent=2,
        ),
        encoding="utf-8",
    )

    print(f"Aggregate: {aggregate}")

    return status


# Function:
# 1. Get all arguments
# 2. Start process
# 3. Get reports
def main() -> int:
    # Get Arguments
    parser = argparse.ArgumentParser(
        description=(
            "Run T3 matching for known benchmark models "
            "or arbitrary valid SU(2) irreps."
        )
    )

    # Mutually exclusive command arguments
    mode = parser.add_mutually_exclusive_group()

    # If we use --smoke we use the 5 T3 models we know
    mode.add_argument(
        "--smoke",
        action="store_true",
        help="five T3 regression models",
    )

    # More running (not used)
    mode.add_argument(
        "--extended",
        action="store_true",
        help="seven historical benchmark points",
    )

    # Dimension input to use a specific diagram
    mode.add_argument(
        "--dims",
        nargs=3,
        type=int,
        metavar=("DS1", "DS2", "DF"),
        help=(
            "run one representation assignment, "
            "e.g. --dims 3 1 2"
        ),
    )

    # We can attach a file for UV scale initial values for our couplings etc to RGE down to EFT
    parser.add_argument(
        "--numerical",
        type=Path,
        default=None,
        help=(
            "optional JSON parameter point for numerical "
            "matched-EFT running E.g python pipeline.py --numerical examples/t3_numerical_example.json"
        ),
    )

    # Adds the alpha hypercharge which is added onto the dims
    parser.add_argument(
        "--alpha",
        type=int,
        default=0,
        help=(
            "T3 hypercharge parameter for --dims mode "
            "(default: 0)"
        ),
    )

    # For logging just in case
    parser.add_argument(
        "--debug-reports",
        action="store_true",
        help=(
            "also keep raw Wolfram logs and generate the full "
            "UV/EFT expression reports"
        ),
    )

    args = parser.parse_args()

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )
    REPORT_OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )
    # ----------------------------------------------------------------------
    # Lagrangian and Weinberg Coefficient pipeline
    # ----------------------------------------------------------------------

    # If the dimensions happen to match T3-A ... T3-E, validate_dimensions()
    # automatically recognises and labels the model appropriately.
    if args.dims:
        d_s1, d_s2, d_f = args.dims

        try:
            record = validate_dimensions(
                d_s1,
                d_s2,
                d_f,
                args.alpha,
                args.debug_reports,
                False,
            )
        except ValueError as exc:
            parser.error(str(exc))

        records = [record]

    else:
        # The benchmark lists still use the familiar A-E notation because
        # it is convenient for regression testing and comparison with the paper.
        if args.smoke:
            mode_name = "smoke"
            points = SMOKE

        elif args.extended:
            mode_name = "extended"
            points = EXTENDED

        else:
            mode_name = "interesting"
            points = INTERESTING

        print(
            f"T3 scan mode: {mode_name}; "
            f"{len(points)} model(s)."
        )

        # This line only exist if we use the special interesting, extended and smoke options where we dont input any dimensions
        # Known A-E models are converted to dimensions first and then sent
        # through exactly the same validate_dimensions() path as generalised models.
        records = [
            obtain_class_dimensions(
                model_class,
                alpha,
                args.debug_reports,
                False,
            )
            for model_class, alpha in points
        ]

    # Organise the matched coefficient before starting the RGE pipeline, so all
    # later stages read C5 from the same final machine-readable location.
    for record in records:
        organise_c5_input(record)

    # ----------------------------------------------------------------------
    # RGE pipeline
    # ----------------------------------------------------------------------

    physics_failed = False

    # For every model, only run the UV RGE stage if the model Lagrangian was built successfully. Also remember whether any model failed.
    # We want UV RGE for comparison with EFT RGE
    for record in records:
        if record.summary.get("BuildStatus") == "Success":
            physics_failed |= not run_uv_rgbeta_stage(record)
        else:
            record.summary["UVRGEStatus"] = "NotRun"

    # After matching has generated the Weinberg coefficient C5, the EFT RGE pipeline starts here. 
    # Each later stage runs only if the previous stage for that model succeeded.
    for record in records:
        summary = record.summary

        # Check if any part works for a model, if a model fails we skip
        if summary.get("BuildStatus") != "Success":
            continue

        if summary.get("MatchingStatus") != "Success":
            continue

        if summary.get("WeinbergExtractionStatus") != "Success":
            continue

        # Evaluate the matched Weinberg coefficient with the general
        # one-generation SMEFT RGE machinery and verify the known SMEFT result.
        if not run_matched_eft_rge_stage(
            record,
            args.debug_reports,
        ):
            physics_failed = True
            continue

        # Get the matched one-generation C5 into the full 3x3 lepton-flavor
        # matrix and calculate its symbolic SMEFT beta matrix.
        if not run_flavor_rge_stage(
            record,
            args.debug_reports,
        ):
            physics_failed = True
            continue

        # Convert the matched symbolic flavor coefficient into the symbolic
        # Majorana neutrino-mass matrix using m_nu = -v^2 C5.
        if not run_symbolic_neutrino_mass_stage(record):
            physics_failed = True
            continue

        # If a numerical parameter point and argument is given, evaluate C5 at the
        # matching scale and numerically integrate the coupled one-loop SMEFT RGEs
        # down to the requested low scale. Then calculate neutrino observables.
        if args.numerical is not None:
            if not run_numerical_rge_stage(
                record,
                args.numerical,
            ):
                physics_failed = True
                continue

    # Reports and final summaries happen only after the full physics pipeline.
    return finish_runs(
        records,
        args.debug_reports,
        physics_failed,
    )


if __name__ == "__main__":
    raise SystemExit(main())
