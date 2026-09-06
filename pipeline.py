# This is the pipeline that starts the process, and gets all the inputs

from __future__ import annotations

import argparse
import json
from pathlib import Path

from RGE.matching.MatchedEFTRGE import run_matched_eft_rge
from RGE.stages.FlavorMatchedRGEStage import run_flavor_matched_rge
from RGE.stages.NeutrinoMassStage import run_neutrino_mass_stage
from RGE.stages.NumericalPipelineStage import run_numerical_pipeline_stage
from RGE.phenomenology.NeutrinoObservables import run_neutrino_observables_stage
from RGE.general.GeneralT3WeinbergRGEStage import (
    run_general_t3_weinberg_rge_stage,
)
from RGE.running.RGBetaT3Running import run_rgbeta_t3

from common.Paths import (
    OUTPUT_DIR,
    PROJECT_ROOT,
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


def run_uv_rgbeta_stage(record: RunRecord) -> bool:
    """Generate and save the one-loop renormalisable UV RGEs for one model."""
    summary = record.summary

    if summary.get("BuildStatus") != "Success":
        summary["UVRGEStatus"] = "NotRun"
        return False

    data_dir = record.output_dir / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    output_path = data_dir / "uv_rgbeta_rge.json"

    print(f"  {record.name}: starting RGBeta UV-RGE stage...", flush=True)

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
    """Move the machine-readable matched coefficient into the data folder."""

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


def finish_runs(
    records: list[RunRecord],
    numerical_config: Path | None = None,
    debug_reports: bool = False,
) -> int:
    """Print the scan summary and generate all Lagrangian reports."""

    # Organise the matched coefficient before printing paths or writing the
    # aggregate summary, so every reported filename points to its final place.
    for record in records:
        organise_c5_input(record)

    # The UV RGE is a separate symbolic stage.  Run it before the summary so
    # its status and output file are included in both terminal and JSON reports.
    uv_rge_failed = False
    for record in records:
        if record.summary.get("BuildStatus") == "Success":
            uv_rge_failed |= not run_uv_rgbeta_stage(record)
        else:
            record.summary["UVRGEStatus"] = "NotRun"

    status = print_summary(records)
    if uv_rge_failed:
        status = 1

    write_reports(records, debug_reports)

    for record in records:
        coefficient_pdf = report_output_dir_for(record) / "c5_coefficient.pdf"
        record.summary["WeinbergCoefficientPDFFile"] = (
            coefficient_pdf.name if coefficient_pdf.exists() else ""
        )

    for record in records:
        summary = record.summary

        if summary.get("T3RGETensorExportStatus") == "Success":
            exchange_file = summary.get("T3RGETensorExchangeFile", "")
            exchange_path = record.output_dir / exchange_file
            general_rge_path = (
                record.output_dir / "data" / "general_t3_weinberg_rge.json"
            )
            print(
                f"  {record.name}: starting representation-generic RGE stage...",
                flush=True,
            )
            try:
                general_rge = run_general_t3_weinberg_rge_stage(
                    exchange_path=exchange_path,
                    output_path=general_rge_path,
                )
            except Exception as exc:
                summary["GeneralT3RGEStatus"] = "Failed"
                summary["GeneralT3RGEError"] = str(exc)
                status = 1
                print(f"  {record.name}: representation-generic RGE failed: {exc}")
            else:
                summary["GeneralT3RGEStatus"] = general_rge["status"]
                summary["GeneralT3RGEFile"] = general_rge_path.relative_to(
                    record.output_dir
                ).as_posix()
                summary["GeneralT3RGEBetaKappaOverKappa"] = general_rge[
                    "beta_kappa_over_kappa"
                ]["complete_one_generation"]["sympy"]
                if general_rge["status"] != "Success":
                    status = 1
                print(
                    f"  {record.name}: representation-generic RGE="
                    f"{general_rge['status']} -> {general_rge_path}"
                )

        if summary.get("BuildStatus") != "Success":
            continue

        if summary.get("MatchingStatus") != "Success":
            continue

        if summary.get("WeinbergExtractionStatus") != "Success":
            continue

        coefficient_file = summary.get("WeinbergCoefficientFile")

        if not coefficient_file:
            continue

        c5_path = record.output_dir / coefficient_file

        if not c5_path.exists():
            continue

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
            status = 1
            print(
                f"  {record.name}: matched-EFT RGE failed: {exc}"
            )
            write_and_compile_rge_report(record)
            continue

        summary.update(rge_summary)

        print(f"  {record.name}: starting symbolic full-flavor RGE stage...", flush=True)
        try:
            flavor_summary = run_flavor_matched_rge(
                c5_path=c5_path,
                output_dir=record.output_dir,
                debug_outputs=debug_reports,
            )
        except Exception as exc:
            summary["FlavorRGEStatus"] = "Failed"
            summary["FlavorRGEError"] = str(exc)
            status = 1
            print(
                f"  {record.name}: full-flavor RGE failed: {exc}"
            )
            write_and_compile_rge_report(record)
            continue

        summary.update(flavor_summary)

        print(f"  {record.name}: starting symbolic neutrino mass stage...", flush=True)
        try:
            mass_summary = run_neutrino_mass_stage(
                c5_path=c5_path,
                output_dir=record.output_dir,
            )
        except Exception as exc:
            summary["NeutrinoMassStatus"] = "Failed"
            summary["NeutrinoMassError"] = str(exc)
            status = 1
            print(
                f"  {record.name}: neutrino mass stage failed: {exc}"
            )
            write_and_compile_rge_report(record)
            continue

        summary.update(mass_summary)

        if numerical_config is not None:
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
                status = 1
                print(
                    f"  {record.name}: numerical RGE failed: {exc}"
                )
                write_and_compile_rge_report(record)
                continue

            summary.update(numerical_summary)

            print(f"  {record.name}: numerical RGE calculation finished.", flush=True)

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
                status = 1
                print(
                    f"  {record.name}: neutrino observables failed: {exc}"
                )
                write_and_compile_rge_report(record)
                continue

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

        print(
            f"  {record.name}: neutrino mass=Success"
            f" -> "
            f"{record.output_dir / mass_summary['NeutrinoMassMatrixFile']}"
        )

        print(
            f"  {record.name}: full-flavor RGE=Success"
            f" -> "
            f"{record.output_dir / flavor_summary['C5FlavorBetaMatrixFile']}"
        )

        print(
            f"  {record.name}: matched-EFT RGE=Success"
            f" -> {record.output_dir / rge_summary['C5BetaFile']}"
        )

        write_and_compile_rge_report(record)

    aggregate = OUTPUT_DIR / "t3_model_comparison.json"

    aggregate.write_text(
        json.dumps(
            [record.summary for record in records],
            indent=2,
        ),
        encoding="utf-8",
    )

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

    # We get to observe Component index information of our objects, getting explicit RGE tensors used for representation RGE
    # We can see the CG coefficients which are allowed for the field combination
    parser.add_argument(
        "--rge-tensors",
        action="store_true",
        help=(
            "export the exact component tensors needed by the "
            "representation-generic RGE engine"
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
                args.rge_tensors,
            )
        except ValueError as exc:
            parser.error(str(exc))

        # This gives us  our reports based on our records
        return finish_runs(
            [record],
            args.numerical,
            args.debug_reports,
        )

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
            args.rge_tensors,
        )
        for model_class, alpha in points
    ]

    # This only gives us our reports created from the records
    return finish_runs(
        records,
        args.numerical,
        args.debug_reports,
    )


if __name__ == "__main__":
    raise SystemExit(main())
