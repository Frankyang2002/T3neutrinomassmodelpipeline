"""Low-energy Weinberg-operator and neutrino calculations.

Physics role
------------
This module starts from the matched dimension-five Weinberg coefficient ``C5``
and performs the EFT-side calculations that follow heavy-field matching:

1. place the matched ``C5`` in its stable machine-readable location;
2. evaluate the one-generation SMEFT Weinberg RGE cross-check used by the
   existing report/output contract;
3. evaluate the symbolic full-flavor Weinberg RGE;
4. construct the symbolic Majorana neutrino-mass matrix;
5. optionally run the numerical EFT trajectory and neutrino observables.

The central ``pipeline.py`` deliberately calls these functions in this physical
order so that it remains the readable backbone of the project.

Conventions retained from the existing pipeline
------------------------------------------------
The numerical and symbolic downstream modules retain the project's current
Weinberg convention, including the existing relation between ``C5`` and the
Majorana neutrino-mass matrix.  This refactor only moves orchestration details;
it does not alter formulas, summary keys, filenames, printed status messages,
or report inputs.
"""

from __future__ import annotations

import json
from pathlib import Path

from common.RunRecords import RunRecord
from RGE.matching.MatchedWeinbergRGE import run_matched_weinberg_rge
from RGE.phenomenology.NeutrinoObservables import run_neutrino_observables_stage
from RGE.running.weinberg.FlavorMatchedWeinbergStage import (
    run_flavor_matched_weinberg_rge,
    run_symbolic_neutrino_mass_stage as _run_symbolic_neutrino_mass_stage,
)
from RGE.running.weinberg.NumericalWeinbergStage import (
    run_numerical_weinberg_stage as _run_numerical_weinberg_stage,
)


def organise_matched_weinberg_coefficient(record: RunRecord) -> Path | None:
    """Move the matched ``C5`` into ``data/`` before EFT-side running.

    The path recorded in ``summary['WeinbergCoefficientFile']`` is part of the
    existing machine-output contract.  The move is therefore intentionally
    conservative: an existing destination is never overwritten.
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


def matched_weinberg_coefficient_path(record: RunRecord) -> Path | None:
    """Return the organised matched ``C5`` used by low-energy EFT stages."""

    coefficient_file = record.summary.get("WeinbergCoefficientFile")
    if not coefficient_file:
        return None

    c5_path = record.output_dir / coefficient_file
    if not c5_path.is_file():
        return None

    return c5_path


def run_smeft_weinberg_rge(
    record: RunRecord,
    debug_reports: bool = False,
) -> bool:
    """Evaluate the matched Weinberg coefficient with the general tensor RGE.

    This preserves the existing one-generation SMEFT reduction check and the
    associated ``RGE...`` summary/output fields.  It remains part of the current
    production/report contract even though the full-flavor calculation below is
    the physically richer downstream description.
    """

    summary = record.summary
    c5_path = matched_weinberg_coefficient_path(record)

    if c5_path is None:
        summary["RGEStatus"] = "NotRun"
        return False

    print(f"  {record.name}: starting matched-EFT RGE stage...", flush=True)

    try:
        rge_summary = run_matched_weinberg_rge(
            c5_path=c5_path,
            output_dir=record.output_dir,
            debug_outputs=debug_reports,
        )
    except Exception as exc:
        summary["RGEStatus"] = "Failed"
        summary["RGEError"] = str(exc)
        print(f"  {record.name}: matched-EFT RGE failed: {exc}")
        return False

    summary.update(rge_summary)

    print(
        f"  {record.name}: matched-EFT RGE=Success"
        f" -> {record.output_dir / rge_summary['C5BetaFile']}"
    )

    return summary.get("RGEStatus") == "Success"


def run_full_flavor_weinberg_rge(
    record: RunRecord,
    debug_reports: bool = False,
) -> bool:
    """Run the symbolic full-flavor SMEFT Weinberg RGE for one model."""

    summary = record.summary
    c5_path = matched_weinberg_coefficient_path(record)

    if c5_path is None:
        summary["FlavorRGEStatus"] = "NotRun"
        return False

    print(
        f"  {record.name}: starting symbolic full-flavor RGE stage...",
        flush=True,
    )

    try:
        flavor_summary = run_flavor_matched_weinberg_rge(
            c5_path=c5_path,
            output_dir=record.output_dir,
            debug_outputs=debug_reports,
        )
    except Exception as exc:
        summary["FlavorRGEStatus"] = "Failed"
        summary["FlavorRGEError"] = str(exc)
        print(f"  {record.name}: full-flavor RGE failed: {exc}")
        return False

    summary.update(flavor_summary)

    print(
        f"  {record.name}: full-flavor RGE=Success"
        f" -> {record.output_dir / flavor_summary['C5FlavorBetaMatrixFile']}"
    )

    return summary.get("FlavorRGEStatus") == "Success"


def build_symbolic_neutrino_mass(record: RunRecord) -> bool:
    """Construct the symbolic Majorana neutrino-mass matrix from matched ``C5``."""

    summary = record.summary
    c5_path = matched_weinberg_coefficient_path(record)

    if c5_path is None:
        summary["NeutrinoMassStatus"] = "NotRun"
        return False

    print(
        f"  {record.name}: starting symbolic neutrino mass stage...",
        flush=True,
    )

    try:
        mass_summary = _run_symbolic_neutrino_mass_stage(
            c5_path=c5_path,
            output_dir=record.output_dir,
        )
    except Exception as exc:
        summary["NeutrinoMassStatus"] = "Failed"
        summary["NeutrinoMassError"] = str(exc)
        print(f"  {record.name}: neutrino mass stage failed: {exc}")
        return False

    summary.update(mass_summary)

    print(
        f"  {record.name}: neutrino mass=Success"
        f" -> {record.output_dir / mass_summary['NeutrinoMassMatrixFile']}"
    )

    return summary.get("NeutrinoMassStatus") == "Success"


def run_numerical_neutrino_observables(
    record: RunRecord,
    numerical_config: Path,
) -> bool:
    """Run numerical SMEFT evolution and calculate neutrino observables.

    ``numerical_config`` retains the existing JSON schema, including the
    optional ``ordering`` entry (default ``"NO"``).  All established summary
    keys and filenames are preserved.
    """

    summary = record.summary
    c5_path = matched_weinberg_coefficient_path(record)

    if c5_path is None:
        summary["NumericalRGEStatus"] = "NotRun"
        return False

    print(f"  {record.name}: starting numerical RGE stage...", flush=True)

    try:
        numerical_summary = _run_numerical_weinberg_stage(
            c5_path=c5_path,
            output_dir=record.output_dir,
            config_path=numerical_config,
        )
    except Exception as exc:
        summary["NumericalRGEStatus"] = "Failed"
        summary["NumericalRGEError"] = str(exc)
        print(f"  {record.name}: numerical RGE failed: {exc}")
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
        print(f"  {record.name}: neutrino observables failed: {exc}")
        return False

    summary.update(observable_summary)

    print(
        f"  {record.name}: neutrino observables=Success"
        f" -> {record.output_dir / observable_summary['NeutrinoObservablesFile']}"
    )
    print(
        f"  {record.name}: numerical RGE=Success"
        f" -> {record.output_dir / numerical_summary['NeutrinoMassMatrixLowScaleFile']}"
    )

    return (
        summary.get("NumericalRGEStatus") == "Success"
        and summary.get("NeutrinoObservableStatus") == "Success"
    )
