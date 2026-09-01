from pathlib import Path

path = Path("pipeline.py")

if not path.exists():
    raise FileNotFoundError(path)

text = path.read_text(encoding="utf-8")

import_anchor = "from NumericalPipelineStage import run_numerical_pipeline_stage\n"
import_line = "from NeutrinoObservables import run_neutrino_observables_stage\n"

if import_line not in text:
    if import_anchor not in text:
        raise RuntimeError(
            "Could not locate NumericalPipelineStage import."
        )

    text = text.replace(
        import_anchor,
        import_anchor + import_line,
        1,
    )

anchor = """            summary.update(numerical_summary)

            print(
"""

replacement = """            summary.update(numerical_summary)

            mass_matrix_path = (
                record.output_dir
                / numerical_summary["NeutrinoMassMatrixLowScaleFile"]
            )

            try:
                observable_summary = run_neutrino_observables_stage(
                    mass_matrix_path=mass_matrix_path,
                    output_dir=record.output_dir,
                )
            except Exception as exc:
                summary["NeutrinoObservableStatus"] = "Failed"
                summary["NeutrinoObservableError"] = str(exc)
                status = 1
                print(
                    f"  {record.name}: neutrino observables failed: {exc}"
                )
                continue

            summary.update(observable_summary)

            print(
                f"  {record.name}: neutrino observables=Success"
                f" -> "
                f"{record.output_dir / observable_summary['NeutrinoObservablesFile']}"
            )

            print(
"""

if "neutrino observables=Success" not in text:
    if anchor not in text:
        raise RuntimeError(
            "Could not locate numerical summary block."
        )

    text = text.replace(
        anchor,
        replacement,
        1,
    )

path.write_text(
    text,
    encoding="utf-8",
)

print("Patched:", path.resolve())
print(
    "Neutrino observable stage present:",
    "neutrino observables=Success" in text,
)
