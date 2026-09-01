from pathlib import Path

path = Path("pipeline.py")

if not path.exists():
    raise FileNotFoundError(path)

text = path.read_text(encoding="utf-8")

import_anchor = "from NeutrinoMassStage import run_neutrino_mass_stage\n"
import_line = "from NumericalPipelineStage import run_numerical_pipeline_stage\n"

if import_line not in text:
    if import_anchor not in text:
        raise RuntimeError(
            "Could not locate NeutrinoMassStage import."
        )

    text = text.replace(
        import_anchor,
        import_anchor + import_line,
        1,
    )

parser_anchor = """    parser.add_argument(
        "--alpha",
"""

numerical_arg = """    parser.add_argument(
        "--numerical",
        type=Path,
        default=None,
        help=(
            "optional JSON parameter point for numerical "
            "matched-EFT running"
        ),
    )

"""

if "--numerical" not in text:
    if parser_anchor not in text:
        raise RuntimeError(
            "Could not locate parser insertion point."
        )

    text = text.replace(
        parser_anchor,
        numerical_arg + parser_anchor,
        1,
    )

finish_signature = "def finish_runs(records: list[RunRecord]) -> int:\n"
new_finish_signature = (
    "def finish_runs(\n"
    "    records: list[RunRecord],\n"
    "    numerical_config: Path | None = None,\n"
    ") -> int:\n"
)

if new_finish_signature not in text:
    if finish_signature not in text:
        raise RuntimeError(
            "Could not locate finish_runs signature."
        )

    text = text.replace(
        finish_signature,
        new_finish_signature,
        1,
    )

insertion_anchor = """        summary.update(mass_summary)

        print(
"""

numerical_block = """        summary.update(mass_summary)

        if numerical_config is not None:
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
                continue

            summary.update(numerical_summary)

            print(
                f"  {record.name}: numerical RGE=Success"
                f" -> "
                f"{record.output_dir / numerical_summary['NeutrinoMassMatrixLowScaleFile']}"
            )

        print(
"""

if "numerical RGE=Success" not in text:
    if insertion_anchor not in text:
        raise RuntimeError(
            "Could not locate neutrino mass summary block."
        )

    text = text.replace(
        insertion_anchor,
        numerical_block,
        1,
    )

text = text.replace(
    "        return finish_runs([record])\n",
    "        return finish_runs([record], args.numerical)\n",
    1,
)

text = text.replace(
    "    return finish_runs(records)\n",
    "    return finish_runs(records, args.numerical)\n",
    1,
)

path.write_text(
    text,
    encoding="utf-8",
)

print("Patched:", path.resolve())
print(
    "Numerical RGE option present:",
    "--numerical" in text,
)
print(
    "Numerical RGE stage present:",
    "numerical RGE=Success" in text,
)
