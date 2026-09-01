from pathlib import Path

path = Path("pipeline.py")

if not path.exists():
    raise FileNotFoundError(path)

text = path.read_text(encoding="utf-8")

import_anchor = "from FlavorMatchedRGEStage import run_flavor_matched_rge\n"
import_line = "from NeutrinoMassStage import run_neutrino_mass_stage\n"

if import_line not in text:
    if import_anchor not in text:
        raise RuntimeError(
            "Could not locate FlavorMatchedRGEStage import."
        )

    text = text.replace(
        import_anchor,
        import_anchor + import_line,
        1,
    )

anchor = """        summary.update(flavor_summary)

        print(
"""

replacement = """        summary.update(flavor_summary)

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
            continue

        summary.update(mass_summary)

        print(
            f"  {record.name}: neutrino mass=Success"
            f" -> "
            f"{record.output_dir / mass_summary['NeutrinoMassMatrixFile']}"
        )

        print(
"""

if anchor not in text:
    raise RuntimeError(
        "Could not locate flavor summary update block. "
        "No changes were written."
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
    "Neutrino mass import present:",
    import_line.strip() in text,
)
print(
    "Neutrino mass stage present:",
    "neutrino mass=Success" in text,
)
