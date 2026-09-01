from pathlib import Path

path = Path("pipeline.py")

if not path.exists():
    raise FileNotFoundError(path)

text = path.read_text(encoding="utf-8")

import_anchor = "from MatchedEFTRGE import run_matched_eft_rge\n"
import_line = "from FlavorMatchedRGEStage import run_flavor_matched_rge\n"

if import_line not in text:
    if import_anchor not in text:
        raise RuntimeError("Could not locate MatchedEFTRGE import.")

    text = text.replace(
        import_anchor,
        import_anchor + import_line,
        1,
    )

anchor = """        summary.update(rge_summary)

        print(
"""

replacement = """        summary.update(rge_summary)

        try:
            flavor_summary = run_flavor_matched_rge(
                c5_path=c5_path,
                output_dir=record.output_dir,
            )
        except Exception as exc:
            summary["FlavorRGEStatus"] = "Failed"
            summary["FlavorRGEError"] = str(exc)
            status = 1
            print(
                f"  {record.name}: full-flavor RGE failed: {exc}"
            )
            continue

        summary.update(flavor_summary)

        print(
            f"  {record.name}: full-flavor RGE=Success"
            f" -> "
            f"{record.output_dir / flavor_summary['C5FlavorBetaMatrixFile']}"
        )

        print(
"""

if anchor not in text:
    raise RuntimeError(
        "Could not locate matched-EFT summary update block. "
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
    "Flavor RGE import present:",
    import_line.strip() in text,
)
print(
    "Flavor RGE stage present:",
    "full-flavor RGE=Success" in text,
)
