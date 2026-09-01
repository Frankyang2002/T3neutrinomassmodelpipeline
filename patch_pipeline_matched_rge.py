from pathlib import Path

path = Path("pipeline.py")

if not path.exists():
    raise FileNotFoundError(path)

text = path.read_text(encoding="utf-8")

import_anchor = "from pathlib import Path\n"
import_line = "from MatchedEFTRGE import run_matched_eft_rge\n"

if import_line not in text:
    if import_anchor not in text:
        raise RuntimeError("Could not locate pathlib import.")
    text = text.replace(
        import_anchor,
        import_anchor + "\n" + import_line,
        1,
    )

old_lines = [
    'def finish_runs(records: list[RunRecord]) -> int:',
    '    """Print the scan summary and generate all Lagrangian reports."""',
    '',
    '    status = print_summary(records)',
    '',
    '    write_reports(records)',
    '',
    '    return status',
]
old = "\n".join(old_lines) + "\n"

new_lines = [
    'def finish_runs(records: list[RunRecord]) -> int:',
    '    """Print the scan summary and generate all Lagrangian reports."""',
    '',
    '    status = print_summary(records)',
    '',
    '    write_reports(records)',
    '',
    '    for record in records:',
    '        summary = record.summary',
    '',
    '        if summary.get("BuildStatus") != "Success":',
    '            continue',
    '',
    '        if summary.get("MatchingStatus") != "Success":',
    '            continue',
    '',
    '        if summary.get("WeinbergCoefficientExtraction") != "Success":',
    '            continue',
    '',
    '        coefficient_file = summary.get("WeinbergCoefficientFile")',
    '',
    '        if not coefficient_file:',
    '            continue',
    '',
    '        c5_path = record.output_dir / coefficient_file',
    '',
    '        if not c5_path.exists():',
    '            continue',
    '',
    '        try:',
    '            rge_summary = run_matched_eft_rge(',
    '                c5_path=c5_path,',
    '                output_dir=record.output_dir,',
    '            )',
    '        except Exception as exc:',
    '            summary["RGEStatus"] = "Failed"',
    '            summary["RGEError"] = str(exc)',
    '            status = 1',
    '            print(',
    '                f"  {record.name}: matched-EFT RGE failed: {exc}"',
    '            )',
    '            continue',
    '',
    '        summary.update(rge_summary)',
    '',
    '        print(',
    '            f"  {record.name}: matched-EFT RGE=Success"',
    '            f" -> {record.output_dir / rge_summary[\'C5BetaFile\']}"',
    '        )',
    '',
    '    aggregate = OUTPUT_DIR / "t3_model_comparison.json"',
    '',
    '    aggregate.write_text(',
    '        json.dumps(',
    '            [record.summary for record in records],',
    '            indent=2,',
    '        ),',
    '        encoding="utf-8",',
    '    )',
    '',
    '    return status',
]
new = "\n".join(new_lines) + "\n"

if old not in text:
    raise RuntimeError(
        "Could not locate finish_runs exactly. No changes were written."
    )

text = text.replace(old, new, 1)
path.write_text(text, encoding="utf-8")

print("Patched:", path.resolve())
print("RGE import present:", import_line.strip() in text)
print("RGE stage present:", "matched-EFT RGE=Success" in text)
