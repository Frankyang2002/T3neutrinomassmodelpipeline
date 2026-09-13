from pathlib import Path

ROOT = Path(__file__).resolve().parent

REPLACEMENTS = {
    "[fresh kernel] EFT1 running heavy insertion loaded. Terms: ":
        "[fresh kernel] EFT1 fixed-order running insertion loaded. Terms: ",
    "[fresh kernel] Added EFT1 leading-log heavy insertion to [C] only.":
        "[fresh kernel] Applied fixed-order EFT1 threshold insertion in [C]; "
        "direct Weinberg running is carried separately.",
}

changed = []

for path in ROOT.rglob("*.wl"):
    text = path.read_text(encoding="utf-8")
    new = text
    for old, replacement in REPLACEMENTS.items():
        new = new.replace(old, replacement)
    if new != text:
        path.write_text(new, encoding="utf-8")
        changed.append(path.relative_to(ROOT))

print("Updated:")
for path in changed:
    print(f"  {path}")

if not changed:
    print("  No .wl file contained the old log messages.")
