from pathlib import Path

path = Path("T3RGECorrectedYukawaPoint.py")
text = path.read_text(encoding="utf-8")

old = '    updated = copy.deepcopy(config)\n    model = updated["t3"]\n'
new = '    updated = copy.deepcopy(config)\n    updated["ordering"] = ordering.upper()\n    model = updated["t3"]\n'

if old not in text:
    raise SystemExit("Patch target not found; no changes made.")

path.write_text(text.replace(old, new, 1), encoding="utf-8")
print("Patched T3RGECorrectedYukawaPoint.py to store ordering in output config.")
