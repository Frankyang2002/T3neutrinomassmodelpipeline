from pathlib import Path

path = Path("pipeline.py")

if not path.exists():
    raise FileNotFoundError(path)

text = path.read_text(encoding="utf-8")

old = '        if summary.get("WeinbergCoefficientExtraction") != "Success":\n'
new = '        if summary.get("WeinbergExtractionStatus") != "Success":\n'

if old not in text:
    raise RuntimeError(
        "Could not find the incorrect WeinbergCoefficientExtraction key."
    )

text = text.replace(old, new, 1)
path.write_text(text, encoding="utf-8")

print("Patched:", path.resolve())
print(
    "Correct extraction-status key present:",
    'summary.get("WeinbergExtractionStatus")' in text,
)
