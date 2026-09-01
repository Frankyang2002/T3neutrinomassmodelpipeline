from pathlib import Path

path = Path("run_t3b_rge_smoke.py")

if not path.exists():
    raise FileNotFoundError(path)

text = path.read_text(encoding="utf-8")

old = '    cleaned = cleaned.replace("Log[", "log(")\n    cleaned = cleaned.replace("]", ")")\n    cleaned = cleaned.replace("^", "**")\n'
new = '    cleaned = cleaned.replace("Conjugate[", "conjugate(")\n    cleaned = cleaned.replace("Log[", "log(")\n    cleaned = cleaned.replace("Sqrt[", "sqrt(")\n    cleaned = cleaned.replace("]", ")")\n    cleaned = cleaned.replace("^", "**")\n'

if old not in text:
    raise RuntimeError(
        "Could not find the current elementary Mathematica conversion block."
    )

text = text.replace(old, new, 1)

old_locals = '                locals={\n                    "log": sp.log,\n                    "conjugate": sp.conjugate,\n                },\n'
new_locals = '                locals={\n                    "log": sp.log,\n                    "sqrt": sp.sqrt,\n                    "conjugate": sp.conjugate,\n                },\n'

if old_locals not in text:
    raise RuntimeError("Could not find SymPy locals block.")

text = text.replace(old_locals, new_locals, 1)

path.write_text(text, encoding="utf-8")

print("Patched:", path.resolve())
print(
    "Generic Conjugate conversion present:",
    'replace("Conjugate[", "conjugate(")' in text,
)
