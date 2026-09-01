from pathlib import Path

path = Path("FlavorMatchedC5.py")
text = path.read_text(encoding="utf-8")

old_factor = "        return sp.factor(value)\n"
new_factor = "        return value\n"

old_matrix = (
    "    return sp.Matrix(\n"
    "        n_lepton,\n"
    "        n_lepton,\n"
    "        entry,\n"
    "    )\n"
)

new_matrix = (
    "    K = sp.MutableDenseMatrix.zeros(n_lepton, n_lepton)\n"
    "\n"
    "    for p in range(n_lepton):\n"
    "        for q in range(p, n_lepton):\n"
    "            value = entry(p, q)\n"
    "            K[p, q] = value\n"
    "            K[q, p] = value\n"
    "\n"
    "    return sp.Matrix(K)\n"
)

if old_factor not in text:
    raise SystemExit("Could not find return sp.factor(value); no changes made.")

if old_matrix not in text:
    raise SystemExit("Could not find matrix-construction block; no changes made.")

text = text.replace(old_factor, new_factor, 1)
text = text.replace(old_matrix, new_matrix, 1)
path.write_text(text, encoding="utf-8")

print("Patched FlavorMatchedC5.py")
print("- removed expensive sp.factor(value)")
print("- computes only the upper triangle")
print("- mirrors it to enforce exact symmetry")
