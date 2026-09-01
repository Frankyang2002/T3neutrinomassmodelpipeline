from pathlib import Path

path = Path("pipeline.py")

if not path.exists():
    raise FileNotFoundError(path)

text = path.read_text(encoding="utf-8")

pairs = []

old = "        if numerical_config is not None:\n            try:\n                numerical_summary = run_numerical_pipeline_stage(\n"
new = "        if numerical_config is not None:\n            print(f\"  {record.name}: starting numerical RGE stage...\", flush=True)\n            try:\n                numerical_summary = run_numerical_pipeline_stage(\n"
pairs.append((old, new))

old = "            mass_matrix_path = (\n                record.output_dir\n                / numerical_summary[\"NeutrinoMassMatrixLowScaleFile\"]\n            )\n\n            try:\n                observable_summary = run_neutrino_observables_stage(\n"
new = "            print(f\"  {record.name}: numerical RGE calculation finished.\", flush=True)\n\n            mass_matrix_path = (\n                record.output_dir\n                / numerical_summary[\"NeutrinoMassMatrixLowScaleFile\"]\n            )\n\n            print(f\"  {record.name}: starting neutrino observables...\", flush=True)\n            try:\n                observable_summary = run_neutrino_observables_stage(\n"
pairs.append((old, new))

old = "        try:\n            mass_summary = run_neutrino_mass_stage(\n"
new = "        print(f\"  {record.name}: starting symbolic neutrino mass stage...\", flush=True)\n        try:\n            mass_summary = run_neutrino_mass_stage(\n"
pairs.append((old, new))

old = "        try:\n            flavor_summary = run_flavor_matched_rge(\n"
new = "        print(f\"  {record.name}: starting symbolic full-flavor RGE stage...\", flush=True)\n        try:\n            flavor_summary = run_flavor_matched_rge(\n"
pairs.append((old, new))

old = "        try:\n            rge_summary = run_matched_eft_rge(\n"
new = "        print(f\"  {record.name}: starting matched-EFT RGE stage...\", flush=True)\n        try:\n            rge_summary = run_matched_eft_rge(\n"
pairs.append((old, new))

changed = 0
for old, new in pairs:
    if old in text and new not in text:
        text = text.replace(old, new, 1)
        changed += 1

path.write_text(text, encoding="utf-8")

print("Patched:", path.resolve())
print("Progress markers added:", changed)
