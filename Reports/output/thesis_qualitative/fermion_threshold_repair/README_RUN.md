# Fermion-threshold runner repair (read-only physics source)

Inspected GitHub `main` at `1b4258b19c8da91afe7f352c2df1b1665de61f26` (2026-11-04). The original `scripts/RunT3GeneralFermionThresholdScan.py` was reconstructed and its Git blob SHA-1 verified as `74dc5fe7ba77e58c3c534d2b3f054986acb62443` **before editing**.

## Scope of replacement

- `scripts/RunT3GeneralFermionThresholdScan.py` — complete replacement script. Supply all four mass-basis tolerance keywords using the verified defaults in `Numerical/running/FinalC5TrajectoryAdapter.extract_final_c5_inputs`; read original or continuation mixed-scalar boundary fractions with consistency and range checks; update `PINNED_HEAD` to current verified main. No physics engines altered.
- `tests/python/test_general_fermion_threshold_scan_repair.py` — new 15-case regression suite for the exact keyword signature, C(-1) continuation provenance, B(-1) original provenance, fail-closed behaviour, and call sites.
- `Reports/output/thesis_qualitative/fermion_threshold_repair/README_RUN.md` — these local instructions.

## PowerShell (local project root)

```powershell
# Check the Git revision; the scan script pins this exact commit.
git rev-parse HEAD
# Expected: 1b4258b19c8da91afe7f352c2df1b1665de61f26

# Copy contents of the archive into your LOCAL project only.
Expand-Archive -LiteralPath "$HOME\Downloads\T3_fermion_threshold_repair.zip" -DestinationPath . -Force

# Test the repair; no matching/RGEs/fits are run by these tests.
python -X utf8 -m pytest -q tests/python/test_general_fermion_threshold_scan_repair.py tests/python/test_general_fermion_threshold_scan.py

# First, new preflight output; do not overwrite the earlier studies.
python -X utf8 -m scripts.RunT3GeneralFermionThresholdScan --preflight --models A_-4,C_-1 --run-id preflight_repaired_01

# Only if preflight is ready, run the TWO-MODEL physics pilot locally.
python -X utf8 -m scripts.RunT3GeneralFermionThresholdScan --models A_-4,C_-1 --run-id pilot_repaired_01

# Inspect statuses; 'ComputedConditional' for all three factors per model
# and saved_baseline_validation.pass=true are prerequisites to all16.
Import-Csv 'output\studies\thesis_qualitative\fermion_threshold\pilot_repaired_01\summary.csv' | Format-Table -AutoSize

# If those pilots succeeded, all 16 with a new output directory.
python -X utf8 -m scripts.RunT3GeneralFermionThresholdScan --models all --run-id all16_repaired_01
```

If your checkout is older than the pinned commit, update it locally only after checking for uncommitted work (`git status --short`), or use `--allow-different-head` deliberately. A commit mismatch is not a physics failure.

The pilot/all16 runs **have not been executed in the file-authoring environment**, because the full local project checkout, symbolic outputs, and Wolfram executable are absent. Only the new 15 synthetic regression tests and Python syntax validation were executed successfully.

## Interpretation

The scan changes physical UV heavy-fermion masses and their common decoupling scale by factors 0.5, 1 and 2, without refitting Yukawa textures. It preserves UV → integrate out all fermions → EFT1 (both scalar multiplets) → integrate out all scalars → SMEFT. Results remain conditional on the one-scale Matchete C5; they are not certified independent sequential hard matching.
