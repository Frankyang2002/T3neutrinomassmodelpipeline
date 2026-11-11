# Worker 1 — A(-4) conditional tensor sensitivity within simultaneous IVP

## Scope / verification

This additional independent batch reuses the **selected** `scripts/RunT3SimultaneousEFT1A4.py` simultaneous DOP853 coupled solver. It does not change production files or any saved historical output. This is a **conditional diagnostic**, not a physically justified competing tensor convention.

- Variant `source`: current source-derived HHSS at all EFT1 RHS evaluations.
- Variant `batch3b`: run with the *same UV inputs* and the *same running renormalisable equations*, but substitute the unverified Batch3B HHSS `CandidateQuarticContext.full_reference(state)[:4,:4,4:,4:]` **only within the EFT1 seven-term Wilson RHS**.
- The source Matchete HHSS evaluator remains in the scalar matching/bubble and other non-HHSS tensor components. Therefore the comparison IS NOT an end-to-end alternate physical potential.
- Exact input hashes for the cached candidate response, quartic seed and RGBeta registry are checked by the existing `CandidateQuarticContext.load` interface. Running scale: MF=1e5 GeV -> MS=1e3 GeV, with both scalars retained; UV start=1e7 GeV and SMEFT low=100 GeV.
- The full sixteen U/E/W/S cases are recomputed in fresh folders separately. No fitted scale/sign is imposed.
- Compares low-energy/surface C5 and complex-proportional shape residual, mass ratios, mixing sin²θ, squared mass splittings, stage-off differences. Stage-off norms are *not additive stage-running fractions*.
- The repository runner itself identifies absent hard one-loop fermion matching, provisional normalization and incomplete higher operator classes.

## Install / tests / preflight / run

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
$zip = "$HOME\Downloads\T3_Worker_1_A4_Conditional_IVP.zip"
Expand-Archive -LiteralPath $zip -DestinationPath . -Force
python -X utf8 -m unittest discover -s tests -p "test_t3_worker1_a4_conditional_ivp.py" -v
python -X utf8 -m scripts.RunT3Worker1A4ConditionalIVP --both --check-only `
    --source-run thesis16_merged_v1
python -X utf8 -m scripts.RunT3Worker1A4ConditionalIVP --both `
    --source-run thesis16_merged_v1 `
    --source-id worker1_A4_source_v1 `
    --candidate-id worker1_A4_Batch3B_conditional_v1 `
    --out-dir "output\audits\worker1_A4_conditional_IVP_v1"
Get-Content "output\audits\worker1_A4_conditional_IVP_v1\REPORT.md" -Encoding utf8
```

If run IDs already exist, choose entirely new names for `--source-id`, `--candidate-id`, and `--out-dir`. Never delete or overwrite historical outputs.

## Project destinations (two new files + documentation)

- `scripts/RunT3Worker1A4ConditionalIVP.py`
- `tests/test_t3_worker1_a4_conditional_ivp.py`
- `docs/WORKER_1_A4_CONDITIONAL_IVP.md`

Generated (new only):

- `output/studies/general_scalar_extended/experimental_simultaneous_eft1_seven_term/worker1_A4_source_v1/A_m4/cases.json`
- `output/studies/general_scalar_extended/experimental_simultaneous_eft1_seven_term/worker1_A4_source_v1/WORKER1_MODE.json`
- `output/studies/general_scalar_extended/experimental_simultaneous_eft1_seven_term/worker1_A4_Batch3B_conditional_v1/A_m4/cases.json`
- `output/studies/general_scalar_extended/experimental_simultaneous_eft1_seven_term/worker1_A4_Batch3B_conditional_v1/WORKER1_MODE.json`
- `output/audits/worker1_A4_conditional_IVP_v1/report.json`
- `output/audits/worker1_A4_conditional_IVP_v1/REPORT.md`

## Required existing local files

- `output/audits/batch3b_quartic_fixed_20261111_234059/A_-4/quartic_response.npz` and `report.json`
- `RGE/general/generated/general_rgbeta_quartic_registry.json`
- Existing A(-4) fitted source run `output/studies/general_scalar_extended/thesis_general_all16/thesis16_merged_v1/A_m4/` (exact location is resolved by `RunT3GeneralAll16FromReference.OUTPUT`)
- Existing original Matchete boundary release `output/general_observable/lo_fermion_boundary/releases/source_021dcc1e8fe8ba0f`, quartic channels and seeds; exact completeness is checked by `--check-only`.

**Acceptance:** run output statuses `Computed`, 16-case input inventory, source provenance match, unchanged UV source ID, valid C5 shape and normalized matrix comparison. Outcome still does *not* certify the physical scalar alternative. If the Stage/Case baseline comparison differs in non-HHSS sectors, no physical conclusion is warranted.
