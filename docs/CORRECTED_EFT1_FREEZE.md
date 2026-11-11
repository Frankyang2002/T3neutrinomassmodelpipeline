# Corrected all-EFT1-frozen calculation — LLSS scope

## Source

The original 16-model generator is `scripts/RunT3GeneralSequentialFitAblation.py`.
Its old implementation uses `C5_W_off = C5_one_scale - C5_leading_log`.
This is a conditional anchored subtraction, not an independently frozen Wilson tensor.

The new local runner is `scripts/RunT3CorrectedEFT1Freeze.py`. It reuses original
Yukawa fits, threshold scales, UV/EFT1 renormalisable RGE code, F-tree LLSS
source tensors, simultaneous two-scalar bubble, and SM+Weinberg low-energy code.
It changes the Wilson matching calculation to

- `W=0`: `C5 = C5_scalar_finite_bubble`;
- `W=1`: `C5 = C5_scalar_finite_bubble + C5_EFT1_leading_log`.

The scalar bubble and log are evaluated by the existing independent
`GeneralT3SequentialLeadingLog.response` using the same F-tree source and
scalar basis. The one-scale Matchete coefficient is recorded as a diagnostic
only; it is never used to create the W=0 result.

**The independent LLSS sector is not the complete physical Weinberg matching:**
fermion hard matching, other intermediate EFT operators, and absolute basis /
spinor normalisation are not certified. This corrected baseline will generally
**differ from the originally locked `1111` benchmark**. It is not permissible to
compare the new EFT1-off results against the old anchored baseline.

## Run on Windows PowerShell

```powershell
python -X utf8 -m scripts.RunT3CorrectedEFT1Freeze --source-run thesis16_merged_v1 --run-id corrected_llss16_v1 --check-only
python -X utf8 -m scripts.RunT3CorrectedEFT1Freeze --source-run thesis16_merged_v1 --run-id corrected_llss16_v1
```

For a fast two-model test first:

```powershell
python -X utf8 -m scripts.RunT3CorrectedEFT1Freeze --source-run thesis16_merged_v1 --run-id corrected_llss_pilot_v1 --models A_-4,D_-2 --check-only
python -X utf8 -m scripts.RunT3CorrectedEFT1Freeze --source-run thesis16_merged_v1 --run-id corrected_llss_pilot_v1 --models A_-4,D_-2
```

New outputs:
`output/studies/general_scalar_extended/corrected_llss_eft1_freeze/<run-id>/`

Cases are 4-bit U/E/W/S. For the intended all-EFT1 running-off comparison:
`1111` versus `1001` **within this new result set**. UV-off `0111` and SMEFT-off `1110` must also come from this new set. Compare exactly the same Yukawas.

The corrected runner refuses `--resume` and will not touch the original locked
`Reports/output/thesis_qualitative/representation_running/` data.

Do not use the old plotting script directly against the corrected `all16_cases.csv`;
it requires a different long CSV format. The corrected results must be reviewed
for finite values, matching-scale sensitivity and Weinberg projection residuals
before plotting for the thesis.

## Plot the entire corrected 16-model result

After all 16 corrected models complete successfully:

```powershell
python -X utf8 -m scripts.PlotT3CorrectedEFT1Freeze --run-id corrected_llss16_v1
```

Graph and CSV output:
`Reports/output/thesis_qualitative/corrected_llss_eft1_freeze/corrected_llss16_v1/`

Outputs: `masssq_all_eft1_frozen.png`, `mixing_all_eft1_frozen.png`,
`running_stage_comparisons.csv`, `provenance.json`. The plotter rejects incomplete
or anchored one-scale runs. These are exploratory LLSS-sector figures only.
