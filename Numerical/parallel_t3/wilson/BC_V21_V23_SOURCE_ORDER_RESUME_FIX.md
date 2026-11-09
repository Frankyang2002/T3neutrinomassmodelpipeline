# Workstream A v21/v22 source-order repair (v23 preflight)

## Observed failure and exact cause

The user's Windows run passed 23/24 tests and failed inside
`bc_mass_eigen_native_bridge_v21.validate_special_matrix()` with
`B_-1: invalid v19 baseline six-of-seven conditional map`.

The v19 6-by-7 matrix uses the *historical Stage22 Matchete* source-coupling
order recorded as `matchete_seven_source_symbols`. The v18 live Matchete
BCEXTRA export independently enumerates by verified scalar-structure order.
These are the same seven named CG operators, but the lists have different
orderings. The old v21 code incorrectly required direct list equality.
This is a representation-of-data issue, not a proven scalar-mass or coupling
normalisation failure.

For B(-1)/C(-1), the historical Stage22 order begins with
`lambdaT3ExtraQ110101I1`, `lambdaT3ExtraQ110101I2`, followed by
`lambdaT3ExtraQ000112I1`, `lambdaT3ExtraQ000202I1`,
`lambdaT3ExtraQ001201I1`, `lambdaT3ExtraQ020020I1`,
`lambdaT3ExtraQ020200I1`. The new native export may list those two 110101
operators after the other five; the v23 preflight reports **the actual**
permutation rather than assuming it.

## Repair

Let `old_names` be the v19 historical source coupling order and `new_names`
the v18 native CG order. Reject missing, duplicated, unfamiliar or invented
couplings; require the *exact same seven-symbol set*. Define

`indices = [old_names.index(symbol) for symbol in new_names]`.

If the v19 map is `r = B_old @ q_old`, the same map in native coordinates is

`B_new = B_old[:, indices]`, `N_new = N_old[indices]`,

where `N_old` is the normalized right null vector of `B_old`. This preserves
`rank(B_new)=6`, `||B_new @ N_new||<1e-9`, `||N_new||=1`, and the source-only
six-of-seven distinction. The v21 operator-transport code uses `N_new` with
the native CG record order. The original v19 witness, matrices, CGs and local
repository physics remain unmodified.

The optional seventh Matchete direction is still tested to give a nonzero
source action after rotation and to preserve its H.c. term. The v22 read-only
census gate remains unchanged. The existing Windows/POSIX SHA256 fix is
preserved, including all seven v19 input-file byte-hash checks.

**No native RGBeta numerical S2 frame, physical scalar matching, complete
physical coupling dictionary, GammaCC-Weyl sign or ADM certification is
claimed. No Wilson running or Match[] is performed.**

## Files in replacement ZIP

- `Numerical/parallel_t3/wilson/bc_mass_eigen_native_bridge_v21.py`
  complete replacement (v19/native source-column matching).
- `tests/test_t3_parallel_a_bc_mass_bridge_v21_v22.py`
  complete replacement with reordered B/C regression, end-to-end synthetic
  v21-to-v22 bridge, and negative controls.
- `scripts/RunT3ParallelABCMassSubspaceOrderV23.py`
  new read-only hash-checked native/historical column-order preflight.
- `scripts/RunT3ParallelABCMassEigenBatchV21V22.ps1`
  complete replacement, runs the new preflight before the v21 audit.
- `Numerical/parallel_t3/wilson/BC_V21_V23_SOURCE_ORDER_RESUME_FIX.md`
  this documentation.

All are within the existing project directory structure. No GitHub write.

## PowerShell

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
$zip = "$env:USERPROFILE\Downloads\t3_workstream_A_bc_v21_v23_source_order_fix.zip"
Expand-Archive -LiteralPath $zip -DestinationPath . -Force
powershell -ExecutionPolicy Bypass -File scripts\RunT3ParallelABCMassEigenBatchV21V22.ps1
if ($LASTEXITCODE -ne 0) { throw "v21/v22 source-order resume failed" }
```

Optional read-only per-model source-axis report:

```powershell
.\.venv\Scripts\python.exe -m scripts.RunT3ParallelABCMassSubspaceOrderV23 --project-root .
```

No v18/v19/v20 JSON regeneration is needed. Existing saved outputs are
revalidated; the v21/v22 scripts refuse altered existing outputs in resume mode.

## Local validation, not installed-Windows certification

- v21/v22 suite: 31 tests executed, 30 passed, 1 skipped (live repository fixture).
- v18/v19/v20 suite: 18 passed.
- Synthetic full v21-to-v22 API test covers a nontrivial seven-column
  permutation with v19 Windows-style backslash SHA256 keys.
- Native Windows v21/v22 mass-action integration: pending this replacement run.

EFT workflow remains UV -> integrate out **ALL fermions** -> integrate out
**ALL heavy scalars** -> SMEFT. No physical Wilson evolution authorized.
