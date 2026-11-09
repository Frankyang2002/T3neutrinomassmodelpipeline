# Workstream A v10-v12 — clean native environment and resumable batch fix

## Verified failure and repair

The first v10-v12 batch successfully executed v11/v12 and their Python
unit tests, then native Wolfram exited **3** before model registration:

```
v10 BLOCKED: optional T3 extensions must be disabled
```

The v10 Wolfram exporter deliberately blocks when any of the following
process environment variables equals the string `1`:

- `T3_EXPERIMENTAL_GENERAL_FERMIONS`
- `T3_EXPERIMENTAL_GENERAL_QUADRATICS`
- `T3_EXPERIMENTAL_GENERAL_QUARTICS`
- `T3_EXPERIMENTAL_FLAVOR_MASSES`

This is **not evidence of a failed `Lam` contraction**: the exporter exited
before loading/building the native model. The actual enabled flag(s) are
unknown from the original console transcript; the revised Wolfram exporter
prints their names if invoked directly with them enabled.

`scripts/RunT3ParallelAScalarBatchV10V12.ps1` now:

1. Runs tests and checks any existing v11 and v12 JSON **read-only**, by
   recomputing both mathematical audits with the current source and comparing
   the full saved records. Missing records are created only once.
2. Captures the values of all four optional flags and clears them **only** in
   this diagnostic PowerShell process before launching each Wolfram child.
3. Preserves the native Wolfram source guard, which still rejects non-baseline
   configurations if run directly with enabled flags.
4. Verifies the native three-model smoke JSON against current source SHA256,
   native model coverage and projector diagonals; validates the all-16 native
   JSON using the existing v10 Python schema/source verifier. If the audited
   result already exists, recomputes and compares it too.
5. Restores **every** original environment variable in a `finally` block,
   including when Wolfram or audit fails. It does not alter user/system environment
   variables, GitHub, or production physics code.

If the native Wolfram check subsequently fails for an **actual** `Lam` or
projector reason, the exporter emits its diagnostic JSON and exits nonzero.
Do not delete or overwrite that evidence; inspect its `models` failures.

## Files

Install these **complete files** into the existing local repository:

1. `scripts/RunT3ParallelAScalarBatchV10V12.ps1` — full replacement.
2. `scripts/RunT3ParallelANativeQuarticProjectorsV10.wl` — full replacement,
   improved explicit environment blocker diagnostics; physics unchanged.
3. `scripts/RunT3ParallelAScalarBatchCheckV10V12.py` — new read-only
   diagnostic checker.
4. `tests/test_t3_parallel_a_scalar_batch_resume_v10_v12.py` — new tests.
5. `Numerical/parallel_t3/wilson/WORKSTREAM_A_V10_V12_RESUME_FIX.md` — this note.

All other v10-v12 modules and prior outputs are retained.

## PowerShell

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
$zip = "$env:USERPROFILE\Downloads\t3_workstream_A_v10_v12_resume_fix.zip"
Expand-Archive -LiteralPath $zip -DestinationPath . -Force

# Recommended after the reported 12/12 v11/v12 tests and saved v11/v12 JSON:
powershell -ExecutionPolicy Bypass -File scripts\RunT3ParallelAScalarBatchV10V12.ps1 -NativeOnly

# Alternatively run the complete, now-resumable batch:
powershell -ExecutionPolicy Bypass -File scripts\RunT3ParallelAScalarBatchV10V12.ps1
```

Expected native output files, created only if missing:

```
output/parallel_t3/wilson/rgbeta_quartic_projectors_native_v10_smoke.json
output/parallel_t3/wilson/rgbeta_quartic_projectors_native_v10.json
output/parallel_t3/wilson/native_quartic_projectors_conditional_v10.json
```

## Certification limits

v11/v12 remained **conditional** and are not altered by this change.
The native v10 symbolic `Lam`/projector evaluation is still **pending** until
run on the user's installed Wolfram kernel. The native component frame,
complex-conjugate projector channels, complete non-portal Matchete coupling
normalization, B/C(-1) scalar mass eigenbases, signed GammaCC→Weyl conversion,
physical ADM, and physical Wilson evolution remain unverified.

EFT order: **UV → integrate out all fermions → integrate out all scalars → SMEFT**.
