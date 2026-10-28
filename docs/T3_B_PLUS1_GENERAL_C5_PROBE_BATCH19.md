# Batch 19 — First non-B(-1)/C(-1) general-scalar C5 source probe: B(+1)

## Scope and verified project interfaces

This experiment targets **B(+1)**, a distinct ordinary T3 model with one
additional complex quartic and **zero additional quadratic mixing operators**
according to `RGE/general/generated/general_rgbeta_quartic_registry.json`:

- RGBeta: `lambdaT3GeneralQ022000I1`;
- Matchete: `lambdaT3ExtraQ022000I1`;
- fields: `Hdag Hdag S1 S1 + h.c.`.

The existing `Lagrangian/interactions/GeneralScalarOperators.wl` enables that
quartic with environment variable `T3_EXPERIMENTAL_GENERAL_QUARTICS=1`, and
`Lagrangian/RunModel.wl` exports an ordered C5 in `c5_coefficient.txt` when
Matchete matching and operator extraction succeed.

This batch executes the **same existing** Matchete driver twice at EFT order 5,
one loop, with **all heavy particles removed at one common step**:

1. Restricted Lagrangian (`T3_EXPERIMENTAL_GENERAL_QUARTICS=0`).
2. Extended Lagrangian (`T3_EXPERIMENTAL_GENERAL_QUARTICS=1`).

It fails closed if the extended interaction is not registered or one of the
matching/extraction statuses is not `Success`. It records source SHA-256 hashes.
The second Wolfram script reads the two actual ordered C5 expressions and
compares them at three *numerical one-generation flavour spurion* points,
including the exact zero-extra-quartic reduction. It probes one nonzero extra
quartic value and records any numerical response. A very large dimensional
regulator variable is used for a finite-limit **sample**, not an MSbar-renormalised
C5 calculation. No neutrino mass fit is performed.

The resulting status distinguishes **source occurrence** from **numeric
sensitivity**. If the extra coupling is absent from the ordered C5 output,
that observation does not imply it has no effect through the UV/EFT1 RGEs.

## Files (new, full files, no replacements)

- `scripts/ProbeT3BPlus1GeneralC5.py`: new-only run coordinator and structural checks.
- `scripts/CompareT3BPlus1OneScaleC5.wl`: native Mathematica numerical ordered-C5 check.
- `tests/python/test_t3_b_plus1_general_c5_probe.py`: offline fail-closed regressions.
- `docs/T3_B_PLUS1_GENERAL_C5_PROBE_BATCH19.md`: this note.

## PowerShell, from project root

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "C:\path\to\T3_B_Plus1_General_C5_Probe_Batch19.zip" -DestinationPath . -Force
python -m pytest -q tests\python\test_t3_b_plus1_general_c5_probe.py
python -m scripts.ProbeT3BPlus1GeneralC5 --run-id b_p1_c5_probe_v1
Get-Content -Encoding UTF8 output\studies\general_scalar_extended\matching_probes\B_p1\b_p1_c5_probe_v1\summary.md
```

The native Matchete calculation can be expensive. If it fails, inspect
`baseline.log`, `general.log` or `numeric_audit.log` in the run directory.
A failed run retains its directory and must **not** be overwritten; rerun
with a new ID after correcting the specific problem.

Output:

`output/studies/general_scalar_extended/matching_probes/B_p1/<run-id>/`

with `baseline/`, `general/`, source certificates, log files, numeric audit,
`report.json` and `summary.md`. Historical data are read-only.

## Scientific limitations

- **Implemented:** symbolic restricted and extended Matchete one-scale matching,
  registered new quartic check, numeric ordered-C5 source comparison.
- **Not implemented:** physical Matchete–RGBeta quartic normalisation, full
  sequential EFT1 Wilson matching, general-scalar mixing for other models,
  new B(+1) UV/EFT1 parameter-to-C5 adapter, oscillation refitting.
- **Local verification pending:** the live Wolfram matching and numeric comparison
  must run successfully in the user's native Windows Matchete installation.

The B(+1) result is the first **source gate** toward extending the original
three numerical studies to the other 14 models, **not yet an accepted new
model-study fit**. The next numerical step depends on whether the new quartic
is present in, or drops out of, the actual one-loop C5 source.

## Batch 19 CLI correction and interrupted-run recovery

`RunModel.wl` identifies B(+1) in `comparison_summary.json` as
`ModelClass: T3-d2-d2-F1`, `Alpha: 1`, not the short CLI label `B`.
The probe now checks the actual representation and hypercharge values.
No physics equations or Mathematica matching source files were modified.

To reuse the **existing baseline match from the original failed run**, replace
the files in this archive and execute:

```powershell
python -m pytest -q tests\python\test_t3_b_plus1_general_c5_probe.py
python -m scripts.ProbeT3BPlus1GeneralC5 --run-id b_p1_c5_probe_v1 --resume
```

`--resume` validates the existing baseline output and only runs missing stages.
The original version did not save source SHA-256s at the start. For that
original interrupted run, the resume checks source modification times and
records `baseline_source_hash_certified=false`; timestamps are **not**
cryptographic proof. New runs save a `source_fingerprints.json` before matching
and verify it on resume. To avoid this limitation altogether, run with a fresh
`--run-id` without `--resume` (which repeats the Matchete baseline).

## Batch 19 correction 2: local Matchete model-class label

The user's live output reports `ModelClass: "T3-B"` while the saved
repository B(+1) reference reports `"T3-d2-d2-F1"`. The prior validator
accepted only the latter and incorrectly rejected two successful baseline
runs. The corrected validator accepts **only those two** class-label spellings,
while requiring **all** of `Alpha=1`, `Scalar1SU2=2`, `Scalar2SU2=2`,
`FermionSU2=1`, `Scalar1Hypercharge="1/2"`,
`Scalar2Hypercharge="3/2"`, and `FermionHypercharge="1"`. A class
label alone is never treated as sufficient evidence of B(+1).

The accepted label is also recorded in `report.json`. No physics code,
matching source, or coupling convention is changed.

Use the latest full-file replacement ZIP and **resume the already completed
baseline in `b_p1_c5_probe_v2`**, which has source fingerprints recorded
before its first native Matchete run:

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath 'C:\path\to\T3_B_Plus1_General_C5_Probe_Batch19_Fix2.zip' -DestinationPath . -Force
python -m pytest -q tests\python\test_t3_b_plus1_general_c5_probe.py
python -m scripts.ProbeT3BPlus1GeneralC5 --run-id b_p1_c5_probe_v2 --resume
Get-Content -Encoding UTF8 output\studies\general_scalar_extended\matching_probes\B_p1\b_p1_c5_probe_v2\summary.md
```

A legacy `b_p1_c5_probe_v1 --resume` is also supported if its timestamp
provenance check passes. Prefer `v2` because it has the full SHA-256 manifest.

The offline tests are mock-based; real Matchete and native numerical C5 checks
remain to be established by this local run.

## Batch 19 correction 3: malformed native Wolfram numeric JSON

The user's B(+1) `b_p1_c5_probe_v2 --resume` completed the restricted and
extended Matchete stages but failed while parsing `native_numeric_audit.json`
with Python `JSONDecodeError` (reported line 41, column 33). Without the actual
file bytes, the precise offending token is **not verified**. High-precision
Wolfram `*^`-style numeric notation is one plausible cause; the updated
comparison writer removes that possible incompatibility without changing any
of the underlying numerical C5 expressions or the high-precision comparisons.

- `scripts/CompareT3BPlus1OneScaleC5.wl` (full replacement): uses finite
  machine-precision real numbers for JSON diagnostics, normalises any Wolfram
  scientific `*^` exponent formatting to JSON `e`, and checks an independent
  Wolfram RawJSON import before writing the result. High-precision arithmetic
  for the three C5 comparisons is unchanged.
- `scripts/ProbeT3BPlus1GeneralC5.py` (full replacement): permits an updated
  **comparison script only** when resuming a source-fingerprinted run. Every
  production matching source must still have the original SHA-256. It reruns
  the *comparison only*, not the expensive Matchete matching stages, and
  archives previous `native_numeric_audit.json` and `numeric_audit.log` as
  `native_numeric_audit.before_json_fix.json` and
  `numeric_audit.before_json_fix.log`. Malformed JSON is never accepted as
  a successful matching result.
- `tests/python/test_t3_b_plus1_general_c5_probe.py`: regressions for corrupted
  numeric JSON, comparison-script-only changes, and rejecting changed matching
  source code. These are mocked tests, not a native Wolfram verification.

Run from the project root, after installing **this** full replacement archive:

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath 'C:\path\to\T3_B_Plus1_General_C5_Probe_Batch19_Fix3.zip' -DestinationPath . -Force
python -m pytest -q tests\python\test_t3_b_plus1_general_c5_probe.py
python -m scripts.ProbeT3BPlus1GeneralC5 --run-id b_p1_c5_probe_v2 --resume
Get-Content -Encoding UTF8 output\studies\general_scalar_extended\matching_probes\B_p1\b_p1_c5_probe_v2\summary.md
```

The new `report.json` distinguishes `baseline_source_hash_certified` (matching
sources), `full_original_source_manifest_unchanged` (false after the comparator
update), and `comparison_script_updated_since_match`. The previous malformed
JSON and log stay available for exact source inspection. If the corrected
report is still malformed, the runner fails and retains the new native audit.
Do not treat this as proof of a new B(+1) physical EFT prediction or a refit.
