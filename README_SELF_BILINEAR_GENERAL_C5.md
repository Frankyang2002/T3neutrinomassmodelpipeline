# General T3 self-bilinear mass-basis matching (local experimental implementation)

## Scope

This package replaces the matching launcher and 16-model fitter **in full**, and
adds a single new Wolfram one-scale matching runner. The objective is to retain
all general Yukawas and scalar quartics while matching A(-2), A(0), D(-2),
D(0), E(-2), E(0) with nonzero scalar self-bilinears.

What is **implemented**:

- Full-general `BuildT3Lagrangian` (all generated extra Yukawa/quartic operators).
- Extract exactly the terms proportional to `lambdaT3ExtraM*` and absorb their
  effect into the split mass eigenvalues; **do not zero the bilinear**.
- Replace the hypercharge-zero complex scalar by **two self-conjugate fields**
  and preserve the conjugate SU(2) index orientation in Yukawa/quartic terms.
- Matchete `CheckLagrangian` after the rotation, then source-derived one-loop
  one-scale matching via `RunT3Matching` and `ExtractWeinbergCoefficient`.
- Return `T31PlusMass/T31MinusMass` or `T32PlusMass/T32MinusMass`, which the
  existing numeric evaluator obtains from the **running** bilinear.
- Distinguish the new matching scheme in provenance; refuse prior self-bilinear
  coefficients and block unsupported complex phases.
- Leave B(-1)/C(-1) to their separate existing mixed-doublet one-scale code.

**Critical limitation:** This is a candidate mass-basis construction. It has
NOT been executed or verified in Matchete in this environment (no Mathematica).
It explicitly fails if the *rotated* Lagrangian does not pass Matchete checks,
if matching fails, or if the expected masses/Yukawa labels are absent from C5.
Real-representation tensor normalization and Matchete-to-RGBeta physical
coupling identification remain provisional. It is NOT full sequential EFT
matching. Do not report the resulting fit as certified general-neutrino physics.

### Convention used

For the positive-real holomorphic coupling in the source sign convention,

    L_mass = -mS^2 S* S + mu (S*S* + SS),  mu > 0
    S = (S_minus + i S_plus)/sqrt(2)
    m_minus^2 = mS^2 - 2 mu
    m_plus^2  = mS^2 + 2 mu

For a real SU(2) irrep, the true conjugation relation carries its invariant
metric. The code retains `Bar[real-field]` on formerly conjugated scalar legs.
Both the SU(2) real-field convention and signs still require live validation.
The current numerical pilot sets the bilinear initially **positive real**.
A genuinely complex running bilinear phase is explicitly rejected until a
phase-correct rotation is implemented.

## Files and intended project paths

- `Lagrangian/RunT3GeneralSelfMassBasis.wl` (NEW)
- `scripts/RunT3GeneralMatcheteMatching.py` (FULL REPLACEMENT)
- `scripts/RunT3GeneralAll16FromReference.py` (FULL REPLACEMENT)
- `Numerical/matching/GeneralT3FlavourC5.py` (FULL REPLACEMENT)
- `tests/python/test_thesis_general16.py` (FULL REPLACEMENT)
- `tests/python/test_thesis_general_c5.py` (UNCHANGED from previous ZIP)
- `tests/python/test_thesis_self_mass_basis.py` (NEW)

The standalone `scripts/RunT3BCGeneralFromReference.py` from the previous
installation remains untouched; B/C direct fit outputs also remain untouched.

## Install and test (PowerShell)

From the project root:

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive "$HOME\Downloads\T3_General_SelfBilinear_MassBasis_C5.zip" -DestinationPath . -Force
python -X utf8 -m pytest -q tests/python/test_thesis_general_c5.py tests/python/test_thesis_general16.py tests/python/test_thesis_self_mass_basis.py
```

The offline tests check model dispatch, real scalar splitting, source guards,
flavour assembly, and campaign configuration. They **do not** validate the
Wolfram field-redefinition algebra or evaluate a live one-loop coefficient.

## Step 1: try ONLY the A(-2) general mass-basis Matchete matching

```powershell
python -X utf8 -m scripts.RunT3GeneralMatcheteMatching --model A_-2
```

Check the log:

```powershell
Get-Content "output\general_matchete_matching\A_m2\general_matching.log" -Tail 130
```

Expected success marker:

    GENERAL SELF-BILINEAR MASS-BASIS C5: SUCCESS

Expected output (ONLY if the live Wolfram calculation passes):

- `output/general_matchete_matching/A_m2/c5_coefficient.txt`
- `output/general_matchete_matching/A_m2/self_mass_basis_metadata.json`
- `output/general_matchete_matching/A_m2/general_matching_status.json`

If the command returns `Failed`, **do not** run optimisation; paste the last
100–150 log lines. We will fix the actual Matchete diagnostic rather than
silently using a restricted coefficient or zeroing the bilinear.

## Step 2: after success, two-point Python/RGE/C5 smoke test

```powershell
python -X utf8 -m scripts.RunT3GeneralAll16FromReference `
    --run-id general_c5_A_m2_massbasis_test `
    --models A_-2 `
    --sobol-points 2 `
    --max-nfev 8 `
    --no-local-refinement
```

This uses the canonical general RGE exports and the saved new general C5
expression. The bilinear is nonzero (default `--mu2-fraction 0.01`).

## Step 3: only after numerical smoke test passes

```powershell
python -X utf8 -m scripts.RunT3GeneralAll16FromReference `
    --run-id thesis16_self_massbasis_v1 `
    --generate-general-c5 `
    --export-bc-betas `
    --reuse-b-run thesis_b_m1_128 `
    --sobol-points 128 `
    --max-nfev 60
```

Use a **new run ID**. The previous `thesis16_v1` and any earlier matching
outputs are not comparable. If reusing C(-1), first ensure its completed direct
run is located where `--reuse-c-run` expects.

All new numerical output stays below:

    output/studies/general_scalar_extended/thesis_general_all16/<run-id>/

## Verification boundary

- Verified offline: syntax-checked Python, plus offline pytest checks.
- Verified from GitHub source: full-general flags, six self-bilinear models,
  Matchete field registration API used by `T3Fields.wl`, and the available
  matching and Weinberg extraction functions.
- NOT verified: Matchete accepts the rotated self-conjugate SU(2) triplet,
  original invariant tensors normalize identically after rotation, full
  one-loop coefficient is generated and evaluated, or optimisation converges.
