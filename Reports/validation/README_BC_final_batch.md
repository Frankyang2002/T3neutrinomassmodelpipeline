# B(-1)/C(-1) focused leading-order mass computation

**Purpose:** Given explicit **same-scale** masses, two 3x3 complex Yukawa matrices, a 2x2 scalar mixing matrix, quartics, and the electroweak VEV, compute the fixed-order, full-theory one-loop neutrino mass matrices using the archived Matchete `A` coefficient validated by Workstream C v22. It does not claim to finish or validate high-scale-to-low-energy RG evolution.

## Required input

File `configs/bc_fixed_order_input.json`, keyed by `models.B_m1` and `models.C_m1`, each containing:

- `MF_GeV`: 3 positive diagonal Majorana masses, all above `MY`.
- `scalar_eigenmasses_GeV`: `[MX, MY]` with `MY>MX>0`.
- `scalar_U`: 2x2 complex unitary matrix in native Matchete convention.
- `y1`, `y2`: 3x3 Yukawa matrices with row lepton flavour, column positive-Majorana mass index.
- `quartics`: `lambdaT3`, `lambdaT3ExtraQ020020I1`, `lambdaT3ExtraQ020200I1`.
- `vev_GeV`: electroweak vacuum expectation value, in the same convention as the Matchete coefficient (usually 246 GeV; verify independently).

All complex entries must be JSON `[real,imag]` pairs, or real numbers. Do not substitute scalar gauge-basis masses for eigenmasses or mix conventions.

The evaluator implements exact source conventions `C5=A+A.T`, `mnu=-vev^2*C5/2`, masses as Takagi singular values. For each fermion column it includes the 480-term direct coefficient with the native C triplet factor, all three quartics, and complex Yukawa contractions. No running, fitting, parameter substitution or approximate finite scalar bubble is added.

## Run

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$HOME\Downloads\T3_BC_fixed_order_mass_batch.zip" -DestinationPath . -Force
python -u -X utf8 scripts\RunT3BCFixedOrderMasses.py --input configs\bc_fixed_order_input.json
Get-Content Reports\validation\bc_fixed_order_masses.json -Raw -Encoding utf8
```

The JSON input is **not supplied**: physical model parameters cannot be inferred from historical summaries without risking mixing different benchmark datasets. Prepare it from your actual chosen common-scale parameter configuration. Existing `v22` native sources and archived v19/v21 report gates are required; the runner fails rather than silently changing source conventions.

## Scope and stopping conditions

The independent direct/sequential v22 comparison has already passed locally at relative `2.8642e-16` (B) and `5.9415e-16` (C). This verifies agreement of two Matchete-based fixed-order routes under tested conditions, not external first-principles correctness. The present script turns the same native source into an explicit-input one-loop calculator, but **has not been executed with user-local input here**.

For the *separate* UV-running-to-100-GeV objective, source `docs/T3_BC_MATCHING_MEGA_BATCH_V1.md` still states grouped scalar threshold is a restricted bubble, lacks certified full one-loop scalar matching and flavour RG normalization. That is a genuine blocker to certifying the previously reported RG-evolved B/C spectra; it is not fixed by this calculator. One-loop RG effects on an already one-loop C5 are formally higher-order in a fixed-order expansion, but evolving high-scale parameters is an RG-improvement choice and requires a specified accuracy convention.

**No production formula was altered.** No complete physical certification is asserted. Only one intermediate EFT is relevant to the separate RG route.

### Immediate, no-input frozen-reference replay

The runner also supports the existing archived v18 frozen benchmark, without an input JSON:

```powershell
python -u -X utf8 scripts\RunT3BCFixedOrderMasses.py --archived-reference --output Reports\validation\bc_fixed_order_archived_reference.json
```

This reuses frozen reference parameters and the direct amplitude. The report explicitly labels them as archived frozen inputs, **not** UV-evolved user inputs or physically certified masses. For new physical parameter choices, use `--input` instead.
