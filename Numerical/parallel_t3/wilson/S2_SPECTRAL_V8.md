# Workstream A v8 — native RGBeta S2 symbolic reduction and basis-invariant portal spectrum

## Source inspection and purpose

Inspected `Frankyang2002/T3neutrinomassmodelpipeline` `main` at
`2bc2e4b516776212645370105cc63495a3c47348` (2026-11-09) and the
upstream `aethomsen/RGBeta` `Package/GroupsAndIndices.m` (blob
`1be60d95d7c6ef6dbde2585c5a811a5284db99d2`). The upstream
`RefineGroupStructures` definitions reduce `tGen[group[S2],A,a,b]` to
**symmetrized fundamental generators**, and `delS2[group,a,i,j]` to
symmetric fundamental delta contractions. These are abstract identities;
there is **no numerical choice of triplet basis** in those rules.

The existing v7 constructs an independently **declared** S2 weight basis but
compares it to **Matchete portal-derived** generators. v7 does not certify a
native RGBeta numerical component intertwiner.

v8 independently checks the **spectrum**, not a fitted intertwiner.
For the Higgs doublet (`j_H=1/2`) and triplet (`j_S=1`), the native
abstract SU(2) generator normalizations imply

```
K = sum_a (sigma_a/2) tensor T_S2^a
J_total^2 = 3/4 + 2 + 2 K
2 tensor 3 = 2 (j=1/2) plus 4 (j=3/2)
eigenvalues(K) = [-1,-1, 1/2,1/2,1/2,1/2]
(K + I_6)(K - I_6/2) = 0
Tr K = 0, Tr K^2 = 3, ||K||_F^2 = 3
```

These checks are invariant under `K -> (I_2 tensor U) K (I_2 tensor U)^dagger`.
The opposite **adjoint sign** has eigenvalues `[-1/2]*4 + [1]*2`, so it
fails the spectral identity even though its Casimir norm is unchanged.

For **each of the 14 triplet sectors**, v8 reads the original
`Matchete\`GetCGTensor` coefficients from the **native v5 export**,
projects the anti-Hermitian `I` phase *inside* `PlusHc` as in v6, removes
the singlet, normalizes to Frobenius norm `sqrt(3)`, and tests **both signs**
against the independent SU(2) total-Casimir eigenvalues. Only *after*
selecting a sign without consulting v6 does it compare to the v6
conditional adjoint matrix. This avoids choosing a Matchete-fitted basis
to define the test. It is nevertheless a basis-independent **relative
portal sign check**, not a complete physical coupling normalization.

## Native RGBeta evidence

The new Wolfram script builds native A(-2) EFT1, obtains the actual loaded
`RGBeta` internal `RefineGroupStructures` function without redefining its
symbols, and evaluates it on **distinct symbolic** S2 generator and `delS2`
indices. It records the actual `InputForm` of the two refinements and whether
they eliminate the original abstract head while retaining fundamental
structures. It does not substitute integers for Einstein indices, which
would risk computing traces instead of components.

The native result can be one of:

- `NativeRGBetaS2SymbolicReductionObserved_NumericalFrameUnspecified`
- `NativeRGBetaS2SymbolicReductionNotExposed_NumericalFrameUnspecified`

**The latter is a diagnostic, not evidence of a broken representation.**
Neither status certifies RGBeta numerical S2 axes. All output flags for
native numeric frame, independent Matchete intertwiner, physical couplings
and Wilson evolution remain `false`.

## Installation and commands (Windows PowerShell, repository root)

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
$zip = "$env:USERPROFILE\Downloads\t3_workstream_A_s2_spectrum_v8.zip"
Expand-Archive -LiteralPath $zip -DestinationPath . -Force

.\.venv\Scripts\python.exe -m unittest discover `
  -s tests -p "test_t3_parallel_a_s2_spectrum_v8.py" -v
if ($LASTEXITCODE -ne 0) { throw "A v8 tests failed" }

wolframscript -file scripts\RunT3ParallelAProbeNativeS2ReductionV8.wl `
  output\parallel_t3\wilson\rgbeta_s2_reduction_native_v8.json
if ($LASTEXITCODE -ne 0) { throw "A v8 native witness failed" }

.\.venv\Scripts\python.exe -m scripts.RunT3ParallelAS2SpectrumV8 `
  --project-root .
if ($LASTEXITCODE -ne 0) { throw "A v8 basis-invariant audit blocked" }
```

New-only, refuse existing output files:

```
output/parallel_t3/wilson/rgbeta_s2_reduction_native_v8.json
output/parallel_t3/wilson/s2_spectrum_conditional_v8.json
```

The v8 Python runner first **revalidates native v5 and v6** through the
existing v7 source-checking audit. It verifies saved-v7 status, each saved
v7 sector's transport residual, and the source hash inventory in the *new*
native Wolfram witness. Any source drift or tampering blocks.

The unit tests contain six standalone numerical scenarios, native witness
status and hash guards, and two all-16 historical registered-CG fixture tests
(requiring the complete repository). They do not simulate a successful
native Wolfram run. The author's container has no installed WolframScript or
complete local project, so the live native witness and complete 14-sector
check must be run on Windows.

## Outputs and status

Expected after a successful live native symbolic reduction:

```
14/14 triplet spectra
S2BasisInvariantSpectrum14of14_NativeRGBetaFormalReductionObserved_NumericalS2FrameStillUnspecifiedConditional
```

The precise maximum numerical residual is measured by the Windows run;
it is **not predefined**. The native symbolic reduction could be
unexposed depending on the installed RGBeta implementation. Both outcomes
are conditional as regards native RGBeta component coordinates.

## Outstanding limitations

- Native RGBeta numerical `delS2` and `tGen[S2]` components have **not**
  been extracted. There is no independently **measured** S2-to-Matchete
  numerical intertwiner.
- Complete scalar-potential action and physical coupling factors uncertified.
- B(-1)/C(-1) require mixed scalar mass diagonalization.
- Signed GammaCC-to-Weyl factor remains separate.
- The Wilson ADM and physical EFT inputs are not certified. Do **not** rerun
  physical Wilson evolution.
- Mandatory EFT order: **UV -> integrate out ALL fermions -> integrate out
  ALL scalars -> SMEFT**.

No existing source model, coupling, RGE, threshold, report or GitHub code is
modified by this batch.

## New files

- `Numerical/parallel_t3/wilson/s2_spectral_v8.py`
- `Numerical/parallel_t3/wilson/S2_SPECTRAL_V8.md`
- `scripts/RunT3ParallelAProbeNativeS2ReductionV8.wl`
- `scripts/RunT3ParallelAS2SpectrumV8.py`
- `tests/test_t3_parallel_a_s2_spectrum_v8.py`
