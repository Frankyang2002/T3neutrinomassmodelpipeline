# All 16 conditional T3 lambdaT3 numerics, including six self-bilinear models (v2)

**Scope:** This batch adds a source-informed **hypothesis test**, not certified physical matching. It reuses the existing, installed `All16PreliminaryNumerics.py` for the ten already supported cases. For A−2, A+0, D−2, D+0, E−2, E+0 it takes the actual source-resolved Matchete fermion-first LLSS boundary, evaluates the independent native `lambdaT3` HHSS source channel, and transports that channel into the source's split-real scalar basis. It computes one fixed-order LLSS-induced EFT1 logarithm and one **joint, same-μS** scalar bubble. It never sequentially integrates scalars.

## Files (new only)

```
RGE/running/intermediate/SixSelfBilinearPreliminary.py
RGE/running/intermediate/README_SixSelfBilinearPreliminary.md
scripts/RunAll16PreliminaryNumericsV2.py
tests/test_six_self_bilinear_preliminary.py
```

Dependencies are the previously installed v1 all-16 runner (`All16PreliminaryNumerics.py`), the B−1 pilot (`B1PreliminaryNumerics.py`), and unchanged repository source matching/quartic/one-scale bubble modules. Nothing in previous workstreams is overwritten. **Do not replace or modify GitHub files.**

## Run from Windows PowerShell in local repository root

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$HOME\Downloads\T3_SixSelf_Conditional_All16_V2.zip" -DestinationPath .
python -m unittest discover -s tests -p 'test_six_self_bilinear_preliminary.py' -v
python -m unittest discover -s tests -p 'test_all16_preliminary_numerics.py' -v
python scripts\RunAll16PreliminaryNumericsV2.py `
  --output-root 'output\general_observable\lo_eft1_mixing\all16_lambdaT3_self_v2'
Get-Content 'output\general_observable\lo_eft1_mixing\all16_lambdaT3_self_v2\summary.json' -Raw
Import-Csv 'output\general_observable\lo_eft1_mixing\all16_lambdaT3_self_v2\summary.csv' | Format-Table -AutoSize
```

Optional: `--model E_-2`, `--check-only` (no writes) or `--require-all` (nonzero exit if even one model is not Weinberg-valid). The output directory must be **new and empty**; choose a fresh name when rerunning.

### Inputs already in the local repository

- F-tree source LLSS for 16 models: `output/general_observable/lo_fermion_boundary/releases/source_021dcc1e8fe8ba0f/<slug>/boundary.json`, with explicit, source-resolved PL channels and uncertified GammaCC sign/factor.
- Native HHSS channels: `output/general_observable/lo_joint_scalar_matching/source_quartic_v4/native_quartic_channels/<slug>.json` plus the source Matchete seeds in `output/full16_symbolic/`.
- Existing illustrative 16-model scalar-only running: `output/general_observable/all16_v1/<slug>/running.json` at μF=100 TeV and μS=1 TeV. No claim these are synchronized to the fermion-first matching benchmark.
- Original F-source split-real metadata: `output/general_observable/fermion_first_self_fix_v2/<slug>/self_mass_basis.json`.
- Native Cartesian triplet rotation: `Lagrangian/RunSixSelfBilinearEigenMatchete.wl` function `T3SBRotate`. The separate `Numerical/experimental/T3SelfBilinearMassBasis.py` uses a **spherical-basis** intertwiner W and must **not** be silently inserted into this native Cartesian coordinate convention.
- The two original mixed-doublet source/mass handoffs B−1/C−1 and other eight preliminary one-channel numerical points are reused from `All16PreliminaryNumerics.py` unchanged.

## Split-real convention and checks

For a native component `S_j` of the self-bilinear multiplet, Matchete's experimental native replacement is

\[
S_j = e^{i\theta}(R_j + i I_j)/\sqrt{2}.
\]

The saved F-source scalar coordinates are ordered `T3Scalar{1,2}Minus`, then `T3Scalar{1,2}Plus`, not the original interleaved complex-field real/imag coordinates. The implementation **hypothesizes** `Minus = I`, `Plus = R` at the F-boundary and enforces a real **positive** RGBeta F-bilinear input so that `theta_F=0`; it never inserts an unsourced phase. The unsplit multiplet retains its canonical original real-component ordering. The resulting real transformation must be exactly orthogonal in dimension 8 (A/D) or 12 (E).

For the *separately declared* potential hypothesis

\[
V_2 = m^2S^\dagger S +\mu_2(S^\dagger)^T S^\dagger +\mu_2^*S^T S,
\]

native Cartesian component mass Hessians are

\[
M^2_{R,I}=\begin{pmatrix}m^2+2\mathrm{Re}\mu_2&2\mathrm{Im}\mu_2\\2\mathrm{Im}\mu_2&m^2-2\mathrm{Re}\mu_2\end{pmatrix},\quad
\mathrm{eig}(M^2)=m^2\pm2|\mu_2|.
\]

**We do not silently equate RGBeta and Matchete bilinear normalizations.** The running's `muT3GeneralS1Sq` or `muT3GeneralS2Sq` is inserted under that explicit hypothesis. This is a parameter-to-native potential mapping that remains unverified. In particular, the real bilinear sign reverses between μF and μS for E−2 and E+0 in the inspected illustrative running outputs. The code **does not swap F-source field labels**; it transforms the full S-scale Hessian and lets the eigenmasses change order correctly.

The original native `lambdaT3` HHSS channel is tensorially rotated with the **same** F-source basis map for both scalar legs, while the other quartic channels remain explicitly excluded (not claimed to vanish). The program verifies the full SU(2) Weinberg tensor in the beta and finite scalar bubble, matching-scale cancellation, and per-source beta recomposition. Only Ward-valid cases receive conditional neutrino-mass singular values in eV per `|N_GammaCC|`. Non-Ward cases retain their diagnostics with `null` masses; missing/invalid provenance is `InputError` and a nonzero CLI exit.

### Outputs

```
output/general_observable/lo_eft1_mixing/all16_lambdaT3_self_v2/
  summary.json
  summary.csv
  A_m4/lambdaT3_one_channel.json
  A_m2/lambdaT3_one_channel.json
  ...
  E_p2/lambdaT3_one_channel.json
```

A complete 16/16 result means **16 conditional, one-channel numerical projections**, **not** 16 certified physical neutrino mass predictions. `physical_neutrino_masses_computed=false`, `source_to_Weyl_factor=null`, and `full_seven_term_EFT1_running_included=false` remain explicit.

## Still genuinely unresolved

1. Signed Matchete `GammaCC` → two-component left-handed Weyl operator normalization.
2. Matchete ↔ RGBeta scalar-quartic `lambdaT3` normalization/sign and native self-bilinear map; source-native scalar component identification has not been independently certified by a matching amplitude.
3. Benchmark synchronization of the 16 fermion boundary, scalar quadratic masses, scalar quartics and Yukawas; the eight ordinary and six self-bilinear RG points are explicitly illustrative.
4. Additional HHSS quartic channels, seven-term complete EFT1 anomalous dimensions, hard F one-loop operators and counterterms.

**EFT sequence is fixed:** UV → integrate out **all fermions** together → EFT1 containing **both scalar multiplets** → integrate out **all scalars together** → SMEFT.
