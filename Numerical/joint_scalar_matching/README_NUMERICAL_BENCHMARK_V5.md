# T3 v5: first source-backed conditional numerical scalar benchmark

## What is actually done

At `mu_F=100000 GeV` all heavy fermions are removed, both heavy scalar
multiplets remain in EFT1, and **both scalar multiplets are integrated at the
same `mu_S=1000 GeV`**. No intermediate scalar EFT is introduced.

This batch reads **existing real repository data**, not a synthetic point:

- `output/general_observable/all16_v1/<model>/running.json`, exactly the
  `eft1_at_scalar_threshold` snapshot at `mu_S=1000 GeV`;
- `output/general_observable/lo_joint_scalar_matching/source_quartic_v4/native_quartic_channels/<model>.json`
  with its source Matchete seed SHA-256 verified;
- `output/general_observable/lo_fermion_boundary/releases/source_021dcc1e8fe8ba0f/<model>/boundary.json`;
- existing source-pivot and full-SU(2) outputs for independent consistency.

All 16 models receive an input audit; **eight ordinary scalar-basis models**
(A_-4, B_-3, B_+1, C_-3, C_+1, D_+2, E_-4, E_+2) can receive a *conditional*
numerical `3x3` complex LLSS-sector `C5/N_GammaCC` using the existing
unbroken-multiplet benchmark `mS1Sq(muS)`, `mS2Sq(muS)` and v4 HHSS tensors.

The remaining **eight special-basis models** (A_-2, A_+0, B_-1, C_-1,
D_-2, D_+0, E_-2, E_+0) **are not evaluated numerically** because their
muS mixed / split-real quadratic masses and independently established
original-to-fermion-boundary field transformations are not supplied. The
renormalisable running of `mS1Sq/mS2Sq` alone is **not** a valid full mass
basis for them.

## CRITICAL: experimental Matchete-to-RGBeta coupling aliases

Matchete's scalar potential uses invariant basis names such as
`lambdaH2Inv2`; RGBeta's numerical EFT1 state has names such as
`lambdaH2Adj`. Their **physical coefficient normalization/relative basis
mapping has not been certified**. Only to obtain rough numerical diagnostics,
the runner offers this explicitly chosen *candidate* map:

| Matchete native source coupling | Provisional RGBeta numeric parameter |
|---|---|
| `lambdaH1Inv1` | `lambdaH1` |
| `lambdaH1Inv2` | `lambdaH1Adj` |
| `lambdaH2Inv1` | `lambdaH2` |
| `lambdaH2Inv2` | `lambdaH2Adj` |
| `lambdaT3` | `lambdaT3` |

**This mapping is an assumption, not a verified physics identity.** Even the
same `lambdaT3` name is not independently a normalization certificate.
No other source coupling is inferred or assigned zero. To make the
assumption visible, calculation requires `--allow-exploratory-map`.

Furthermore, the `all16_v1` scalar runs are themselves explicitly marked
`illustrative_benchmark=true`, `physical_mass_parameters_available=false`,
and `physical_matchete_rgbeta_mapping_certified=false`. We preserve all three
caveats. This batch **does not** turn the benchmark into a physical scalar
potential or certify a neutrino mass prediction.

## Numerical output

For ordinary-basis models we construct, in the *original source real* basis
(which equals the fermion boundary's basis in exactly these eight cases),

\[
M_S^2(\mu_S)=\operatorname{diag}(m_{S1}^2\mathbf1_{2d_{S1}},
                                   m_{S2}^2\mathbf1_{2d_{S2}}),
\]

and evaluate all exact v4 source channels
`lambda_HHSS_channels_muS` (shape `k,4,4,n,n`) using the chosen
candidate numerical couplings. All channel signs, CG coefficients and
factorials come from the v4 Matchete source tensors; no new portal potential
is invented. The source-resolved fermion LLSS coefficient remains per unknown
signed normalization `N_GammaCC`.

The existing v3 simultaneous scalar matcher evaluates, for this one-loop
LLSS sector only,

\[
 C_{5,\mathrm{LLSS}}^{\mathrm{conditional}}(\mu_S)/N_{\Gamma CC}
 =\frac{\ln(\mu_S/\mu_F)}{16\pi^2}\beta_{C_5}^{\mathrm{formal}}
   + \Delta C_5^{\mathrm{joint\ scalar,formal}}(\mu_S).
\]

The result has unknown overall sign/factor, uncertain quartic aliases, and
**does not include** certified heavy-fermion hard one-loop matching,
additional EFT1 operator sectors or complete scalar-potential matching.
Never use it as a physical neutrino mass matrix.

The runner verifies **full Weinberg weak-SU(2) Ward residual**, source-pivot
parity, full-SU2 export identity, independent contraction and explicit
matching-scale cancellation. If they fail it **reports an error**, not a
numerical neutrino mass.

## PowerShell (no placeholder paths)

Download the v5 ZIP to Downloads, then run from the project root:

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
$zip = Get-ChildItem "$HOME\Downloads" -File -Filter 'T3_joint_scalar_numerical_benchmark_v5*.zip' |
       Sort-Object LastWriteTime -Descending | Select-Object -First 1
if (-not $zip) { throw 'Download the v5 ZIP first.' }
Expand-Archive -LiteralPath $zip.FullName -DestinationPath . -Force
python -m unittest discover -s tests -p 'test_t3_joint_scalar_*py' -v
python -m scripts.RunT3JointScalarNumericalBenchmarkV5 `
    --allow-exploratory-map `
    --output-dir 'output/general_observable/lo_joint_scalar_matching/exploratory_benchmark_v5'
```

Output is *only* under the chosen `lo_joint_scalar_matching` directory:

- `summary.json`: all16 model input status, per-model blockers, source SHA,
  candidate alias table and benchmark scalar masses;
- `numerical_backgrounds/A_m4.json` and `.npz` etc. for eligible models;
- `C5_per_unknown_NGamma_exploratory/A_m4.json` etc.: **conditional** matrices
  (never a certified physical C5).

Expected, **not yet independently verified on your Windows checkout**:
8 exploratory numerical ordinary models and 8 correctly blocked special
models. If the actual source tensors/coupling values violate the expected
Ward projection or any hash/scale/source condition, the runner will list
`ERROR_source_or_trajectory_validation` and return exit code 2. Do not turn
errors into zeros or bypass checks. Use a *fresh* output directory if re-running.

The v5 code adds no modifications to any existing module or upstream
F-tree/RGE source. GitHub is read-only; unpack locally.
