# T3 joint scalar matching — 16-model integrated batch v2

## Files (new; no existing files replaced)

- `Numerical/joint_scalar_matching/large_batch.py`: independently audits the source-resolved all-16 `lambda_[H0R,H0R,a,b]` pivot kernels and computes optional conditional fixed-order leading-log plus joint scalar-bubble response.
- `scripts/RunT3JointScalarIntegratedBatch.py`: all-16 single-command runner.
- `tests/test_t3_joint_scalar_integrated_batch.py`: 10 synthetic/algebraic tests, including real component layouts for all 16 models, exact inter-source cancellations, matrix basis changes, double counting checks, and missing-data blocking.
- `Numerical/joint_scalar_matching/README_INTEGRATED_BATCH_V2.md`: this document.

**The program only writes inside** `output/general_observable/lo_joint_scalar_matching/`. Its inputs are read-only; no GitHub operations occur.

## Run now: audit real published source kernels for all 16

PowerShell in the local repo root:

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
python -m unittest discover -s tests -p "test_t3_joint_scalar_integrated_batch.py" -v
python -m unittest discover -s tests -p "test_t3_joint_scalar_*py" -v
python -m scripts.RunT3JointScalarIntegratedBatch --output-dir "output\general_observable\lo_joint_scalar_matching\integrated_batch_v2"
```

Inputs are the **existing local** files:

- `output/general_observable/lo_fermion_boundary/releases/source_021dcc1e8fe8ba0f/{model}/boundary.json`
- `output/general_observable/lo_eft1_mixing/source_kernel_v1/{model}/quartic_pivot_kernel.json`

The CLI checks the exact documented format, common source boundary hash, all Matchete PL source identifiers, the heavy-scalar real basis, every nonzero scalar pair in each source, the *source-specific* `4*(2-delta_ab)*C_formal[L_SU2=1,L_SU2=1,a,b]` pivot coefficient, and the correctly summed source kernels. It preserves exact zero-benchmark source terms and allows a summed scalar pair to be omitted **only if its source coefficients genuinely cancel** (important in E-class models). It does not invent physical quartics, spinor signs, masses, or EFT1 running results.

Expected outputs:

```
output/general_observable/lo_joint_scalar_matching/integrated_batch_v2/
  summary.json
  scalar_background_contracts/
    A_m4.json
    ...
    E_p2.json
```

`scalar_background_contracts/` files are **explicit incomplete specifications**, NOT ready-to-run numerical endpoints. Their mass/scale/normalization values are left null, and no NPZ is synthesized.

If the default output already exists, use a fresh `--output-dir`, e.g. `integrated_batch_v2b`. `--model B_-1` is available for a focused run.

## Optional numerical computation: one-loop *conditional* response

The single same-`mu_S` scalar background input for each model is a JSON metadata file beside one `.npz` file. This is a deliberately **smaller interface** than the earlier *fully normalized* EFT1 endpoint: it requires only independent quartics, the scalar mass matrix, and the exact source-to-matching real orthogonal basis map. It produces a numerical response per **unit unknown signed normalization**. It does **NOT** waive the physical normalization requirement of the earlier certified endpoint adapter.

Example: `backgrounds/B_m1.json`:

```json
{
  "format": "t3-joint-scalar-conditional-background-v1",
  "model": "B_-1",
  "source_boundary_sha256": "<exact SHA256 from the source boundary>",
  "source_scalar_real_labels": ["<all n ordered labels from its boundary file>"],
  "renormalization_scheme": "MSbar",
  "tensor_convention": "T3_REAL_LLSS_HHSS_WEINBERG_v1",
  "mu_f_gev": 100000.0,
  "mu_s_gev": 1000.0,
  "reference_scale_gev": 1000.0,
  "scalar_mass_origin": "unbroken_electroweak",
  "same_muS_quartic_and_mass": true,
  "no_separate_scalar_threshold": true,
  "tensor_npz": "B_m1.npz",
  "quartic_channel_names": ["<independent source quartic 1>", "<independent source quartic 2>"],
  "provenance": {
    "source": "<actual, identified, independently constructed scalar EFT1 background>",
    "quartic_potential_normalization": "V=lambda_abcd*phi_a*phi_b*phi_c*phi_d/4!"
  }
}
```

The NPZ must contain these four arrays for `n=2dS1+2dS2` real scalars, in matching-component order:

| Array | Shape | Interpretation |
|---|---|---|
| `source_to_matching_real_U` | `(n,n)` | real orthogonal, `s_source = U @ s_matching` |
| `mass_squared_gev2_muS` | `(n,n)` | positive symmetric heavy-scalar mass squared at common `mu_S` |
| `lambda_HHSS_muS` | `(4,4,n,n)` | complete symmetric **real** quartic tensor, including every source invariant, `V = lambda phi^4/4!` |
| `lambda_HHSS_channels_muS` | `(k,4,4,n,n)` | all source quartic contributions, whose sum equals the total quartic |

The scalar source order is checked byte-for-byte against the source boundary's component labels. For B(-1), C(-1) the source labels are mixed X/Y; **do not silently treat them as unmixed S1/S2**. For self-bilinear cases, source labels are already split real fields; **do not double the dimension**. A positive-definite mass matrix and the full SU(2)/Higgs Weinberg projection are mandatory. Exact EWSB-resummed masses can be declared as `EWSB_resummed`, but higher-dimension contamination may cause the Ward test to fail; this does not license a weakened physical statement.

Run when actual backgrounds are available:

```powershell
python -m scripts.RunT3JointScalarIntegratedBatch `
    --background-dir "output\general_observable\lo_joint_scalar_matching\scalar_backgrounds" `
    --output-dir "output\general_observable\lo_joint_scalar_matching\integrated_batch_v2_with_backgrounds"
```

The `--background-dir` is **optional** for the immediately runnable pivot audit. Without it, the program never reports an absent scalar mass/quartic as zero. With it, every selected model must have a complete `A_m4.json` + `A_m4.npz` etc.; missing/invalid models produce explicit blockers and exit code `2`.

## Numerics, scale and operator accounting

With the **formal** source LLSS polynomial `C^formal` and an unknown signed source-to-master/Weyl normalization `N`, define `C_LLSS = N*C^formal` (the numeric output below takes `N = 1` only as a coefficient of that unknown). There is no sourced physical value for `N`.

In the **specific tensor convention implemented by the scalar-bubble kernel**, the leading-order, fixed-coupling EFT1 insertion is

\[
B^{fgij}_{AB}=\sum_{a,b}C^{fgij,\mathrm{formal}}_{ab}\lambda_{ABab},\quad
\beta_{C_5}^{\mathrm{LLSS,formal}}=\mathrm{Project}_{\mathrm{Weinberg}} B.
\]

For `C5(muF)` without a documented independent hard matching contribution, the **LLSS-induced** part alone is

\[
C_{5,\mathrm{log}}^{\mathrm{formal}}(\mu_S)=
\frac{\ln(\mu_S/\mu_F)}{16\pi^2}\beta_{C_5}^{\mathrm{LLSS,formal}},
\]

\[
C_{5,\mathrm{bubble}}^{\mathrm{formal}}(\mu_S)=
-\frac{1}{2(16\pi^2)}\,
\mathrm{Project}_{\mathrm{Weinberg}}\sum_{a,b}
\widehat C^{\mathrm{formal}}_{ab}\widehat\lambda_{ABab}
\overline B_0(0;m_a^2,m_b^2;\mu_S^2),
\]

\[
C_{5,\mathrm{LLSS}\text{-}\mathrm{sector}}^{\mathrm{formal}}
=C_{5,\mathrm{log}}^{\mathrm{formal}}+C_{5,\mathrm{bubble}}^{\mathrm{formal}}.
\]

The actual conditional LLSS-sector result is `N * C5_formal`. There is **no** claim that `C5(muF)` vanishes physically, or that a full-theory hard loop has been included. In this one-loop LLSS sector, changing `C_LLSS` by its own one-loop running inside a one-loop scalar bubble enters only at **two-loop** order; using `C_LLSS(muF)` in both terms avoids that form of double counting. Running in the quartic/mass inputs is also higher-order for this specific one-loop contraction when evaluated consistently.

The program checks, at fixed couplings and masses:

\[
\frac{d}{d\ln\mu_S} C_{5,\mathrm{log}}^{\mathrm{formal}}
=+\frac{\beta_{C_5}}{16\pi^2},\quad
\frac{d}{d\ln\mu_S} C_{5,\mathrm{bubble}}^{\mathrm{formal}}
=-\frac{\beta_{C_5}}{16\pi^2}.
\]

Numerical outputs include full complex symmetric `3×3` matrices for total, running and finite bubble; matrices per independent PL source and quartic channel; scalar mass eigenvalues; full Higgs-SU(2) Ward residuals; and a matching-scale cancellation residual. The result is named `C5_*_per_N_gev_inverse` and always has `physical_C5_certified=false`, `source_to_eft_weyl_normalization_factor=null`.

**Important normalization caveat:** the existing independent source-pivot worker reports its neutral-component source-polynomial kernel as `4*(2-delta_ab)*C_formal[H0R,H0R,a,b]`. That normalization is audited exactly as *source-polynomial data*. The new complete tensor projection uses the declared scalar-bubble adapter conventions and does not silently infer the relative source-to-master factor from a **single Higgs component**. A signed Matchete `GammaCC` and independent operator-normalization comparison remain required. Comparing the two without operator translations can introduce a factor-of-two error.

## Verified vs outstanding

**Verified by test:** all 16 synthetic layouts, all 16 synthetic source/pivot channel inventories, cancellation of multiple source terms, actual `3×3` complex symmetry and phase, complete four-index quartic reconstruction, positive scalar eigenmasses, mixed basis rotation, per-source and per-quartic accounting, Ward projection, fixed-order matching-scale-log cancellation, missing-input failure, and safe output routing.

**Verified by reading actual GitHub data:** the source boundary and conditional pivot kernels are available for all 16 models; the source code and kernel metadata explicitly mark their physical normalization as unresolved. The pivot formula was independently checked against all 16 existing per-source numerical kernel exports with no coefficient discrepancy (read-only inspection). The user must run the CLI to verify their *local* checkout's source and kernel files.

**Not implemented / certified:** unambiguous signed `GammaCC` -> left-handed Weyl mapping, automatic correct physical `lambda_HHSS` and scalar masses at muS from all general T3 scalar interactions, independent hard-F finite matching, and other EFT1 operator/finiteness counterterms. The batch never promotes its conditional result to a complete one-loop physical neutrino-mass calculation. The existing fully gated certified bridge (`BridgeT3JointScalarEFT1.py`) remains unchanged and is still required before certifying a physical endpoint.
