# T3 joint scalar matching — all-16 full-SU(2) readiness batch v3

## Installation and paths

This is **additive** to joint scalar matching v1 + fermion bridge v1 + integrated batch v2. It creates four **complete new files**:

- `Numerical/joint_scalar_matching/source_su2_batch.py`
- `scripts/RunT3JointScalarBatchV3.py`
- `tests/test_t3_joint_scalar_batch_v3.py`
- `Numerical/joint_scalar_matching/README_PHYSICAL_READINESS_BATCH_V3.md`

It does **not** modify any existing shared code or outputs, or GitHub. It reads existing repository artifacts from the fixed paths below. All new output is restricted to `output/general_observable/lo_joint_scalar_matching/`.

## Run on Windows PowerShell (from repository root)

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
$zip = Get-ChildItem "$HOME\Downloads" -File -Filter 'T3_joint_scalar_batch_v3*.zip' |
    Sort-Object LastWriteTime -Descending | Select-Object -First 1
if (-not $zip) { throw 'Download T3_joint_scalar_batch_v3.zip to your Downloads folder first.' }
Expand-Archive -LiteralPath $zip.FullName -DestinationPath . -Force
Test-Path 'scripts\RunT3JointScalarBatchV3.py'
python -m unittest discover -s tests -p 'test_t3_joint_scalar_*py' -v
python -m scripts.RunT3JointScalarBatchV3 `
    --output-dir 'output\general_observable\lo_joint_scalar_matching\physical_readiness_v3'
```

Expected: `Test-Path` returns `True`, 35 joint-scalar tests pass, and `summary.json` reports 16 independently checked source/pivot/full-SU(2) source records if your local data correspond to inspected GitHub `main` 2026-11-02. Failures remain explicit **per model** and return exit code 2. A missing or incompatible all-16 summary fails at startup. The runner never overwrites an existing `summary.json`: choose `physical_readiness_v3b` for another run.

### Inputs (read-only, already present in inspected GitHub main)

1. Fermion matching source release, 16 real-component formal PL tensors at `mu_F = 100000 GeV`:
   `output/general_observable/lo_fermion_boundary/releases/source_021dcc1e8fe8ba0f/{model_slug}/boundary.json` plus `summary.json`.
2. Quartic pivot kernel at `mu_F`:
   `output/general_observable/lo_eft1_mixing/source_kernel_v1/{model_slug}/quartic_pivot_kernel.json`.
3. New independent full-SU(2) factorised source response:
   `output/general_observable/lo_eft1_mixing/source_su2_response_v1/{model_slug}/full_su2_response.json` and `summary.json`.
4. Existing **non-certified** signed normalization evidence:
   `output/general_observable/lo_fermion_boundary/normalization_probes/signed_zero_momentum_d734e270ed0138fa.json`, `.../native_feynman_a_m4_b669a99ea16d439b8c50e81b16dbd1c8.json`.

Every source SHA, model, source index, scalar pair, scalar-pair multiplicity and **all** `(p,q,i,j)` leptonic SU(2)/three-flavour coordinates are independently checked, for every individual PL source, against the original fermion boundary. A zero-at-benchmark source is retained as an empty channel; exact PL-source cancellations are retained without fake nonzero totals. The pivot audit additionally checks the independent neutral-Higgs conditional kernels. Nothing assigns an absolute `GammaCC` sign or propagates muF LLSS to muS implicitly.

The native Feynman-rule diagnostic currently reports **`Blocked`, `blocked_phase="EFT1 Feynman rules"`, `reason="No nonempty native EFT1 LLSS contact vertex"`**. The separate signed Grassmann zero-momentum calculation proves identities *for a declared two-component action*, not Matchete's installed native `GammaCC` convention. Thus the signed factor is **not** certified and every C5 result from this module, if any, is conditional per unit unknown factor.

## Optional conditional scalar response: full numerics as soon as source-hashed scalar tensors exist

The previous v2 runner already outputs the required all-16 empty scalar-background **contracts**:

`output/general_observable/lo_joint_scalar_matching/integrated_batch_v2/scalar_background_contracts/{slug}.json`.

To supply independently verified scalar backgrounds, create an input folder **of your choosing** (e.g. `output/general_observable/lo_joint_scalar_matching/backgrounds_at_muS`) with `A_m4.json/.npz`, etc. Use the exact `t3-joint-scalar-conditional-background-v1` metadata and sibling NPZ schema documented in `README_INTEGRATED_BATCH_V2.md`:

- `source_to_matching_real_U` real orthogonal `(n,n)`;
- `mass_squared_gev2_muS` positive symmetric `(n,n)`;
- `lambda_HHSS_muS` real fully symmetric appropriate Higgs/scalar indices `(4,4,n,n)`;
- `lambda_HHSS_channels_muS` `(k,4,4,n,n)` and all channel names;
- source SHA, real-component labels, common `mu_S_gev`, `mu_f_gev`, MSbar and unbroken-EW or declared EWSB-resummed mass origin, and explicit independent provenance.

Use **one** common scalar threshold for both scalar multiplets, including split-real/self-bilinear and B/C mixed cases. Do not use the muF boundary directly as a muS-renormalized endpoint or claim that a symbolic scalar-mass-basis placeholder is a numerical matrix.

Then run:

```powershell
python -m scripts.RunT3JointScalarBatchV3 `
    --background-dir 'output\general_observable\lo_joint_scalar_matching\backgrounds_at_muS' `
    --output-dir 'output\general_observable\lo_joint_scalar_matching\physical_readiness_v3_with_backgrounds'
```

This additionally performs an independent contraction over **all** SU(2) indices and unordered scalar pairs in the original source basis, and compares its projected beta to the rotated matching-basis beta. It then delegates to v2's tested `evaluate_conditional()` for the fixed-order leading-log EFT1 mixing and simultaneous finite scalar bubble, with source-resolved and quartic-resolved contributions:

\[
 C_5^{\mathrm{LLSS,1loop}}\big/ N_{\Gamma CC}
 = \frac{\ln(\mu_S/\mu_F)}{16\pi^2}\,\Pi_W[\lambda(\mu_S) C_{LLSS}^{(0)}(\mu_F)]
   + \Delta C_5^{\mathrm{scalar\ bubble}}(\mu_S)/N_{\Gamma CC}.
\]

`N_GammaCC_to_Weyl` is **unknown**; the expression is **per unit source-to-master factor**, not a physical C5. At this order the running of the LLSS coefficient inside an already one-loop scalar bubble would be two-loop and is not double-counted. The two terms cancel the artificial matching-scale logarithm at fixed input masses/couplings; this is *not* a proof of complete renormalization-group invariance with running scalar masses/couplings.

Without `--background-dir` the default run deliberately performs **no physical C5 computation**, and reports the exact missing inputs in `summary.json` for each model. With an incomplete background directory, each missing or incompatible model is explicitly blocked rather than generating a zero C5.

## Verified versus outstanding

**Verified here in synthetic regression (35/35 tests including older suites):** all-16 source indices/layouts, four-coordinate SU(2) flavour parity, split-real and B/C cases, pair multiplicity, missing/duplicate pair blocking, exact inter-source cancellations, independent beta contraction under nontrivial orthogonal rotations, missing/inconsistent evidence handling, explicit no-background status, and fail-closed CLI. The earlier eight scalar-matching tests and ten integrated-batch tests remain included.

**Verified from GitHub source inventory:** latest `main` contains the 16 full-SU(2) source response files, 16 pivot kernels, signed Grassmann diagnostic and currently blocked native Feynman rule report, as listed above.

**Not verified:** actual Windows live 16-model run of this new batch (must execute locally); full numerical scalar quartic and mass matrices at one muS; native signed four-point amplitude and master spinor normalization; hard one-loop F matching and any other operator sectors/finite counterterms. No physical neutrino mass is predicted here.

## Outputs

Only `output/general_observable/lo_joint_scalar_matching/physical_readiness_v3/summary.json` on the default run; optional `conditional_C5_per_GammaCC/{slug}.json` for valid user-provided backgrounds. No shared worker output is rewritten.
