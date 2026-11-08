# T3 all-16 full-tensor *exploratory* Wilson trial

This is a **new numerical trial** intended to extend the A(-4) calculation to all 16 source models using existing local data. It is **not** the fully mapped 17-beta corrected A(-4) renormalizable Wilson background applied to the other 15 models. The v6 saved scalar background supplies `lambda_HHSS(muS)` but **does not supply the complete `lambda_abcd(mu)` for the other models**, so the trial uses an explicit reduced scalar potential.

## Exact input sources

- `output/general_observable/lo_fermion_boundary/releases/source_021dcc1e8fe8ba0f/{slug}/boundary.json`: all-F-first formal 3-family `LLSS` tensor.
- `output/general_observable/lo_joint_scalar_matching/exploratory_all16_v6_native_batch1_local/numeric_backgrounds/{slug}.json` and `.npz`: complete model-specific heavy scalar quadratic matrix and `HHSS` quartic, including special mass matrices/basis transformations.
- `output/general_observable/lo_joint_scalar_matching/exploratory_all16_v6_native_batch1_local/C5_per_unknown_NGamma_exploratory/{slug}.json`: original fixed-order control.
- `output/general_observable/all16_v1/{slug}/running.json`: actual two renormalizable EFT1 scale snapshots at 100 TeV and 1 TeV.

## Exact numerical procedure

1. Apply explicitly declared `C_master(muF) = 0.5 * C_old_formal(muF)`, with a tested two-projector convention identity. The factor 0.5 is **not** an independently measured native `GammaCC` factor. The master `C5` convention includes `1/2` in its operator; the old joint-scalar projector does not. Test that the new zero-running log and bubble exactly reproduce *both* old all-16 terms, before adding RG effects. No invented Wilson boundary support or guessed signs.
2. Embed the source in all **45 surviving SM Weyl coordinates** and the full real scalar basis. For A/B/C/D, `n_s=12` and the number of independent coefficients is **80,730**; for E, `n_s=16`, it is **140,760**.
3. Build the actual SU2/U1/SU3 generators using `T3Hypercharges(alpha) = {S1:alpha/2, S2:(alpha+2)/2}` from `T3ModelCatalog.wl`, and rotate scalar generators into the special mixed/self-split basis using the existing `special_mass_v6.build_special_basis` map. All heavy fermions have already been removed.
4. Reuse `GeneralEFT1SixTermReferenceCandidate.six_term_beta` and evolve the **complete complex symmetric 45x45xnsxns** Wilson tensor with explicit log-mu midpoint RK2. Interpolate **actual** F- and S-endpoint SM gauge/Yukawa couplings and Higgs self quartic; keep v6 `HHSS(muS)` fixed along the log interval. **Other scalar quartic channels (heavy-only and extra independent operators) are not available and set to zero explicitly.** This is a reduced-potential trial, not a complete model-specific quartic RGE.
5. Project full generated LLHH and perform **one simultaneous scalar threshold** with the existing `scalar_bubble_matching`. Do not integrate split/mixed scalars separately. Check both full Weinberg projections and output conditional 3x3 mass matrices and singular values.
6. Report failure per model and **never create zero fallback masses**. Each file carries `full_model_specific_quartic_tensor_included=false`.

## PowerShell

From project root, install the ZIP, then:

```powershell
python -m unittest discover -s tests -p "test_all16_conditional_full_wilson_trial.py" -v

# Check A(-4) before running 16 models:
python scripts\RunT3All16ConditionalFullWilsonTrial.py --model A_-4 --steps 2 --output-root "output\general_observable\all16_full_wilson_A_m4_smoke_v1"

# Full trial, lower precision first pass:
python scripts\RunT3All16ConditionalFullWilsonTrial.py --steps 2 --output-root "output\general_observable\all16_full_wilson_reduced_potential_v1"
```

Add `--save-tensors` to save the full boundary/endpoint C arrays. Use `--steps 4` in a fresh output directory to check integration-step sensitivity. A single whole-batch summary/CSV is produced even if some models fail, with explicit per-model error messages and exit code 2.

## Not the same as the existing full A(-4) calculation

A(-4)'s recently completed high-quality benchmark used all 17 native RGBeta beta channels, the complete known A(-4) quartic tensor, and 7 running samples. The new all-16 trial necessarily uses only the scalar HHSS/Higgs quartic and 2 snapshots, and its boundary is normalized to the **older all-16 fixed-order benchmark**. Do not directly interpret the difference as an A(-4) physics correction.

`C5` and masses are **conditional, rough numerical proxies**, not certified physical predictions, and no SMEFT running to v is performed. The fixed-order control must agree algebraically with the saved all-16 result before any Wilson trial is accepted.

No existing production source is changed. All calculations stay in local project files; GitHub is read-only.
