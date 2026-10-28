# Reduced T3 scalar-sector final numerical summary

## Scope and accuracy

Primary data: **two frozen-UV benchmarks**, 16 models each, exact zero additional scalar boundary,
18 real fitted Yukawa parameters per model, candidate one-loop scalar EFT1 running,
model-specific existing leading fixed-order one-loop C5 matching. Frozen boundary is imposed at 10^5 GeV.

Reference: **newly optimised scenario-specific A(-2) Yukawa matrices**, transferred unchanged
to each model. Other model-specific inputs and group-theory factors are retained.

UV-running information is **exploratory**: existing high-scale inputs are evolved, not refitted,
and mass/hierarchy failures are not interpreted as neutrino-data exclusions.

## Primary frozen comparison

| Scenario | Accepted individual fits | Acceptable shared A(−2) evaluations | Total |
|---|---:|---:|---:|
| Small topology coupling | 16 | 5 | 16 |
| Large topology coupling | 16 | 5 | 16 |

Acceptance is defined by the five-observable Gaussianised NuFIT target χ² ≤ 5.86.
The 18 real Yukawa fitting degrees of freedom exceed the five fitted observables: excellent χ²
is fit flexibility, not a prediction or a guarantee of uniqueness.

| Model | small: individual | small: A(−2) shared | large: individual | large: A(−2) shared |
|---|---:|---:|---:|---:|
| A_+0 | 1.121489e-11 | 0.00817089 | 0.005162948 | 0.007585785 |
| A_-2 | 5.301052e-17 | 5.301052e-17 | 3.539857e-11 | 3.539857e-11 |
| A_-4 | 2.081946e-16 | 11527.68 | 2.860277e-14 | 11464.09 |
| B_+1 | 1.198898e-21 | 1.473904 | 3.446356e-11 | 2.187401 |
| B_-1 | 8.937932e-18 | 37.39544 | 1.067709e-08 | 35.45188 |
| B_-3 | 8.297949e-15 | 12177.75 | 1.482725e-09 | 12124.08 |
| C_+1 | 6.602229e-20 | 9273.324 | 3.484557e-20 | 9256.323 |
| C_-1 | 2.286049e-10 | 9749.82 | 2.453295e-18 | 9739.222 |
| C_-3 | 8.040784e-11 | 17818.81 | 8.57308e-15 | 17797.39 |
| D_+0 | 2.575589e-12 | 0.0002023746 | 1.046566e-08 | 0.0002090227 |
| D_+2 | 4.125362e-12 | 232.8243 | 8.591204e-11 | 238.6801 |
| D_-2 | 1.2035e-10 | 0.01288426 | 3.4551e-11 | 0.01225028 |
| E_+0 | 1.368353e-09 | 3922.104 | 4.127976e-11 | 4106.154 |
| E_+2 | 2.435921e-17 | 2932.508 | 5.58612e-11 | 3106.756 |
| E_-2 | 1.185518e-09 | 3915.549 | 3.026984e-15 | 4099.545 |
| E_-4 | 8.270065e-12 | 15394.62 | 8.651976e-11 | 15437.31 |

## Running-UV diagnostic at fixed inputs

Historical running-UV Yukawas and scalar input boundaries at 10^7 GeV are left unchanged.
The nominal 1 TeV EFT1 endpoint masses are a **formal diagnostic**, not physically matched
threshold masses when the scalars become heavy or reverse the decoupling hierarchy.
For all 16 historical seeds the original 0.5–5 TeV scalar window fails.

| Model | formal mS1 at 1 TeV [TeV] | formal mS2 at 1 TeV [TeV] | rhoF at 100 TeV |
|---|---:|---:|---:|
| A_+0 | 40.9155 | 42.3856 | 0.4110 |
| A_-2 | 39.8900 | 40.5393 | 0.3805 |
| A_-4 | 35.1197 | 59.4371 | 0.5250 |
| B_+1 | 30.4184 | 37.8432 | 0.3593 |
| B_-1 | 72.8071 | 100.0762 | 1.0704 |
| B_-3 | 31.3093 | 36.4894 | 0.3562 |
| C_+1 | 63.5017 | 66.1264 | 0.5868 |
| C_-1 | 118.5830 | 169.1006 | 1.6976 |
| C_-3 | 63.5579 | 65.0881 | 0.5903 |
| D_+0 | 34.1921 | 51.7206 | 0.4955 |
| D_+2 | 45.4834 | 46.3077 | 0.4075 |
| D_-2 | 34.7952 | 50.8104 | 0.4925 |
| E_+0 | 46.4420 | 64.0919 | 0.6441 |
| E_+2 | 43.4201 | 44.4576 | 0.3729 |
| E_-2 | 45.3040 | 56.5035 | 0.5255 |
| E_-4 | 45.7779 | 62.6794 | 0.5309 |

UV trajectory replay: **candidate renormalisable UV 10^7→10^5 GeV replays completed**.

## Interpretation and limitations

1. Model-specific Yukawa refitting reproduces oscillation data in both frozen scenarios.
   The shared A(−2) reference is often not transferable; the comparison does not exclude models.
2. The selected high-scale boundary yields heavy scalar running masses and violates the
   *chosen benchmark window*, not a general theoretical requirement.
3. One-loop RG-improved couplings are used with leading fixed-order C5; this is not a uniform
   two-loop calculation. Radiatively generated extra B/C quartics are not inserted into C5.
4. Eight additional holomorphic quadratics have zero input boundary and are not evolved;
   the zero-boundary relevance argument is conditional, not an exact symmetry.
5. Four extreme-hypercharge models use nonuniform legacy one-scale C5 artifacts.
   Full Matchete/RGBeta component normalization and general sequential Wilson matching
   have not been certified. No vacuum stability, confidence intervals or predictivity claim.

## Provenance

This report was generated from preexisting local result files; no fit values were modified.
- `RGE/general/generated/general_rgbeta_quartic_registry.json`: `86a096410384641efa8f751b2aa98e8fee4bb6d3cccb4cfd116c559b71182eea`
- `configs/nufit_6_1_ic24_no.json`: `27bd4d4297fddb3ef71b5b804afd8e31d44d302b569f7cf8c34d37d54005e15e`
- `output/general_t3_stage22/stage30_A_p0_full_uv_betas.json`: `6f712e8e0406b68c0387923b0f374f99474707c6ad32c662d0abc50c081f4342`
- `output/general_t3_stage22/stage31_A_m2_full_uv_betas.json`: `3f48e43fbe82e59f51ccbde12f9148d28e40e6d646028355de461147a07d2018`
- `output/general_t3_stage22/stage31_A_m4_full_uv_betas.json`: `fa828e761e7e96da76c82f0c2d61c1c1c832a014d47864219f62bbcd7cbeeee1`
- `output/general_t3_stage22/stage31_B_m1_full_uv_betas.json`: `a810fa96d6794177a129732320b446cbd7d4e94a260ee3eda10044b0274ce8f6`
- `output/general_t3_stage22/stage31_B_m3_full_uv_betas.json`: `6ea2d9faaa8773e5aba7842eb5bd95307cb452cb77066139be69f071b8edf9b0`
- `output/general_t3_stage22/stage31_B_p1_full_uv_betas.json`: `dbb6e6124c13130496dd9e2f78ba3c9730010a6c799b50d25c2f40b99f235540`
- `output/general_t3_stage22/stage31_C_m1_full_uv_betas.json`: `1e7c063a5c17ce907b1dd0dae299881b8a66f2d9e1a47dedd447d7588e91ac5d`
- `output/general_t3_stage22/stage31_C_m3_full_uv_betas.json`: `759dcbe185b940dd133db492160cdad0792e81b4495da081d945a18c3e8fccb2`
- `output/general_t3_stage22/stage31_C_p1_full_uv_betas.json`: `87cf0df4ae6eebf249f123965e3f25a7e5a3b92aa93cc0f69e3ab3ce86c1aa2c`
- `output/general_t3_stage22/stage31_D_m2_full_uv_betas.json`: `d8058f1c24f6be0b92b35504cd67a8ea858c1c73ce3664dab73b77a3f20e13c6`
- `output/general_t3_stage22/stage31_D_p0_full_uv_betas.json`: `025b7fe84571a4e323d44c55c8c10c2ec8c737c12e9efc23a9f039c3717ca5e8`
- `output/general_t3_stage22/stage31_D_p2_full_uv_betas.json`: `63f271e895be6074945b844538437d883c0533aea547a0bff067d02ea0dc3f17`
- `output/general_t3_stage22/stage31_E_m2_full_uv_betas.json`: `5effb3d2d5819299d96c6b21acd519f0072e42daaaeb2f937432e1dda04ad379`
- `output/general_t3_stage22/stage31_E_m4_full_uv_betas.json`: `2762cd535c49b771a65d0cb2ce96260332c92ea480f249a858a387174e6e84d8`
- `output/general_t3_stage22/stage31_E_p0_full_uv_betas.json`: `2bf3a887921540a8465a8c6dfb92a359754058ac2069d84e06b9139a34e2e1c4`
- `output/general_t3_stage22/stage31_E_p2_full_uv_betas.json`: `74dc4dfe4ab621abb92a00a98f73ccbef76ffc8d6006c6d27fe937d6b18816c8`
- `output/studies/lambda_5e-9_refit/fits/T3_dS1_1_dS2_3_dF_2_alpha_m2.json`: `f63b873e004261613312975908587a12d6a8965cefb17d6521a4db88bf1acab4`
- `output/studies/lambda_5e-9_refit/fits/T3_dS1_1_dS2_3_dF_2_alpha_m4.json`: `dd1fb84061c46c224409bf823a0efcc17b711912ca49a63e903c0c572a683ef1`
- `output/studies/lambda_5e-9_refit/fits/T3_dS1_1_dS2_3_dF_2_alpha_p0.json`: `3816b7028577add29225eb2d53ba1b0e39d13a8e166a56e8381d9debafd4d7da`
- `output/studies/lambda_5e-9_refit/fits/T3_dS1_2_dS2_2_dF_1_alpha_m1.json`: `bd8076ed9ea5fdb2cffcbd09918d2110a8398544be78d2e8e8d7defe3723ea20`
- `output/studies/lambda_5e-9_refit/fits/T3_dS1_2_dS2_2_dF_1_alpha_m3.json`: `cdddc842340ca41f69733dd3dcda863c98f7e88b7a29bb494765bb424b41cee5`
- `output/studies/lambda_5e-9_refit/fits/T3_dS1_2_dS2_2_dF_1_alpha_p1.json`: `4adac4a15776b94e11129642bce362961468edb9eecedf0e5923aadc5fa1e58c`
- `output/studies/lambda_5e-9_refit/fits/T3_dS1_2_dS2_2_dF_3_alpha_m1.json`: `10ed72e8f0872c6429b218014033ec606e8da9cbc13e5f8e1a90b62e84fce91c`
- `output/studies/lambda_5e-9_refit/fits/T3_dS1_2_dS2_2_dF_3_alpha_m3.json`: `563a49a2af858b02a562cf254cc9025eb7f0537a4f15abc20ee3b9fd0a31fcfd`
- `output/studies/lambda_5e-9_refit/fits/T3_dS1_2_dS2_2_dF_3_alpha_p1.json`: `f7f92c9daa572dccc5124688f7b2403a748bfc2bc603c046425f44ef34ec12f3`
- `output/studies/lambda_5e-9_refit/fits/T3_dS1_3_dS2_1_dF_2_alpha_m2.json`: `b2a624e1962fa5f70f8ab93663ded7854fd8eb6d783c3e12e1f9fde8943af2d1`
- `output/studies/lambda_5e-9_refit/fits/T3_dS1_3_dS2_1_dF_2_alpha_p0.json`: `4d45c4cd3e93c206e00f06b876ec68a35d7b82c661fd734c35cec2d9b0226062`
- `output/studies/lambda_5e-9_refit/fits/T3_dS1_3_dS2_1_dF_2_alpha_p2.json`: `c6c3857150d3b89c2b3626e266afb383049b660e1b00456b646fff8eaa5e4130`
- `output/studies/lambda_5e-9_refit/fits/T3_dS1_3_dS2_3_dF_2_alpha_m2.json`: `c370bc3cb88b8a05a94e76df87da920d09e685532816b90e54a5843e43369dca`
- `output/studies/lambda_5e-9_refit/fits/T3_dS1_3_dS2_3_dF_2_alpha_m4.json`: `e22d7b209ea98b531060e1f847ceb05558d21a19f2304a269beb4d5d79272b98`
- `output/studies/lambda_5e-9_refit/fits/T3_dS1_3_dS2_3_dF_2_alpha_p0.json`: `fa7cafcd9cd3cff0bb5b35e39f811942674d2f1da5e26a400a695356337d9935`
- `output/studies/lambda_5e-9_refit/fits/T3_dS1_3_dS2_3_dF_2_alpha_p2.json`: `399e4aa8a3166de911ccf6528d1b638f5e4edd5324879c872d88869d1904b5e1`
- `output/studies/lambda_5e-9_refit/fitted/T3_dS1_1_dS2_3_dF_2_alpha_m2.json`: `446f9d288662c4d155f5af9a4989de3d49deadc840eb8a90b3b68baf8a13c4fb`
- `output/studies/lambda_5e-9_refit/fitted/T3_dS1_1_dS2_3_dF_2_alpha_m4.json`: `a5d2fe3d10b2a1b5af834e1cabd2bc98f27175d9826e7dc197fc2c2abdd4105f`
- `output/studies/lambda_5e-9_refit/fitted/T3_dS1_1_dS2_3_dF_2_alpha_p0.json`: `81b1ace1a22e0a0bb1cbdcebad30c9e9c80678ffe5d4705a0b9b15a6e4c28add`
- `output/studies/lambda_5e-9_refit/fitted/T3_dS1_2_dS2_2_dF_1_alpha_m1.json`: `e8bb30476227feda2bc4ed101417a9c0fded1adc3cbbbb5a0b9e3fdc9f949e85`
- `output/studies/lambda_5e-9_refit/fitted/T3_dS1_2_dS2_2_dF_1_alpha_m3.json`: `37e0e930357a68d1d7b6170479cbecc91e96d3b440f3d141e2ab6e4ed7568859`
- `output/studies/lambda_5e-9_refit/fitted/T3_dS1_2_dS2_2_dF_1_alpha_p1.json`: `f06dfafef056aeec8c5ec8afb43f5503d6dce73e38067be64ca9e0edd4ffb1ef`
- `output/studies/lambda_5e-9_refit/fitted/T3_dS1_2_dS2_2_dF_3_alpha_m1.json`: `b68ab5884cd4f5600d808f598b9d7691331dc60e9f617646c73ee3e28dae9f20`
- `output/studies/lambda_5e-9_refit/fitted/T3_dS1_2_dS2_2_dF_3_alpha_m3.json`: `7cb37e139de86feba959fb8eeaeb606691f1ae12f2286d764e13abf0ec58b856`
- `output/studies/lambda_5e-9_refit/fitted/T3_dS1_2_dS2_2_dF_3_alpha_p1.json`: `e2044fde78aef1b36c13d8a35d5cf72ebc9f8482e0eb1c618c781bf17b4665f4`
- `output/studies/lambda_5e-9_refit/fitted/T3_dS1_3_dS2_1_dF_2_alpha_m2.json`: `5b8a1e3547a0ba0b3460b08620dba5ec26479ced24fa044223cb5eeebf6eb190`
- `output/studies/lambda_5e-9_refit/fitted/T3_dS1_3_dS2_1_dF_2_alpha_p0.json`: `36c1fd0a571356cf7ba3e044d5ad273bee882b77997a63a3c5432caa204b3582`
- `output/studies/lambda_5e-9_refit/fitted/T3_dS1_3_dS2_1_dF_2_alpha_p2.json`: `4e117c21977fa3600b772d084d1a79314902456c34d2bfa0f4c5143ec9b6d8d9`
- `output/studies/lambda_5e-9_refit/fitted/T3_dS1_3_dS2_3_dF_2_alpha_m2.json`: `de37c0041ddf7932af2f4e8f03a8830e1d1273032ef372888f9d22bd54cab856`
- `output/studies/lambda_5e-9_refit/fitted/T3_dS1_3_dS2_3_dF_2_alpha_m4.json`: `61b91654ffb1b12ea0e60ca709b3aa2335d162c8f2af28e8bd21b8eca0e2c0d9`
- `output/studies/lambda_5e-9_refit/fitted/T3_dS1_3_dS2_3_dF_2_alpha_p0.json`: `772279807f975d17be94521724511cd72767720014bc7de6a4fdc5f01bf674b6`
- `output/studies/lambda_5e-9_refit/fitted/T3_dS1_3_dS2_3_dF_2_alpha_p2.json`: `f03f0bbf0e3d44c9e8933c9f529b95b184b3aded8f665f2214c9047df41e170e`
- `output/studies/reduced_scalar_sobol/reduced_48_v1/small_uv_running_5e-9/A_m2/fit.json`: `d999bdfa05b42d578eba04a9a82176d483920f4219f137d77177ecd25266803b`
- `output/studies/reduced_scalar_sobol/reduced_48_v1/small_uv_running_5e-9/A_m4/fit.json`: `6f9d04e9a3addcfea97e80921d4442b82f95137146e44e96e792f3698bf9256a`
- `output/studies/reduced_scalar_sobol/reduced_48_v1/small_uv_running_5e-9/A_p0/fit.json`: `c81a4ae725b044d5261cd481a1700bfc34e2465793b5d1e479827d0a5aec125b`
- `output/studies/reduced_scalar_sobol/reduced_48_v1/small_uv_running_5e-9/B_m1/fit.json`: `45a11ca26998f24501e573846672000d43c1197dd65d67a13a6d869ae548dc19`
- `output/studies/reduced_scalar_sobol/reduced_48_v1/small_uv_running_5e-9/B_m3/fit.json`: `b24b2297a3cd93d5c81475cca46b4438658aa76ae42528f48243088647001046`
- `output/studies/reduced_scalar_sobol/reduced_48_v1/small_uv_running_5e-9/B_p1/fit.json`: `bfe7e3058f67de033fbfbc7b451ae1f2246f0f05807b84d94cd364957da4b9db`
- `output/studies/reduced_scalar_sobol/reduced_48_v1/small_uv_running_5e-9/C_m1/fit.json`: `264e80c42c9195cae4c4356dcda8df9879d1da69fb62111a50fe9f33a9f6d0f0`
- `output/studies/reduced_scalar_sobol/reduced_48_v1/small_uv_running_5e-9/C_m3/fit.json`: `41d8de54a5fd15dddff0558d2810305ff51a979bb578d8982dd9c9f70670e819`
- `output/studies/reduced_scalar_sobol/reduced_48_v1/small_uv_running_5e-9/C_p1/fit.json`: `22af3de1cf14db38d60a563282b197a8dc8a401e9a8e53c3b45f274daed44f06`
- `output/studies/reduced_scalar_sobol/reduced_48_v1/small_uv_running_5e-9/D_m2/fit.json`: `687ef89a636448050fd115c75aff06141f46bdb04a0b481c11d62273fa72d585`
- `output/studies/reduced_scalar_sobol/reduced_48_v1/small_uv_running_5e-9/D_p0/fit.json`: `4a2e07cb6a26f450bfc0f14ffc0e9a66220647fea396b5d539176324858748b9`
- `output/studies/reduced_scalar_sobol/reduced_48_v1/small_uv_running_5e-9/D_p2/fit.json`: `e2ffe53a416417d11dd0392a039717cd1e4600aceff98375648a7621471d8417`
- `output/studies/reduced_scalar_sobol/reduced_48_v1/small_uv_running_5e-9/E_m2/fit.json`: `23e445ccceedda06a1d2ec7926fa8b7296a1ceb756590a38c6b5869297eb8ad9`
- `output/studies/reduced_scalar_sobol/reduced_48_v1/small_uv_running_5e-9/E_m4/fit.json`: `2b032fe1ecccd4f652228a5b65b550e7c76cd53076729e8a7aaa9e20177fa91b`
- `output/studies/reduced_scalar_sobol/reduced_48_v1/small_uv_running_5e-9/E_p0/fit.json`: `cb89c8e676054529620ad67d890fe244ec89c0d95254df549e19d81f801ca5c0`
- `output/studies/reduced_scalar_sobol/reduced_48_v1/small_uv_running_5e-9/E_p2/fit.json`: `41aff00975fa2bf6d2da144ac56d5796ef97dbe7942d0e6686636516e20097d1`
- `output/studies/reduced_scalar_sobol/reduced_frozen_overnight_v2/comparison.csv`: `2c6e1acaceb7ea7e798bfe7f8166985c46c805b12daf45fdc16c93049fa910c2`
- `output/studies/reduced_scalar_sobol/reduced_frozen_overnight_v2/comparison.json`: `c79341b59f31b2591ebd1e7b0161f6b187be2baf4f8b7a236081b6b62316c1ed`
- `output/studies/reduced_scalar_sobol/reduced_frozen_overnight_v2/large_uv_frozen_0p5/A_m2/fit.json`: `1c7ba27dba5d1d1075ec6e2863367c2c2a11610612d66e8985b54b83868a9afe`
- `output/studies/reduced_scalar_sobol/reduced_frozen_overnight_v2/large_uv_frozen_0p5/A_m4/fit.json`: `31b31d9ae4b57c167acc6ff4b1579186234e79a1ed0e2a121d58fa4439f80991`
- `output/studies/reduced_scalar_sobol/reduced_frozen_overnight_v2/large_uv_frozen_0p5/A_p0/fit.json`: `5b16bc7519e6d81a0823daecf87c098ca67d6673a60f108c1168cc7ae5c3689a`
- `output/studies/reduced_scalar_sobol/reduced_frozen_overnight_v2/large_uv_frozen_0p5/B_m1/fit.json`: `2bc55237f9e36235cf04b10da218ec020f509f1eababd2da11b2aac02e1a4ea2`
- `output/studies/reduced_scalar_sobol/reduced_frozen_overnight_v2/large_uv_frozen_0p5/B_m3/fit.json`: `5b459a7eeb09ad0e0103ced674f7f234ec9e13d788e93bab71471283d4ab708d`
- `output/studies/reduced_scalar_sobol/reduced_frozen_overnight_v2/large_uv_frozen_0p5/B_p1/fit.json`: `d5102581e2b1766923a0861b8470dddf36de92e309f2249e28a4c58e81b13a83`
- `output/studies/reduced_scalar_sobol/reduced_frozen_overnight_v2/large_uv_frozen_0p5/C_m1/fit.json`: `60cd3577ca66f062d82efdee7dc7b4507e849563587907d82a5d1fa8865d3a83`
- `output/studies/reduced_scalar_sobol/reduced_frozen_overnight_v2/large_uv_frozen_0p5/C_m3/fit.json`: `d2fc3adaa966136e2b7e1eaf1a86d5e4b33c6fd8c4aec10f6a459c86fed9a494`
- `output/studies/reduced_scalar_sobol/reduced_frozen_overnight_v2/large_uv_frozen_0p5/C_p1/fit.json`: `69d090efd43ade0334954dfd6720a88450f56c477f7c5ce69cafc8edc2ccb62c`
- `output/studies/reduced_scalar_sobol/reduced_frozen_overnight_v2/large_uv_frozen_0p5/D_m2/fit.json`: `43afdb8cb420d5f16472795a1568180dcf9c6ca9aba127df484dadbfdaa7187e`
- `output/studies/reduced_scalar_sobol/reduced_frozen_overnight_v2/large_uv_frozen_0p5/D_p0/fit.json`: `c2ccc46c6b5c5a14f44ccf974303201b93e97697d68fed62f921c56a2b178052`
- `output/studies/reduced_scalar_sobol/reduced_frozen_overnight_v2/large_uv_frozen_0p5/D_p2/fit.json`: `0732767779dd68631dab3cf34d0e2ec47c5a1e16f7890e182d9ec9f3edd248c6`
- `output/studies/reduced_scalar_sobol/reduced_frozen_overnight_v2/large_uv_frozen_0p5/E_m2/fit.json`: `e0696992200e52867c415ae9cdccc857a28ffc51f7a4fca759044dadf9e6012e`
- `output/studies/reduced_scalar_sobol/reduced_frozen_overnight_v2/large_uv_frozen_0p5/E_m4/fit.json`: `c775ad326f58923c90f3f502a72b60922c9d28e1b9a612887ed4870d4cdad25c`
- `output/studies/reduced_scalar_sobol/reduced_frozen_overnight_v2/large_uv_frozen_0p5/E_p0/fit.json`: `8d52a30fa87983966b64c57a866767ef23e9658565cb5dee379907e86dd902b7`
- `output/studies/reduced_scalar_sobol/reduced_frozen_overnight_v2/large_uv_frozen_0p5/E_p2/fit.json`: `be8d5eed3f8f4c4563ea003808b1f34258b8edb76d73ae692b53cdac032322e4`
- `output/studies/reduced_scalar_sobol/reduced_frozen_overnight_v2/small_uv_frozen_5e-9/A_m2/fit.json`: `c90dbce59fcf98d9e6d9f7df34ea8c4b20e87ffd39d33c73bc68fb0f0a08a91c`
- `output/studies/reduced_scalar_sobol/reduced_frozen_overnight_v2/small_uv_frozen_5e-9/A_m4/fit.json`: `218b817a85c425eb7d0d62887f6d75e8ac535c49e51869d926c1b86237ea93c6`
- `output/studies/reduced_scalar_sobol/reduced_frozen_overnight_v2/small_uv_frozen_5e-9/A_p0/fit.json`: `ca112e1879f4b492d92bb8af604c6989a8b3d3438453f06b5a8aa939cc81ea1a`
- `output/studies/reduced_scalar_sobol/reduced_frozen_overnight_v2/small_uv_frozen_5e-9/B_m1/fit.json`: `f55128a39f49f9c01dad67bfe99b7ae218e0aa7d06a22cfe326eea377be7ff7e`
- `output/studies/reduced_scalar_sobol/reduced_frozen_overnight_v2/small_uv_frozen_5e-9/B_m3/fit.json`: `e8bf5c2ceb43010427875d7c25c65e9205cdb6c66125e1193b5168d4587363ba`
- `output/studies/reduced_scalar_sobol/reduced_frozen_overnight_v2/small_uv_frozen_5e-9/B_p1/fit.json`: `b7d0d7ee68c345c12ef2dc06b99c16b7b8174d14e97ca46051c2f76f95672f8a`
- `output/studies/reduced_scalar_sobol/reduced_frozen_overnight_v2/small_uv_frozen_5e-9/C_m1/fit.json`: `1deddc775da17e91f3ae6984770260be95d683a6b71c647025636b55eee3dc89`
- `output/studies/reduced_scalar_sobol/reduced_frozen_overnight_v2/small_uv_frozen_5e-9/C_m3/fit.json`: `6020b1f05c56a1dd4872248111f9862c3b873933a92bac2b0f11b57ba702aa90`
- `output/studies/reduced_scalar_sobol/reduced_frozen_overnight_v2/small_uv_frozen_5e-9/C_p1/fit.json`: `2178f6be48f188bf567675a57a0f47019ec984869013d41f0a0a0a6b8bd68bbb`
- `output/studies/reduced_scalar_sobol/reduced_frozen_overnight_v2/small_uv_frozen_5e-9/D_m2/fit.json`: `302550dc9b81df3f2d3c7f2eb2ac2d03a3fd9b542827c794c8eaea38e753ab00`
- `output/studies/reduced_scalar_sobol/reduced_frozen_overnight_v2/small_uv_frozen_5e-9/D_p0/fit.json`: `f857a7e89d2ebd31d33584f94b6e828b6107ee41d812932985b95225c95518b1`
- `output/studies/reduced_scalar_sobol/reduced_frozen_overnight_v2/small_uv_frozen_5e-9/D_p2/fit.json`: `747d7d1f71ce378f112aaa01fc18895fcdecac8f11d62d6316d99835f0fb737a`
- `output/studies/reduced_scalar_sobol/reduced_frozen_overnight_v2/small_uv_frozen_5e-9/E_m2/fit.json`: `3bd6ec78e4b75aa2450cf1a29797c80595018673a1ae64183163d6577b700275`
- `output/studies/reduced_scalar_sobol/reduced_frozen_overnight_v2/small_uv_frozen_5e-9/E_m4/fit.json`: `a0ecfd5570bc45b9d543d9c78322983956f44fb1d60de74d4f39f6fbc7965660`
- `output/studies/reduced_scalar_sobol/reduced_frozen_overnight_v2/small_uv_frozen_5e-9/E_p0/fit.json`: `3878489f24805f975de26cf2516c6d7f9193bbc2e1fab7ff126f9a1b0abf2ed1`
- `output/studies/reduced_scalar_sobol/reduced_frozen_overnight_v2/small_uv_frozen_5e-9/E_p2/fit.json`: `f15502531d6f266ef1d87bb1b6b13e0120f764295bb881276220bd7d96e564f6`
- `scripts/BuildT3ReducedScalarFinalResults.py`: `2013353affa1d57ae2f11e6e797cce26eeb1b78e02d5839c4ff6c3da39809309`
