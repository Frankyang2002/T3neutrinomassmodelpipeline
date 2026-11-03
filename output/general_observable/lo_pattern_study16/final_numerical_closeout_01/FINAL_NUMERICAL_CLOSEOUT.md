# Final numerical closeout — T3, 16 provisional one-loop benchmarks

**Numerical phase: CLOSED FOR SOURCE-CONDITIONAL RESULTS ONLY.**

**Not a physical neutrino-mass certificate:** the Matchete/RGBeta →
left-Weyl signed action normalisation and basis/phase conversion remain
unresolved. Physical certification is **0/16**. No EFT, RG running,
threshold matching, oscillation fit or experimentally validated masses
are inferred by this closeout.

## Frozen results (all masses in eV; source-conditional)

| Model | F | m1 | m2 | m3 | m1/m3 | A1 | C1 | Mechanism | B6 baseline eps |
|---|---|---:|---:|---:|---:|---:|---:|---|---:|
| A_-4 | Dirac | 0.000474483 | 0.00758273 | 0.0436006 | 0.0108825 | 0.102495 | 0.057116 | mixed_or_neither_dominant | 6.45e-16 |
| A_-2 | Dirac | 0.00721475 | 0.0110085 | 0.0485602 | 0.148573 | 0.168426 | 0.480645 | mode_geometry_dominant | 0 |
| A_+0 | Dirac | 0.00721427 | 0.011 | 0.048499 | 0.148751 | 0.168679 | 0.480074 | mode_geometry_dominant | 6.36e-16 |
| B_-3 | Dirac | 0.00721823 | 0.0109814 | 0.0483324 | 0.149346 | 0.168877 | 0.480121 | mode_geometry_dominant | 4.53e-10 |
| B_-1 | Majorana | 0.00452249 | 0.00980338 | 0.0504413 | 0.0896585 | 0.157008 | 0.208537 | mixed_or_neither_dominant | 4.53e-10 |
| B_+1 | Dirac | 0.00721583 | 0.0109461 | 0.0480743 | 0.150097 | 0.159493 | 0.505219 | mode_geometry_dominant | 4.53e-10 |
| C_-3 | Dirac | 0.0269469 | 0.0281923 | 0.0550407 | 0.48958 | 0.268784 | 1 | mode_geometry_dominant | 4.53e-10 |
| C_-1 | Majorana | 0.0374973 | 0.0385889 | 0.062285 | 0.602027 | 0.5691 | 0.285778 | mixed_or_neither_dominant | 4.53e-10 |
| C_+1 | Dirac | 0.000365488 | 0.00826871 | 0.0477861 | 0.00764841 | 0.319197 | 0.0131221 | heavy_eigenspace_cancellation_dominant | 4.53e-10 |
| D_-2 | Dirac | 0.00721417 | 0.0109978 | 0.0484824 | 0.1488 | 0.169299 | 0.478505 | mode_geometry_dominant | 5.69e-16 |
| D_+0 | Dirac | 0.00721467 | 0.011007 | 0.0485499 | 0.148603 | 0.168458 | 0.480577 | mode_geometry_dominant | 5.13e-16 |
| D_+2 | Dirac | 0.000474733 | 0.00758174 | 0.0435948 | 0.0108897 | 0.10255 | 0.0571179 | mixed_or_neither_dominant | 1.49e-15 |
| E_-4 | Dirac | 0.0033194 | 0.00741694 | 0.0384885 | 0.086244 | 0.157262 | 0.214 | mixed_or_neither_dominant | 3.93e-10 |
| E_-2 | Dirac | 0.00471419 | 0.00867883 | 0.0424695 | 0.111002 | 0.0389577 | 1 | mode_geometry_dominant | 3.93e-10 |
| E_+0 | Dirac | 0.00464358 | 0.00864859 | 0.0425154 | 0.109221 | 0.0387518 | 1 | mode_geometry_dominant | 3.93e-10 |
| E_+2 | Dirac | 0.00332597 | 0.00744248 | 0.0386336 | 0.0860901 | 0.157661 | 0.213222 | mixed_or_neither_dominant | 3.93e-10 |

Values are the **original frozen** candidate left-Weyl
one-loop results. `A1`, `C1` and the mechanism labels are the
Batch 5 source-conditional mode diagnostics; `eps` is the
Batch 6 Yukawa-independent coefficient residual at a common
|lambdaT3|=0.10 comparison (not the frozen-model physical error).

## Findings justified by the stored numerics

- **16/16** frozen candidate mass matrices evaluated;
  min m1/m3 = 0.00764841 (C_+1);
  max m1/m3 = 0.602027 (C_-1).
- **C_+1:** m1/m3 = 0.00764841;
  heavy-eigenspace coherence C1 = 0.0131221.
  The saved phasor analysis attributes its small lightest
  singular value to destructive heavy-eigenspace interference
  within the candidate Weyl convention.
- **288/288** recorded scan points passed: 16 models ×
  [5 y1 + 5 y2 + 5 lambdaT3 + 3 heavy-mass points].
  All 48 y1/y2/lambdaT3 zero points vanish;
  all 64 t=1 points reproduce the frozen
  values. Max positive Yukawa-rescaling relative mass error:
  4.09e-14.
- **64** common-mixing Yukawa-texture comparisons are preserved
  separately in `controlled_texture_comparisons_64.csv`.
  These are *counterfactual* substitutions, not new fitted models;
  they isolate which texture effects persist for fixed candidate
  representation/fermion constructions.
- **80** Batch 6 coefficient comparisons and 80 independent
  tensor-reconstruction checks passed the recorded grid.
  Maximum reconstruction residual:
  4.79e-16.

## Batch 6 coefficient-comparison scenarios

| Scenario | Max eps | Max heavy-group spread | Interpretation |
|---|---:|---:|---|
| baseline | 4.53e-10 | 5.952e-10 | Conditional fixed-source comparison |
| scalar_common_shift | 4.359e-10 | 5.705e-10 | Conditional fixed-source comparison |
| heavy_generation_tilt | 2.365e-09 | 4.282e-09 | Conditional fixed-source comparison |
| large_scalar_mixing | 1.818e-07 | 2.382e-07 | Conditional fixed-source comparison |
| charge_split_diagnostic | 0.0002235 | 0.0002929 | Diagnostic only; SU(2)-split masses outside restricted ansatz |

Near-zero coefficient residual at the chosen masses establishes
numerical proportionality of the stored candidate coefficient
tensors **at those parameters**; it is not proof for arbitrary
masses or a physical cross-package normalisation.

## Source and evidence boundaries

The run compares the frozen summary JSON/CSV, all 288 scan CSV
rows, Batch 5 mass/texture diagnostics, and Batch 6 coefficient
JSON/CSV and reconstruction grid. It checks the recorded
study→Batch 5→Batch 6 SHA256 lineage, including Git CRLF-only
checkout normalization when applicable. It does **not** remeasure
native Matchete/RGBeta source actions or claim current-HEAD
equivalence to the frozen source revision.

| Input/evidence | SHA256 |
|---|---|
| Frozen study JSON | `05017be5a8ee07711ef2a2f7ca11f1b2ecabadbfefc4f47563a77177b15a065e` |
| Batch 5 JSON | `5839573f1088ae7c33024f08bb81b79b8f2bf296c8b23f253a8c433e8ee09534` |
| Batch 5 manifest | `9d624bf26fd66890d07365b851fbe3160b6685b026ec6efe4b0ad4c525037fdb` |
| Batch 6 JSON | `94a3e72b598a877fe90e32b5de878acfb30fc5f692ed3ecf91c26dc824344587` |
| Batch 6 manifest | `1a893d819f16c862bfd546b3fdea26f89aece30273ad496af90906b0c3cab477` |
| Frozen scan CSV | `d39cf6ffcef91ac5b42b3f780882d25eefcdd8c6668b1cbae41c5fe7bbeff3b0` |

Original frozen-study recorded commit: `0b210f18f6a2cfb10fff03059e2ec4a73d2195a4`.
Its earlier revision is retained as provenance, not imposed as the
current local checkout HEAD.

## Stop rule

**The requested rough numerical milestone is complete** if
this closeout validates against the Windows checkout. Do not
continue SU(2) phase, physical normalisation or EFT work as part
of this milestone. Reopen those only if independently required
by the thesis scope or by a specific physical claim.
