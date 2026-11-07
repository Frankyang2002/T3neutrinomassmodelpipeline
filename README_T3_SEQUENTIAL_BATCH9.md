# Batch 9 v11 — ordered Matchete A5 versus half-normalized C5 correction

## Why mixed_inherited_10 failed

Both B(-1) and C(-1) read the archived native F-first one-loop LLHH sector
successfully (68 pole, 68 MF-log, 68 finite Matchete terms) but reported a
relative pole residual of exactly 0.500. The two common-scale checks also
failed. The inherited one-loop source was not missing or malformed.

**Verified source-level normalization mismatch**:
`Lagrangian/RunMatching.wl` defines
`ExtractWeinbergCoefficient[...]["Coefficient"]` as the **ordered** Matchete
coefficient `A5`. Its own `WeinbergOperatorConvention` states explicitly

```
L_EFT contains (1/2) C5^(ij) O_ij + h.c.
C5^(ij) = A5^(ij) + A5^(ji)
```

The previously executed Batch 8 numerical routine
`Numerical/running/T3MixedBatch6Numerics.py` symmetrized this ordered scalar
coefficient as `A_sym = (A5 + A5.T)/2`, then *labelled* its matrix `C5`. The
independent F-first one-loop component projector in
`BuildT3BCYFiniteBubbleInterface.experimental.py` and Batch 9 v10 already
returned **half-normalized symmetric C5**, from `A5 + A5.T`.

Therefore the source-dictated scalar contribution must be converted as

```
C5_scalar_pole   = 2 * old_Batch8_ordered_symmetric_pole
C5_scalar_finite = 2 * old_Batch8_ordered_symmetric_finite
C5_total         = C5_scalar + C5_F_first_inherited
```

This is **not** an empirical best-fit factor. The v11 runner reads and checks
`Lagrangian/RunMatching.wl` at runtime, and rejects an extractor with
unsupported semantics. Both pole and finite pieces receive the same factor.
The F-first one-loop contribution is NOT rescaled. The saved Batch 8 and v10
outputs are NOT rewritten.

The expected diagnostic if F and scalar residues differ by precisely the
factor-two convention is: v10 unconverted relative residual 0.5, followed
by a v11 converted pole residual near numerical zero. The converted scalar
and inherited finite terms must **also** have a negligible change under a
common matching-scale shift. If either condition fails, no combined finite
C5 is exported. The observed 0.5 alone is not a pole-cancellation proof.

## Source evidence and certification limitations

- Reuses `mixed_hard_08` native scalar Match and `mixed_fixed_08` conditional
  three-heavy-flavour lift, including their exact recorded source fingerprints.
- Reuses the archived F-first one-loop source under
  `output/general_observable/lo_certification_spinor/native_hard_e515d82d118284ba_native01/`
  when the legacy `output/general_t3_mass_eigen_matchete/eft1_fermion_first_oneloop/`
  source is missing. The fallback is allowed only with the committed native
  report, inspection report, 204-term inventory, and source SHA256 fingerprints.
- The F source must pass 68/68/68 pole/log/constant split, 34 PL/34 PR source
  field orientation, exact pole/log equality and source CG projection.
- The scalar source is independently recomputed from the unchanged native
  `hard_c5_candidate.txt` and verified against the prior Batch 8 result.
- The factor-of-two converter is used only after checking the **actual** local
  `Lagrangian/RunMatching.wl` source for the ordered-A5 return convention.
- The scalar pole and finite matrices now consistently use half-normalized
  symmetric C5 and are reported alongside the unconverted ordered-symmetric
  matrices, so the change is explicitly auditable.

A result of `CommonCouplingFInheritedAndScalarRawPolesCancelled_Conditional`,
if obtained, means only that these **projected raw** F-first and scalar
sectors cancel at common frozen couplings. It does **not** prove EFT1 RGE,
independent UV/EFT renormalization counterterms, full signed Majorana/SMEFT
matching, independently Matchete-matched three-MF first threshold, or an
unconditional Matchete-to-RGBeta quartic map. In particular, the archived
one-loop F expression has not been certified as a completely subtracted
physical hard threshold. Physical `C5` and physical neutrino masses remain
uncertified, regardless of a passing raw-pole check.

## Files to install

Complete files in the ZIP, preserving the following project paths:

- `scripts/RunT3MixedInheritedPoleBatch9.py` — replace v10, adds source-gated
  `A5 -> C5` normalization, separate original/corrected pole matrices and
  corrected common-scale finite test.
- `tests/python/test_t3_mixed_inherited_batch9.py` — replace tests, includes
  the exact factor-0.5 regression and a source-convention tamper rejection.
- `README_T3_SEQUENTIAL_BATCH9.md` — replace these instructions.

No native Matchete source, numerical Batch 8 implementation, earlier result,
or GitHub file is modified.

## PowerShell

```powershell
Set-Location 'C:\thesistemp\T3neutrinomassmodelpipeline'

Expand-Archive `
  -LiteralPath "$HOME\Downloads\T3_Batch9_A5C5Normalization_v11.zip" `
  -DestinationPath . -Force

python -X utf8 -m pytest -q tests/python/test_t3_mixed_inherited_batch9.py

python -X utf8 -m scripts.RunT3MixedInheritedPoleBatch9 `
  --hard-run mixed_hard_08 `
  --fixed-run mixed_fixed_08 `
  --check-only

python -X utf8 -m scripts.RunT3MixedInheritedPoleBatch9 `
  --hard-run mixed_hard_08 `
  --fixed-run mixed_fixed_08 `
  --run-id mixed_inherited_11
```

## Result files

```
output/general_sequential_pilot/mixed_inherited_pole_batch9/
  mixed_inherited_11/
    B_m1/result.json
    C_m1/result.json
    summary.json
```

Key JSON fields: `scalar_ordered_to_symmetric_C5_convention`,
`source_pole_check.ordered_batch8_unconverted_relative_residual`,
`source_pole_check.common_pole_residual_relative`,
`common_scale_finite_matrix_relative_change`,
`hard_scalar_ordered_A5_sym_pole_residue_gev_inverse`,
`hard_scalar_pole_residue_c5_gev_inverse`.

Only when the converted pole and scale checks pass are
`conditional_common_scale_sum_c5_gev_inverse` and
`conditional_common_scale_sum_singular_values_ev` written.

Required EFT ordering: UV → integrate out **all F** → EFT1(SM+X+Y) → integrate
out **both X,Y** → SMEFT. Nothing here changes that ordering.
