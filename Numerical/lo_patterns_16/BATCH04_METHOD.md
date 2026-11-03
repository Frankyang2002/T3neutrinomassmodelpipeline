# Chat A, Batch 4 — analytical flavour geometry of all 16 T3 spectra

## Source and status

This is a **read-only postprocessor** of frozen, restricted source-conditional outputs:

1. `provisional_01/all16_summary.json`, `provisional_01/models/*.json`;
2. `analysis_batch02_01/batch02_analysis.json`;
3. `analysis_batch03_01/batch03_analytic_explanation.json`, its model files and `SOURCE_MANIFEST.json`.

It reuses Batch 3's independently reconstructed heavy-eigenspace matrices; it does **not** call `core.evaluate`, regenerate Weyl inputs, refit Yukawas, match EFTs or change physical-source status. **0/16 physically certified.** It checks the original executable-source hashes, CG audit hashes and benchmark/Weyl-input hashes via the existing Batch-2 guard, and verifies all source summaries and matrices again.

## Exact heavy-eigenspace identity

For each *distinct* heavy singular mass, Batch 3 constructs a spectral projector of the mass-squared matrix, then evaluates the exact scalar charge kernels to obtain a symmetric complex matrix `H_g` in eV:

\[
M_\nu=\sum_g H_g,\qquad
H_g=\frac{1}{16\pi^2}\sum_q\sum_{k\in g} K_q(m_k^2)\,
\mathrm{sym}\!\left[Z_{s_q}(M_F^\dagger P_k)Y_{t_q}^{T}\right].
\]

A degenerate heavy eigenspace is treated as **one** projector sum. A decomposition into individual degenerate Takagi columns would not be basis invariant. The source scalar Hessian has no `B` mixing under the restricted benchmark assumptions; the scalar loop kernel from Batch 3 is exact for its isolated 2×2 charge blocks.

For a complex matrix SVD `Mν = U diag(m1,m2,m3) V†`, with ascending singular values and the corresponding left/right singular vectors `u_i,v_i`, define:

\[
a_{gi}=u_i^\dagger H_g v_i,\qquad \sum_g a_{gi}=m_i,\qquad
G=\sum_g\lVert H_g\rVert_F.
\]

Thus the **exact** identity for each neutrino singular mode is

\[
\boxed{\frac{m_i}{G}=\underbrace{\frac{\sum_g|a_{gi}|}{G}}_{A_i:\text{ flavour projection}}
\underbrace{\frac{m_i}{\sum_g|a_{gi}|}}_{C_i:\text{ coherent sum}}.}
\]

`A_i` is bounded by unity via the spectral/Frobenius norm inequality. Small `A_i` means even the separate heavy-group terms project weakly onto the final singular mode (flavour alignment or nearly null direction). Small `C_i` means their complex amplitudes substantially cancel in that direction. This is **not** a model-independent distinction between dynamical mechanisms; it is an invariant split for the specified heavy-eigenspace decomposition.

For numerical description only:
- `A_i < C_i/2`: `mode_geometry_dominant`;
- `C_i < A_i/2`: `heavy_eigenspace_cancellation_dominant`;
- otherwise: `mixed_or_neither_dominant`.

These thresholds are heuristic descriptive labels, **not** physical significance tests. If the projected amplitude sum is smaller than a floating-point threshold, the coherence is unresolved.

## Other exact and conditional identities

The matrix Frobenius group coherence is `||Mν||F / Σ_g ||H_g||F`, which is distinct from both mode coherence and E-family charge-channel alignment. Pairwise overlaps are `Re tr(H_g† H_h) / (||H_g||F ||H_h||F)`. The determinant identity `|det Mν| = m1 m2 m3` characterises the hierarchy through `|det Mν| / ||Mν||F³`. None of these are obtained by summing eigenmasses of separate matrices.

For every group, deleting `H_g` and recomputing all singular values diagnoses nonlinear sensitivity of the *fixed matrix sum*. This operation does **not** represent a consistent UV model in which one physical fermion was removed, because it does not update couplings, gauge representations, masses or counterterms.

## Questions addressed

- **C_-1 versus C_+1:** Both have a single active topology charge, yet one has a much less suppressed `m1/m3`. Batch 4 measures whether this originates from small `A_1`, small `C_1`, or both. The frozen Yukawas and heavy-field construction differ (C_-1 Majorana, C_+1 Dirac), so the result **does not** establish a one-variable causal comparison.
- **C_-3 versus C_-1:** Diagnose whether the second large-lightest-mass C spectrum shares the same mode-level suppression pattern.
- **A/D pairs:** Compare mode geometry for cases already near identical in the frozen spectra; any equalities are benchmark-specific unless independently proven symbolically.
- **E channels:** The two scalar-charge contributions were seen to align in Batch 3. Batch 4 separately tests heavy-generation cancellation, which is logically compatible with perfect charge-channel alignment.

The generated `BATCH04_FLAVOUR_REPORT.md` gives the actual 16-model quantitative mechanism classification. The statements above define what will be tested; they do not preclaim unmeasured mode-level results.

## Outputs (all newly generated under `analysis_batch04_01`)

- `BATCH04_FLAVOUR_REPORT.md`: central all-16 analytical interpretation, targeted C comparison and per-model explanations.
- `batch04_flavour_analysis.json`: all models and singular-mode diagnostics.
- `batch04_model_mechanisms.csv`: selected model parameters, spectra and group-coherence summaries.
- `batch04_mode_projections.csv`: real and imaginary complex amplitudes for each mass mode and heavy group.
- `batch04_heavy_groups.csv`: projector multiplicities, group matrix Frobenius norms and ranks.
- `batch04_heavy_group_overlaps.csv`: complex pairwise matrix overlaps.
- `batch04_leave_one_group_out.csv`: explicitly diagnostic omission experiment.
- `SOURCE_MANIFEST.json`, `HANDOFF.json`: source hashes, status and non-certification guard.

The script fails closed on missing cases, changed source bytes, spurious certification, mismatched matrices, non-symmetric groups, changed heavy spectra/multiplicities and existing nonempty output directories. It writes only after validating all 16.
