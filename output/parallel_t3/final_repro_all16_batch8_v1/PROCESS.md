# Final T3 numerical benchmark process (conditional prescription)

1. Saved UV/fermion matching inputs: completed Batch 1 source-resolved heavy-F triangle.
2. All three heavy fermions are integrated out at the common MF = 100000 GeV.
3. EFT1 retains both heavy scalar multiplets; strict one-loop Wilson channels are evolved to MS = 1000 GeV.
4. Both scalars are integrated out at MS using the existing full-tensor scalar bubble threshold.
5. The resulting SMEFT C5 evolves to the saved low-energy scale and gives mnu = -v^2 C5/2.
6. Batch 7 retains named native Laurent poles, explicit EFT1 scalar counterterms, source matrices, and three conditional matching prescriptions over three subtraction scales.
7. This report validates all 16 model outputs, 14 saved prerequisite source hashes per model, per-model data hashes, nominal threshold identities and masses derived again from stored C5(low). No source is recalculated or altered during collection.

**Status:** Complete and reproducible *conditional numerical benchmark*. No independent bare-UV amplitude, fitted missing counterterm, full physical MSbar matching, or experimental fit is asserted.

**Scenarios:** `analytic_triangle` (the primary conditional source-triangle result); `native_formal` (saved Matchete pole-dropped prescription); `native_as_UV_minus_soft_hypothesis` (comparison hypothesis). Subtraction-scale variation uses MF-frozen couplings and is *not* a new physical threshold.

**Eigenvalues:** m1 <= m2 <= m3 are ascending Takagi singular values; they are not a normal/inverted-ordering identification. Neutrino masses are in eV, C5 in GeV^-1.

**Outputs:** `all16_scenarios.csv` (144 rows), `all16_nominal.csv` (48 rows, nominal subtraction scale), `all16_overview.csv` (16 primary nominal rows), `all16_wilson_matrices.npz` (complex 3x3 Wilson matrices), `manifest.json` (provenance and hashes).
