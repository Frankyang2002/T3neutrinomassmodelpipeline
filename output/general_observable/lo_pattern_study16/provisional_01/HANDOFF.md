# Chat A handoff — restricted all-16 direct UV

- Reviewed GitHub main: `0b210f18f6a2cfb10fff03059e2ec4a73d2195a4`; actual local HEAD: `0b210f18f6a2cfb10fff03059e2ec4a73d2195a4`; dirty worktree: `True`.
- This run: `16/16` conditional UV evaluated, `0/16` blocked, `0/16` physically certified.
- Input provider: `C:\thesistemp\T3neutrinomassmodelpipeline\output\general_observable\lo_direct_one_loop\provisional_inputs`; frozen input and audit SHA256 are recorded in per-model JSON.
- Exact conditional computation: existing `core.evaluate` in canonical left-Weyl/complex-scalar basis; independent reference uses `M^dagger f(M M^dagger)` and scalar Cartesian Hessian.
- Source-CG rule (as audited): `T3Y1CG * T3Y2CG * conjugate(T3MixCG)` in common source-weight indices; NOT certified cross-package normalisation.
- Yukawa scans: conditional linearity for topology-only cross terms, tested numerically; lambda scan varies only `A[S1,S2]` off-diagonal; mass scan scales heavy `M`.
- Critical unresolved: native four-spinor chirality/phase/normalisation, RGBeta/Matchete scalar portals/T3 sign and component-axis mapping, physical fermion multiplicities.
- This package **is not safe for asserting physical neutrino masses**. It is safe for numerical regression and source-conditional structure. No EFT/running/threshold matching.
- Plot warnings: [].
- No previously existing output directory overwritten.
