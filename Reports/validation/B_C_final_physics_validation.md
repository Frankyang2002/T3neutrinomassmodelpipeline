# B(-1)/C(-1) final physics validation — source audit

Date: 2026-11-01
Repository: Frankyang2002/T3neutrinomassmodelpipeline
GitHub tree inspected: b7fe51ae12a68cd8b2e4eac7ca7923af33dff88f
Scope: one intermediate EFT; no repository source files changed.

## Conclusion
FAIL / NOT CERTIFIED: Direct one-loop UV amplitude equality with fully
renormalized two-threshold sequential EFT amplitude has not been demonstrated.
This is an evidence/completeness failure, NOT a demonstrated numerical
disagreement between two complete amplitudes.

## Source-backed observations
- `Lagrangian/RunGeneralT3EFT1ScalarThresholdProbe.wl` has
  `hard_pole_counterterms_combined -> False`,
  `hard_pole_rge_certified -> False`,
  `physical_sequential_C5_computed -> False`.
- `Lagrangian/RunT3BCFullEFT1Residual.experimental.wl` reports
  `CompleteSymbolicResidualExportedEquivalenceUnproven`,
  `physical_matching_certified -> False`,
  `renormalized_finite_C5_computed -> False`.
- `Numerical/diagnostics/T3BEndToEndControl.py` explicitly records
  `full_flavour_c12_c5_coupled_rge: NOT_IMPLEMENTED` and
  `scalar_threshold_resummed_matching: NOT_IMPLEMENTED`.
- `Numerical/diagnostics/T3BPerturbativeAccuracy.py` states that its
  perturbative bookkeeping is not complete NLO matching and sets
  `validated_physical_prediction_available: False`.
- `Lagrangian/RunT3UVFinalEvidence.py` labels its checks
  `ConditionalUVOperatorEvidencePassedNotFullCompleteness` and
  does not certify full Matchete/RGBeta operator normalization.

## Exact outstanding mathematical task
Choose a fixed order (leading nonzero one-loop amplitude, plus explicitly
specified RG improvement if desired) and one operator basis/scheme.

For each B(-1), C(-1), independently compute the symmetric 3x3 amplitude
A_full[LLHH] at vanishing external momenta (or a specified equivalent
projection) from the complete UV Lagrangian, for the same complex Yukawa and
scalar couplings, running masses and fermion mass basis used by the EFT.

Evaluate:
1. UV -> EFT1(F removed): tree C_LLSS and all hard one-loop
   contributions at the retained order; include pole subtraction and
   renormalized finite parts, with unambiguous Majorana symmetry factors.
2. EFT1: relevant one-loop operator counterterms/anomalous dimensions,
   including any mixing feeding C5 at the retained order, with no
   double counting of logs included at the thresholds.
3. EFT1 -> SM+C5(S1,S2 removed): full one-loop finite threshold projection
   including inherited operators, scalar interactions and required
   counterterms, not only a restricted scalar-bubble term.
4. Compare the full complex symmetric C5 matrices, not just masses or
   Frobenius norms, at identical parameters and scale. Verify SU(2)
   contractions, charged-lepton flavour basis, conjugations and
   `m_nu`/C5 vev normalization separately.

Define R = ||C5_direct-C5_seq||_F /
           max(||C5_direct||_F,||C5_seq||_F,1e-30 GeV^-1).
A proposed numerical acceptance threshold is R <= 1e-5
AFTER projecting both calculations onto the same EFT power and loop order,
and after establishing that errors from omitted powers and higher loops
are below that threshold; otherwise use an explicit theory-error budget.
Also require per-entry residuals for non-negligible elements, symmetry
residuals <= 1e-10, and tests with nontrivial complex inputs for BOTH
B(-1) and C(-1). These are proposed criteria, NOT tests executed here.

## Current benchmark classification
B(-1): [0.03930781275397567, 0.01578031440186623,
        0.0005343778626890098] eV
C(-1): [0.047061201781537984, 0.04676977395460441,
        0.022337982786261482] eV
Status: conditional pipeline outputs; NOT certified physical predictions.
Formal F-hard finite norms 1.529332607e-16 and 2.662294178e-16 GeV^-1
do not determine the missing full amplitude.

No verified source-level bug was repaired; modifying production formulas
without the counterterm calculation would be unjustified.
No independent quantitative direct-versus-sequential comparison was run.
The local Windows directory is not accessible from this environment.
