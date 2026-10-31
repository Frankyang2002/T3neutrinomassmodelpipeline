
## Experimental B/C scalar-bubble finite-matching diagnostic (2026-10-31)

* New files: `Numerical/experimental/T3BCScalarBubbleThreshold.py`, `scripts/RunT3BCScalarBubbleCandidate.py`, `tests/test_t3_bc_scalar_bubble_threshold.py` and the technical note `docs/T3_BC_SCALAR_BUBBLE_MATCHING_DIAGNOSTIC.md`.
* Extends the already completed conditional EFT2 sequential numerical stage with a source-normalized one-loop scalar-pair bubble candidate at the Y crossing. Provides full-internal-mass and massless-X hard-limit variants, while preserving the full lepton flavour pair array. Both candidates are then propagated through the existing EFT2 and SM+Weinberg numerical solvers.
* Internal checks: heavy-leg quartic beta partition, `d B0bar/d ln(mu)=2`, exact input provenance, EFT2 quartic fit parity, no-correction full numerical replay, symmetry and Weinberg projection.
* The Y and X physical finite matching, complete counterterms, scalar renormalizable threshold matching, higher-dimensional operators, complete three-generation Yukawa Wilson beta and all RGBeta/Matchete quartic convention certification remain **unresolved**. No physical neutrino mass claim is made.
* Previous study outputs and production files are unchanged. This milestone must not change `full_general_scalar_matching_certified` or `physical_neutrino_mass_prediction_certified` to True.
