# A(-4) production source verification: two regression repairs and full-beta audit

The 2026-11-08 production run regenerated 17/17 RGBeta EFT1 channels, matched all nine complete native quartic betas to the previous *both fresh* experiment at three nonzero gauge/Yukawa states, passed the independent 55 scalar-only checks, and integrated seven scale points. The accompanying structural test file had two incorrectly formulated assertions; neither indicates failed code physics.

**Replacement**: `tests/test_a_minus4_fixed_source_batch.py` (whole file). It uses full-file round-trip reconstruction of the unchanged original from the new source by reversing exactly two anchored invariant-body regions, rather than expecting unmodified declarations in diff output. It inspects the parsed Python `ValueError` message for the physical Wilson guard, rather than looking for two adjacent source string literals as one contiguous text sequence.

**New script**: `scripts/InspectAminus4FixedFullBetaDifferences.py`. It uses existing local JSON outputs, requires no Wolfram calls and no ODE rerun, and creates one diagnostic report. At three scale points and nine quartic beta channels, the script evaluates both the permanent-fix beta and earlier scalar-only-diagnostic beta on *identical numerical coupling inputs*. It uses four configurations (all gauge/Yukawa present; gauge only; Yukawa only; neither) to calculate

- `scalar_only = Δβ(0,0)`;
- `gauge_dependent = Δβ(g,0) - Δβ(0,0)`;
- `yukawa_dependent = Δβ(0,y) - Δβ(0,0)`;
- `joint_gauge_yukawa_interference = Δβ(g,y) - Δβ(g,0) - Δβ(0,y) + Δβ(0,0)`.

The sum equals the complete fixed-point difference. The script also compares endpoint values from the three already-solved trajectories. This distinguishes beta-expression changes from the accumulated nonlinear ODE response. It does **not** establish that RGBeta and Matchete couplings are physically identical, does not solve a Wilson RGE, and does not match the scalars.

**New tests**: `tests/test_a_minus4_fixed_beta_origin.py` (synthetic algebra-only checks of inclusion-exclusion).

All produced outputs are read-only diagnostics; the original corrected RGBeta source and its source hash are unchanged. This batch is intentionally restricted to closing the failed regression and explaining the unexpectedly large complete-beta vs scalar-only-diagnostic difference.
