T3 GENERAL SCALAR THRESHOLD HANDOFF (conditional)

Copy both files into the project repository's scripts/ directory:
  scripts/RunT3GeneralScalarThreshold.py
  scripts/RunT3GeneralSMEFTBridge.py

This connects:
  saved general scalar UV/EFT1 numerical endpoint
  + separately available FULL canonical EFT1 Wilson tensor at MS
  + separately available FULL canonical quartic tensor at MS
  -> existing scalar-bubble matcher -> existing SMEFT runner -> PMNS.

The tensor inputs MUST be in identical canonical real-interleaved
Weyl45/real scalar coordinates. A native Matchete tensor, RGBeta
coefficient list, or source-basis tensor cannot be passed as-is.
This script does not convert a noncanonical source frame.

The scalar mass matrix is rebuilt from the general trajectory's
mS1Sq, mS2Sq, complex extra mu2, using the existing general
quadratic implementation, then permuted from grouped (all R, all I)
to canonical interleaved (R1,I1,R2,I2,...).
U(1)_Y and SU(2) Ward residuals are checked against the actual
canonical representation generators. Do not use the opt-out flags
for physical/theoretical claims.

PowerShell, from project root:
  python -X utf8 -m scripts.RunT3GeneralScalarThreshold --help

Example with REAL, COMPATIBLE PRODUCER ARTIFACTS (replace paths):
  python -X utf8 -m scripts.RunT3GeneralScalarThreshold `
    --trajectory "output/studies/general_scalar_extended/continuous_full_scalar/RUN/A_m4/continuous_trajectory.json" `
    --wilson "PATH/TO/full_wilson_at_MS.npz" `
    --quartic "PATH/TO/full_quartic_at_MS.npz" `
    --assume-same-canonical-frame `
    --output "output/studies/general_scalar_extended/connected/A_m4"

The "same canonical" flag is an explicit assumption; it does not
automatically establish tensor provenance. The Wilson must already
have been evolved after integrating out all fermions; the quartic
must be an actual full component tensor, NOT reconstructed by
incorrectly treating raw RG beta coordinates as tensor components.
If the matching tensor sources are not yet generated, this command
must fail rather than invent them.

Outputs:
  output/manifest.json
  output/general_scalar_threshold_c5.npz
  output/smeft/summary.json
  output/smeft/smeft_low.npz

Limitations: first threshold and EFT1 Wilson-running inputs are external;
scalar threshold is conditional bubble only; extra finite contributions
are not supplied; exact native coefficient mapping is not assumed.
