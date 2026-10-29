T3 Batch 1 readiness gate (additive, non-destructive)

Place scripts/AuditFullScalarReadiness.py in the existing project at:
  T3neutrinomassmodelpipeline/scripts/AuditFullScalarReadiness.py

Run from repository root (PowerShell):
  python scripts/AuditFullScalarReadiness.py

Report:
  output/audits/full_scalar_readiness.json

This tool checks the model inventory and extra-quartic registry consistency.
It does NOT claim full RG/matching coverage, and it does NOT change physics.
A full Batch 1 implementation requires a local copy of the repository source
and an installed Wolfram + RGBeta/Matchete runtime for end-to-end validation.
