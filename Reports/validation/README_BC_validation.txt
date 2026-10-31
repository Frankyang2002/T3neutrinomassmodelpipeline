Run from repository root:
  python -u -X utf8 scripts\RunT3BCLeadingOrderValidation.py

Report: Reports\validation\bc_leading_order_validation_result.json
Native report: output\audits\t3_bc_direct_vs_sequential_full_mass_v22\summary.json

This runs the repository's existing Workstream C v22 comparison and inspects
both model passes. It intentionally never declares the RG-evolved benchmark
masses physically certified. If archived Matchete evidence is missing it fails
and records the missing input instead of fabricating an answer.

The v22 computation is the available frozen-parameter direct vs sequential
one-loop comparison; it does NOT establish physical correctness of the five-
mass RG-evolved B/C spectra by itself. The independent Lagrangian amplitude,
inter-package coupling map and complete physical threshold prescription remain
to be checked before those numbers are certified.
