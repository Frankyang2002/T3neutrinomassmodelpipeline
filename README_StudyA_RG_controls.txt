Study A fixed-controls replacement (threshold scale-label correction)

Project path:
  Numerical/orchestration/RunCommonScalarRGControls.py

Replace the existing file at this path in your LOCAL project.
No files on GitHub are modified.

Cause of failure:
  build_intermediate_scalar_boundary() requires UV state.mu_gev to equal the
  fermion matching scale. The frozen control passed the original UV-scale
  state, so all 80 comparisons failed before the matching formula ran.

Correction:
  The frozen control now uses dataclasses.replace() to set ONLY mu_gev to the
  fermion scale; projects through the existing fermion boundary; relabels the
  intermediate state to the scalar matching scale without integrating RGEs;
  then calls the existing authoritative final-C5 adapter.
  Model couplings, masses, gauge couplings and charged-lepton Yukawas remain
  frozen at their input values. The one-loop final-C5 matching remains active.
  No new matching formula is introduced.

This is a diagnostic frozen-parameter control, NOT a full fixed-order
amplitude and NOT a guarantee that all points will pass numerical validation.

From the project root in PowerShell:

  python -m Numerical.orchestration.RunCommonScalarRGControls --dry-run
  python -m Numerical.orchestration.RunCommonScalarRGControls --run --limit 2
  python -m Numerical.orchestration.RunCommonScalarRGControls --run

Result directory:
  Reports/output/full/common_scalar_rg_controls/smallY_smallL/

The program rewrites the diagnostic results directory for the selected models;
original common-scalar running scan data remain untouched.

Validation performed on delivered file:
  Python syntax compilation and ZIP content integrity. Full integration with
  local project artifacts has NOT been run in this environment.
