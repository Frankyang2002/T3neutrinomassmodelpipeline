# Experimental Stage 7: fixed formal C5 atom/linearity audit

## Why Stage 6 stopped

Your log shows 37 Python checks passed, then the formal response audit returned
`ReviewRequired` for B(-1) and C(-1). Importantly, the previous algebraic proof
already established `C5_extended|extra=0 == C5_baseline` in both cases; a false
linearity audit does **not** overturn that identity. The Stage 6 report also
identified the expected `y1`/`y2` terms.

The old algebraic extraction used `Expand[C5_extended-C5_baseline]` to prove
linearity even though `C5_extended|extra=0-C5_baseline` was proved zero through
potentially stronger simplification. It also used a recursive `Cases` pattern
that counted an underlying direct `Coupling[...]` *inside* `Bar[Coupling[...]]`
as if it were an independently present atom, producing empty derivative rows.

## Changes

1. `CompareGeneralT3C5.wl`: still proves the complete zero-extra limit against
   the independent baseline; extracts actual response separately as
   `C5_extended-C5_extended|extra=0`.
2. Conjugated coupling atoms are masked before searching for direct couplings,
   preventing duplicated ghost atoms.
3. Checks polynomiality, coefficient independence, vanishing zero-response,
   and linear reconstruction residual using exact symbolic reduction. Filters
   mathematically zero coefficients out of the physical response report.
4. Exports residual diagnostics when checks fail rather than silently passing.
5. `AssessGeneralT3C5Response.py`: demands the new proof fields. Old Stage 6
   JSON files *cannot* be misinterpreted as a successful Stage 7 audit.
6. `RunGeneralT3Stage6.ps1`: executes a standalone Wolfram atom self-test and
   regenerates formal reports from the *existing* matching outputs by default.
   Does **not** rerun any 16-model matching.

## Windows PowerShell commands (project root)

```powershell
Expand-Archive -LiteralPath .\T3_General_Scalar_Batch7.zip -DestinationPath . -Force
python -m pytest -q tests\python\test_general_scalar_trial.py tests\python\test_scalar_quadratic_prototype.py tests\python\test_general_scalar_c5_proof.py tests\python\test_general_rgbeta_coverage.py tests\python\test_general_t3_stage6.py
powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralT3Stage6.ps1
```

Artifacts:
- `output/general_t3_trial/c5_algebraic_report.json` (refreshed)
- `output/general_t3_trial/B_m1/general_c5_algebraic/delta_c5_atoms.json` (refreshed)
- `output/general_t3_trial/C_m1/general_c5_algebraic/delta_c5_atoms.json` (refreshed)
- `output/general_c5_response.json` (refreshed)
- `output/general_rgbeta_source_proof.json` (if C5 proof passes)
- `output/general_t3_trial/matchete_active_quartic_tensors.json` (if all prior stages pass)

Expected (not yet verified on your Matchete artifacts): `VerifiedFormalResponse`,
with two active quartic names for B/C(-1). If that still fails, the exported
`response_reconstruction_residual_inputform` and `response_at_zero_inputform`
fields will identify whether the symbolic coefficient is truly nonlinear or
has an algebraic-form issue; do not treat such a failure as a model result.

### Important scope

This remains a **scalar-quartic-only** symbolic extension. UV scalar quadratic
mixing, full RGBeta coupling coverage and numerical propagation have not been
implemented or verified. In particular, Matchete/RGBeta normalisation factors
must be determined before using the additional symbolic couplings numerically.
