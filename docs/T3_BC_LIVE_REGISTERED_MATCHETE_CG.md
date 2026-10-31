# B/C(−1) `110101`: independently read the registered Matchete CG tensors

## Why this audit is the last convention-specific gate

The preceding local reports established: 128/128 live RGBeta `Lam` component
coefficients with zero residual; all four live Matchete `110101` interactions
validated by `CheckLagrangian`, plus both combined interactions; and 128/128
exact signed action coefficient identities in a *saved exported* SU(2)
component basis (`bc_live_action_bridge_v1`). The remaining missing
component evidence is the tensor currently **registered under the actual
`CG` symbols** used in each Matchete `PlusHc` interaction.

This probe calls Matchete's public `GetCGTensor[cg]` on
`T3ExtraCGQ110101I1` and `T3ExtraCGQ110101I2` immediately after
`T3ExtendedBuildExtraQuarticCandidates` registers them, separately for
B(−1) and C(−1). It verifies each symbol occurs in the selected live
interaction. **No Stage22 source basis is substituted into the Wolfram
probe.** It records the 16 exact InputForm components and their dimensions.
Any API failure or unexpected shape is **blocked**, with a JSON diagnostic.
The full build and CG extraction stay inside `UsingFrontEnd`.

The Python verifier then checks each of the four live tensors against the
saved exact native Matchete Stage22 records, checks the preceding signed
bridge report's input hashes, and independently repeats the signed
component identity using the live Matchete components:

\[
(R_1,R_2)=(M_1,M_2)
\begin{pmatrix}\tfrac12&1\\-\tfrac{\sqrt3}2&0\end{pmatrix},
\qquad h=-Ag,
\qquad \mathcal{L}_M=-V_{\mathrm{RGBeta}}.
\]

It compares 32 exact registered CG components per model against Stage22,
and 64 direct + conjugate action coefficients per model. No tolerance is
used for the component identities; negative controls reject an overall sign
error and a spurious `/4`.

## Files (all local only)

- `scripts/ProbeT3BCMatchete110101RegisteredCG.wl` (new)
- `scripts/VerifyT3BCMatchete110101RegisteredCG.py` (new)
- `tests/test_t3_bc_live_matchete_110101_registered_cg.py` (new)
- `docs/T3_BC_LIVE_REGISTERED_MATCHETE_CG.md` (new)
- `docs/T3_GENERAL_SCALAR_MATCHING_RESOLUTION_LOG.md` (full updated replacement)

Requires previous audit scripts and existing output files from your passing
B/C tests. No production files or saved reports are modified.

## PowerShell (repository root)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "C:\path\to\T3_BC_live_registered_matchete_cg_v1.zip" -DestinationPath . -Force
python -X utf8 -m pytest -q tests/test_t3_bc_live_matchete_110101_registered_cg.py
wolframscript -file scripts/ProbeT3BCMatchete110101RegisteredCG.wl output/audits/t3_bc_matchete_registered_cg/live_registered_cg_v1.json
if ($LASTEXITCODE -ne 0) {
    if (Test-Path output/audits/t3_bc_matchete_registered_cg/live_registered_cg_v1.json) {
        Get-Content output/audits/t3_bc_matchete_registered_cg/live_registered_cg_v1.json
    }
    throw "Registered Matchete CG extraction failed"
}
python -u -X utf8 -m scripts.VerifyT3BCMatchete110101RegisteredCG --input output/audits/t3_bc_matchete_registered_cg/live_registered_cg_v1.json --run-id bc_live_registered_cg_v1
```

New output: `output/audits/t3_bc_matchete_registered_cg/bc_live_registered_cg_v1/summary.json`. The probe refuses to overwrite its input JSON; the verifier refuses to overwrite its summary.

## What this does *not* establish

This test directly verifies the **registered Matchete CG tensors** against
Stage22 and connects them to the already validated RGBeta component identity.
It does not subtract the two full evaluated runtime Lagrangians in one
Mathematica session or establish the Yukawa component conventions. It does
not insert the mixed scalar bilinear into Matchete `HeavyMassBasis`, invoke
`Match[]`, do full intermediate Wilson running or scalar threshold matching.
Consequently production `physical_coupling_mapping_certified`,
`mass_basis_matching_certified`, and `sequential_eft_matching_certified`
remain **false**. After a passing registered CG audit, stop expanding the
`110101` quartic checks and prioritize B/C mass-eigenstate Yukawa and
mixed-heavy matching.
