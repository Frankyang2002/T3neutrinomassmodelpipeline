# B/C(-1): independent live Matchete construction of 110101 quartic terms

## Purpose and current evidence

The user ran the `native_lam_coeff_v4` probe against installed RGBeta v1.2.0,
then `scripts.VerifyT3BCNativeLam110101Coefficients` successfully:

| Model | Direct/h.c. SU(2) components | Mismatches | Maximum residual |
| --- | ---: | ---: | ---: |
| B_-1 | 64/64 | 0 | 0 |
| C_-1 | 64/64 | 0 | 0 |

The live RGBeta tensor entries are therefore verified at coefficient level
for the selected two `110101` couplings, with **ordered real-tensor
coefficient `1/4`**. This does not alone certify Matchete's action sign,
absolute coupling map, the mixed scalar mass-basis Lagrangian, or matching.

The next independent gate evaluates the Matchete side **using the existing
T3 Lagrangian builder**. The new script:

1. Loads the installed Matchete package within `UsingFrontEnd`.
2. Builds the project's B_-1 and C_-1 representations/fields/couplings.
3. Calls `T3ExtendedBuildExtraQuarticCandidates[model]` in the original
   scalar basis; insists on 7 extra operators and exactly 2 with counts
   `{1,1,0,1,0,1}`.
4. Captures the actual Matchete action expressions for the two couplings:
   `lambdaT3ExtraQ110101I1`, `lambdaT3ExtraQ110101I2`.
5. Calls live `CheckLagrangian` on `LSM + LFree +` each selected interaction,
   and the two interactions together.
6. Writes a diagnostic JSON, **including blocked reports**, and refuses to
   certify cross-tool equality or physical matching.

The Matchete source at
`Lagrangian/interactions/GeneralScalarOperators.wl` defines each extra
interaction with `contribution = coupling[] product CG[...]` and
`"Expression" -> PlusHc[contribution]`. The companion Python verifier
checks these source lines remain unchanged; it does not infer runtime
Hermitian-conjugation behavior solely from a source string.

**Important:** This diagnostic does **not** insert the symmetry-allowed
holomorphic B/C scalar mass bilinear and does **not** check the off-diagonal
heavy scalar mass basis. Its live Matchete `CheckLagrangian` is a quartic
vertex check only. It is not `HeavyMassBasis`, `Match`, a numerical EFT
threshold correction, nor the final cross-tool polynomial identity.

## Install and execute from PowerShell

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "C:\path\to\T3_BC_live_matchete_110101_action_v1.zip" -DestinationPath . -Force
python -X utf8 -m pytest -q tests/test_t3_bc_matchete_110101_action.py
wolframscript -file scripts/ProbeT3BCMatchete110101Action.wl output/audits/t3_bc_matchete_action/matchete_110101_v1.json
if ($LASTEXITCODE -ne 0) {
    if (Test-Path output/audits/t3_bc_matchete_action/matchete_110101_v1.json) {
        Get-Content output/audits/t3_bc_matchete_action/matchete_110101_v1.json
    }
    throw "Matchete 110101 action probe failed; inspect JSON"
}
python -u -X utf8 -m scripts.VerifyT3BCMatchete110101Action `
    --input output/audits/t3_bc_matchete_action/matchete_110101_v1.json `
    --output output/audits/t3_bc_matchete_action/matchete_110101_v1_verified.json
```

Choose a fresh suffix if an output file already exists; the scripts refuse
to overwrite prior evidence. No GitHub or production code is changed.

## After the live probe

Compare **both** Matchete evaluated interaction tensors and Hermitian
conjugates against the verified native RGBeta coefficients in the same
complex component convention, testing directly that

\[
  \mathcal L_{M,110101} = -V_{R,110101},\qquad
  \begin{pmatrix}h_1\\h_2\end{pmatrix}
  =-\begin{pmatrix}1/2&1\\-\sqrt3/2&0\end{pmatrix}
  \begin{pmatrix}g_1\\g_2\end{pmatrix}.
\]

The old candidate identity is **not** being promoted here to a certified
physical one-loop coupling map. Even with a cross-tool scalar potential
identity, Yukawa CG conventions, mixed scalar mass-basis action validation,
threshold matching and the rest of the EFT calculation remain outstanding.

## Verification limits

- Local Python tests are synthetic schema/source-gate tests only, **not**
  Matchete runtime tests.
- The standalone Wolfram source has not been executed in the assistant's
  environment. The user must run the above commands locally.
- If the Matchete probe fails, the diagnostic JSON and console output should
  be inspected before proposing another code change; the failure must not
  be inferred to be a physics mismatch.

## v1 Windows FrontEnd failure and v2 correction (2026-10-31)

The user ran v1 on the installed Matchete. Matchete itself loaded, but both
B(-1) and C(-1) emitted `FrontEndObject::notavail`, then the script aborted
with `Could not construct Matchete action diagnostic`. No valid JSON or
Matchete action certificate was produced. The Python verifier was not reached.

Verified script defect: v1 used `UsingFrontEnd[Needs["Matchete`"]]` for
package loading but performed `ResetAll`, `LoadModel`, representation/CG
registration and `CheckLagrangian` in `T3BCMOne` **after** the front-end scope
had ended. The project’s existing
`scripts/ExportGeneralT3ScalarTensorMatchete.wl` encloses the full native
computation inside `UsingFrontEnd[ ... ]`. v2 follows that existing pattern,
without changing the two quartic terms, any coupling normalization or any
physical matching flag.

Use a new output file; do not reuse the failed v1 path:

```powershell
python -X utf8 -m pytest -q tests/test_t3_bc_matchete_110101_action.py
wolframscript -file scripts/ProbeT3BCMatchete110101Action.wl output/audits/t3_bc_matchete_action/matchete_110101_v2.json
if ($LASTEXITCODE -ne 0) { throw 'Matchete front-end action probe failed' }
python -u -X utf8 -m scripts.VerifyT3BCMatchete110101Action `
  --input output/audits/t3_bc_matchete_action/matchete_110101_v2.json `
  --output output/audits/t3_bc_matchete_action/matchete_110101_v2_verified.json
```

The corrected v2 Wolfram runtime has **not yet been executed** in this
environment. A passing Python unit test checks source structure and evidence
schema only, not Matchete front-end functionality. If v2 still emits
`FrontEndObject::notavail`, inspect the running Mathematica installation's
front-end availability; do not infer a quartic inconsistency from it.
