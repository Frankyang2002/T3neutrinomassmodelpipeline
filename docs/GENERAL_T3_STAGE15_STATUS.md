# Stage 15 — quartic Lagrangian sign convention and source-mapped B/C(-1) response

Stage 14 successfully calculated 16 hypothetical complex-coupling mappings per
model. For its toy benchmark the induced extra-only neutrino-mass norms span
B(-1) `5.112920e-14`–`3.103345e-13` eV and C(-1)
`4.542737e-14`–`2.757276e-13` eV. These were **not physical predictions**.

Stage 15 derives the source-level coupling convention for the *two C5-active*
quartics using three independently checked facts:

1. `Lagrangian/interactions/GeneralScalarOperators.wl` inserts the new scalar
   interaction as `PlusHc[lambdaMat[] * fields * CG]` into the UV Lagrangian,
   **with a plus sign and without a 1/2 vertex factor**.
2. `RGE/running/rgbeta/T3RGBetaModel.wl` registers the identical field
   monomials and the `delS2 delS2` invariant using `AddQuartic[..., SelfConjugate -> False]`:
   `Bar[H] Bar[H] S2 S2` and `Bar[H] Bar[H] Bar[S1] Bar[S1]`.
   `output/general_t3_stage13/doublet_tensor_map.json` already certifies
   their RGBeta/Matchete group-tensor factor is **+1**.
3. **RGBeta's published model-builder Lagrangian** includes scalar quartics
   with the overall sign **minus**, Eq. (2.3) in A. E. Thomsen,
   *Introducing RGBeta*, EPJC **81** (2021) 408,
   <https://doi.org/10.1140/epjc/s10052-021-09142-4>.
   The package's `AddQuartic` interface does **not** insert optional user-level
   symmetry normalisation factors automatically, as described in Section 4.

In the **same named scalar-field basis**, equating the interaction terms gives

\[
    \lambda_{\text{Matchete},i}=-\lambda_{\text{RGBeta},i},\qquad
    i\in\{020020,020200\}.
\]

**No complex conjugation** is needed in that source-level identification;
conjugated couplings inside the exported Matchete C5 response remain
conjugated when numerically evaluated.

This is a *source-convention derivation*, not a claim of a full end-to-end
RGE/matching validation. The scripts fail if the relevant sources change, the
previous tensor report is incomplete, or the new numerical evaluation
fails to reproduce Stage 14's `negative/negative` scenario.

## Run from project root

```powershell
Expand-Archive -LiteralPath .\T3_General_Scalar_Batch15.zip -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_stage15.py
powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralT3Stage15.ps1
```

To run only the source convention check:

```powershell
powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralT3Stage15.ps1 -ConventionOnly
```

Read-only inputs: Stage 11, 13, 14 report JSONs, experimental C5 kernel JSONs,
current Matchete and RGBeta source files, and the toy parameters.

Outputs (all under `output/general_t3_stage15/`):

- `source_lagrangian_convention.json`
- `B_m1_source_mapped_extra_c5.json`
- `C_m1_source_mapped_extra_c5.json`
- `stage15_summary.json`, `stage15_summary.csv`

## Important limits

- Stage 13's running stops at 3.5 TeV above the toy fermion masses, while the
  Matchete hard responses are evaluated at `MF = 1.5, 2.1, 2.8 TeV`. The
  parameters are not evolved between these scales or matched across thresholds.
- This is an extra-only response; the original C5 and scalar/intermediate EFT
  evolution are not included.
- Quartics induced from zero by one-loop running and inserted into the one-loop
  C5 response contribute at formal order \(\hbar^2\). This is **not** a complete
  two-loop neutrino-mass calculation.
- A field-basis rotation or changed coupling normalisation would require a
  new independent vertex-level check; source signatures are pinned to prevent
  silently applying this map elsewhere.
- None of the production numerical, RGBeta or matching modules are modified.
