# Stage 20 — full ordinary-T3 scalar interaction programme, first generalisation batch

## Objective

Extend the scalar basis across **all 16 ordinary T3 assignments**, rather than
continuing B/C(-1)-only numerical neutrino-mass studies. Everything here is
**experimental** and is kept separate from the production baseline.

## Source-verified inventory

Input: `docs/general_scalar_trial_operator_manifest.json` and the Matchete
`Lagrangian/interactions/GeneralScalarVerifiedQuartics.wl` table.

- 16 independent model assignments.
- 50 additional **complex** independent quartic coefficients summed *across*
  the assignments, NOT 50 coefficients within each individual model.
- 12 of these 50 have existing RGBeta **structural** analogues (six for B(-1),
  six for C(-1)); **38 new RGBeta candidate vertices** are generated.
- Eight additional complex *quadratic* coefficients exist, one each in
  A(-2), A(0), B(-1), C(-1), D(-2), D(0), E(-2), E(0).
- Additional Z2-even, gauge-invariant renormalisable cubic interactions do
  not arise for these assignments (odd BSM fields forbidden, and one Higgs
  plus two BSM scalars has half-integer hypercharge).

## Implemented in this batch

1. `scripts/GenerateGeneralT3RGBetaBasis.py`: deterministic field-monomial
   generator using RGBeta's native `eps`, `del`, `delS2` and triplet coupling
   structures, with Bose-symmetrised rank selection and a fail-closed check
   against all 50 independently enumerated structures.
2. `RGE/general/generated/*_extra_quartics.wl`: 16 model-specific source files,
   covering 38 previously unrepresented quartics. B/C(-1) reuse six existing
   RGBeta couplings **without creating duplicates**; the second independent
   H H† S1† S2† structure gets a new vertex.
3. `RGE/running/rgbeta/RunGeneralT3RGBeta.wl`: **experimental** UV and
   intermediate scalar-EFT beta exporter. Loads the unchanged production
   `T3RGBetaModel.wl`, adds generated vertices, requests the new `BetaTerm`s,
   and includes them in `report_betas`. `-CheckGauge` runs an independent, explicit U(1)Y/Z2/field-type/SU(2)
   candidate-rank audit. RGBeta v1.2.0 rejects `CheckInvariance -> True`;
   therefore its native invariant Ward check is **not** invoked.
4. `Numerical/experimental/GeneralT3ScalarRGState.py`: a general
   per-representation, manifest-checked state extension and generic
   UV/intermediate beta evaluator/ODE runner. It retains all generated
   extra-quartic components across a tree-level F threshold; no silent ignored
   beta keys and no claim of complete matching.
5. `Numerical/experimental/GeneralT3ScalarQuadraticMasses.py`: **standalone,
   spherical SU(2) basis prototype** of mass-squared Hessians for all eight
   additional bilinears. It verifies a zero-coupling diagonal limit and tests
   the B/C(-1) two-multiplet analytic spectrum. Not incorporated into the
   physical Matchete mass basis, RGBeta mass running, or scalar thresholds.
6. New tests, scripts and machine-readable coverage reports.

## Important verification boundaries

The **38 generated extra quartic vertices are *candidate* RGBeta bases**.
The Python rank tests confirm their span in the generator's index convention;
they DO NOT demonstrate equivalent invariant normalisation, triplet basis
rotation, phase, or the individual coefficient map against Matchete.
`-ExportBetas -CheckGauge` performs a separate structural charge/parity/
index-basis audit and exports the candidate betas under the locally installed
Wolfram kernel. **It does not run RGBeta's native `CheckInvariance` Ward check**:
that option is rejected by the installed RGBeta v1.2.0. In addition, the
B/C(-1) new delta/epsilon contraction was corrected to respect dual SU(2)
index orientations; the pre-existing B/C mixed contraction remains subject to
a separate native tensor validation. No Wolfram kernel is present
in the package-construction environment, so this local run is **mandatory**
before claiming RGBeta success. The independent structural audit still does not establish a full
native SU(2) Ward identity or Matchete/RGBeta coupling equality.

No existing production matching, RGBeta source, UV numerical states or their
existing test semantics have been changed. The eight quadratic terms are **not
physically installed**; the most general T3 neutrino-mass pipeline is NOT yet
complete. The new generic RG numerical code is intentionally an experimental
candidate until exports and coupling maps are verified.

## Commands (PowerShell, project root)

```powershell
Expand-Archive -LiteralPath .\T3_General_Scalar_Batch20.zip -DestinationPath . -Force

python -m pytest -q tests\python\test_general_t3_stage20.py tests\python\test_general_t3_stage20_quadratics.py

powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralT3Stage20.ps1
python -m scripts.RunGeneralT3QuadraticAudit
```

To export candidate beta functions, start with **two high-value controls**
and request the independent structural charge/parity/rank check:

```powershell
powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralT3Stage20.ps1 -ExportBetas -CheckGauge -Models B_-1,C_-1
```

After that succeeds, export/verify both UV and scalar-EFT beta functions for
**all 16 ordinary model assignments**:

```powershell
powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralT3Stage20.ps1 -ExportBetas -CheckGauge
```

Partial runs may select `-Models A_-2,E_-2` and/or `-Stage UV` / `-Stage EFT1`.
Model-level JSON and logs are written to `output/general_t3_stage20/`.

## Next required work (not claimed as complete)

1. For every model, independently validate native RGBeta SU(2) Ward identities, beta closure and its
   exact Matchete tensor change-of-basis/normalisation (especially complex
   conjugations for triplets).
2. Incorporate the eight quadratic mass mixings in **Matchete**, RGBeta,
   scalar mass diagonalisation and threshold matching (no diagonal-scalar
   approximation where mass mixing is nonzero).
3. Use the generic extended RGBeta numerical state to run the full coupling
   sets and revalidate zero-extra-interaction baseline recovery for all 16.
4. Re-do multi-threshold physical neutrino masses after the general scalar
   mass-basis and extra EFT1 Wilson sectors have been completed.
