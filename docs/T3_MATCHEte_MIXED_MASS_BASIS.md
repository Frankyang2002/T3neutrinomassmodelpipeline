# T3 Matchete mixed-heavy-mass audit fix

## Observed failure (B_-1)

The original experiment generated an independent gauge-invariant bilinear
`S1dag S2dag + h.c.` and added it to a *heavy-field* Lagrangian passed to
`CheckLagrangian`. Matchete reported:

`HeavyMassBasisQ::mixedfields` / `CheckLagrangian::HeavyMassBasis`.

This is a **mass-basis obstruction**. It does not by itself establish a
violation of U(1)Y, SU(2)L, or Z2. The same obstruction applies to C_-1.

The model hypercharges are

    Y(S1) = alpha/2 = -1/2,
    Y(S2) = (alpha+2)/2 = +1/2.

Their conjugate bilinear is neutral and their SU(2) doublets admit one
antisymmetric-singlet contraction. A nonzero coefficient mixes scalar
states. To use Matchete for heavy-scalar matching, one must diagonalize the
quadratic Lagrangian and consistently rotate all Yukawa and quartic vertices.

## What this batch changes

**Complete replacement** of:

- `Lagrangian/RunGeneralQuadraticMatcheteAudit.wl`
- `scripts/RunGeneralT3QuadraticMatcheteAudit.py`

**New**:

- `tests/python/test_general_t3_matchete_mass_basis.py`
- `docs/T3_MATCHEte_MIXED_MASS_BASIS.md`

For all 16 models the audit builds and checks the original heavy-field
Lagrangian. It generates all additional candidate quartics and quadratic
bilinears in the existing Matchete generator and checks the expanded
Lagrangian.

For the 14 models other than B_-1/C_-1, it validates in the original
heavy-field basis. For B_-1/C_-1, it uses the *existing* repository helper
`ReclassifyT3Fields[model, {"F"}, {"F", "S1", "S2"}]`, which classifies S1/S2
as light but keeps F heavy. It rebuilds the free Lagrangian and checks the
entire coupled interaction in the **pre-diagonalization scalar basis**.
This deliberately bypasses only the requirement that external heavy fields
enter Matchete already in their matching mass basis. It does *not* perform
matching or provide a physical mass-basis certificate.

The Python wrapper checks all JSON flags and refuses to label the two
pre-diagonalization checks as heavy mass-basis passes. There must be 14
heavy-basis plus 2 pre-diagonalization audits and 50 quartic/8 quadratic
model-specific slots. An unsuccessful light-basis Matchete run remains an
explicit failure, not a reinterpreted success.

## Commands (PowerShell, project root)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath 'C:\path\to\T3_Matchete_Mixed_Mass_Basis_Fix.zip' -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_matchete_mass_basis.py
python -m scripts.RunGeneralT3QuadraticMatcheteAudit --model B_-1
python -m scripts.RunGeneralT3QuadraticMatcheteAudit --model C_-1
python -m scripts.RunGeneralT3QuadraticMatcheteAudit --all
```

Results: `output/general_t3_quadratic_eft1/matchete/`; combined audit
summary: `matchete_16_model_summary.json`.

## Certainty/status

- **Verified from source**: two doublet mixed-field mass structures; production
  `T3Fields.wl` supports heavy/light reclassification; the original audit
  unconditionally added a mixed mass to a heavy-field Lagrangian.
- **Verified offline**: this replacement's Python syntax and report basis
  flags, model totals, CLI arguments, ZIP structure and tests.
- **Must be tested on your machine**: the live Matchete `UsingFrontEnd` behavior
  after reclassifying S1 and S2 as light. This environment has no Wolfram
  front end and cannot assert that Matchete accepts the reclassified model.
- **Not established**: physical Matchete/RGBeta coupling mapping,
  diagonalized heavy-scalar states, mass-matrix RG conventions,
  generalized scalar matching, or neutrino-mass predictions.

The two models' pre-diagonalization operator test **is not** evidence that
Matchete supports mixed heavy states during matching. The next physics task
is scalar mass-basis transformations (including Yukawa/quartic rotation) and
mass-eigenstate matching.
