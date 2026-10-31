# T3 Workstream A — UV-to-extended-C5 conditional numerical bridge

## Purpose and exact physics scope

This stack **reuses, without replacing**, the successfully run 16-model coupled
RGBeta UV solver and the six source-verified Matchete one-step extended-C5
reductions. The bridge runs **six** affected models (`A_p0`, `A_m2`, `D_p0`,
`D_m2`, `E_p0`, `E_m2`) and creates a numerical UV endpoint that stays **above**
all relevant heavy thresholds, then evaluates the source-derived **one-step**
matched Matchete-basis C5 with an **explicitly unverified coupling dictionary**.

It does **not** perform EFT1 intermediate running, sequential fermion/scalar
matching, conversion to certified physical mass eigenstates, electroweak scalar
mass mixing or low-energy C5 evolution. Do not present these numbers as
physical neutrino masses or use them for a final fit.

The explicit hypothesis file sets `y_M = -y_R`, including the independent
conjugate-scalar Yukawa, and `lambda_M = -lambda_R`. Those overall signs are
**choices**, not an established cross-package vertex identity. The A(0)
Weyl-action analysis in the existing source supports one *conditional* such
choice, but physical mass/phase transport and the doublet/triplet CG coordinate
identification are not certified across all six models. The extra quartic is
varied by an additional phase 0, +pi/2, pi, -pi/2 so the sensitivity to its
unknown orientation is visible, not hidden.

The relationship between the Matchete source quartic name
`lambdaT3ExtraQ...` and RGBeta source name `lambdaT3GeneralQ...` is checked
against the actual six-model AddQuartic field monomials. **Field content alone
is not an equality of the SU(2) tensor conventions**. The native Matchete C5
source itself is re-verified from the 24 native files before evaluation.

## Implemented equations

The verified Matchete ordered one-step expression (after explicitly
identifying the unindexed MF in display metadata with MF[r]) is

    A5_pq = k/(16*pi^2) sum_r M_r I(M_r^2,MS1^2,MS2^2) * [
      lambdaT3_M * conj(y1_M[p,r]) * conj(y2_M[q,r])
      + conj(lambdaExtra_M) * conj(ya_M[p,r]) * conj(yb_M[q,r]) ]

where `k=-2/sqrt(3)` for A/D and `k=2*sqrt(2/3)` for E; the extra pair is
`(y1Conj,y2)` for alpha=0 or `(y1,y2Conj)` for alpha=-2.

`C5=A5+A5.T` and `mnu=-v^2*C5/2` in the previous code's conventions. The
new explicitly **hypothetical** C5 evaluator uses this formula and the existing
neutrino-observable functions. A cross-check against the existing Matchete-only
evaluator verifies equal numerical kernels for identical inputs, but the
Matchete-only evaluator is not used to bypass its strict no-transport guard.

UV mass matrices are biunitarily diagonalized, using the same **hypothetical
Yukawa heavy-slot rotation** as the current production
`FinalC5TrajectoryAdapter.py`: `y1`/`y1Conj` on the conjugate-left SVD
factor, `y2`/`y2Conj` on the right factor. Off-diagonal MF is handled rather
than being treated as diagonal. The scalar masses are sqrt of the running
**bare** mS1Sq,mS2Sq: EWSB scalar mixing and the scalar mass matrix are **not**
claimed to be included.

For each extra-quartic phase scenario, output includes the baseline, direct
extra C5, the seeded-minus-zero-extra-seed running effect, and conditional
neutrino masses/PMNS (using the running charged-lepton Yukawa). Zero-seed
control is integrated independently: radiative regeneration is not assumed to
vanish.

## Files in ZIP (new files only)

- `Numerical/experimental/T3RGAssistedExtendedC5.py`
- `Numerical/configs/t3_rg_assisted_c5_convention_hypothesis.json`
- `scripts/RunT3RGAssistedExtendedC5.py`
- `scripts/TestT3RGAssistedExtendedC5.py`
- `README_T3RGAssistedExtendedC5.md`

No existing production source, UV output, C5 matching JSON or full numerical
scan file is overwritten. All outputs go to
`Reports/output/workstream_a/rg_assisted_extended_c5/`.

## Local PowerShell

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
$zip = "$HOME\Downloads\T3RGAssistedExtendedC5_Local.zip"
Expand-Archive -LiteralPath $zip -DestinationPath . -Force

python scripts\RunT3RGAssistedExtendedC5.py --self-test
python scripts\TestT3RGAssistedExtendedC5.py
python scripts\RunT3RGAssistedExtendedC5.py --pilot
if ($LASTEXITCODE -eq 0) {
    python scripts\RunT3RGAssistedExtendedC5.py --all
}
$LASTEXITCODE
```

If you only want a model: `python scripts\RunT3RGAssistedExtendedC5.py --model E_m2`.
To use a different expressly declared map, copy/edit
`Numerical/configs/t3_rg_assisted_c5_convention_hypothesis.json` and supply
`--scenario <path>`. These edits change numerical hypotheses, **not** the
verified RGBeta or Matchete equations.

The output `summary.json` and `models/*.json` always carry false physical
certification flags. A successful process exit means **conditional UV-to-one-
step numerical evaluation succeeded**, not certified hierarchical matching.

## Tests performed before release

- Python compile and isolated six-model synthetic C5 tests using the existing
  Matchete-basis evaluator.
- Exact reduction validation against all six uploaded native Matchete reports.
- Source monomial/tamper checks and configuration-negative tests.
- Zero-beta SciPy integration exercised through the new orchestration layer.
- No actual Windows Wolfram or local 16-model RGBeta run can be executed here:
  the decisive real numerical run is the command above.

## What genuinely remains

1. Establish native absolute Matchete/RGBeta amplitude-level coupling maps
   for the charged fermion, both Yukawa chiralities and the triplet tensor
   bases, including extra-quartic CG normalizations/phases.
2. Include allowed scalar bilinears and proper mass eigenstates where they
   mix. Derive and run sequential hard threshold matching/EFT1 RGEs and
   one-loop counterterm-consistent final C5; then run SM+Weinberg RGEs.
3. Only after those two gates may physical mass and PMNS fits use transported
   RGBeta parameter values.
