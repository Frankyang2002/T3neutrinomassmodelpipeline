# T3 extended C5 — coordinated large-batch low-energy study

## Scope and physical status

This is a **local-only**, opt-in study built on the user's already working 16-model RGBeta / 6-model extended-C5 code. It **does not replace production files**, edit GitHub, or claim a certified physical prediction.

The six applicable models are `A_p0, A_m2, D_p0, D_m2, E_p0, E_m2`. The runner:

1. Re-verifies the six native Matchete one-loop ordered-coefficient reductions and the currently installed RGBeta extra quartic source field words, saving SHA-256 source hashes and a *blocked* physical-convention ledger.
2. Re-integrates the exact full coupled UV RGBeta state, **including** independent `y1Conj` / `y2Conj` matrices and every general scalar quartic, for both seeded and zero-initial-extra-Yukawa cases.
3. Repeats these UV evolutions for three common **above-all-heavy-mass** endpoint scales, specified as multipliers `[1.5, 1.8, 2.0]` of the existing fermion threshold/reference scale. It refuses to continue across an evolving heavy mass.
4. Evaluates the source-checked one-step, symmetric, full-flavour Matchete `C5` (baseline plus explicit extra Yukawa term) at those endpoints using the **already declared, unverified** RGBeta-to-Matchete conversion scenario. It retains the four independent relative extra-quartic phase hypotheses from the previous package.
5. As a **sudden-decoupling approximation**, constructs a SM+Weinberg EFT at that single common endpoint, and uses the existing project's 76-real-variable `Numerical.running.SMWeinbergEvolution.evolve_sm_weinberg` to integrate to `mu_low_gev = 91.1876`.
6. Computes low-energy neutrino mass eigenvalues, squared mass differences and `|PMNS|`, plus seeded-vs-zero-extra response and direct-extra contribution using the actual existing Takagi routines; checks the numerical **linearity of the Weinberg RGE in `C5`** for a fixed SM running background.

**Important physics limitation:** The one-step decoupling is performed at a scale *above* every heavy mass proxy; all heavy fields are removed simultaneously there, and only the SM runs below. The program does **not** calculate the intermediate scalar-only EFT, matching counterterms, mass-mixing effects, sequential fermion and scalar matching, or a physical absolute Matchete/RGBeta component-basis coupling dictionary. Its matching-scale dependence is a **diagnostic of this approximation**, not a certified error bar. The study outputs must not be presented as validated physical neutrino-mass predictions or physical fits.

A direct extra Yukawa contribution from the native reduced Matchete expressions is already source-verified for all six models; the new step here is the coupled numerical **SM+Weinberg continuation**, not a new certified matching calculation.

## Install (Windows PowerShell)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline

$zip = "$HOME\Downloads\T3ExtendedLowEnergyBatch_Local.zip"
Expand-Archive -LiteralPath $zip -DestinationPath . -Force

python scripts\RunT3ExtendedLowEnergyBatch.py --self-test
python scripts\TestT3ExtendedLowEnergyBatch.py
python scripts\RunT3ExtendedLowEnergyBatch.py --preflight

# First test A(0), one endpoint and four phases
python scripts\RunT3ExtendedLowEnergyBatch.py --pilot
$pilotExit = $LASTEXITCODE

# Then all six models, three endpoint scales and four phases
if ($pilotExit -eq 0) {
    python scripts\RunT3ExtendedLowEnergyBatch.py --all
}
$LASTEXITCODE
```

All commands require the preceding local projects/stacks already installed (you have run them), particularly:

- `Numerical/experimental/T3GeneralYukawaCoupledUV.py`
- `Numerical/experimental/T3RGAssistedExtendedC5.py`
- `Numerical/experimental/T3ExtendedC5ExpressionReduction.py`
- `Numerical/configs/t3_rg_assisted_c5_convention_hypothesis.json`
- Prior 6 native one-step Matchete C5 JSON and InputForm files, or `--native-source PATH_TO_ExtendedC5Expressions.zip`.
- Real project `Numerical/running/SMWeinbergEvolution.py`, `RGE/running/weinberg/WeinbergRGE.py`, SciPy, NumPy and SymPy.

The options `--scenario`, `--study`, `--native-source`, and `--output-dir` are available. A valid scenario is still **explicitly hypothetical**. Do not flip either certification flag to `true` to remove the gate; it is intentionally rejected.

## New files / intended installation paths

- `Numerical/experimental/T3ExtendedLowEnergyStudy.py` — full conditional UV / one-step / SM+Weinberg study.
- `Numerical/experimental/T3ExtendedConventionLedger.py` — source-checked names/monomials, native proof manifest and physical normalization blockers.
- `Numerical/configs/t3_extended_low_energy_study.json` — opt-in scale-study settings.
- `scripts/RunT3ExtendedLowEnergyBatch.py` — batched CLI, logging/JSON/CSV.
- `scripts/TestT3ExtendedLowEnergyBatch.py` — negative controls, synthetic integrated tests and production SM+Weinberg smoke when available.
- `README_T3ExtendedLowEnergyBatch.md` — this document.

**No existing file is replaced.** The previous 16-model UV and native C5 reports remain unchanged.

## Outputs

`Reports/output/workstream_a/extended_c5_low_energy/`

- `physical_dictionary_gate.json`: a **blocked** sign/phase/basis ledger with SHA-256 evidence.
- `models/A_p0.json` ... `models/E_m2.json`: full phase- and scale-resolved initial and evolved complex `C5` matrices and low-energy observables.
- `low_energy_scan.csv`: single comparison table of the six models × three endpoint scales × four phases (72 rows if all pass).
- `summary.json`: status, settings, source, model-by-model data and physical certification flags.

The CSV can be used for **rough conditional numerical comparisons**, not authoritative physical fits. These are not guaranteed to reproduce an experimental oscillation dataset without fitting.

## Status distinctions

**Verified previously on the user's Windows system:** 16/16 real coupled UV RGBeta reports, 6/6 native Matchete extended one-step `C5` extractions, 6/6 exact symbolic C5 reductions, 6/6 conditional one-step UV-to-C5 evaluation.

**Tested in this package's isolated environment:** Python compilation, synthetic SciPy UV integrations (two scale points), all six C5 model paths, source-word tampering, invalid threshold rejection, rejected false certification, broken C5 linearity negative control, report-code imports.

**Requires the user's live local execution:** actual 76-variable SM+Weinberg continuation for the six real UV running endpoints. The test file includes a production SM+Weinberg smoke test for that environment.

**Still pending physical theory:** absolute CG/Weyl/Dirac convention map and full sequential heavy-fermion/scalar matching with extra operators and threshold corrections. These cannot be inferred merely by connecting two working numerical packages.
