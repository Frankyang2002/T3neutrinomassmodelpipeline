# A(-4): exploratory coupled seven-term dimension-five EFT1 evolution

This package supplies **one complete replacement Python file** without changing the existing UV/EFT1 numerical pipeline, previous coupled quartic runner, saved fits, or saved results.

## Files and project paths

- `scripts/RunT3CoupledEFT1SevenTerm.py`: complete new seven-term experimental runner.
- `Reports/notes/SevenTermEFT1Experiment.md`: this documentation.

**Dependency:** the previous experiment's `scripts/RunT3CoupledEFT1Quartic.py` must already exist in the local project. It is imported solely to compute the scalar-quartic-only comparison. The user has already run this previous script successfully.

## Implemented comparison

1. Fixed LLSS leading log: existing `Numerical/running/GeneralT3SequentialLeadingLog.response`.
2. Coupled LLSS↔LLHH **quartic-only** evolution: existing experimental `RunT3CoupledEFT1Quartic.integrate_coupled_quartic`.
3. Coupled **all-seven-term** master \(\psi^2\phi^2\) anomalous-dimension evolution using **the actual repository function** `RGE/running/intermediate/GeneralEFT1C5Mixing.calculate_master_beta`, integrated with `solve_ivp(method='DOP853')`.

For (3), all 12 real scalar components are included: four Higgs components, two components of \(S_1\), and six of \(S_2\) in A(-4). The Wilson tensor contains all initially-zero and generated components in the **nine-Weyl-field lepton sector**, not just LLSS and LLHH. These nine fields are the six components of three lepton doublets and the three charged-conjugate right-handed leptons. The SM quark Yukawas contribute to the Higgs scalar anomalous dimension through their full three-generation, three-colour trace. With only SM Higgs Yukawas in this fermion-decoupled EFT, they do not source a quark-external Wilson tensor from an initial LLSS lepton tensor. This reduction was cross-checked against the unreduced 45-field beta function in a separate local numerical test at the initial boundary (all seven sectors agreed and quark-external beta entries remained zero).

The anomalous dimension uses:
- Yukawa wavefunction term;
- scalar quartic and gauge-exchange term;
- mixed Yukawa/gauge term;
- crossed Yukawa term;
- conjugate-Wilson-coefficient Yukawa term;
- scalar anomalous dimension term;
- fermion anomalous dimension term.

**The external renormalisable coupling values are frozen to the EFT1 \(\mu_S\) endpoint** during the Wilson ODE. The existing *renormalisable* EFT1 trajectory was still computed by the original pipeline. This is an intermediate full-seven-term **Wilson-ODE** experiment, not yet a simultaneously evolved \((g,Y,\lambda,C)\) integration.

A(-4) has no additional quartic candidates in `RGE/general/generated/general_rgbeta_quartic_registry.json`; the HHSS tensor is taken *unchanged* from the authoritative LLSS scalar-tensor calculation. The remaining HHHH and SSSS blocks are assembled using the following **provisional baseline potential convention**:

\[
V_4\supset \frac{\lambda_H}{2}(H^\dagger H)^2 + \frac{\lambda_{S_1}}{2}(S_1^\dagger S_1)^2 + \frac{\lambda_{S_2}}{2}(S_2^\dagger S_2)^2 + \lambda_{12}(S_1^\dagger S_1)(S_2^\dagger S_2) +\frac{\lambda_{S_2\mathrm{Adj}}}{2}\sum_a (S_2^\dagger T^a S_2)^2.
\]

This source-to-master convention, the left-Weyl Yukawa embedding and Matchete GammaCC normalisation are **not independently certified**. We are deliberately testing the equation numerically, not claiming a verified physical prediction.

The final scalar matching uses the *same finite scalar bubble as the fixed-order run* so that differences isolate the newly integrated Wilson RGE. No independent hard fermion threshold, extra dimension-five operator classes or new scalar-loop finite diagrams are added.

## Fix for the first A(-4) local run (2026-11-11)

The previous `seven_A_m4_v1` run stopped with:

```text
A_-4 Blocked TypeError: must be real number, not dict
```

Cause: `Numerical.orchestration.GeneralRunningBenchmark._state_values` encodes
full SM Yukawa matrices as nested 3x3 lists of `{"real":...,"imag":...}`
complex entries. The first experimental runner called the LLSS `decode()`
function, which expects a *whole-matrix* `real`/`imag` dictionary; NumPy then
attempted to convert a per-entry dictionary into a number.

Fix: `_finite_complex_matrix` now calls the repository's established
`scripts.RunT3GeneralSequentialLeadingLog._complex_matrix` decoder, already
used for SM+Weinberg evolution. This supports both encoded formats. A
regression test has been added to `--self-test`.

No physics equations, Wilson integrator, scalar tensor definitions, or ablation
logic were changed in this repair. The schema version was advanced to v2.

## Install the complete replacement (PowerShell)

Project paths included in this ZIP:

- `scripts/RunT3CoupledEFT1SevenTerm.py` — **replace existing file**
- `Reports/notes/SevenTermEFT1Experiment.md` — updated instructions

Place the ZIP in `Downloads`, then:

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive `
  -Path "$env:USERPROFILE\Downloads\T3_Seven_Term_EFT1_ParserFix.zip" `
  -DestinationPath . -Force
```

## Rerun using a NEW run ID

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline

python -X utf8 -m scripts.RunT3CoupledEFT1SevenTerm --self-test

python -X utf8 -m scripts.RunT3CoupledEFT1SevenTerm `
  --source-run thesis16_merged_v1 `
  --run-id seven_A_m4_v2 `
  --models A_-4 `
  --check-only

python -X utf8 -m scripts.RunT3CoupledEFT1SevenTerm `
  --source-run thesis16_merged_v1 `
  --run-id seven_A_m4_v2 `
  --models A_-4
```

New results (if completed):

`output/studies/general_scalar_extended/experimental_coupled_eft1_seven_term/seven_A_m4_v2/`

Inspect `summary.json`, `A_m4/cases.json` and, if the model is blocked,
`A_m4/error.json` for the full Python traceback. The previous failed
`seven_A_m4_v1` results are left intact.

## Verification of this repair

- Full Python replacement compiles successfully.
- Both source-compatible Yukawa encoding formats decode correctly using the
  verified existing decoder logic.
- Full A(-4) physics calculation **not rerun here**; it needs the local fitted
  project data.
- The original seven-term analytic test was reported as passed in the user's
  local run before the input-decoder error.

The underlying physical limitations of the experimental seven-term solver
remain unchanged.
