# General T3: status and first direct one-loop Weinberg probe

## What is already verified

- All 16 ordinary models: general quartic candidate inventories and **experimental** quadratic+quartic EFT1 numerical runs have completed for the user's benchmarks. These are scalar-sector running results, not final neutrino-mass results.
- 14 ordinary models had diagonal Matchete heavy-basis checks; `B_-1` and `C_-1` now also pass native `HeavyMassBasis` with two mixed scalar doublet eigenfields. The user's local live runs reported the two pairs of scalar masses near (1035.58,1188.53) and (1032.61,1186.00) GeV.
- The B/C eigenbasis has 2 Yukawa source couplings and 17 scalar quartic source slots. The Matchete SU(2) tensor recoupling and numeric scalar-mass diagonalization were separately checked. Physical RGBeta-to-Matchete coupling normalization is not certified.

## What is still missing for general physical results

1. Matching to a physical Weinberg operator after the heavy states are removed, preserving the coefficient convention and flavour structure.
2. Correct numerical identification of Matchete couplings with the evolved RGBeta/EFT1 couplings and the benchmark Yukawa matrices, including relative CG signs and finite terms.
3. Consistent fermion-first and split-scalar matching at the appropriate scales, rather than comparing quantities at mismatched threshold scales.
4. Run the final SM+Weinberg RGEs and obtain symmetric 3x3 C5 matrices, neutrino masses and residual checks. Extend/check the matching and quadratic mass treatment across all 16 models; B/C do not establish this for all models.

## This experimental batch

Two **new** project files introduce a symbolic exploratory diagnostic while leaving `RunGeneralT3MassEigenMatchete.wl`, the validated heavy-basis checker, unchanged.

- `Lagrangian/RunGeneralT3MassEigenOneLoop.wl`: forks the verified eigenfield construction in a fresh Wolfram process, rechecks `HeavyMassBasis`, calls the existing `Lagrangian/RunMatching.wl` entry point `RunT3Matching[lag,5,1]`, exports the matched EFT in raw `InputForm`, and invokes `ExtractWeinbergCoefficient[matchedEFT]`. If a holomorphic dimension-five Weinberg term is found, saves the **ordered**, raw Matchete coefficient in `InputForm`. The built-in optional display canonicalizer assumes `MS1/MS2`, so it is not certified for `MSX/MSY` and is deliberately not used for a numerical C5.
- `scripts/RunGeneralT3MassEigenOneLoop.py`: requires the previous successful local native heavy-basis result, rechecks the numerical mixing input/provenance, executes the symbolic probe, and validates the status and exported artifacts. Failures report the `*.progress.json` stage and full `*.log`.
- `tests/python/test_general_t3_mass_eigen_one_loop.py`: offline interface and fail-closed checks. They do not constitute a live Matchete one-loop matching test.

The probe integrates out all heavy species **in one step at a common symbolic matching point**; it does **not** follow the project's production UV -> F -> EFT1 -> S -> Weinberg threshold sequence, does not choose a numerical matching scale and does not evaluate a numerical 3x3 C5 matrix. It therefore **does not establish physical neutrino masses or finish the general 16-model calculation.**

### Install (PowerShell, repository root)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath 'C:\path\to\T3_General_One_Loop_Matching_Probe.zip' -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_mass_eigen_one_loop.py tests\python\test_general_t3_mass_eigen_matchete.py tests\python\test_general_t3_mixed_scalar_mass_basis.py
python -m scripts.RunGeneralT3MassEigenOneLoop --model B_-1
python -m scripts.RunGeneralT3MassEigenOneLoop --model C_-1
```

### Output files

`output/general_t3_mass_eigen_matchete/one_loop_probe/` contains model-specific `*_one_loop_probe.json`, `*_one_loop_probe.progress.json`, `*_one_loop_probe.log`, and—on successful matching—`*_matched_eft.inputform.txt` and `*_weinberg_ordered_coefficient.inputform.txt`.

Next task, once the symbolic coefficient is confirmed by a **live** run: compare with the project matching convention, independently verify the loop and flavour normalization and supply benchmark inputs to evaluate C5 and `m_nu=-v^2 C5/2` at a fixed scale. Then connect to the fermion-first EFT1 evolution for final general results.
