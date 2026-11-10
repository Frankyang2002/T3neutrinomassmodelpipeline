# T3 Batch 10 — A(−4)/D(−2) native hard source and general HHSS completion

## Scope

This batch continues the **successful Batch 9** computations. It is a new, local-only source-backed numerical stage. It **does not** redo Batch 9, does not modify GitHub, and does not overwrite earlier outputs.

The EFT decoupling sequence stays fixed:

\[
UV\;(10^7\,\mathrm{GeV})\longrightarrow F\text{-only threshold}\;(10^5\,\mathrm{GeV})\longrightarrow EFT1(S_1,S_2)\longrightarrow (S_1,S_2)\text{ jointly out}\;(10^3\,\mathrm{GeV})\longrightarrow SMEFT\;(100\,\mathrm{GeV}).
\]

**New calculations:**

1. Read the current thesis Batch 2 **full Matchete operator ledger**, not a reduced analytic-triangle substitute. Source-classify the complete 6-term `A_-4` / 12-term `D_-2` LLHH hard sector, sum the actual three-flavour UV endpoint, and export its ordered holomorphic **regular and pole** matrices independently. The code reuses the existing repository `hard_llhh_source.py`, which blocks a nonzero `y2Conj` in D if its native extra-quartic normalization has not been established. In the saved thesis benchmark, `y2Conj` was zero at the F endpoint; changed inputs are rechecked.
2. **D(−2) scalar-potential extension omitted by the archived Batch 9 HHSS seed.** The general scalar registry includes five extra quartics, **two directly of HHSS type**: `lambdaT3GeneralQ020101I1` and `lambdaT3GeneralQ110002I1`. Read their existing exact algebraic SU(2) reference tensors, evaluate the *nonzero* running couplings at `muS`, rotate to the named real `T3Scalar2Minus/T3Scalar2Plus` fields, and compute both the scalar-pair leading logarithm and the simultaneous full-mass scalar B0 loop. The other three extra quartics have no direct HHSS legs in the present tree-LLSS scalar-pair vertex; they are not reset in the existing renormalizable RG history.
3. Assemble **only a formal, convention-dependent** `hard regular + scalar log + scalar B0` source-sum, without adding Batch 9 scalar contributions twice. Compare the raw hard pole with the explicit candidate scalar-bubble pole **as a discrepancy**, never artificially cancelling it.
4. For **A(−4)**, call the repository's full 45-Weyl/12-real-scalar candidate tensor RGE at the F endpoint, with all seven source-named beta channels and the `fonseca_eq22_six_term` variant. Compare its scalar-pair projection against the existing Matchete-basis logarithm. Since the scalar CG and spinor dictionary is still unverified, **the difference is diagnostic only**. If the candidate beta contracts fail, the failure is recorded without deleting the independently sourced hard or scalar results.
5. Optional `--illustrative-smeft`: pass the formal source-sum through the existing SM+Weinberg evolution, Takagi masses, PMNS and five-Gaussian-residual numerical machinery with an **explicitly assumed, unverified signed factor of +1**. These masses and chi-square are **illustrative and not physical predictions**. The normal run does not make an unsupported mass claim.

The D(−2) extra-tensor **reference SU(2) basis is not a verified numeric Matchete–RGBeta coupling-coordinate identification**. The two additional HHSS pieces are therefore conditional reference-basis sensitivities, not source-verified physical amplitudes. The raw Matchete F hard term is **unsubtracted**; no MSbar finite boundary has been established. A good numerical Ward residual, a small pole mismatch, or a low illustrative chi-square cannot turn this result into a physical C5. The model output records these limits explicitly.

## Files in ZIP and intended project paths

| Project-relative replacement/new file | Role |
|---|---|
| `Numerical/parallel_t3/complete/two_model_hard_general_batch10.py` | Source-based hard matrices, D extra HHSS tensors, threshold B0/log, seven-channel A diagnostic and optional low-energy calculation |
| `scripts/RunT3HardGeneralBatch10.py` | One-command two-model runner and preflight; never overwrites earlier outputs |
| `tests/test_t3_hard_general_batch10.py` | Synthetic exact-tensor, flavour symmetry, missing-channel and counterterm-isolation tests |
| `README_T3_HARD_GENERAL_HHSS_BATCH10.md` | Scope, install and physical qualifications |

Batch 9 source modules and its existing result files are **not replaced**.

## Windows PowerShell

Run from the local repository:

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline

Expand-Archive `
  -LiteralPath "$HOME\Downloads\t3_hard_general_hhss_batch10.zip" `
  -DestinationPath . -Force

.\.venv\Scripts\python.exe -m pytest -q `
  tests\test_t3_hard_general_batch10.py

.\.venv\Scripts\python.exe -X utf8 `
  -m scripts.RunT3HardGeneralBatch10 --all --preflight

.\.venv\Scripts\python.exe -X utf8 `
  -m scripts.RunT3HardGeneralBatch10 --all --debug
```

The default command runs the **A full seven-channel candidate contraction** (which can be the slowest substep). To skip that unverified but useful comparison and still calculate both actual hard/scalar-source stages, use a new output path:

```powershell
.\.venv\Scripts\python.exe -X utf8 `
  -m scripts.RunT3HardGeneralBatch10 --all --skip-a4-seven-term `
  --output output/parallel_t3/thesis_hard_general_hhss_batch10_quick
```

For explicitly *illustrative* neutrino masses (not physical MSbar matching), use **another new** run directory:

```powershell
.\.venv\Scripts\python.exe -X utf8 `
  -m scripts.RunT3HardGeneralBatch10 --all --illustrative-smeft `
  --output output/parallel_t3/thesis_hard_general_hhss_batch10_illustrative
```

## Outputs

Default output directory:

`output/parallel_t3/thesis_hard_general_hhss_batch10_v1/`

- `A_-4/result.json`: formal Matchete hard regular/pole, Batch 9 baseline response, optional A full-master seven-channel candidate diagnostic, explicit missing physical steps.
- `D_-2/result.json`: same hard and baseline matrices, two individually resolved added general-HHSS leading logs and simultaneous finite bubbles in the split-real mass basis; non-HHSS quartics retained as inventory.
- `D_-2/D_m2_extra_HHSS_reference_channels.npz`: two numerical 4×4×8×8 tensors in the named split-real basis (conditional coefficient assignment).
- `summary.json`: statuses and matrix Frobenius norms.

A successful output reports `HardAndGeneralScalarSourceComponentsComputed_ConditionalNotPhysicalC5` for each model. Output existence does not assert that a fully renormalized physical C5 has been computed.

### Verification performed here

- Current GitHub source and saved model-specific Batch1/2 artifacts were inspected **read-only**.
- The code includes three self-contained unit tests for exact SU(2) component symmetries, both mandatory D HHSS channels, and non-cancellation of formal poles; these pass in the build environment.
- The complete physical pipeline, Windows Matchete runtime, and your saved Batch9 numerical outputs are **not available here**. Actual two-model integration must be checked by the PowerShell command above; this is not a claim of local physics execution.

## Remaining implementation (not another audit project)

The physical next stages are **(i)** derive and implement the finite UV-hard-minus-EFT1-soft matching/CT with a common native-to-Weyl and Matchete-to-RGBeta coefficient dictionary, including D's extra Yukawa when nonzero; and **(ii)** complete the joint scalar one-loop operator inventory and connect the resulting genuine source-driven C5 evaluator to the existing full 18-real-UV-Yukawa fit. Roughly two more **large** batches are plausible if the native subtraction is tractable; a third may be needed if the source expressions require additional dedicated Matchete extraction.


## Windows/JSON hotfix for Batch 10 (10 November 2026)

The original Batch 10 runner blocked for **both** models with
`Joint-scalar fermion bridge: kernel complex matrix must have real/imag fields`.
The F-hard source (from `general_batch1_api.encode_native`) emits a **3x3
list of per-entry `{real, imag}` objects**. Batch 9's scalar response emits a
**single `{real: 3x3, imag: 3x3}` object**. The original helper only accepted
the second form. The replacement `two_model_hard_general_batch10.py` accepts
both, strictly checks shape and finiteness, and does not change any physics.

The original unsuccessful run has already created
`output/parallel_t3/thesis_hard_general_hhss_batch10_v1`; it is deliberately
not overwritten. Run again with a new destination:

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$HOME\Downloads\t3_batch10_json_schema_fix.zip" -DestinationPath . -Force
.\.venv\Scripts\python.exe -m pytest -q tests\test_t3_hard_general_batch10.py
.\.venv\Scripts\python.exe -X utf8 -m scripts.RunT3HardGeneralBatch10 --all --output output/parallel_t3/thesis_hard_general_hhss_batch10_fix_v2 --debug
```

The tests have five cases, including regression coverage for both actual
complex-matrix JSON shapes and malformed entries. The standalone tests passed
locally; the full integration run still requires your Windows checkout and data.
All physical certification disclaimers elsewhere in this README remain unchanged.
