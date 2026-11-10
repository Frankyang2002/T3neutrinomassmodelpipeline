# Batch 13 — A(−4), D(−2) strict one-loop master EFT1 Wilson transport

**Purpose:** Use the successful local Batches 9 and 11 to calculate the full
`C[i,j,a,b]` (45 Weyl × 12 real scalars) **candidate** one-loop anomalous
dimension, including the seven named RGE channels, and connect its direct
`LLSS -> LLHH` projection to the existing 100 GeV neutrino pipeline.

This is a **large numerical implementation batch**, not a physical
renormalization or normalization certificate. The generic seven-term ADM
remains unverified for the full gauge-complete general complex-Weyl basis.

## What is implemented

- Rebuild source-verified **three-flavour** tree `LLSS` tensors at `mu_F=1e5`
  from the saved native first fermion matching. All fermions have already
  been integrated out; neither scalar is integrated out until `mu_S=1e3`.
- Map the six lepton Weyl components into the actual SM 45-component Weyl
  tensor basis. Keep all 12 real scalar components and both separate exchange
  symmetries. For D(−2), reuse the existing split-real scalar rotation.
- Build full SM Yukawa and SU(3) × SU(2) × U(1) generators in that basis.
- Construct the RG-closed baseline scalar quartic tensors from the existing
  RGBeta potential definitions, retaining the source-derived topology quartic
  with native Clebsch–Gordan tensors. This includes portal and scalar-self
  structures **absent** from the archived Matchete seed. All five D(−2) extra
  quartics enter: two HHSS exact SU(2) reference tensors and three explicit
  scalar-only reference polynomials. Their unverified RGBeta-to-Matchete
  normalization is recorded, not suppressed.
- Calculate **all seven** candidate anomalous-dimension terms on the *full*
  `45 × 45 × 12 × 12` boundary tensor at the fermion threshold.
- Use the first-order expansion
  `C(mu_S) = C_tree(mu_F) + log(mu_S/mu_F) beta[C_tree(mu_F)]/(16 pi^2)`.
  This is strict one-loop power counting: do **not** evolve one-loop hard
  matching again or add a one-loop correction to the one-loop scalar bubble.
- Compare the `LLHH` projection of this tensor result against Batch 11's
  *independent* Laurent-pole-derived EFT1 scale log. **Replace**, never add,
  the Batch 11 log with the candidate RGE log; save both for comparison.
- Run the existing one-loop SMEFT evolution, Takagi/PMNS solver, and NuFIT
  Gaussian objective at 100 GeV, **without rerunning the Yukawa fit**.

## What this does NOT establish

The signed Matchete GammaCC → two-component Weyl normalization is still
unresolved; the master-to-native source coefficient factor is an **explicit
conditional input**, not a derived constant. The unverified Ward/gauge
completion of the general seven-term anomalous dimension and the full
one-loop native finite operator/counterterm inventory remain open. A full
physical `C5(mu_S)` and physical NuFIT fit cannot yet be claimed. The
source-generated RGE and native pole-derived logarithms must be compared;
any discrepancy is diagnostic and must **not** be numerically fitted away.

### Install on Windows PowerShell

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive `
  -LiteralPath "$HOME\Downloads\t3_master_wilson_batch13.zip" `
  -DestinationPath . -Force

.\.venv\Scripts\python.exe -m pytest -q tests\test_t3_master_wilson_batch13.py
.\.venv\Scripts\python.exe -X utf8 -m scripts.RunT3MasterWilsonBatch13 --all --preflight
.\.venv\Scripts\python.exe -X utf8 -m scripts.RunT3MasterWilsonBatch13 --all --debug
```

Defaults are `--master-to-native-factor 1 --signed-factor 1` (explicit
**illustrative** factors, NOT verified). If diagnosing only the EFT1 tensor
without SMEFT, pass `--skip-smeft`. If you want the unverified crossed-
Yukawa candidate rather than the reference six-term mode, use
`--beta-mode ward_covariant_candidate`. Neither is a certified ADM.

### Files installed (project-relative paths)

1. `Numerical/parallel_t3/complete/two_model_master_wilson_batch13.py`
2. `scripts/RunT3MasterWilsonBatch13.py`
3. `tests/test_t3_master_wilson_batch13.py`
4. `README_T3_MASTER_WILSON_BATCH13.md`

Results: `output/parallel_t3/thesis_master_wilson_batch13_v1/`, with
`A_-4/result.json`, `D_-2/result.json`, their `full_tensor_candidate.npz`
files, and `summary.json`.

**No GitHub changes.** Existing Batches 9–12 are used read-only.
