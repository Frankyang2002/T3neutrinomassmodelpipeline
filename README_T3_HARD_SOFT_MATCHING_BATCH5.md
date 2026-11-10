# T3 Batch 5 — hard triangle / EFT1 soft subtraction and native pole ledger

## Installation

Complete new files, relative to the local Windows project root
`C:\thesistemp\T3neutrinomassmodelpipeline`:

1. `Numerical/parallel_t3/wilson/matching_soft_subtraction_batch5.py`
2. `scripts/RunT3HardSoftSubtractionBatch5.py`
3. `scripts/RunT3HardSoftSubtractionBatch5.ps1`
4. `tests/test_t3_matching_soft_subtraction_batch5.py`
5. `README_T3_HARD_SOFT_MATCHING_BATCH5.md`

All are additions, not replacements. No GitHub modifications. The implementation
was developed against public main `4bd1230256701ac3801af4e83e8c73e5d448ee1e`.

## What is implemented (and what is not)

**Actual algebraically complete subtraction for ONE triangle topology:** Use
`Numerical/lo_certification/completeness/hard_matching.py` to evaluate every
mass-eigenstate pair in the existing one-loop heavy-fermion/two-scalar triangle:

\[
 C_{5,\mathrm{full\ UV\ triangle}}=
 C_{5,\mathrm{hard\ dim5}}(\mu)
 + C_{5,\mathrm{EFT1\ scalar\ soft}}(\mu)
 + \Delta C_{5,\mathrm{power\ suppressed}}.
\]

In that original source convention, for each heavy eigenstate \(M_k\) and two
scalar eigenstates \(m_a,m_b\), with common source tensor weight \(w_{kab}\),

\[
 K_{\rm hard}=1+\log(\mu^2/M_k^2),\qquad
 K_{\rm soft}=-\bar B_0(0;m_a^2,m_b^2;\mu),\qquad
 K_{\rm UV}=\operatorname{triangle\_ratio}(M_k^2,m_a^2,m_b^2).
\]

Here `w` includes the loop factor `1/(2*16π²)` and complete real-scalar
contractions. The library's Weyl projector is converted to the established
`C5` convention by a factor **2**. All **three fermion eigenstates and all**
scalar real eigenstate pairs are included. The \(\mu\) dependence cancels
numerically in `C5_hard_dim5+C5_EFT1_scalar_soft` with couplings frozen, for
\(\mu/M_F = 1/2,1,2\). The full UV triangle is explicitly \(\mu\)-independent.
The difference is a **one-loop dimension-7-and-higher** contribution, not a
two-loop term; it is never inserted into the dimension-5 SMEFT coefficient.

**Native Matchete pole/constant/log inventory:** The saved Batch-4 numerical
`C5_native_by_kernel_and_heavy` is divided by its verified
`1+log(MF²/Mk²)` factor to reconstruct the source's actual 88 PL Laurent
kernel residues. The exact native pole is compared numerically with the EFT1
scalar soft pole from the above triangle contraction. An **unaccounted pole
residual is reported unchanged**, because UV counterterms and the precise
meaning of the saved Matchete `MatchedEFT` output have not been established.

**A separate hypothetical sensitivity branch** treats the native pole-dropped
finite output *as though it were the unsubtracted UV amplitude*, and subtracts
an explicit renormalized EFT1 scalar soft contribution at \(M_F\):

\[
 C_5^{\rm native-as-UV\ hypothesis}(M_F)
 = C_5^{\rm native\ formal}(M_F)-C_5^{\rm EFT1\ soft,finite}(M_F).
\]

This assumption is **NOT established**. If native `MatchedEFT` already had a
soft subtraction, this scenario subtracts it twice and must not be used as a
physical result. It exists only to quantify the size and effect of the
unresolved matching distinction. **This is not a completed full one-loop
MSbar Wilson boundary**.

The hypothetical coefficient **replaces** (does not add to) the Batch-4
formal native C5 boundary. The existing Batch-2 strict EFT1 Wilson running
with Weinberg mixing, Workstream-C common two-scalar threshold, and SM+
Weinberg running are reused as-is. No additional matching contributions,
spinor signs, wavefunction counterterms or beta-function terms are fabricated.

## Windows PowerShell commands

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive `
  -LiteralPath "$HOME\Downloads\t3_hard_soft_matching_batch5_all16.zip" `
  -DestinationPath . -Force

.\.venv\Scripts\python.exe -m pytest -q tests\test_t3_matching_soft_subtraction_batch5.py

.\.venv\Scripts\python.exe scripts\RunT3HardSoftSubtractionBatch5.py
.\.venv\Scripts\python.exe scripts\RunT3HardSoftSubtractionBatch5.py --show
```

**Resume blocked models only:**

```powershell
.\.venv\Scripts\python.exe scripts\RunT3HardSoftSubtractionBatch5.py --resume
```

**Recompute all Batch 5 outputs (does not modify Batch 1–4):**

```powershell
.\.venv\Scripts\python.exe scripts\RunT3HardSoftSubtractionBatch5.py --overwrite
```

**Audit all models without SMEFT running**, using a distinct output folder:

```powershell
.\.venv\Scripts\python.exe scripts\RunT3HardSoftSubtractionBatch5.py `
  --audit-only --output output\parallel_t3\hard_soft_matching_audit_batch5_v1
```

### Input prerequisites

Require successfully completed Batch 1, Batch 2 and Batch 4 (all 16), and the
unchanged source tensors they use. The runner reads, but never edits:

- `output/parallel_t3/fhard_triangle_all16_batch1_v1/`
- `output/parallel_t3/strict_one_loop_all16_batch2_v1/`
- `output/parallel_t3/native_finite_scenario_all16_batch4_v1/`

The same `MF=1e5 GeV`, `MS=1e3 GeV`, `mu_low=173 GeV` and conditional
scalar/fermion bases are used. The **saved MS scalar mass matrix** is used as
the triangle's IR mass regulator, matching Batch 1's already-implemented
one-loop kernel. The MF quartic is from the same saved Batch-1 input.

### Output

`output/parallel_t3/hard_soft_matching_all16_batch5_v1/`

- `summary.json`: complete/blocked counts and provenance hashes
- `inputs/<model>/hard_soft_and_native_Laurent.npz`: at three renormalization
  scales: local hard, EFT1 soft, full massive UV, power-suppressed remainder;
  plus native pole, soft pole, unaccounted residue, hypothetical new C5
- `inputs/<model>/hard_soft_and_native_Laurent.json`: matching cancellation
  residuals, source hashes, explicit physical limitations
- `inputs/<model>/heavy_scalar_pair_kernels.json`: ALL fermion/scalar-pair
  integral kernels and the common-MF source weights, without truncation
- `inputs/<model>/wilson.json` / `wilson_at_MS.npz` and scalar background:
  strict hypothetical branch input in the same existing Workstream-C contract
- `results/<model>/C5_smeft_at_MS.npz`, `C5_smeft_low.npz`,
  `observables.json`, `PMNS_matrices.npz`: hypothetical SMEFT scenario
- `batch5_analytic_vs_native_vs_native_as_UV_all16.csv`: 16-model mass
  comparison across Batches 2, 4, 5 and the pole/cancellation diagnostics

## Remaining physical implementation

The existing native source states there are **no independently identified
soft-region subtractions or counterterms**, and a native GammaCC-to-Weyl
normalization remains conditional. This batch does not infer either. To finish
true physical one-loop matching, calculate and renormalize the actual native
UV/EFT1 amplitudes and operator counterterms in a common scheme, then match
*all* operator classes and upgrade the candidate EFT1 anomalous dimension.

Do not mistake the existence of 16/16 numerical masses in a conditional branch
for having completed physical one-loop matching.
