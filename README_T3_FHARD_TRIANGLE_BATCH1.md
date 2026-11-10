# Batch 1 / 5 — Heavy-fermion hard one-loop triangle in the all-16 T3 chain

## Scope and status

This is a **large, separately runnable Batch 1**, not a replacement of the completed 16/16 v2 benchmark. The sources were read from GitHub `main` at commit `4bd1230256701ac3801af4e83e8c73e5d448ee1e` (2026-11-09 23:30 UTC). All GitHub access was read-only. No original Workstream B input or v2 output is changed.

The batch computes **one explicit one-loop fermion-hard matching topology** and includes it at the correct first EFT boundary, without introducing an extra leading-log term:

`UV -> [ALL three F integrated at MF=100,000 GeV] -> EFT1 with BOTH S1,S2 and hard C5 -> complete C_ijab one-loop candidate RGE -> [BOTH scalars at common MS=1,000 GeV, bubble matching] -> SMEFT Weinberg RGE -> mnu`.

For mass eigenstate `k`, the stored source-exact LLSS coefficient is reconstructed as

`C_k ~ y_(i k) y_(j k)/M_k` with exact Matchete operator coefficients, Clebsches and real scalar field maps; no substitution of a generic Yukawa proxy or sum-only C is made.

The existing `Numerical/lo_certification/completeness/hard_matching.py:evaluate_hard_scalar_channel` implements the hard local triangle kernel

`K_hard(MF,Mk) = 1 + ln(MF^2 / Mk^2)`.

At the common heavy-fermion threshold the computed contribution is

`C5_hard,triangle(MF) = (in project conventions) sum_{k,a,b} [C_k,LLss * lambda_HHss(MF)] * K_hard / (2*16*pi^2)`

with the original hard-matching function's amplitude-normalized Weinberg projector converted by an explicit factor **2** to the canonical full-sum convention used by Workstream A/C (`T_hard = 2*T_master`). This factor is tested in `tests/test_t3_fermion_hard_triangle_batch1.py`.

The MF quartic is obtained from the saved RGBeta background, rebuilt with the v2 source-invariant reduced-potential conventions. The source CG and Yukawa matrices are reconstructed from the **same-run** Workstream B UV endpoint. The native LLSS contribution of all three heavy mass eigenstates must sum back to Workstream B's saved formal LLSS coefficient; otherwise this case is blocked with its exact error. In mixed B(-1)/C(-1), use the saved MF X/Y scalar mixing `U` and the source's epsilon map. In self-bilinear triplet cases use the verified source CG rotation/permutation. The scalar eigenmass Hessian for full-triangle diagnostics is the corrected MS Hessian.

### What Batch 1 does NOT compute

- The **remaining** one-loop heavy-fermion hard topologies, gauge/counterterm and wavefunction/EOM pieces. No complete UV-to-EFT1 one-loop claim is made.
- The complete finite scalar matching at MS; the existing scalar bubble stays unchanged.
- A certified signed `GammaCC -> Weyl` factor, full RGBeta/Matchete scalar coupling normalizations, or a certified full EFT1 anomalous-dimension matrix. The selected `fonseca_eq22_six_term` candidate is unchanged.
- Strict order-by-order one-loop truncation: this is **RGE-improved conditional** evolution, not the strict-one-loop mode.

**Expected milestone:** 16 numerical candidates with the **computed F-hard triangle added** and side-by-side v2 comparison. This batch is not guaranteed to complete until run against the local saved B data. Any source mismatch fails rather than silently inserting zero.

## Complete project files contained in the ZIP

| Intended project path | Purpose |
|---|---|
| `Numerical/parallel_t3/wilson/fermion_hard_triangle_batch1.py` | New source-resolved per-heavy matching adapter, hard kernel, boundary injector |
| `scripts/RunT3FHardTriangleBatch1.py` | New full all-16 runner; source-v2 workflow with MF hard boundary injection |
| `scripts/RunT3FHardTriangleBatch1.ps1` | PowerShell all-16 resume/fresh/test entry point |
| `tests/test_t3_fermion_hard_triangle_batch1.py` | Normalization, frame transport and saved endpoint tests |
| `README_T3_FHARD_TRIANGLE_BATCH1.md` | This document |

These are **new project paths**. Existing six v2 files remain intact. No GitHub write operations are performed.

## PowerShell install and run

```powershell
Set-Location 'C:\thesistemp\T3neutrinomassmodelpipeline'
Expand-Archive -LiteralPath "$HOME\Downloads\t3_fhard_triangle_batch1_all16.zip" -DestinationPath . -Force

.\.venv\Scripts\python.exe -m pytest -q tests\test_t3_fermion_hard_triangle_batch1.py

# First run, or resume any blocked/missing models:
.\.venv\Scripts\python.exe scripts\RunT3FHardTriangleBatch1.py --resume
.\.venv\Scripts\python.exe scripts\RunT3FHardTriangleBatch1.py --show

# Full deterministic recomputation in this batch's OWN result tree:
.\.venv\Scripts\python.exe scripts\RunT3FHardTriangleBatch1.py --overwrite

# One model first (separate output from default all-16, optional):
.\.venv\Scripts\python.exe scripts\RunT3FHardTriangleBatch1.py --model A_-4 --output output\parallel_t3\fhard_triangle_A4_smoke --overwrite
```

Or run the PowerShell wrapper:

```powershell
powershell -ExecutionPolicy Bypass -File scripts\RunT3FHardTriangleBatch1.ps1 -WithTests
powershell -ExecutionPolicy Bypass -File scripts\RunT3FHardTriangleBatch1.ps1 -Fresh
```

## Output paths

- `output/parallel_t3/fhard_triangle_all16_batch1_v1/summary.json`: all 16 statuses and masses, hard C5 norms and blocker details.
- `output/parallel_t3/fhard_triangle_all16_batch1_v1/hard_fermion/<slug>/fermion_hard_triangle.json`: per-model hashes, source/residuals, eigenmasses, explicit limitations.
- `output/parallel_t3/fhard_triangle_all16_batch1_v1/hard_fermion/<slug>/fermion_hard_triangle.npz`: explicit `C5_hard_F_triangle_MF`, source- and canonical-basis `C_LLSS_per_fermion`, scalar mass and quartic, and **diagnostic-only** full UV triangle, soft bubble and EFT1 scalar log.
- `output/parallel_t3/fhard_triangle_all16_batch1_v1/wilson/<slug>/`: full coupled Wilson ODE results including the injected hard MF C5.
- `output/parallel_t3/fhard_triangle_all16_batch1_v1/results/<slug>/`: final scalar matching, SMEFT running and conditional observables.
- `output/parallel_t3/rough_inputs/fhard_triangle_all16_batch1_v1/`: new isolated A/C stream inputs.

The previous output `output/parallel_t3/source_audited_all16_v1/` remains unchanged. The program does not rewrite `Numerical/lo_certification/completeness/hard_matching.py` or any other prior source file.

## Numerical interpretation

The new hard coefficient should be **nonzero**, but it is **not necessarily a large correction** to existing rough masses. Do not interpret a small difference as a missing calculation until inspecting `fermion_hard_triangle.npz`. The source-resolved MF hard term is injected into `C_full`, then participates in the **same coupled EFT1 RGE as every other operator component**. At MS the original C runner still supplies the scalar bubble and evolves the final C5 with the SM + Weinberg RGE.

The next batches will implement the remaining F hard topologies, full EFT1 operator closure, remaining scalar threshold terms and strict numerical aggregation, without erasing this independent Batch 1 comparison.

### Automatic comparison with the completed v2 benchmark

If the completed v2 `source_audited_all16_v1/summary.json` exists, the Batch 1
runner also writes `mass_comparison_vs_source_audited_v2.csv` under its OWN
output directory. Each model has the three original masses, the three
Batch 1 masses, their absolute and relative differences, and the MF hard
C5 Frobenius norm. Missing or blocked models are not fabricated. These
comparisons measure the implemented matching correction in the chosen
conditional benchmark, not a physical uncertainty or certified sensitivity.

```powershell
Import-Csv 'output\parallel_t3\fhard_triangle_all16_batch1_v1\mass_comparison_vs_source_audited_v2.csv' |
    Format-Table model,m1_delta_eV,m2_delta_eV,m3_delta_eV -AutoSize
```
