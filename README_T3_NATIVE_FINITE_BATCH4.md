# T3 Batch 4 — conditional all-16 native Matchete LLHH finite scenario

**Install** these complete new files under `C:\thesistemp\T3neutrinomassmodelpipeline`:

- `Numerical/parallel_t3/wilson/native_finite_batch4.py`
- `scripts/RunT3NativeFiniteBatch4.py`
- `scripts/RunT3NativeFiniteBatch4.ps1`
- `tests/test_t3_native_finite_batch4.py`
- `README_T3_NATIVE_FINITE_BATCH4.md`

## Objective and what is actually calculated

Batch 3 demonstrated all-16 one-loop Coleman–Weinberg Higgs-potential
thresholds, but left the native Matchete first-threshold one-loop LLHH sector
**unanalyzed numerically**. Batch 4 numerically evaluates the **entire saved
native PL finite-kernel inventory** in the tracked all-16 pole/log bridge,
including all **88** PL source kernels, conditional additional quartics,
three evolved heavy-fermion flavour eigenstates, and the B(-1)/C(-1)
X/Y mixing matrix. PR is the conjugate and is excluded from the PL C5.

The original raw Matchete F-only expression is, term by term,

```
K_r * (1 + 1/epsilon + Log[mubar2 / MF_k**2])
```

The new **conditional pole-dropped comparison scenario** uses

```
C5_native_formal(MF) = signed_factor * Sum_{r,k} K_r^{PL}(MF,k)
                               * (1 + log(MF**2 / M_k**2))
```

where each `K_r` already includes exactly one `hbar=1/(16*pi**2)`. It preserves
source signs, Yukawa bars/conjugations, scalar quartic phases, native
`T3MEUij` spurion powers, and full three-by-three complex flavour.
The selected Higgs/spinor SU2 template is the native double-epsilon LLHH
structure; `--native-spinor-factor 1` is **only the previous conditional
normalization**, not an independently established GammaCC-Weyl factor.

Source conditional `r_RGBeta = A @ m_Matchete` is solved with the entire
V3 complex matrix, including non-diagonal B/C and E extra-quartic mixing;
no zipped/diagonal mapping is assumed. The original running B files are
only **read**, never changed.

**Important limitation:** The native raw Matchete release states explicitly
that its F-only one-loop `MatchedEFT` has **not** been independently shown to
be `UV hard - EFT1 soft + counterterms`, and its 1/epsilon poles are not
certified as MSbar counterterms. Therefore, removing those poles produces
an *exploratory, reproducible finite scenario*, **not the physical one-loop
hard matching**. This package does NOT make a claim of full one-loop matching.

## Strict matching and no double counting

The existing strict Batch 2 computation is

```
C5_B2_SMEFT = C5_analytic_triangle(MF) + delta_C5_EFT1_running(tree0)
                                     + delta_C5_scalar_pair(MS)
```

Batch 4 runs **an alternative**:

```
C5_B4_SMEFT = C5_native_formal(MF) + delta_C5_EFT1_running(tree0)
                                     + delta_C5_scalar_pair(MS)
```

The previously included analytic fermion triangle is **removed, not added**.
The strict full Wilson coefficient is updated at MF and propagated directly
to the strict MS endpoint without re-evolving an O(hbar) hard term. The
full EFT1 candidate anomalous-dimension **Weinberg mixing** remains the
same as Batch 2. Both scalar multiplets are integrated out simultaneously
at common MS=1000 GeV. The existing SM + Weinberg solver runs below MS.
Batch 3 Higgs-mass/quartic CW matching is unchanged and **not silently fed
back into strict C5**.

## PowerShell

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive `
  -LiteralPath "$HOME\Downloads\t3_native_finite_batch4_all16.zip" `
  -DestinationPath . -Force

.\.venv\Scripts\python.exe -m pytest -q tests\test_t3_native_finite_batch4.py

# Full 16-model F native conditional -> strict EFT1 -> scalar matching -> SMEFT:
.\.venv\Scripts\python.exe scripts\RunT3NativeFiniteBatch4.py

.\.venv\Scripts\python.exe scripts\RunT3NativeFiniteBatch4.py --show

# Resume only blocked models:
.\.venv\Scripts\python.exe scripts\RunT3NativeFiniteBatch4.py --resume

# Audit all 16 native source/mapping and preflight, no matching/SMEFT run:
.\.venv\Scripts\python.exe scripts\RunT3NativeFiniteBatch4.py --audit-only `
   --output output\parallel_t3\native_finite_audit_only_batch4_v1

# Recompute existing Batch 4 output, leaving Batches 1-3 untouched:
.\.venv\Scripts\python.exe scripts\RunT3NativeFiniteBatch4.py --overwrite
```

Optional explicit conditional comparison: `--native-spinor-factor -1` changes
the assumed signed source-to-Weyl conversion; choose a distinct `--output`
for this sensitivity scenario. No phase/sign is claimed to be certified.

## Output artifacts

`output/parallel_t3/native_finite_scenario_all16_batch4_v1/`

- `summary.json`: all-16 source coverage, assumptions, results and blockers
- `inputs/<slug>/native_MF_LLHH_finite_source_candidate.npz`: native C5,
  88 (distributed per model) source-kernel × 3 fermion matrices, old triangle,
  alternative strict MS C5, unchanged EFT1 running delta
- `inputs/<slug>/native_MF_LLHH_finite_source_candidate.json`: all kernel and
  quartic contributions, source SHA256 and unrenormalized-status flags
- `inputs/<slug>/strict_hard_replacement_ledger.json`: additive closure
- `inputs/<slug>/wilson.json`, `wilson_at_MS.npz`, `background.json`,
  `background.npz`: source-scoped Workstream C strict pair
- `results/<slug>/C5_smeft_at_MS.npz` and `C5_smeft_low.npz`: existing C solver
- `results/<slug>/observables.json`: conditional masses and PMNS data
- `batch4_native_vs_batch2_strict_all16.csv`: masses and MF hard changes

The runner does not overwrite any Batch 1, 2, or 3 output.

## Input dependency list

- Completed Batch-2 output (`strict_one_loop_all16_batch2_v1`)
- Original Workstream B (`B_all16_qualitative_v1`)
- Native frozen Matchete Laurent trios (`pole_log_e515d82d118284ba_bridge01`)
- Original conditional V3 quartic dictionary
- Batch-1 `fermion_hard_triangle_batch1.py` for reconstruction of the same
  evolved heavy-flavour mass basis

If a source, numerical mass basis, or CG structure changes, the model is
marked BLOCKED with an explicit reason. The full actual Windows 16/16 run
must be executed before assuming coverage and numerical success.

## Remaining genuine physics gap

The next physical batch still requires a true `UV hard - EFT1 soft + CT`
calculation in the **same** subtraction scheme, wavefunction/operator
renormalization and counterterms, plus an independently fixed signed
GammaCC-to-Weyl factor. The current strict anomalous dimension remains a
candidate. Even if Batch 4 completes 16/16, those gaps are not solved.
