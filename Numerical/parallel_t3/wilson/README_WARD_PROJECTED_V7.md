# A(-4) Ward-singlet boundary projection — controlled diagnostic, v7

**Purpose:** Test the numerical EFT1→simultaneous scalar threshold→SMEFT
pipeline with an *explicitly altered* A(-4) `LLSS` boundary. This is **not**
a correction to Matchete's source-to-canonical tensor dictionary. It does not
certify neutrino masses, Wilson normalisations, native CG phases, ADM closure
or complete finite matching.

## Evidence from local v6 run

- `A_-4` initial full Wilson at `MF=100000 GeV`: `SU2 Ward ≈1.509`.
- A(-4) Yukawa and full quartic at MF: Ward residual ≈0.
- A(-4) evolved MS LLHH SU2 Ward ≈1.535.
- Most other models also suffer noninvariant source quartic tensors, so **v7
  deliberately does not attempt to repair or rerun the other 15 models**.

The Stream1 source's formal-to-canonical heavy-field transformation was
explicitly selected as identity and is not independently certified. The
unphysical original tensor is never overwritten or silently accepted.

## Mathematical operation

Take the full canonical `C_ijab(MF)` from Stream1's existing A(-4) input.
Only select `i,j` in the three lepton-doublet Weyl pairs and `a,b` in BOTH
heavy scalar multiplets. For each symmetric flavour pair `(p,q)` construct
an explicitly gauge-invariant tensor `T_pq` using the repository's
`synthetic_invariants.synthetic_llss` **solely as an algebraic basis**.

The orthogonal Frobenius projection is

```
C_LLSS^(proj) = sum_p<=q T_pq * <T_pq,C_LLSS>/<T_pq,T_pq>.
```

The six flavour templates are orthogonal because they occupy disjoint
flavour-index blocks. The code validates all gauge Ward identities in the
full 45-Weyl tensor after replacement; all non-LLSS entries, including any
independent initial hard C5, are retained bitwise. The report includes the
absolute and relative size of the discarded non-invariant component and the
original source hashes.

No claim is made that this projected component is the actual Matchete
matching coefficient. Its phase and normalisation **can be different** from
the true native-to-canonical map. If almost no invariant component survives,
the run is blocked instead of producing trivial neutrino masses.

## Numerical evolution and gate

1. Read Stream1 canonical background and boundary with original SHA checks.
2. Check the unmodified MF Yukawa and full quartic gauge Ward identities.
3. Project the A(-4) boundary to the invariant template span; record purity.
4. Check the **existing** candidate ADM beta at MF preserves Ward identities.
5. Run **existing** `workflow.run_wilson` DOP853 MF→MS with both scalars.
6. Require the evolved complete `C_ijab(MS)` to satisfy Ward identities.
7. Write an isolated Stream2 diagnostic endpoint, run unmodified strict C
   preflight and existing C simultaneous-MS scalar matching plus SM RGE.
8. Before publishing any C output, relabel both observable and manifest to
   `WardProjectedBoundaryDiagnosticOnly`, never `RoughNumericalBenchmark`.

A failure at any stage produces `output/.../ward_projected_diagnostic_v7.blocked/`
with a JSON blocker. For a successful run, results are written under
`output/parallel_t3/ward_projected_diagnostic_v7/`:

- `ward_projection_summary.json`: input hashes, original/restored Ward residuals,
  fraction of LLSS source retained, and **diagnostic-only** mass outputs.
- `wilson/A_m4/wilson.json` and `wilson_at_MS.npz`: isolated projected boundary
  evolution, source explicitly marked `ward_projected_boundary_diagnostic`.
- `adapted/A_m4/*`: C-owned adapter inputs.
- `results/A_m4/manifest.json`, `observables.json`, plus existing C output NPZs:
  conditionally generated, all **diagnostic-only**.

## PowerShell (repo root)

```powershell
$zip = Join-Path $env:USERPROFILE 'Downloads\t3_rough_A4_ward_singlet_v7.zip'
Expand-Archive -LiteralPath $zip -DestinationPath . -Force

.\.venv\Scripts\python.exe -m pytest -q tests\test_t3_ward_singlet_v7.py

# Fast first step: check the A4 source projection AND candidate beta at MF.
# No RGE or scalar threshold is integrated.
.\.venv\Scripts\python.exe scripts\RunT3WardProjectedA4V7.py `
    --project-only --overwrite

# If the above passes, evolve EFT1 only (no scalar matching).
.\.venv\Scripts\python.exe scripts\RunT3WardProjectedA4V7.py `
    --skip-c --overwrite

# Re-run full EFT1 + C + SMEFT in a separate, diagnostic-only directory.
.\.venv\Scripts\python.exe scripts\RunT3WardProjectedA4V7.py --overwrite

$d = Get-Content output\parallel_t3\ward_projected_diagnostic_v7\ward_projection_summary.json -Raw | ConvertFrom-Json
$d.boundary_projection | Select-Object kept_LLSS_frobenius_fraction,discarded_LLSS_frobenius_fraction,original_Ward,projected_Ward
$d.MF_projected_boundary_beta_Ward
$d.MS_full_Ward
$d.masses_eV_DIAGNOSTIC_ONLY
```

If C blocks due to a genuine gauge or operator-contract mismatch, retain that
blocker: **do not** use `--allow-projection-residual`. Note the native scalar
map and matching diagram coefficients are still unverified even if this run
returns three numerical masses. Do not mix any of these results into the
original `rough_results/rough_fermion_first_v1` genuine output tree.

## Source files

`Numerical/parallel_t3/wilson/ward_singlet_v7.py`

`scripts/RunT3WardProjectedA4V7.py`

`tests/test_t3_ward_singlet_v7.py`

`Numerical/parallel_t3/wilson/README_WARD_PROJECTED_V7.md`
