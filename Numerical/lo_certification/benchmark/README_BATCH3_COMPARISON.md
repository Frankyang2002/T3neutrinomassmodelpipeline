# Batch 3 — A_-4 and B_-1 conditional EFT versus independent direct UV

**Status:** comparison only, not an EFT matching calculation or physical neutrino mass prediction.

The current GitHub repository (`main` read-only inspection on 2026-11-03,
commit `ab01e9eccbbfe79c1277f27ee3da0e992f01fd6b`) already reports an
independent full-mass direct-UV numerical reference and a separate conditional
all-16 EFT1/SMEFT scalar-matching release. This package **does not replace**
either solver or claim native Matchete GammaCC/Weyl normalization.

## Source inputs

In `--source-mode github-reference` (default), use the exact 3x3 numerical
matrices and recorded mass singular values copied from these inspected GitHub
reports into `batch3_github_reference.json`:

- `output/general_observable/lo_certification_direct_uv/native_uv_reference_01/models/A_-4.json`
- `output/general_observable/lo_certification_direct_uv/native_uv_reference_01/models/B_-1.json`
- `output/general_observable/lo_joint_scalar_matching/conditional_all16_batch2_v1/summary.json`
- `output/general_observable/lo_fermion_boundary/releases/source_021dcc1e8fe8ba0f/{A_m4,B_m1}/boundary.json` (flavour support)

The snapshot includes the Git blob IDs and source commit; it is a portable
reference if the corresponding source releases have not yet been checked out
locally. It does not claim to be a new numerical computation.

In `--source-mode local`, read and validate those **existing** exact project
files in the user's Windows checkout. The validator recomputes three Takagi
singular *magnitudes* for each complex symmetric neutrino mass matrix using
NumPy SVD (no Takagi basis reconstruction), compares them to the recorded
values, checks stated conditional status, and checks source scales.

## Meaning of the comparison

Direct UV masses use frozen input couplings at `mu_UV=10^7 GeV`. The
conditional EFT response uses an **illustrative** running scalar endpoint,
with all fermions integrated together at `mu_F=10^5 GeV`, both scalars active
in EFT1, and both scalars integrated together at `mu_S=10^3 GeV`.

The direct-UV result is expressed in eV and the EFT result in eV per
unknown signed factor `N_GammaCC`.  For a pure global multiplicative
convention difference, the matrices in one **common flavour basis** would
be related by one complex scale `z`. We compute:

```text
z = vdot(M_EFT, M_UV) / vdot(M_EFT, M_EFT)
residual = ||M_UV - z M_EFT||_F / ||M_UV||_F
```

**The diagnostic is not a determination of `N_GammaCC`.** Different input
scales, Yukawa running, scalar masses, relative phases, incomplete hard-loop
matching and potential flavour basis differences prevent assigning this best
scale any physical meaning. The reported UV off-diagonal Frobenius fraction
and the post-fit residual identify where further source alignment is needed.
The original provisional LLSS source tensors have **only p=q lepton-flavour coordinates**:
60/60 for A_-4, 108/108 for B_-1. Their 3x3 EFT matrices are flavour
diagonal in the reported axes, unlike the independent UV references. This
is an unresolved source/input/flavour-basis provenance discrepancy; the
current evidence does not establish which calculation (if either) is wrong.

## Installation / run in PowerShell

All paths in this ZIP are **complete new files** relative to the root:

- `Numerical/lo_certification/benchmark/uv_eft_compare.py`
- `Numerical/lo_certification/benchmark/batch3_github_reference.json`
- `Numerical/lo_certification/benchmark/README_BATCH3_COMPARISON.md`
- `scripts/RunLOBenchmarkBatch3Comparison.py`
- `tests/test_lo_benchmark_batch3_comparison.py`

```powershell
$repo = 'C:\thesistemp\T3neutrinomassmodelpipeline'
$zip = Join-Path $HOME 'Downloads\T3_Batch3_UV_EFT_Comparison.zip'
if (-not (Test-Path -LiteralPath $zip)) { throw "Missing ZIP: $zip" }
Expand-Archive -LiteralPath $zip -DestinationPath $repo -Force
Set-Location $repo
$env:PYTHONPATH = (Get-Location).Path
python -X utf8 -m unittest discover -s tests -p 'test_lo_benchmark_batch3_comparison.py' -v
python -X utf8 scripts\RunLOBenchmarkBatch3Comparison.py --check-only
python -X utf8 scripts\RunLOBenchmarkBatch3Comparison.py --run-id source01
```

If the *full* Batch 2 release and direct UV results are already present in
this checkout, compare **local** files in a second, separately owned run:

```powershell
python -X utf8 scripts\RunLOBenchmarkBatch3Comparison.py --source-mode local --check-only
python -X utf8 scripts\RunLOBenchmarkBatch3Comparison.py --source-mode local --run-id local01
```

The runner writes only under
`output/general_observable/lo_certification_benchmark/batch3_uv_eft_<run-id>/`
(`summary.json` and `comparison.md`). It never overwrites an existing run.
No source code, previous outputs, or GitHub repository is modified.

## Gates remaining

- Native signed four-point amplitude / GammaCC-to-Weyl and master operator
  prefactor; Dirac/Majorana crossing with original CG factors.
- UV and EFT1 couplings/fermion mass/flavour axes consistently at `mu_F`,
  and scalar inputs transported to `mu_S` with matched conventions.
- Full hard one-loop heavy-fermion matching and EFT1 operator running/soft
  subtraction at the requested perturbative order.

The numerical comparison does not attempt any of these gates and cannot
promote conditional numbers into physical masses.
