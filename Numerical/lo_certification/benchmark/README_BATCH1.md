# Restricted T3 Batch 1 — A_-4 benchmark and signed spinor audit

**Ownership:** exclusively `Numerical/lo_certification/benchmark/`,
`scripts/RunLOBenchmarkBatch1.py`, `tests/test_lo_benchmark_batch1.py`,
and `output/general_observable/lo_certification_benchmark/`.

## What is implemented

- `A_m4_input.json` is a fully explicit, source-pinned restricted numerical
  benchmark contract. Its Yukawa matrices and scalar/fermion input masses are
  the pre-existing **frozen UV** values, not a newly fitted point or values
  already RG evolved to `mu_F`.
- `batch1.py` checks these data against the original frozen benchmark's
  SHA-256 and exact parameter values. It checks SU(2) dimensions,
  hypercharges, positive unbroken scalar quadratic masses, perturbativity,
  Dirac fermion masses, and the hierarchy of both grouped thresholds.
- It reuses the repository's `spinor.boundary.validate_model` and
  `verify_reference_hashes` for the original source-resolved Matchete
  F-tree LLSS tensor. Every original source/reference SHA-256 is checked.
  It checks the V6 A_-4 **zero additional HHSS** result.
- It reuses `spinor.spinors.convention_probe` and `exchange_probe`, and
  independently tests both allowed signs `C` and `C^-1` in a declared Weyl
  basis. Both pass the abstract `C^2=-I`/antisymmetry identities but differ
  by a sign on the Weyl `epsilon` bilinear. Thus these identities alone
  **cannot** fix the native `GammaCC` factor.
- It compares pinned Git blob SHA-1 values against actual local files. A
  local difference is **reported** as source drift; it is never silently
  represented as identical to the GitHub version.
- It writes an immutable JSON numerical/source summary and a physical
  normalization report. It neither recomputes an LLSS tensor from guessed
  couplings nor touches another workstream.

## Gauge representation and selected frozen parameters

For `A_-4`, `T3ModelCatalog.wl` gives SU(2) dimensions `(1,3,2)` for
`(S1,S2,F)` and hypercharges `(-2,-1,-3/2)` in the repository's usual
`Q=T3+Y` convention. The complex fermion doublet is Dirac. The source
contains one native left-chiral Matchete F-tree operator, with `T3MixCG`,
60 independent numerical nonzero LLSS coordinates, and no counted PR
Hermitian conjugate. The previously validated symmetric orbit factors must
not be applied again.

At the level of **source operator families**, the restricted interaction
terms have the schematic form (native SU(2) contractions suppressed)

\[
\mathcal L_Y\supset (y_1)_{\alpha k}\,\overline L_\alpha F_k^c S_1
 +(y_2)_{\alpha k}\,\overline L_\alpha F_k S_2^\dagger
 +\mathrm{h.c.},\qquad
V\supset m_{S_1}^2 S_1^\dagger S_1+m_{S_2}^2 S_2^\dagger S_2
 +\lambda_{T3}\,\mathrm{CG}(H H S_1 S_2^\dagger)+\mathrm{h.c.}
\]

The full contracted potential also retains the listed baseline quartics.
The displayed signs/prefactors are not an identification with the project's
signed **two-component** master action; the original native input remains
authoritative for that identification. The other independent couplings
listed as zero in the JSON are restricted **at mu_UV**.

The frozen benchmark uses `mu_UV = 10^7 GeV`, `mu_F=10^5 GeV`,
`mu_S=10^3 GeV` and `mu_low = 100 GeV`, with Dirac masses
`(100,110,120) TeV`, **unbroken** scalar quadratic mass parameters
`m_S1 = 1142.857142857... GeV`, `m_S2 = 1000 GeV`, and
`lambdaT3 = 0.1` real. The minimum fermion-to-maximum scalar mass
ratio is **87.5**. Full real 3x3 `y1,y2` matrices are in the JSON contract.
All other source-independent extra couplings are fixed to zero **at the
source UV input scale**, not assumed RG-invariant.

This establishes a useful heavy-fermion expansion hierarchy at the frozen
source parameter point. It does not establish physically normalized scalar
mass eigenstates after electroweak symmetry breaking or prove matching
between the UV input at 10^7 GeV and EFT1 at 10^5 GeV.

## The precise unresolved sign/prefactor

In the **declared** Weyl representation with `epsilon=i*sigma2`,
`C=i*gamma2*gamma0` and `PL=(1-gamma5)/2`:

- `Psi^T*C*PL*Chi = +psi^T*epsilon*chi`;
- `Psi^T*C^-1*PL*Chi = -psi^T*epsilon*chi`.

Both candidates obey the native abstract symbolic facts `C^2 = -1` and
`C^T=-C`. These facts do **not** select which numerical matrix (or which
external-leg prefactor) the installed Matchete `GammaCC` represents.
The conditional Dirac exchange factor `1/M` follows from the previously
checked declared UV Weyl action, not an independently normalized native
Matchete connected amplitude.

**The required missing object** is a signed native connected
`l l S1 S2` amplitude for A_-4 in the installed Matchete basis, with
Dirac crossings, `T3MixCG`, external Grassmann order and the precise master
operator prefactor. One must project that amplitude onto a specified Weyl
`epsilon_12=+1` convention. A bare choice of `GammaCC=C` would not be
certification. The saved native one-loop matching also requires a check of
Matchete's hard/soft subtraction semantics before declaring a complete
leading one-loop F threshold.

## Physical status

**Blocked**, by design: `physical_C_LLSS=null`, `physical_C5=null`,
`m_nu=null`. The scalar threshold is **not** executed. Do not pass the
provisional formal polynomial to a physical EFT1 running or scalar-matching
worker. This is a deliberately restricted benchmark, not the most general
scalar potential and not a full neutrino-oscillation fit.

## Install and execute (PowerShell)

From a directory containing the downloaded ZIP:

```powershell
$repo = 'C:\thesistemp\T3neutrinomassmodelpipeline'
Expand-Archive -LiteralPath '.\T3_Batch1_A_m4.zip' -DestinationPath $repo -Force
Set-Location $repo
$env:PYTHONPATH = (Get-Location).Path
python -X utf8 -m unittest discover -s tests -p 'test_lo_benchmark_*.py' -v
python -X utf8 scripts\RunLOBenchmarkBatch1.py --check-only
python -X utf8 scripts\RunLOBenchmarkBatch1.py --run-id source01
```

The last command creates two files:

```
output/general_observable/lo_certification_benchmark/batch1_A_m4_source01/summary.json
output/general_observable/lo_certification_benchmark/batch1_A_m4_source01/physical_normalization_report.md
```

If the local checkout differs from the inspected commit, the summary lists
all differing source blob hashes. If a required frozen source, original
provenance reference or original Matchete boundary is absent/changed, the
runner fails closed with exit code **2**. It neither fabricates a boundary
nor restores the missing file. `--require-physical` is deliberately expected
to return **2** until the missing amplitude and matching evidence exist.

The Python runner does **not** run Wolfram, Matchete, an EFT RGE or a scalar
matching calculation. It will not modify GitHub or any pre-existing result.

## Inspected read-only source and what was actually executed

Remote GitHub main commit inspected on 2026-11-03:
`e581a434c4114ea1b100c50378c085529fb85be7`.
The `GIT_BLOBS` dictionary in `batch1.py` gives individual inspected Git
blob SHA-1 values. The original frozen input also has a pinned SHA-256:
`3cb8dc3e94f2bd6f60571a7a85477638e8fc02109bccfcec0e04055852d81247`.

The archive's isolated `unittest` tests were run in the packaging container:
**10 passed, 1 skipped**. The skipped test requires the user's original
repository source/reference files and needs execution in the Windows checkout.
No Windows execution, native Matchete calculation, signed physical LLSS,
C5 or neutrino mass is claimed.
