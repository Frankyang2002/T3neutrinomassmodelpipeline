# Workstream 4: all-16 direct UV LO one-loop certification (additive)

## Exact scope and physical status

This layer **does not** modify `Numerical/direct_one_loop/core.py`,
`Numerical/direct_one_loop/comparison.py`, any pre-existing benchmark,
Matchete or RGBeta code, or any earlier UV/EFT result. All new code and
reported products stay within the owned paths.

It evaluates the complete loop **for all component interactions actually
supplied in a source-derived input** using the original exact `p=0` UV engine:

\[
\mathcal{L}_{m}=-\tfrac12\psi^T M_F\psi+\mathrm{h.c.},\quad
\mathcal{L}_{Y}=-\nu_a\psi_i(Y_{ais}\phi_s+Z_{ais}\phi_s^*)+\mathrm{h.c.},
\]
\[
V_2=\phi^\dagger A\phi+\tfrac12(\phi^T B\phi+\mathrm{h.c.}),\quad
\phi=(x+iy)/\sqrt 2,
\]
\[
\mathcal{H}_{\mathrm{real}}=\begin{pmatrix}
\Re(A+B)&-\Im(A+B)\\
\Im(A-B)&\Re(A-B)
\end{pmatrix},
\]
\[
(m_\nu)_{ab}=\frac{1}{16\pi^2}\sum_{k,r}M_k\Gamma_{akr}\Gamma_{bkr}
\,B_0^{\overline{\mathrm{MS}}}(0;M_k^2,m_r^2;\mu).
\]

No heavy-mass expansion is performed. The exact engine performs fermion Takagi
and scalar real-Hessian diagonalization. The normalisation of the actual
Matchete four-component `GammaCC` convention is **not inferred** from the
above independently defined two-component convention. The existing source
adapter remains *provisional*.

The **EFT** being benchmarked has precisely two simultaneous heavy-particle
thresholds: all heavy BSM fermions together at `mu_F` to EFT1 with both BSM
scalars, followed by both BSM scalars together at `mu_S` to SMEFT. SM fermions
are retained. The direct loop has neither threshold nor EFT expansion.

## Source inspection (read-only GitHub main)

Commit: `2c29e5ba2c0350b7e9be292e02ceb8bc7c52c0d3`.

- `Lagrangian/model/T3ModelCatalog.wl`: A (1,3,2), B (2,2,1), C (2,2,3),
  D (3,1,2), E (3,3,2); total 16 allowed alpha cases.
- `Lagrangian/interactions/T3Topology.wl`: ordinary vertices
  `Bar[L] F^c S1`, `Bar[L] F S2^dag`; mixing
  `H H S1 S2^dag + h.c.`; opt-in independent opposite-scalar extra Yukawa
  for A/D/E with alpha=-2 or 0. No frozen extra Yukawa matrices were supplied.
- `Lagrangian/interactions/ScalarPotential.wl`: distinct native
  `lambdaHiInvk` CG portal families, with `PlusHc` convention.
- `Numerical/direct_one_loop/core.py`: exact finite MSbar `B0`, Weyl Takagi,
  Hermitian A/symmetric B, 3x3 complex symmetric neutrino matrix, retained
  vertex-pair and `(fermion,real scalar)` contributions, cancellation and pole
  diagnostics.
- `scripts/BuildT3DirectOneLoopInputs.py`: restricted candidate translation;
  it *explicitly* does not certify spinor phases, Majorana charge pairing,
  native portal dictionary or all general couplings.
- `output/general_observable/lo_direct_one_loop/convention_batch1/`
  `integrated_convention_audit.json`: 16 algebraic candidates, **0 physical**
  source certifications.
- `output/general_observable/lo_direct_one_loop/source_audits/`
  `native_portal_identity_decomposition.json`: 32/32 tensor geometries in
  native span, **no physical coupling map**.

These findings come from the repository's existing evidence; a successful
synthetic test or a zero UV pole does not upgrade the physical status.

## New checks, without phase invention

For all 16 models the runner:

1. checks the exact frozen benchmark SHA-256, representation and UV base scale;
2. checks active source guard SHA-256 values against the previous audit and
   reports stale/missing evidence (never silently promotes old results);
3. constructs the candidate electric charges of **every** Weyl and scalar
   component, validating every nonzero `M_F`, `A`, `B`, `Y` and `Z` entry;
4. checks Hermiticity, symmetric tensors and an independent complex-field
   witness for the complete real quadratic Hessian, including possible complex B;
5. checks existing source charge-channel evidence **and candidate complex
   group-contraction weights** when available; neither establishes physical
   loop multiplicity or spinor normalisation;
6. evaluates the **unchanged** exact direct UV engine, checks UV-pole
   cancellation and complex symmetry, and preserves its full individual terms;
7. independently decomposes the exact mass into unordered pairs of **source
   scalar charges** before mass-basis rotations. It retains charge interference,
   per-heavy-fermion, per-real-scalar and per-vertex-pair matrix contributions.
   It checks that their coherent sum reproduces the unchanged engine;
8. marks weak channels relative to the largest norm **without discarding
   them**; records strong cancellations separately;
9. conditionally compares to an EFT physical `3x3` matrix only where the
   benchmark SHA, scale and same-renormalized-parameters flag agree.

A Dirac Takagi mass eigenvector is *not* assigned a unique electric charge:
charge projectors act on source component vertices **before** diagonalization.
Charge decomposition of finite MSbar matrix contributions is bookkeeping,
not a statement that every component is separately physical or scale invariant.

**No automatic physical certification is produced** by this layer because
independent Matchete-to-Weyl phases, exact portal coupling normalisation, mass
conventions and optional interaction inputs are absent. The code deliberately
sets `physically_certified=false`; this is not a newly discovered failing
numerical calculation.

## Input, EFT and output contracts

Default input: existing
`output/general_observable/lo_direct_one_loop/provisional_inputs/A_-4.json`
etc. Alternative: pass `--input-dir input/lo_direct_one_loop` after running the
source adapter. Inputs remain read-only and source provenance is verified.
Each input supplies `model`, `fermion_mass_gev`, `scalar_A_gev2`,
`scalar_B_gev2`, `vertices` and `renormalization_scale_gev`. Complex tensors
use `{"real":...,"imag":...}`. An input may contain arbitrary additional
verified vertices and scalar bilinears; they are not fabricated from missing
benchmark couplings.

EFT reference: see `eft_reference_schema.json` here. The existing
`compare_same_point()` requires `physical_dimension="mass_matrix_gev"`,
`benchmark_sha256`, `renormalization_scale_gev`, `scheme`,
`same_renormalized_uv_parameters=true`, and a complex `matrix_gev`.
At frozen UV base `mu=1e7 GeV`, the benchmark's two EFT thresholds are
`mu_F=1e5 GeV`, `mu_S=1e3 GeV`; an EFT mass at `mu_S` is *not* the same
input/scale point as the frozen UV mass at `1e7 GeV` unless complete running
has been provided. Exact UV/EFT equality is not required outside the expansion.

The runner creates, only inside
`output/general_observable/lo_certification_direct_uv/`:

- `all_16_certification_status.json`, `all_16_summary.csv`;
- `models/A_-4.json` etc, with full source audits, exact 3x3 complex matrix,
  Takagi masses, **all** original `(pair,fermion,scalar)` terms and charge
  decomposition matrices, suppression/cancellation metrics and EFT comparison.

It refuses to overwrite existing report files. To repeat, specify a **new**
`--output-dir` subdirectory in the owned output root. The package includes
read-only, existing GitHub provisional matrix snapshot and status files
independent of the new local runner outputs.

## PowerShell commands

From the existing **local** repository root (GitHub remains read-only):

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
# Refuse to collide with any previous worker-owned file before extraction.
$zip = "$HOME\Downloads\T3_LO_direct_UV_all16_batch.zip"
$existing = (tar -tf $zip) | Where-Object { Test-Path -LiteralPath $_ }
if ($existing) { throw "Package path already exists: $($existing -join ', ')" }
Expand-Archive -LiteralPath $zip -DestinationPath .
python -m pip install numpy sympy
python -m unittest discover -s tests -p 'test_lo_direct_uv_*.py' -v
python scripts/RunLODirectUVCertification.py --project-root .

# A further run must use a NEW, unused output path:
python scripts/RunLODirectUVCertification.py --project-root . `
  --output-dir output/general_observable/lo_certification_direct_uv/recheck_01

# Optional strict same-input/same-scale EFT reference:
python scripts/RunLODirectUVCertification.py --project-root . `
  --output-dir output/general_observable/lo_certification_direct_uv/eft_compare_01 `
  --eft-reference-dir input/eft_same_uv_point
```

Do not interpret `physical_source_certified=0` as 16 failed loop integrals.
A file can be exactly evaluated and still be provisional because its source
normalisation was not independently certified.

## Verification limits and outstanding tasks

- Source-allowed general Yukawas at A/D/E alpha=-2,0 have no frozen numeric
  matrices. General scalar bilinear/portal coefficients are not globally fixed.
- Native `GammaCC` phase and four-spinor-to-Weyl Feynman-rule normalisation,
  and RGBeta-to-Matchete coupling signs (including lambdaT3), remain blocking.
- The frozen Yukawas are real; the 16 existing numeric matrices have vanishing
  imaginary parts. Complex tests exercise the kernel, but do **not** prove
  all-16 physical CP phases.
- Strong cancellations can make the unchanged `core.takagi` relative symmetry
  gate sensitive to sub-machine-level asymmetry of complex toy inputs. No
  change was made to that separately owned engine. This is a future numerical
  stability check, not grounds to manufacture physical phases.
- The Windows local project tree is not mounted into the artifact-building
  environment. **No claim is made to have re-executed real source-derived
  all-16 numerical inputs locally here.** The 16 matrix snapshots are explicitly
  the existing GitHub provisional results at the recorded commit. The runner
  is intended for the user's local files.
