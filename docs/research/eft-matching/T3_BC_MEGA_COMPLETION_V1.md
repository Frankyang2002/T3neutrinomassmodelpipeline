# B/C(-1) Mega-batch v1 — concrete numerical integration and physics gates

## Scope

One integrated local campaign, reusing the already successful **all-sixteen
UV**, B/C **five-real-mass**, B/C **running Matchete 480-term** and **mixed-quadratic
resummation** results. The main physics goal is to determine whether the
heavy-spinor transformation required by the native B/C source is compatible
with RGBeta's matrix beta functions **at complex heavy flavour**. It also
runs the existing one-loop SM+Weinberg solver on the same 1e5 GeV candidate
C5 and running SM inputs, producing rough reference masses at 100 GeV.

No other workstream source is modified and no connection to GitHub is required
at runtime. In particular, this package does NOT modify `T3Topology.wl`,
RGBeta, the native Matchete matching scripts, or the accepted results.

## Complete files to install

- `Numerical/experimental/T3BCMegaPhysicalBridge.py` — explicit complex
  conjugation and spinor-action redefinition, tagged source provenance,
  validated same-scale C5, numerical SM+Weinberg continuation.
- `Numerical/experimental/T3BCComplexUVStress.py` — full coupled UV B/C
  RGBeta beta-function basis covariance with complex 3x3 flavour matrices,
  including two separate real ODE integrations per test case. Checks the
  `y1@conj(U), y2@U, U.T@MF@U` convention against native betas and compares
  the deliberately wrong `y1@U` case; returns actual numerical discrepancies.
- `scripts/RunT3BCMegaCompletion.py` — checkpointed B/C campaign, genuine
  source hash gates for the saved 480-term direct Matchete source, checks
  all-seven native quartic tensors and mixed kernel via v2 source-provenance,
  full 3-generation low-energy SM+Weinberg proxy with two tolerances.
- `scripts/RunT3BCMegaCompletion.ps1` — Windows tests + campaign.
- `tests/test_t3_bc_mega_completion_v1.py` — 15 regression tests, including
  live source-dependent checks where the project files are installed.

## Installation

From `C:\thesistemp\T3neutrinomassmodelpipeline`:

```powershell
Expand-Archive -LiteralPath "$HOME\Downloads\t3_bc_mega_completion_v1.zip" -DestinationPath . -Force
& .\scripts\RunT3BCMegaCompletion.ps1
```

To resume independently successful B and C studies:

```powershell
& .\scripts\RunT3BCMegaCompletion.ps1 -Resume
```

If the native complex UV beta run is impossible because the latest source
exports changed, **do not** label the remaining results certified. The
per-model result records a blocked complex-UV stage and preserves the other
calculations. A fresh output directory can be selected with:

```powershell
& .\scripts\RunT3BCMegaCompletion.ps1 -OutputId t3_bc_mega_completion_second_run
```

For a shorter initial diagnostic only (no direct complex RGE covariance):

```powershell
& .\scripts\RunT3BCMegaCompletion.ps1 -OutputId t3_bc_mega_fast -SkipComplexUV
```

## Exact existing inputs

The runner loads, hashes and checks:

- `output/studies/t3_full_stack_overnight_v1/bc_full_mass/{B_m1,C_m1}_uv_to_eft1_coupled_five_mass_v1.json`
- `output/studies/t3_bc_running_one_scale_triangle_v1/{B_m1,C_m1}_running_one_scale.json`
- `output/studies/t3_bc_mixed_quadratic_kernel_v2/{B_m1,C_m1}_mixed_kernel.json`
- `output/general_t3_stage22/matchete_tensors.json`
- `output/general_t3_mass_eigen_matchete/one_loop_probe/{B_m1,C_m1}_final_one_loop_probe_weinberg_ordered_coefficient.inputform.txt`
- Corresponding direct native Matchete matching JSON.
- Full native RGBeta B/C UV dimensionless and five-real-mass beta exports,
  via the existing `RunT3BCCoupledUVToEFT1` source contracts.

The `T3Topology.wl` source head is checked directly. It generates `y1` with
`CConj[F]` and `y2` with `F`. For a **fixed complex flavour field change**
`F = U N`, the algebraic identities are

\[
CConj[F]=U^* CConj[N],\qquad y'_1=y_1U^*,\qquad y'_2=y_2U,\qquad M'_F=U^TM_FU.
\]

This exact field redefinition is necessary, but does not by itself certify
Matchete's full 4-component Majorana spinor/CG convention. The numerical
RGE covariance study directly evaluates whether the native beta functions
obey the corresponding transformation.

## Outputs

`output/studies/t3_bc_mega_completion_v1/`

- `B_m1_mega.json`, `C_m1_mega.json` — full comparison records.
- `summary.json` — per-model progress, explicit errors and physical gates.

Each model includes (1) file-hashed same-UV-scale coefficient and Matchete
source verification, (2) mixed-quadratic/full-vs-strict matching errors and
seven native tensor checks from the already completed v2 calculation,
(3) algebraic complex spinor action tests, (4) true native UV complex-flavour
beta and finite-interval evolution covariance comparison, and (5) a SM-only
`C5(1e5 GeV) -> C5(100 GeV)` solver convergence test and mass spectrum.

## Strict physics interpretation

1. **The SM-only continuation is a computational proxy, NOT a valid
   sequence of threshold matchings.** It treats F, Y, X as absent even while
   those particles would be active between `1e5 GeV` and their respective
   masses. It must NOT be combined with the prior *different-trajectory*
   EFT1/EFT2 Y/X output.
2. **Complex RGE covariance can fail.** That is a result identifying
   incompatible field/coupling transformations, rather than a numerical
   license to choose the visually closer answer. The report retains both
   same-U and conjugated-U variants; a source convention is not certified
   solely by matching a fixed benchmark.
3. The Matchete one-loop source still represents a **restricted order
   `y^2 lambda`** UV model. The independent mixed quadratic modifies the
   propagator and masses, but its full native Matchete mass-insertion
   treatment and any higher-order interaction vertices have NOT been
   independently matched.
4. The **finite F/Y/X counterterms, original-vs-canonical EFT EOM/field
   dictionary and physical EFT Wilson RG handoff** are not created by this
   package. Physical C5 and neutrino masses remain uncertified.
5. The full 16-model additional Yukawa physical matching is not performed.

This batch is intentionally larger than a single kernel audit; it consolidates
five physics/numerical checks in one pass and produces explicit reference
masses from a single consistent running UV state. It does **not** pretend
that the still-missing physical renormalized matching equations were derived.
