# Workstream A v9 — 58 native Higgs-portal real-scalar actions (conditional)

## Scope

**Do not mistake the constructed real `Lam` for RGBeta's native `Lam` evaluator.**
This is a single, bounded Workstream A follow-up to the successful v8 native
S2 symbolic reduction and 14/14 basis-invariant spectral sign audit. It
independently expands the **existing native v5 registered portal operators**
into real scalar components, using RGBeta's source-defined `6*Matrix[coupling]`
coefficient, `SelfConjugate=True` metadata and 24 unique arrangements.

Read-only public GitHub source reviewed at
`Frankyang2002/T3neutrinomassmodelpipeline`
`2bc2e4b516776212645370105cc63495a3c47348` (2026-11-09):

- `Lagrangian/interactions/ScalarPotential.wl`
- `RGE/running/rgbeta/T3RGBetaModel.wl`
- `Numerical/parallel_t3/wilson/portal_tensor_audit.py` (v5)
- `Numerical/parallel_t3/wilson/portal_orientation_v6.py` (v6)
- `Numerical/parallel_t3/wilson/s2_intertwiner_v7.py` (v7)
- `Numerical/parallel_t3/wilson/s2_spectral_v8.py` (v8)
- `Numerical/parallel_t3/wilson/scalar_action_audit.py` (v4)

The package adds only a new v9 module, test file, launcher and document. No
production source, prior audit output, GitHub file or physical Wilson RGE is
modified. All input source hashes are revalidated via the existing v8 chain.

## Numerical assertion

For each `Hbar H Sbar S` native Matchete portal tensor `T_k`, with real
coupling `c_k` and native `phase_k` in `{1,i}`, the **Matchete** interaction
polynomial is

\[
V_{M,k} = \frac{c_k}{2}\left[p_k T_{k,ijuv}
  H_i^* H_j S_u^* S_v + \mathrm{h.c.}\right].
\]

For all 14 triplet sectors, the adjoint sign is selected from the independent
v8 spectrum

\[
\operatorname{spec}\!\left(\sum_A t_H^A\otimes t_{S2}^A\right)
=\{-1,-1,\tfrac12,\tfrac12,\tfrac12,\tfrac12\},
\]

**before** comparison against the saved v6 coupling matrix. For doublets the
v6 native RGBeta fundamental Fierz identity fixes the generator-product
contraction. Scalar singlets use the identity tensor. These checks are in the
identified Matchete tensor coordinates and do *not* expose native numeric
RGBeta S2 axes.

Write each complex field `z=(R+iI)/sqrt(2)` and let `P` be the unnormalised
real-to-complex projector (`P x= sqrt(2) z`). v9 constructs the quartic
real-coordinate Lam tensor

\[
\Lambda_{abcd}^{(R)} = 6\,\mathrm{Sym}_4\left[
  \sum_k \lambda_k T_{k,ijuv}
  P^{H*}_{ia} P^H_{jb} P^{S*}_{uc} P^S_{vd}\right],
\quad
V^{(R)} = \frac{\Lambda_{abcd}^{(R)}}{4!}x_a x_b x_c x_d.
\]

Here `lambda_k` is either the original Matchete coupling (after its `1/2
PlusHc` phase projection) or the **conditional** singlet/adjoint coefficient
vector. No second conjugate is inserted for the self-conjugate portal.
The identity `6*2^2/24=1` follows from RGBeta's complex-scalar projectors
and explains the source prefactor. Random *real* scalar inputs and *real*
native Matchete coefficients compare four independently evaluated forms:

1. Matchete `1/2 PlusHc` on the raw native `GetCGTensor` components;
2. Matchete projected complex scalar action;
3. Symmetrized real Lam action with Matchete projected tensors;
4. Symmetrized real Lam action with the conditional RGBeta invariant tensors.

The report requires **16 models, 32 portal sectors, 58 Matchete slots,
14 triplet spectral signs and 14 anti-Hermitian `I` phases**. Missing source,
non-Hermitian terms, wrong `6*Matrix` coefficient, wrong 24-arrangement
metadata, wrong v6 map, wrong v8 sign, stale source, or altered v8 saved
records are rejected. Expected numerical residual: of order machine epsilon;
actual 58-slot Windows result must be measured, never pre-filled.

## What v9 does *not* certify

- The installed RGBeta kernel's native `Lam` / `sDelS` / projector operation is
  **not executed**. Constructed Lam equals the source algebra in the chosen
  coordinates; it is not a native projector observation.
- Native numerical RGBeta S2 axes / cross-package field map remain unspecified.
- `lambda12Cross`, non-portal independent quartics, complete scalar masses,
  exact physical coupling normalization and entire scalar potential remain
  **uncertified**. Prior v4 50 extra quartic + 96 norm quartic action checks
  are *separate conditional results* and do not close this gap.
- B(-1)/C(-1) require mixed-scalar mass diagonalisation before mass-basis matching.
- Signed GammaCC-to-Weyl factors, physical fermion matching, and the complete
  one-loop Wilson ADM are separate gates.
- **Physical Wilson evolution remains blocked.** The mandatory EFT order is
  **UV -> integrate out ALL fermions -> integrate out ALL scalars -> SMEFT**.

## Windows installation and invocation

The ZIP contains full new files at these project-relative paths:

- `Numerical/parallel_t3/wilson/portal_real_action_v9.py`
- `Numerical/parallel_t3/wilson/PORTAL_REAL_ACTION_V9.md`
- `scripts/RunT3ParallelAPortalRealActionV9.py`
- `tests/test_t3_parallel_a_portal_real_action_v9.py`

```powershell
Set-Location 'C:\thesistemp\T3neutrinomassmodelpipeline'
$zip = "$env:USERPROFILE\Downloads\t3_workstream_A_portal_real_action_v9.zip"
Expand-Archive -LiteralPath $zip -DestinationPath . -Force

.\.venv\Scripts\python.exe -m unittest discover `
  -s tests -p "test_t3_parallel_a_portal_real_action_v9.py" -v
if ($LASTEXITCODE -ne 0) { throw 'v9 tests failed' }

.\.venv\Scripts\python.exe -m scripts.RunT3ParallelAPortalRealActionV9 `
  --project-root .
if ($LASTEXITCODE -ne 0) { throw 'v9 real-action audit blocked' }
```

Inputs: the same **immutable** v5–v8 native JSON already generated, including
`rgbeta_s2_reduction_native_v8.json` and
`s2_spectrum_conditional_v8.json`. The script does not call Wolfram; v8's
native witness is read and source-revalidated. The full-repository tests and
launcher require the user's local project. The standalone tests deliberately
do not masquerade as live native RGBeta output.

New output (refuses overwrite):

`output/parallel_t3/wilson/portal_real_action_conditional_v9.json`

Expected conditional status after a successful local run:

`NativePortal58RealActionAlgebraVerified_NativeRGBetaLamProjectorsAndFullScalarNormalizationUnverifiedConditional`

## Remaining blocker order (recommended)

1. Native RGBeta `Lam` / scalar projector action and complete independent
   scalar quartic/mass dictionary (not another S2 numerical basis search).
2. Resolve B(-1)/C(-1) mixed scalar mass diagonalisation and physical scalar
   field coordinates; supply the normalized B workstream background.
3. Native signed GammaCC-to-Weyl conversion, physical fermion boundary and
   one-loop ADM counterterm/crossed-channel certification.
4. Only with certified inputs: evolve the full 45-Weyl real-scalar Wilson
   tensor from `M_F` to `M_S`, perform **one** all-scalar threshold into
   SMEFT, then perform SMEFT running/observables.

This v9 batch does not claim to resolve items 1–4.
