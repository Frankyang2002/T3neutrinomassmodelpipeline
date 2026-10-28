# B(-1)/C(-1): conditional scalar-normalisation and mixed-mass C5 audit (Batch 10)

## Scope

**New independent report only.** Read saved Batch 9 candidate RGBeta UV/EFT1 trajectories for the three historical benchmark cases. Re-evaluate the existing one-scale Matchete *ordered* C5 expression with the actual 1 TeV scalar mass matrix, alongside two explicit **conditional** Matchete/RGBeta coupling identifications. All Yukawa and heavy-fermion matrices remain frozen at the saved historical fitted input. No RGE is rerun; no fit, Matchete matching, historical result, or Batch 9 result is changed.

The script uses **both** `pilot_q0p1_c5_unmixed` and optionally the prior `pilot_q0p1_mu0p01`, with their SHA-256 provenance. The outputs are **not** a certified physical general-scalar neutrino prediction or a sequential two-threshold EFT result.

## Source-verified normalisation issue

From the actual Stage 23 report `output/general_t3_stage22/stage23_quartic_prefactor_candidates.json`:

- For both B(-1) and C(-1), `020020` identifies `lambdaHHdagS2S2` with `lambdaT3ExtraQ020020I1` at tensor ratio `R/M = 1`.
- For both, `020200` identifies `lambdaHHdagS1barS1bar` with `lambdaT3ExtraQ020200I1` at the **corrected-covariant candidate** tensor ratio `R/M = 1`.
- RGBeta's `AddQuartic` factor `6` is cancelled by its real-field/tensor symmetrisation and complex-field convention; it is **not** an extra physical scalar-potential factor.
- Stage 26 `output/general_t3_stage22/stage26_conditional_quartic_action.json` checks the source Lagrangian signs: Matchete uses `+m O + h.c.`, RGBeta uses `-r O + h.c.` under a common field convention.

Thus the source-supported **conditional** map for *these two one-dimensional blocks* is

\[
  m_1 = -r_1,\qquad m_2=-r_2
\]

rather than the earlier Batch 9 identity assumption `m_i = r_i`. Stage 23 and Stage 26 explicitly say that a **physical component-basis map is not certified**. This script fails if their stated source signs, relevant block matrices, total inventories, or uncertified status change. It does not invent a map for the other five B/C additional quartics, whose RG effects are retained but whose direct C5 insertions are not present in the saved coefficient.

## Scalar mass basis

At the 1 TeV RG endpoint, diagonalise

\[
 A=\begin{pmatrix}m_{S1}^2&\mu_2^2\\(\mu_2^2)^*&m_{S2}^2\end{pmatrix},\quad
 \Phi=(S_1,\epsilon S_2^*),\quad \Phi=U(X,Y),\quad U^\dagger A U=\operatorname{diag}(M_X^2,M_Y^2).
\]

Reuse `Numerical.experimental.GeneralT3MixedScalarMassBasis.MixedDoubletMassBasis.build` and feed **its** complex `U`, positive eigenmasses, the saved historical 3×3 Yukawas, `MF`, and the two quartic values into the **existing** ordered `B_m1_final` / `C_m1_final` Matchete one-loop expression. Symmetrise as `C5 = A5 + A5.T` (no complex conjugation), with `hbar=1/(16π²)`.

Unlike Batch 9's original preview, the candidate `C5_joint` here uses the **actual nonzero RG-generated mixing**. The script also provides unmixed controlled projections to isolate mixing effects. This re-evaluation is **not a new matching derivation**, and success does not establish the Matchete native mixed-field `CheckLagrangian`, `HeavyMassBasis`, one-loop EFT1 scalar threshold matching, or the physical signed map.

## Exact telescoping decomposition

For each of the two coupling maps, define, with `Y1,Y2,MF` held at their historical fitted input:

- `C0`: zero-extra endpoint diagonal masses and lambdaT3; zero direct quartics and zero mixing;
- `Cd`: general-quartic endpoint diagonal masses and lambdaT3; zero direct quartics and **mixing projected to zero**;
- `Cm`: same general-quartic endpoint with its actual mixing; zero direct quartics;
- `Cq`: same general-quartic endpoint, nonzero mapped direct quartics, **mixing projected to zero**;
- `Cj`: same general-quartic endpoint, nonzero mapped direct quartics and its actual mixing.

Then *exactly as complex 3×3 matrices*:

\[
 C_j-C_0=(C_d-C_0)+(C_m-C_d)+(C_q-C_d)+(C_j-C_m-C_q+C_d).
\]

The terms are, respectively, diagonal-background/`lambdaT3` effect, mixed-quadratic effect at zero direct quartics, direct-quartic effect at unmixed masses, and quartic–mixing interference. The displayed values are **individual Frobenius norms divided by** `||C0||_F`. Norms are *not additive*, although the complex matrix decomposition is. The report checks numerical closure.

A separate calculation at `C0`'s diagonal masses and `lambdaT3`, with positive identity-mapped quartic insertion and zero mixing, must reproduce the saved Batch 9 preview **exactly within numerical tolerance**. This guards against changing the historical definition of that diagnostic.

An optional `--mixing-run-id pilot_q0p1_mu0p01` compares the **separately evolved** explicit-mixing input scenario to the quartics-only scenario at the same mapping sign; this includes any accompanying RG mass/quartic changes, so it is *not* a pure fixed-background derivative with respect to `mu2`.

The script also records *singular values* of `-v² C5 / 2` in eV, for illustrative one-scale comparison only. They are neither low-energy RG-evolved fitted neutrino masses nor independent oscillation predictions.

## Installation paths

All new files:

- `scripts/AnalyzeT3BCMatchingConventions.py`
- `tests/python/test_t3_bc_matching_conventions.py`
- `docs/T3_BC_CONDITIONAL_MIXED_C5_CONVENTION_AUDIT.md`

All new results go to
`output/studies/general_scalar_extended/<historical_case>/<new_run_id>/summary.{json,md}`.

The three case directories are `small_uv_running_5e-9`, `small_uv_frozen_5e-9` and `large_uv_frozen_0p5`. Existing run-ID directories are never overwritten. The original `output/studies/lambda_5e-9_refit/` and `output/studies/fermion_start_lambda_comparison/` stay read-only.

## PowerShell — repository root

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline

Expand-Archive `
  -LiteralPath "C:\path\to\T3_BC_Conditional_Mixed_C5_Conventions_Batch10.zip" `
  -DestinationPath . -Force

python -m pytest -q tests\python\test_t3_bc_matching_conventions.py

# Re-evaluate existing quartic-only candidate RG endpoint, with its induced mu2:
python -m scripts.AnalyzeT3BCMatchingConventions `
  --all --source-run-id pilot_q0p1_c5_unmixed `
  --run-id matching_conventions_v1

# Optional: use the existing nonzero-input mu2 scenario as well (new run ID):
python -m scripts.AnalyzeT3BCMatchingConventions `
  --all --source-run-id pilot_q0p1_c5_unmixed `
  --mixing-run-id pilot_q0p1_mu0p01 `
  --run-id matching_with_mu_v1

Get-ChildItem output\studies\general_scalar_extended\*\matching_conventions_v1\summary.md |
  ForEach-Object { Get-Content -Encoding UTF8 $_.FullName }
```

If a stage-23/stage-26 JSON, fitted input, beta file, or saved Matchete ordered coefficient has been regenerated, the script **fails closed** instead of substituting a different reference. In that case review source compatibility before rerunning Batch 9 with a *new* run ID. If `matching_conventions_v1` already exists, change `--run-id` to another unused value; do **not** delete prior results.

## Result interpretation

The current calculations can quantify **conditional sensitivity** to an explicit normalisation convention and a nonzero scalar mass-mixing matrix. They **cannot** certify a unique physical nonzero-general-scalar `C5`. The native Matchete mass-eigenstate Lagrangian and the Matchete/RGBeta physical field identification remain unverified, and a complete sequential EFT1 matching/evolution of general-scalar operators has not been performed. Do not use the resulting percentage changes as experimental uncertainties or fitted neutrino predictions.
