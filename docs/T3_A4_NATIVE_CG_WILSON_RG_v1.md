# A(-4) source-CG Wilson basis and EFT1 RG (v1)

Inspected GitHub `main` at `75e471df5522e6a852c3cca6811aa546a51a7ae5`.
**Local-only files; no GitHub changes.** This is an algebraic SU(2) CG coordinate
map plus a formal dimension-five master-tensor RG transport. It is not a
complete renormalised matching or a physical neutrino mass calculation.

## Implemented

1. Read the exact six nonzero Cartesian `T3MixCG` components and representation
   list in the **fresh** `A_-4/data/eft1_after_F_wilson_seed.json`. Changes to
   its tensor or orientation cause a hard stop rather than an inferred map.
2. Construct the **SU(2)-invariant** spherical `LL S2` contraction. It is
   *not* the same as the equivariant `HH -> 3` coupling used in the reference
   `lambdaT3 HH S1 S2dag` vertex. The triplet invariant metric is
   `J=[[0,0,1],[0,-1,0],[1,0,0]]`. Check all three generator Ward identities.
3. Solve `C_Matchete[i,j,k] = b sum_s U[k,s] C_invariant[i,j,s]`, with
   `b=sqrt(2)/3^(1/4)`, `U` unitary. On the actual PL operator,
   `S2_Matchete = conjugate(U) S2_RGBeta`. Export the associated **orthogonal
   6x6 real-scalar transformation** for three complex S2 components.
4. Form `Wbar[p,q] = sum_r conj(y1[p,r]) conj(y2[q,r])/MF[r]` for all three
   nearly diagonal fermions, with a strict `MF` non-diagonality cutoff of
   `1e-8`. The full non-diagonal chiral rotation is not approximated beyond
   that cutoff.
5. Expand `S1bar*S2` into eight real components and place *all* six lepton
   SU(2) components, all flavour pairs, and both real-scalar orders into the
   **45 x 45 x 12 x 12** symmetric master `C_ijab` polynomial. A factor of
   `1/4` arises from distributing the two interchangeable index pairs; this
   is algebraically tested against the original ordered scalar polynomial.
6. Reuse `Batch 1/trajectories.npz` (already integrated RGBeta renormalisable
   EFT1, F -> S). Reconstruct `Y_ija`, `lambda_abcd`, and 3+1+8 gauge
   generator matrices from the current `Aminus4EFT1Couplings.py` without a
   fresh Matchete or RGBeta run.
7. Call the existing *seven-term* `calculate_master_beta` implementation.
   The default fixed-order mode evaluates

       C(muS) = C(muF) + ln(muS/muF)/(16 pi^2) beta_C^(1)(muF)

   consistently at one-loop order. Optional `--mode full-ode` calls the
   existing full tensor ODE solver, reconstructs running Y/lambda/gauge
   couplings from the Batch 1 saved trajectory, and resums the corresponding
   leading logs in the specified master operator class.
8. Save both the formal projected LLHH coefficient and the **full master
   tensor**, with Ward and flavour residuals, each of seven beta contributions,
   source-file SHA256, coordinate map and explicit non-certification flags.

## Unresolved physics (do not bypass)

- The signed four-component Matchete `GammaCC P_L` -> canonical two-component
  Weyl bilinear factor is still **not independently certified**. The source
  tensor is built with a *formal unit* factor; no physical Wilson normalization
  is asserted.
- Although the **triplet CG coordinate transform** is verified, the complete
  Matchete-to-RGBeta scalar **quartic and Yukawa coupling dictionary** is not.
  The RG uses the named RGBeta parameter values as a provisional input, not
  a proved identification of the Matchete coupling normalisations.
- Full dimension-five operator closure, F-hard counterterms and finite
  UV-hard minus EFT1-soft subtraction remain unverified.
- The joint scalar matching with inherited hard operators, the physical
  symmetric `C5(MS)`, SMEFT RG and neutrino fits are **not** produced.
- The direct hard/soft Laurent result from Patch 2A is **not added here**;
  adding it before operator-normalisation and RG matching is certified
  would risk double counting or an incorrect sign.
- Running the one-loop hard Wilson boundary through an additional one-loop
  RGE would contribute at two loops; this patch correctly starts RG from
  the complete two-tree-term LLSS source within the inspected master class.

## Project paths in the ZIP

- `Numerical/running/T3A4NativeWilsonBridge.py` — new
- `scripts/RunT3A4WilsonMixing.py` — new
- `tests/test_t3_a4_native_wilson_bridge.py` — new
- `docs/T3_A4_NATIVE_CG_WILSON_RG_v1.md` — new

## PowerShell

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$HOME\Downloads\T3_A4_NativeCG_WilsonRG_v1.zip" -DestinationPath . -Force
python -X utf8 -m unittest discover -s tests -p test_t3_a4_native_wilson_bridge.py -v

# Fast fixed-order one-loop calculation (7 RGE terms, complete master tensor):
python -X utf8 -m scripts.RunT3A4WilsonMixing --mode fixed-order

# Optional genuine coupled master ODE with running renormalisable inputs.
# This operates on all 45^2 x 12^2 complex components and is more costly.
python -X utf8 -m scripts.RunT3A4WilsonMixing --mode full-ode `
  --out output\parallel_t3\a4_native_wilson_rg_v1_full_ode\A_-4
```

Requires installed `numpy`, `scipy` and existing unchanged Batch 1 files,
including `trajectories.npz`, plus the already committed A(-4) fresh matching
output. Both modes fail closed if expected source conventions change. Output
folders are never overwritten.

Inspect fixed-order output:

```powershell
Get-Content output\parallel_t3\a4_native_wilson_rg_v1\A_-4\A4_formal_tree_master_RG.json
```

The solver also writes `A4_formal_master_tensor_muS.npz` with the complete
complex `(45,45,12,12)` tensor. **Never copy this formal result into**
`output/parallel_t3/full_general_complete_input/A_-4/renormalized_c5_at_MS.json`.

## Verification scope

Offline source-level tests cover the exact Cartesian tensor, SU(2) Ward
identities, unitary CG map, real orthogonal scalar map, three-flavour heavy
contractions, operator-index symmetrisation and reconstruction, input
fail-closed checks, and the Batch 1 trajectory adapter using a mocked saved
trajectory. They do not reproduce a real Mathematica matching session or
run the full coupled ODE on the user's real project data. That execution
must happen locally after installation.
