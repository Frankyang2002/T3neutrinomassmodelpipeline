# Workstream C v6 — EFT1/EFT2 light-XX scalar-bubble overlap

## Why this is the next check

The user's v5 **local** result is successful for both B(-1) and C(-1): from independently reconstructed source vertices, `18 YY + 32 XY` complex-scalar kernel structures reproduce the native Y/LLHH UV pole, and the associated scalar-bubble MS-bar counterterm cancels that pole. **That is only the Y-containing bubble subset.** It is not a full EFT1/EFT2 calculation.

The immediate missing contribution on the *EFT2 side* is the light scalar `X` loop. Its coefficient must be reconstructed from the **saved Matchete Y tree output**, not copied from the EFT1 source by assumption. This script independently extracts the `LLXX` and `HHXX` tree vertices from:

- Original EFT1 saved F-first sources:
  - `output/general_t3_mass_eigen_matchete/eft1_fermion_first_tree/{B_m1,C_m1}_final_eft1_fermion_first_tree_llss_tree.inputform.txt`
  - `output/general_t3_mass_eigen_matchete/eft1_fermion_first_tree/{B_m1,C_m1}_final_eft1_fermion_first_tree_eft1_transition_eom.inputform.txt`
- Independent **Y-tree-matched EFT2 source**:
  - `output/general_t3_mass_eigen_matchete/eft1_scalar_Y_only_canonical_probe/{B_m1,C_m1}_final_y_only_canonical_match_sm_tree.inputform.txt`

The repository source was inspected read-only. This last EFT2 file exists for both models and contains **no Y field**. For B it has six additive `LLXX` terms and fourteen `HHXX` terms; for C the same. In the original EFT1 EOMEFT the `HHXX` source has twelve additive terms, and the `LLXX` polynomial may be grouped differently. These **source-term counts cannot be compared as coefficients**; the script instead expands their exact native SU(2) components and monomials to test equality.

Pre-existing local prerequisites: v2, v3, v4 and v5 scripts, the v3 native-CG JSON exports for both models, the Y-pole ledger, the successful v5 `summary.json`, and original tree files. **No Matchete rerun is required.**

## Restricted one-loop calculation

In the ordered Matchete neutral component (`lepton SU2 index=1`, `H SU2 index=2`, `GammaCC/Proj[±1]`), define source-derived contraction kernels, with the shared `hbar=1/(16π²)` and `1/MF` removed:

\[
K_{XX}^{(1)} = (LLXX)_{\mathrm{EFT1}} \star (HHXX)_{\mathrm{EFT1}},
\qquad
K_{XX}^{(2)} = (LLXX)_{\mathrm{EFT2}} \star (HHXX)_{\mathrm{EFT2}}.
\]

The star means the *same exact, source-verified* native complex-scalar Wick contraction, SU(2) CG tensor evaluation, external component, and chirality projection as in the v3/v5 audit. The new audit compares `LLXX`, `HHXX` **and** `KXX` separately between EFT1 and EFT2. It fails if a source tensor or coefficient differs.

For the scalar-quartic bubble **subset only**, the UV residues on the two sides are

\[
P^{(1)}_{\rm bubble}=K_{XX}^{(1)}+K_{XY}+K_{YY},
\quad P^{(2)}_{\rm bubble}=K_{XX}^{(2)}.
\]

For each bubble, the logarithmic UV residue is one in the earlier convention (proved by the Euclidean integral in v5). The **restricted**, derived MS-bar counterterm difference is therefore

\[
\delta C^{(1)}_{\rm bubble}-\delta C^{(2)}_{\rm bubble}
=-\frac{\mathrm{hbar}}{M_F\bar\epsilon}
\bigl(K_{XY}+K_{YY}\bigr)
\]

**if and only if** the source-derived `KXX` coefficients agree across EFT1 and EFT2. The script tests this as an identity and compares it independently with the native Y/LLHH pole and the successful v5 coefficient map. It does not use the native pole to assign either counterterm.

The finite light-only bubble at the actual nonzero mass `M_X`, evaluated at `μ=M_Y`, has the common coefficient

\[
B_{XX}^{\overline{\rm MS},\,\rm fin}(0;M_X,M_X;M_Y)
=\ln(M_Y^2/M_X^2),
\]

which multiplies the same `KXX` on **both** sides. The script tests that these finite terms cancel algebraically. Do **not** use `mX=0` and blindly set the dimensionally regulated `XX` integral to zero: the scaleless integral has compensating UV and IR poles, and the scheme comparison still requires proper operator renormalization.

## Installation and exact PowerShell commands

Extract `T3_WorkstreamC_YEFT2LightXXOverlap_v6.zip` to the root of the local project, adding exactly two new files:

1. `scripts/CheckT3BCYEFT1EFT2LightXXOverlap.experimental.py` (complete Python script)
2. `docs/WorkstreamC_YEFT2LightXXOverlap_README.md` (this file)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline

Expand-Archive `
  -LiteralPath "$HOME\Downloads\T3_WorkstreamC_YEFT2LightXXOverlap_v6.zip" `
  -DestinationPath . -Force

python scripts\CheckT3BCYEFT1EFT2LightXXOverlap.experimental.py --self-test

python scripts\CheckT3BCYEFT1EFT2LightXXOverlap.experimental.py `
  --output output\audits\t3_bc_y_eft2_light_xx_overlap_v1\summary.json
```

Output: `output/audits/t3_bc_y_eft2_light_xx_overlap_v1/summary.json`. It includes complete individual `XX` and EFT1/EFT2 counterterm difference maps, SHA256 provenance, exact source equality and mismatch details, and explicit anti-certification flags. On a B/C mismatch the command exits with code 1 and writes the detailed report whenever possible. No GitHub writes and no production-file changes.

**The actual B/C EFT2 XX equality has not been executed in the assistant environment.** Only synthetic v6 checks passed before packaging; the original local Matchete output and exported CG JSON are needed for the genuine run.

## Remaining after v6

Even a passing v6 would establish **only** the scalar-quartic/`LLSS` beta-function piece of the EFT1/EFT2 counterterm difference. It does **not** establish complete operator mixing (gauge/Yukawa/field terms), the full mixed-XY soft-region subtraction, higher-dimension terms at `MX^2/MY^2 ~ 0.76`, matched scale-dependent couplings, or the ordered `GammaCC` to Weyl flavour-symmetric Weinberg `C5` mapping. Neither a physical Y threshold coefficient nor a neutrino mass prediction follows from v6 alone. The immediate next full-theory task is to derive the **complete EFT1/EFT2 renormalization/matching scheme** consistently, then normalize and propagate the physical coefficient to the X threshold.
