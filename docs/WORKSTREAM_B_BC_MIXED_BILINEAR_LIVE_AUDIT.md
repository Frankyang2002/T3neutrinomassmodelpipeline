# Workstream B — B(-1)/C(-1) live mixed-bilinear convention gate

## Scope and evidence

The existing `Lagrangian/interactions/GeneralScalarOperators.wl` generates the complex mixed quadratic monomial `S1dag S2dag` using a **Matchete** `PlusHc[contribution]` with no additional explicit `1/2`. The relevant registered coupling is `lambdaT3ExtraM000101I1` and registered CG is `T3ExtraCGM000101I1`.

`RGE/running/rgbeta/RunGeneralT3QuadraticRGBeta.wl` registers `muT3GeneralS1S2` with `{Bar @ S1, Bar @ S2}`, `eps[SU2L @ fund,#1,#2]`, and `SelfConjugate -> False` in the **RGBeta** scalar *potential* convention. In the upstream RGBeta `AddScalarMass`, the internal two-complex-leg normalization is `24*(1/Sqrt[2])^2=12`. The four-leg massive `Lam` uses the registered couplings and includes the conjugate term; **12 is NOT an extra physical mass-potential factor**.

These two independent live probes avoid treating previous status labels as proof:

- `scripts/ProbeT3BCMixedBilinearMatchete.wl` builds only the quadratic candidate in the native B/C model and obtains the **actual registered** CG via `GetCGTensor`. It fails if the current coupling/CG is absent.
- `scripts/ProbeT3BCMixedBilinearRGBeta.wl` builds the scalar EFT1, declares the native mixed bilinear, reads its registered fields/normalization/permutations, and requires its actual installed `RGBeta\`PackageScope\`Lam[...,True]` to contain mass entries.
- `scripts/VerifyT3BCMixedBilinearMap.py` combines these outputs, checks the observed Matchete matrix is proportional to the **declared** epsilon, derives the algebraic map and compares complex polynomials, real Hessians, mixed mass eigenvalues and negative controls. It checks source contracts at run time.

## Conditional convention equation

Let `C_M` be the **live** Matchete registered 2-by-2 tensor and define `epsilon = [[0,1],[-1,0]]`. If the live test gives `C_M = k epsilon`, then for identical labelled complex fields and the same holomorphic monomial orientation:

\[
\mathcal L_M = +[\mu_M^2 (S_1^*)^T C_M S_2^*+\mathrm{h.c.}],\qquad
V_{RGB} = +[\mu_R^2(S_1^*)^T\epsilon S_2^*+\mathrm{h.c.}]
\]

Action equality `L_M = -V_RGB` requires

\[
\boxed{\mu_R^2=-k\mu_M^2},\qquad \boxed{\mu_M^2=-\mu_R^2/k}.
\]

There is **no complex conjugation** under the stated *same monomial orientation* assumption. A different cross-tool barred-doublet intertwiner can change this map and must be established independently; the test intentionally does not set `physical_mapping_certified=true`.

For the RGBeta convention, a mass matrix can be represented in the fundamental-doublet reference basis `Phi=(S1,epsilon conjugate(S2))`, `Phi = U (X,Y)`, with `A=[[mS1Sq,muR],[conjugate(muR),mS2Sq]]` and `U^dag A U=diag`. No scalar-threshold matching is carried out here.

## Installation and exact Windows commands

Extract the ZIP **at the repository root**, preserving these 5 paths:

```
scripts/ProbeT3BCMixedBilinearMatchete.wl
scripts/ProbeT3BCMixedBilinearRGBeta.wl
scripts/VerifyT3BCMixedBilinearMap.py
tests/test_t3_bc_mixed_bilinear_map.py
docs/WORKSTREAM_B_BC_MIXED_BILINEAR_LIVE_AUDIT.md
```

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath 'C:\path\to\Workstream_B_BC_MixedBilinear_Live_Audit.zip' -DestinationPath . -Force
python -X utf8 -m pytest -q tests\test_t3_bc_mixed_bilinear_map.py

wolframscript -file scripts\ProbeT3BCMixedBilinearMatchete.wl output\audits\workstream_b_bc_mixed_matchete_v1.json
wolframscript -file scripts\ProbeT3BCMixedBilinearRGBeta.wl output\audits\workstream_b_bc_mixed_rgbeta_v1.json

python -X utf8 -m scripts.VerifyT3BCMixedBilinearMap `
  --matchete output\audits\workstream_b_bc_mixed_matchete_v1.json `
  --rgbeta output\audits\workstream_b_bc_mixed_rgbeta_v1.json `
  --output output\audits\workstream_b_bc_mixed_bilinear_map_v1.json
```

**Important:** choose new JSON filenames if these outputs already exist. The probes refuse to overwrite prior evidence and the Python verifier fails closed if native extraction, scalar mass-registration or expected package implementation changed. No production files or reports are modified.

## Certification boundaries

A passing result certifies **the captured installed-package source registration**, the native Matchete `CG` components, conditional same-component action equality and its Hessian. It does *not* independently certify the installed RGBeta epsilon component phase, the global Matchete/RGBeta kinetic/generator map, the physical coupling map, Workstream A's final operator basis, or finite sequential matching. In particular, the RGBeta native massive `Lam` test is a registration witness, not an independent fully contracted complex scalar action evaluation.

Next gate: extract the actual RGBeta fundamental representation epsilon component phase and kinetic/generator intertwiner, then perform a **live** cross-tool action-vertex subtraction with the final common field-component map. Only after this gate can the coupling map become a numerical-workstream-certified conversion.
