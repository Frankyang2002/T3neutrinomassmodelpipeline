# Batch 3 — UV vs conditional EFT matrix comparison

**DIAGNOSTIC ONLY — NOT PHYSICAL MATCHING.**

Source mode: `local`. The two calculations use different source scales/coupling states and unmatched operator content.

| Model | UV masses (eV) | EFT masses (eV / abs N) | UV off-diagonal fraction | Best-scale residual |
|---|---|---|---:|---:|
| `A_-4` | 0.0436006, 0.00758273, 0.000474483 | 0.0259021, 0.0215851, 0.0194265 | 0.803112 | 0.908928 |
| `B_-1` | 0.0504413, 0.00980338, 0.00452249 | 0.0236911, 0.0197426, 0.0177682 | 0.794740 | 0.886970 |

The residual is `||M_UV-z M_EFT||_F / ||M_UV||_F`, where `z=vdot(M_EFT,M_UV)/vdot(M_EFT,M_EFT)` and both matrices are expressed in eV. It is evaluated in the **reported** lepton-flavour axes; there is no basis optimisation.

The inspected source-resolved F-tree LLSS tensors contain only lepton flavour pairs `p=q` for both models, while the reported direct-UV matrices have nonzero off-diagonal flavour entries. This is an **unresolved source/basis/input alignment diagnostic**; it is not proof the tensors should agree before common-scale matching.

An overall spinor phase or factor alone cannot align the matrices in these axes. This is NOT evidence of a physical matching failure until source couplings, scales, flavour basis, CG normalisations, and missing hard/soft matching are reconciled.

Both heavy scalar multiplets are integrated together at `muS`. Both direct-UV and EFT masses are conditional; no physical `GammaCC`, `C5` or `mnu` is assigned.
