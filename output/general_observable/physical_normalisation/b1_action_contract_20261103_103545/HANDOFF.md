# B1 next-batch: native source coefficient and chirality witness

- Local Git revision: `e2a5def5326cc452f668ef593bd81ae13fdff30b`. Source byte equivalence checked against pinned Git blobs.
- Native symbolic actions A_+0, B_-1, C_-1: **3/3 structurally checked** (original and HC for both Yukawas, spinor/projector orientation, mass and symmetry prefactors).
- Extracted mass coefficient: Dirac `-1`, Majorana `-1/2` **once**; Yukawa originals and HC `+1` for each.
- Prior conditional 16/16 maximum relative residuals: mass `2.182e-16`, Y `1.664e-16`, Z `0`.
- A native GammaCC probe, if supplied, verifies symbolic Matchete identities only; an overall numerical sign relative to the declared Weyl representation remains **unmeasured**.
- Relative CG phases/normalisations from B2 remain **unmeasured** in this B1 audit. No source-certified neutrino mass coefficient.
- Common neutrino sign reversal of both Y and Z is a field rephasing; not a physical factor-of-two discrepancy.
- Physical certifications: **0/16**, even after all checks pass.
- Safe to integrate into production: **NO**. Native mass/Yukawa structural coefficients and chosen-basis identities may be used as source-grounded inputs to the later joint physical-dictionary step.
- No production files, frozen inputs, GitHub, EFT, or RGE changed.
