# B2 — independent gauge/coupling track, initial evidence handoff

- **Source inspected:** public GitHub `main` at `0b210f18f6a2cfb10fff03059e2ec4a73d2195a4` (3 November 2026, Australia/Sydney). The user's local Windows revision was **not observable**. Source Git blob hashes are listed in `cross_package_dictionary.json` and checked by the local runner.
- **Shipped report:** exact charge skeletons for 16/16 and conditional Weyl action rephasing; **no local native CG JSON or Wolfram execution** was available here. The `conditional_native_CG_channels` are deliberately null. The script will compute them against source-guarded local audits if available.
- **Basis:** U1=(1), U2=identity, U3 Cartesian-to-weight in derivation; doublet epsilon antisymmetric, triplet mass metric U3.T U3. RGBeta eps/del/delS2 numeric axes remain unknown.
- **Operator orientation:** Matchete y1 uses `CConj[F]`/S1, y2 F/Bar[S2], mixing H H S1 Bar[S2] + HC. RGBeta y1 Bar[S1],l,F; y2 S2,l,Fc (Dirac) or F (Majorana); mass F,Fc epsilon for dF2 and F,F for Majorana.
- **Missing action factors:** physical RGBeta real-scalar/Yukawa/MF/quartic coefficients, numeric native SU2 tensor axes, B1 GammaCC spinor phase `eta`. Source-to-source signed conversion cannot be proved from registry metadata.
- **Candidate spurion:** `y1*y2*conjugate(MF)*conjugate(lambdaT3)` conditional on B1 Weyl monomials; wrong lambda HC/same F/Fc fail negative tests.
- **Portals:** previous native report for 10/10 has Hermitian null sectors A/S2, D/S1, E/S1, E/S2. The runner imports **all exact nullspace vectors** only when prior proof JSON and all source SHA256 match; no unique map claimed for those four sectors.
- **Tests here:** `12 passed` in Python; no Wolfram native tests. Includes source drift and false certificate controls.
- **Safe to integrate:** analytic charge skeletons and explicitly conditional spurion. **Not safe:** absolute physical `y1,y2,MF,lambdaT3` ratios, relative phases or neutrino mass normalisation.
- **Flags:** `physical_crosspackage_complete=false`; `physical_source_certified=false`; `native_physical_source_certified_models=0/16`.
- **Next evidence:** execute local Python audit and native RGBeta formal definition export, then obtain explicit source-aligned numeric tensor/action evidence and B1 spinor convention. Do not use CG diagnostic products as physical loop factors.
