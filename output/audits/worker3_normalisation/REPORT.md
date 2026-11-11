# Worker 3 convention and normalisation audit

**Result:** Common **real** Wilson scaling protects neutrino mixing angles, mass ratios, and appropriately defined stage-running ratios **if the full coupled Wilson system and all source/matching terms scale together**. Masses scale as `|r|`, squared differences as `r²`. A source-certified common Matchete-to-Weyl factor is **not established**.

## Evidence ledger

Read-only GitHub repository `Frankyang2002/T3neutrinomassmodelpipeline`, `main`, inspected 2026-11-12. Exact blobs of eight files are recorded in `source_manifest.json`. The observed implementation specifies a canonical symmetric master Wilson tensor and its Weinberg projector; it explicitly blocks legacy Matchete tensors without nine certificates. The existing B1 chosen-gamma-basis work is **conditional** on a native Matchete gamma/epsilon and SU(2) dictionary. Tree Matchete LLSS in the fermion-first intermediate EFT is distinct from any direct LLHH Weinberg source. See `convention_table.csv` and `scripts/worker3_normalisation/README.md`.

Source references:

- `RGE/running/intermediate/GeneralEFT1C5Mixing.py`: module convention header; `CERTIFICATES`; `_read_wilson`; `_weinberg_template`; `project_c5`; `calculate_master_beta` with conjugate coefficient term.
- `RGE/general/WilsonTensorRGE.py`: seven master beta structures and symmetry permutations.
- `Numerical/running/GeneralT3WilsonNormalizationAudit.py`: existing native PL/PR and polynomial diagnostics, explicit spinor normalisation exclusion.
- `Numerical/lo_certification/physical_normalisation/B1_MATCHETE_WEYL_DERIVATION.md`: conditional 4-spinor/Weyl derivation, SU(2) triplet mass metric, outstanding native equivalence.
- `Lagrangian/RunGeneralT3FermionFirstEFT1Tree.wl`: tree LLSS, X/Y heavy scalar retained, physical C5 false.
- `Lagrangian/RunT3EFT1WilsonHierarchyNative.wl`: independent W[p,q] EFT1 coefficients and finite source status.

For the master polynomial `C_ijab psi_i psi_j phi_a phi_b+h.c.` and symmetric Weinberg polynomial `(1/2) C5_pq (L_p.H)(L_q.H)+h.c.`, the Higgs neutral real component yields `C_nunu,hR,hR=C5/4`. This **projector identity is internal to explicitly declared conventions** and cannot determine the factor between an arbitrary Matchete ordered LLSS coefficient and canonical `C_ijab`.

## Numerical results produced here

- 32 complex-symmetric nondegenerate matrix cases: passed common real scaling and isolated complex global mass-matrix scaling, mass-squared differences, ratios, PMNS modulus/angles and synthetic stage ratios (`regression_cases.json`).
- 16 independently solved two-stage coupled-IVP cases, including Wilson-conjugate terms, coupled cross-channel evolution and a matching injection: passed common real equivariance, mass-ratio/angle and stage-ratio invariance; nonreal scaling does not commute with conjugate beta (`coupled_ivp_regression.json`).
- 11 unit tests: 10 passed, 1 production-import projector test skipped in this isolated environment (local source absent). An additional coupled-IVP unit test is included for local execution.

**Important:** No production T3 coupled-IVP solve, exact native Matchete action, or 16-model amplitude comparison was run. These numerical tests are diagnostics and proofs of scaling *conditional on homogeneity*, not a new physical source certificate. Output files do not supersede your existing benchmark calculations.

## Unresolved absolute conversion

There is no justified global number from the inspected source. For each independent `LLSS_XX`, `LLSS_XY`, `LLSS_YY`, chirality-conjugate and direct `LLHH` channel, certify the native 4-spinor charge-conjugation phase, field/CG basis, identical-index orbit multiplicity and functional derivative amplitudes. Repeat with nontrivial complex flavour tensors across representation classes A–E; only then test whether conversion ratios coincide. Without this, normalization may be channel-, representation- or flavour-map-dependent. Treat any rough numerical fit as internally defined, not physically source certified.
