# Track C — B1 independent B/C charged-Dirac UV calculation

**Fixed scope:** `B_-3`, `B_+1`, `C_-3`, `C_+1` only. Do **not** recalculate `B_-1` or `C_-1`. One batch, distinct filenames, no writes to production physics, other tracks, frozen data, RGE, EFT or GitHub. This roadmap is immutable after initial creation.

**Source inspected before roadmap creation (read-only GitHub main, 2026-11-03):** `Lagrangian/model/{T3ModelCatalog,T3Fields}.wl`, `Lagrangian/interactions/{T3Topology,ScalarPotential}.wl`, existing `physical_normalisation` B1 / BM1 / CM1 derivations and evidence; four each of `output/general_observable/lo_direct_one_loop/source_{cg,portal_cg}_all16/...` registries; `Numerical/benchmarks/fitted_reference/T3_dS1_2_dS2_2_dF_{1,3}_alpha_{m3,p1}.json`. Inspected Git tree SHA `d2163e6bb95109c0fac419f7e41f335f0d347b62`; not a claim about the user's local working tree. The local Windows repository and Matchete kernel are not mounted in this environment.

## C0 — Provenance and four-model inventory

Verify representations, hypercharges, benchmark inputs, current source action, and available native tensors. Record per-file Git blob SHA (or local hash) and distinguish archived *native-produced* registry from a newly executed native action. `B: (dS1,dS2,dF)=(2,2,1)`; `C: (2,2,3)`; `alpha=-3` gives `(YS1,YS2,YF)=(-3/2,-1/2,-1)` and `alpha=+1` gives `(+1/2,+3/2,+1)`. F is charged, complex and Dirac in all four. Benchmark 3-flavour diagonal F masses are 100, 110, 120 TeV; native default `MF` is one common mass unless separately enabled and verified.

**Gate:** evidence for *each* of the four, no undocumented assumptions. Failure must be reported separately for that model.

## C1 — Native scalar/Yukawa/fermion contractions

Establish charge-selected CG components and all contributing internal fermion/scalar channels directly from actual oriented Matchete tensor axes. Check **original** `y1: Bar[L] CConj[F] S1` and `y2: Bar[L] F Bar[S2]`, both H.c.s, native scalar-mixing invariant `H H S1 Bar[S2]`, and `LFree` Dirac mass coefficient. Distinguish triplet **Cartesian real SU(2) representation** from a self-conjugate fermion: since YF=±1, Dirac mass is Hermitian pairing `Bar[F]F`; no Majorana `1/2` and no weight-space anti-diagonal Majorana mass metric are permitted. Verify full four-index Inv1/Inv2 portal identity and all charge-sector scalar quadratic blocks, not just a neutral 2x2 truncation; enforce eigenvalue positivity. Source-native tensor registries may be reused read-only as archived evidence but **fresh installed Matchete action extraction is a distinct pending gate until run**.

**Gate:** derive native charge-selected vertices, action terms, scalar mass matrices and fermion normalisation without reference to frozen matrices.

## C2 — Independent complete four-model loop calculation

Use native CG tensors and benchmark matrices in a single generalised Python evaluator, with a fresh numerical Feynman-parameter integration of the finite zero-momentum mixed-scalar propagator. Do not import `Numerical/direct_one_loop`, provisional input builders, existing frozen masses, or frozen C5. Include all allowed internal electric charge blocks and full scalar mixing. Derive coefficients from action/tensors only; give each of the four 3x3 complex symmetric matrices and Takagi masses for **zero portals** and separately **selected identity portals**, with common-MF versus three-diagonal-MF explicitly labelled.

**Gate:** independent numerical matrix for each model, source provenance and scoped conditionals retained.

## C3 — Numerical and physical cross-checks

Verify scalar mixing eigenvalues/positivity for all charge sectors, cancellation of the B0 ultraviolet pole in the mixed finite difference, numerical integration against its analytic finite function, matrix symmetry, Takagi diagonalisation, phase-rephasing covariance on genuinely complex Yukawa/portal data, and negative controls (wrong CG/charge, omitted mixing, wrong Yukawa conjugation, etc). Only **after** the independent numerical matrix exists may the frozen reference be opened to measure per-model relative Frobenius matrix residual. Never use agreement to select a normalisation coefficient. Track A retains certification of the shared absolute `GammaCC` spinor convention; any conditional Weyl prescription is labelled as such.

**Gate:** tests and residuals actually executed; otherwise mark pending, not passed.

## C4 — One final batch package

Bundle all Track-C-only full replacement scripts, the immutable roadmap, source manifests, native-output review and concise mathematical derivation, JSON/CSV four-row results, exact safe Windows PowerShell commands and remaining-blocker summary into **one ZIP** retaining project-relative paths. Write numerical outputs only into newly timestamped `output/general_observable/physical_normalisation/parallel_C_BC_<timestamp>/` directories. Never overwrite an existing file. Report exactly four rows with native evidence, scalar and fermion normalisation, independent mass matrix/Takagi spectrum, frozen comparison residual and blockers. Stop after delivery; do not expand to other classes or EFT.

**Certification rule:** partial evidence or source-only calculations do not count as native executed or as physically certified. An unresolved model may be blocked without stopping independent work on the other three.
