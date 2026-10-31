# Workstream A — Live baseline Matchete CG bridge

**Status:** `LiveMatcheteBaselineCGsPassedRGBetaPhysicalBridgeUncertified`

- Live registered baseline CGs: **48/48**.
- Symmetric Higgs-Higgs topology CGs: **16/16**.
- Exact projective cross-model comparisons: **36**.
- Independent source: installed Matchete `GetCGTensor` (not generated tensor candidates).
- Native RGBeta baseline vertex definitions: checked against `T3RGBetaModel.wl`.

| Model | CGs | Components | HH symmetric |
|---|---:|---:|---|
| A_m4 | 3 | 28 | PASS |
| A_m2 | 3 | 28 | PASS |
| A_p0 | 3 | 28 | PASS |
| B_m3 | 3 | 24 | PASS |
| B_m1 | 3 | 24 | PASS |
| B_p1 | 3 | 24 | PASS |
| C_m3 | 3 | 40 | PASS |
| C_m1 | 3 | 40 | PASS |
| C_p1 | 3 | 40 | PASS |
| D_m2 | 3 | 28 | PASS |
| D_p0 | 3 | 28 | PASS |
| D_p2 | 3 | 28 | PASS |
| E_m4 | 3 | 60 | PASS |
| E_m2 | 3 | 60 | PASS |
| E_p0 | 3 | 60 | PASS |
| E_p2 | 3 | 60 | PASS |

## What is not certified

Matchete registered CGs now have a direct live extraction and exact internal
consistency checks. This is **not** an independent RGBeta component-level
comparison or a physical signed coupling dictionary. The prior RGBeta native
`S2` generator evidence remains symbolic (no numerical triplet matrices),
and the required conjugate-representation/Bar field intertwiners have not
been identified in a common basis. Therefore **cross-package normalization
and full physical UV equivalence remain unresolved**.

The 46 Stage23 / 50-slot quartic map is conditional. This stage neither
calls `Match[]` nor registers an exact discrete `Z2` symmetry in Matchete.
It leaves Workstreams B, C and D and existing output files unchanged.
