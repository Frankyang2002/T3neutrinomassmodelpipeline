# Workstream A — native triplet coordinate gate

**Status:** `LiveMatcheteTripletRepresentationVerifiedRGBetaCoordinatesUnavailable`

- Native Matchete HH-triplet CG representation: **6/6**, exact SU(2) tests.
- Previous baseline native Matchete CGs: **48/48**.
- E-family cross-representation covariance probes passed: **4/4** (diagnostic).
- RGBeta `delS2` exact numerical components exposed: **False**.
- RGBeta `tGen` exact numerical components exposed: **False**.
- Installed RGBeta is probed directly; **no numerical component assignments injected**.

## Verified/conditional distinction

The Matchete triplet representation, SU(2) Casimir, generator algebra, conjugate
bilinear metric and exact unitary map to the declared symmetric-doublet reference
are established from the **live registered** A/D topology Clebsch–Gordan tensors.
This fixes a Matchete-side reference convention, **not an RGBeta physical basis**.
Literal integer indices are only a diagnostic probe of RGBeta. A failure to
resolve numeric components does not invalidate its abstract group algebra.
Even if the diagnostic exports arrays, the candidate is not a certified
cross-package full Lagrangian identity, signed quartic/Yukawa dictionary,
or Matchete matching. The source Stage23 maps remain conditional.

| Model | Native CG norm squared | Metric scalar identity | SU2 exact |
|---|---|---|---|
| A_m4 | `2*sqrt(3)/3` | True | PASS |
| A_m2 | `2*sqrt(3)/3` | True | PASS |
| A_p0 | `2*sqrt(3)/3` | True | PASS |
| D_m2 | `2*sqrt(3)/3` | True | PASS |
| D_p0 | `2*sqrt(3)/3` | True | PASS |
| D_p2 | `2*sqrt(3)/3` | True | PASS |

## E-family cross-representation diagnostic

This separately checks HH S1 S2dag against the A-derived S2dag and
D-derived S1 triplet generators. Any nonzero residual is retained
rather than fixed with an assumed phase or basis rotation.

| Model | Exact generator covariance | Nonzero residual counts |
|---|---|---|
| E_m4 | True | `[0, 0, 0]` |
| E_m2 | True | `[0, 0, 0]` |
| E_p0 | True | `[0, 0, 0]` |
| E_p2 | True | `[0, 0, 0]` |

No production files, RGE beta functions or EFT outputs were modified.
