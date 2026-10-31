# Workstream A — UV closure stage

Outcome: **StructuralAndNumericalStagePassed**. This is not a full UV certificate.

| Model | Extra quartics | Extra bilinear | Scalar structure | Actual Z2 terms | Eigenbasis |
|---|---:|---:|---|---|---|
| A_m4 | 0 | 0 | PASS | PASS | PASS |
| A_m2 | 5 | 1 | PASS | PASS | PASS |
| A_p0 | 5 | 1 | PASS | PASS | PASS |
| B_m3 | 1 | 0 | PASS | PASS | PASS |
| B_m1 | 7 | 1 | PASS | PASS | PASS |
| B_p1 | 1 | 0 | PASS | PASS | PASS |
| C_m3 | 1 | 0 | PASS | PASS | PASS |
| C_m1 | 7 | 1 | PASS | PASS | PASS |
| C_p1 | 1 | 0 | PASS | PASS | PASS |
| D_m2 | 5 | 1 | PASS | PASS | PASS |
| D_p0 | 5 | 1 | PASS | PASS | PASS |
| D_p2 | 0 | 0 | PASS | PASS | PASS |
| E_m4 | 0 | 0 | PASS | PASS | PASS |
| E_m2 | 6 | 1 | PASS | PASS | PASS |
| E_p0 | 6 | 1 | PASS | PASS | PASS |
| E_p2 | 0 | 0 | PASS | PASS | PASS |

## Scope and boundaries

Independent scalar counting uses Bose-symmetrized SU(2) Sym^n weight multiplicities, U(1)Y, and the assumed exact Z2.
Actual-terms parity uses Matchete-built BSM interactions and free terms with extra bilinears disabled, so B/C can be built in their original basis.
Eight native mass eigenbases are certified only by the earlier user-generated Matchete reports; this stage adds numerical potential and arbitrary-tensor transport checks.
The six self-scalar Matchete eigenmass spurions remain independent symbolic masses; the numerical relation is conditional on the source potential sign and native factor.

## Not certified

- Symbolic original-to-eigenmass and phase parameter relation in native Matchete
- Exact independent full interaction CG normalization and Matchete/RGBeta map
- Explicitly registered Z2 selection rule, not just parity of assembled terms
- Independent exhaustive Lorentz/fermion operator completeness and flavour tensors
- Physical EFT matching / threshold normalization (outside Workstream A)
