# A4 — Pass/fail evidence report (Track A, 2026-11-03)

**Overall: NOT PHYSICALLY CERTIFIED — 0/16.** The mathematical convention tests pass, but installed-native Matchete-to-Weyl relative phase is unresolved. This package does not make production changes, use frozen UV mass matrices, or count already completed native tests as new results.

## Actual public-source and installed-native evidence

- Latest inspected GitHub `main` commit: `d2163e6bb95109c0fac419f7e41f335f0d347b62`, read-only. Source file blob SHA identifiers are in `EVIDENCE_PROVENANCE.json`.
- Matchete `T3Topology.wl`: source `+y1 Bar[L] CConj[F] S1`, `+y2 Bar[L] F Bar[S2]`, `+lambdaT3 HH S1 Bar[S2]`, each with H.c. `PlusHc`. Restricted flags exclude extended interactions.
- Previously obtained **installed** action witness, inspected through current repository JSON: `y1` and `y2` original/H.c. factors `+1`, Dirac mass `-1`, Majorana mass `-1/2`, `Proj[1]` originals, and native `GammaCC` in y1/H.c. and Majorana free mass. It reports `explicit_numeric_4x4_exposed=false`, `installed_absolute_phase_established=false`.
- Previously obtained **installed** CG witnesses provide neutral components `(C1,C2,Cm)` for A/B/C. Reused as **audited constants** only for this execution; not a newly rerun Matchete source extraction. Source and witness SHA256 are verified only by the `--mode native` code on the user machine.
- `FeynmanRules` was exposed by the earlier installed native probe but **not called** there. Matchete v0.6.0 publicly documents that API; the new `native_rule_probe.wl` calls it and uses the wrong-H.c. operator control. It was **not executed here** because `wolframscript` is absent. No native amplitude/operator coefficient residual is reported.

## Declared-Weyl independent complex numerical experiment

All three models were evaluated once in a single substantive batch with deterministic independent complex `3x3` `y1,y2` draws (seeds 704, 805, 906), `MF=(100,111,123) TeV`, `lambdaT3=0.13 exp(0.29i)`, `v=246.22 GeV`, and deliberately chosen **zero portals** `a=1306122.448979592 GeV²`, `b=1000000 GeV²`. Not fitted values and not the frozen UV amplitude. Heavy masses are real and positive in the base point; controlled complex-M mass rephasings test the conjugate numerator. Eigenvalues of the complex scalar squared-mass matrix are used, with `A21=-lambda v² Cm/2`.

| Model | symmetry residual | Feynman parameter vs analytic `B0` | full field rephase | source global-phase test | wrong `y1` conjugation difference | source certificate |
|---|---:|---:|---:|---:|---:|---|
| `A_+0` | `0.000e+00` | `0.000e+00` | `3.106e-16` | `2.969e-16` | `1.273e+00` | **conditional only** |
| `B_-1` | `0.000e+00` | `0.000e+00` | `1.209e-16` | `2.006e-16` | `1.322e+00` | **conditional only** |
| `C_-1` | `0.000e+00` | `0.000e+00` | `1.460e-16` | `1.836e-16` | `1.408e+00` | **conditional only** |

All displayed residuals are relative Frobenius norms. Numerical `0.000e+00` means the independent high-precision calculations agree to **float64 output precision**, not an exact zero identity for the loop.

Additional group-algebra results: doublet epsilon invariance residual `2.665e-17`; triplet `G=U^T U` residual `1.822e-16`; triplet charge-intertwiner residual `6.328e-17`; wrong triplet identity charge violation `2.828e+00`.

### Numerical masses, CP and negative controls

| Model | singular masses (eV, descending) | CP-odd invariant | wrong lambda conjugate | wrong complex-M conjugate | extra 1/2 | lambda=0 |
|---|---|---:|---:|---:|---:|---:|
| `A_+0` | 0.25353716, 0.09342800, 0.02293859 | `-7.651368e-02` | `0.5719` | `0.7563` | `0.5000` | `0.0` |
| `B_-1` | 0.18911527, 0.11124788, 0.03610871 | `1.356376e-01` | `0.5719` | `0.8772` | `0.5000` | `0.0` |
| `C_-1` | 0.15720697, 0.07192280, 0.02407056 | `5.731506e-02` | `0.5719` | `0.4804` | `0.5000` | `0.0` |

Negative controls are **changed relative matrices under deliberately inconsistent choices**, not new installed-native measurements. These differences do not automatically imply physical observability of an overall matrix sign. The `mnu→-mnu` test leaves singular masses unchanged.

## Pass/fail gates

| Gate | Result | Evidence / obstruction |
|---|---|---|
| A0 read-only public source audit | PASS | GitHub commit and blob SHA in provenance |
| A0 Windows working-tree/source-hash equality | NOT TESTED | Local `C:\...` repository unavailable in this container |
| A1 declared gamma/C and exterior Grassmann identities | PASS | `evidence/exact_weyl_algebra.json` |
| A1 physical native-to-Weyl relative sign | **BLOCKED** | Installed numerical `GammaCC` orientation unmeasured |
| A2 installed rule-level FeynmanRules + wrong-H.c. control | **NOT RUN HERE** | Standalone executable supplied; requires user Matchete runtime |
| A3 symmetric matrix, rephasings, independent 1-loop | PASS (CONDITIONAL) | `evidence/audited_constants_three_models.json` |
| A3 source-native witness-gated complex mode | **NOT RUN HERE** | Requires local native JSONs/source hashes and successful A2 witness |
| Physical source certification | **FAIL CLOSED: 0/16** | Cannot infer native relative spinor/CG phase from conditional covariance |

## Precise unresolved native gate

To upgrade any T3 source certificate, extract a **relative**, convention-sensitive `y1`–`y2`–mass insertion amplitude or operator coefficient from the **installed** Matchete rules in a **common numerical Weyl basis** with the actual source SU(2) epsilon, triplet mass metric and scalar phase rotations. An isolated global `C` sign or full neutrino mass sign does **not** fix a physically measurable quantity. `native_rule_probe.wl` records Feynman rules and intentionally detects the unphysical wrong-H.c. action; it does not pretend a raw symbolic `GammaCC` has an installed 4×4 value.

**Track B/C handoff:** Use the declaration in `DERIVATION.md` only as a conditional reference; the relative Matchete↔Weyl charge-conjugation and numerical SU(2) orientation factors remain unresolved. No physical certification or production integration is authorised.
