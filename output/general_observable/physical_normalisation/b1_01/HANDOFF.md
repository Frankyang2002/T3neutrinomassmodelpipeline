# B1 handoff: Matchete → Weyl, first restricted batch

Reference GitHub `main`: `0b210f18f6a2cfb10fff03059e2ec4a73d2195a4` (read-only, verified 2026-11-03). No local Windows revision observed. Existing source files were not changed.

## Results actually established

1. The source action includes `y1 Bar[L] CConj[F] S1`, `y2 Bar[L] F S2†` and Hermitian conjugates; actual `FreeLag` is assembled separately. `T3Y1CG`, `T3Y2CG`, and `T3MixCG` are separately registered by the source. The native result **has not been executed in this environment**.
2. Explicit choice `C=+i gamma²gamma⁰`, `epsilon12=+1`, `F=(chi,bar xi)` yields `BarL CConj F=(l chi)†`, `BarL F=(l xi)†`, `BarF F=chi xi+h.c.` in exact exterior algebra.
3. In `core.py`'s `-nu psi(Y phi+Z phi*)` convention, `Z_y1=-conj(y1*C1)` and `Y_y2=-conj(y2*C2)`, with Dirac `psi=(xi,chi)`. Existing builder uses simultaneous positive signs; **both vertices negated together are a common neutrino field rephasing**, not a physical discrepancy.
4. Majorana free action uses the identical-Weyl `1/2` exactly once. Triplet Majorana metric in weight basis is symmetric, unitary `G=U.T@U`; its nonzero entries pair opposite weights. `I3` is incorrect in a standard weight basis.
5. Canonical real-scalar two-component self-energy contains `1/(16 pi²)` with no additional global Majorana half/two; for `phi=(R+iI)/sqrt2` the single-scalar Yukawa normalization is `1/sqrt2`.
6. **13 local mathematical/synthetic tests passed** on the generated files. No Matchete/RGBeta or all-16 benchmark run was possible here.

## Inputs safe to integrate?

**No physical integration approval.** The expressions may be used as a **conditional candidate** in an explicitly chosen basis. Neither the source's concrete `GammaCC` phase nor the Matchete/RGBeta SU(2) gauge-index map (B2) has been measured. The conditional reconstruction of the existing provisional candidate cannot be upgraded to `physical_source_certified`.

## Missing evidence

* User's latest local Git SHA and actual on-disk source files (audit runner fails closed on drift from reference main); repository's current `Reports/output/full/multiplicity/` native tensor audit and 16 benchmark JSONs.
* Fresh native action files for `A_p0`, `B_m1`, `C_m1` from `scripts/RunPhysicalB1MatcheteAction.wl`; check every original/H.c. and heavy mass expression, not just status.
* B2's concrete physical cross-package SU(2) tensor/metric/epsilon/coupling dictionary and relative phases.

## Output use

`b1_01/` is this chat's immutable first-batch mathematical handoff; `all16_component_check.json` explicitly records sixteen unexecuted cases with `null` residuals. On the user's local repo, the independent audit defaults to fresh `b1_02/` and computes new residuals without overwriting this report. All new code paths are listed in the delivery instructions. No shared production module, benchmark, RGE, or EFT file is modified.
