# Stage 26 — Source sign and conditional quartic action identities

## Scope and result

This stage is **local/read-only** and does **not** change RGBeta production definitions or promote physical-coupling certification. It uses the user's **fresh post-Stage24** Stage22/23 evidence and successful **live Stage25** native assembly report, then checks the quartic *polynomials*, including repeated identical scalar fields and Hermitian-conjugate contributions, in the **provisional common component labels** used by Stage22.

The result is stronger than matching the ordered invariant tensors alone. It is still **not** an independently established global Matchete/RGBeta Cartesian real-field intertwiner; the six `CalibratedTripletCGCandidate` cases require that separate verification.

## Relative scalar-action sign: published and project source

RGBeta's formal convention is explicit in A. E. Thomsen, *Introducing RGBeta*, EPJC **81** (2021) 408 (DOI [10.1140/epjc/s10052-021-09142-4](https://doi.org/10.1140/epjc/s10052-021-09142-4)):

- Eq. (2.3): the model-builder Lagrangian contains `- Σ lambda (phi phi phi phi)`.
- Eq. (2.9): the complex field is `phi = Delta.Phi/sqrt(2)`.
- Eq. (2.13): the quartic piece in the generic real-field Lagrangian has `- 1/4 (lambda Delta^4 Phi^4)`.
- Eq. (2.15): the real symmetric tensor uses `1/4` times its 24-permutation sum.

RGBeta's `AddQuartic` native implementation stores, for four complex scalar legs,

\[
  \Lambda\text{-normalisation} = 4!\,(1/\sqrt2)^4 = 6.
\]

Combining the generic `1/4!` quartic action, the native factor and the complex-field projection gives

\[
  \frac{6}{4!}(\sqrt{2})^4 = \frac14\times4 = 1.
\]

Thus the **RGBeta** contribution in model-builder complex fields is **`-lambda * invariant * fields`**, not `-6*lambda` and not `-lambda/4`. The conjugate interaction is present as a separate term for `SelfConjugate -> False`.

In this project's `Lagrangian/model/LagrangianBuilder.wl`, the BSM interactions enter the action through `LInt = Total[Lookup[valid, "Expression", {}]]`. `scripts/ExportGeneralT3Stage22Matchete.wl` records these new quartics as `+lambda[] times contracted fields in PlusHc[contribution]` with prefactor 1. Stage26 verifies these exact source strings locally, and **fails closed** if they change. This is a **source-level relative-sign argument**, not a new live Matchete-vs-RGBeta action evaluation.

Consequently, under the **same physical field-component basis** and identical holomorphic orientation, a Stage22 tensor block `R=M A` requires

\[
  r = -A^{-1}m,
\]

where `m` and `r` are Matchete and RGBeta coupling column vectors. Do not insert this map in production until the global component basis is established.

## Stage26 implementation

`python -m scripts.Stage26QuarticActionIdentity`:

1. Recomputes Stage23 and Stage25 diagnostics from their native inputs, compares them with stored JSON, and fails on stale evidence. In particular it needs the **live Stage25** RGBeta assembly file and its 50-slot/46-block result.
2. Checks the local Matchete builder and Stage22 exporter source conventions.
3. Constructs an **exact commuting polynomial** from every Matchete component tensor and an **independently reconstructed** RGBeta invariant tensor in the Stage22 provisional component basis. All repeated-field monomials are added, retaining exact coefficients (`sympy`, no floats).
4. Explicitly adds Hermitian conjugates using the formal species swaps `H↔Hdag`, `S1↔S1dag`, `S2↔S2dag` and conjugates coefficients.
5. Substitutes the **conditional** signed Stage23 coupling matrices and checks the *entire* polynomial identity `L_M + V_R = 0` exactly, plus negative controls for a wrong sign or a spurious factor `1/4`.
6. Requires exactly **50 slots / 46 blocks** across 16 models and refuses non-post-Stage24 reindexing classes.

**Important limitation:** Formal conjugation in the provisional component labels is not proof that the underlying native triplet components are physically identical or correctly normalised. The source sign is well documented; the physical coupling map is **still UNCERTIFIED**. No JSON flag in this stage claims otherwise.

## Files

- `scripts/Stage26QuarticActionIdentity.py` — complete new diagnostic
- `tests/python/test_general_t3_stage26_quartic_action.py` — regression tests and optional live integration
- `docs/GENERAL_T3_STAGE26_ACTION_IDENTITY.md` — this document

The output report goes to:

`output/general_t3_stage22/stage26_conditional_quartic_action.json`

## PowerShell commands (project root)

```powershell
Expand-Archive -LiteralPath .\T3_General_Scalar_Stage26_ActionIdentity.zip `
    -DestinationPath . -Force

python -m pytest -q tests\python\test_general_t3_stage26_quartic_action.py

python -m scripts.Stage26QuarticActionIdentity
```

Requires successful output from the previous commands:

- `python -m scripts.Stage22RGBetaAlignedGeometry`
- `python -m scripts.Stage23QuarticPrefactorAudit`
- `python -m scripts.Stage25QuarticActionAudit`

**Expected output after the user's successful Stage25 live run:** 16 models, **50 quartic slots**, **46 complete conditional action blocks**, with source-supported sign and **physical coupling mapping UNCERTIFIED**.

## Next gate

Independent **physical component-basis** verification of the six calibrated triplet structures. In particular the calibrated `sqrt(2)/3^(1/4)` tensor coefficient must not be mistaken for a proven, kinetic-term-preserving unitary field redefinition. Only after the global field identification passes can the signed coupling matrices be used in the numerical RGBeta pipeline.
