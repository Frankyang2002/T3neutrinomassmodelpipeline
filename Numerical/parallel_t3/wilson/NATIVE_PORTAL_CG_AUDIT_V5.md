# Workstream A v5: corrected native Higgs-portal CG audit

## Scope and sources

This batch adds **six full new files** alongside the v1–v4 A modules. It does
not modify `ScalarPotential.wl`, `T3RGBetaModel.wl`, RGBeta's installation,
existing outputs, any B/C-owned file, or GitHub.

Inspected public repository HEAD: `a944e0030b8ef24b23342d99e0d8f33cef3cffd5`.
Source contracts:

* `Lagrangian/interactions/ScalarPotential.wl` defines real
  `lambdaH{1,2}Inv{k}` and `1/2 PlusHc[phase coupling[] ... CG]` with
  `T3HiggsPortalHermitianPhase`. For an anti-Hermitian CG, `phase = I` **inside**
  `PlusHc`; otherwise the real-coupling portal vanishes.
* `Lagrangian/model/LagrangianBuilder.wl` registers CGs and assembled
  interactions. `Matchete`'s public `GetCGTensor` reads the registered CG.
* `RGE/running/rgbeta/T3RGBetaModel.wl` registers source baseline
  `lambdaH1`, `lambdaH2`, and, for non-singlet BSM scalars, independent
  `lambdaH1Adj`, `lambdaH2Adj`. These are *symbolic* `tGen` contractions.
* RGBeta's `AddQuartic` uses a factor 6 for four complex fields before
  the `Lam/4!` action; `GroupInvariant` and permutation data are retained.
  `scripts/RGBetaNativeAccess.wl` handles package-context keys correctly.

Do not reuse the older five-class restricted native CG report as a replacement
for a fresh export. Its Matchete source SHA predated the anti-Hermitian-portal
repair, so treating its `1/2 PlusHc` as unphased would artificially reduce the
real operator rank.

## Native exports and rank logic

1. `RunT3ParallelAExportPortalMatchete.wl` builds and validates **all 16**
   native Matchete Lagrangians; extracts each registered Higgs-portal CG,
   exact components and source-defined phase; verifies the `CG` and coupling
   occur in `LInt`. It fails if `A_-4/lambdaH2Inv1` is not `I` phased.
2. `RunT3ParallelAExportPortalRGBeta.wl` builds native **EFT1** RGBeta
   (SM + S1 + S2) for all 16; records actual registered scalar fields,
   symbolic singlet/adjoint invariants, permutation counts, real-coupling
   metadata and the quartic prefactor. It does **not** pretend abstract
   `tGen[SU2L[S2],A,i,j]` is a numerical component matrix.
3. `RunT3ParallelAPortalAudit.py` consumes both exports, validates SHA256
   against the current local source bytes, and constructs the complete
   rank-4 `Hbar H Sbar S` tensor with singleton axes for singlets.

For each native Matchete CG tensor \(T_k\), let

\[
 (T^{\rm HC})_{ijab}=T^*_{jiba}, \qquad
 H_k=\tfrac12(p_k T_k+p_k^*T_k^{\rm HC}),\quad p_k\in\{1,i\}.
\]

This is exactly the source `1/2 PlusHc` for **real** `lambdaH{which}Inv{k}`.
The audit requires each native CG to be purely Hermitian or anti-Hermitian,
checks `phase` against its actual numerical tensor, and demands full real
Hermitian rank: one portal direction for a scalar singlet, two for a doublet
or triplet. Missing, rank-deficient or stale source evidence fails closed.

The singlet/identity tensor is

\[
 I_{ijab}=\delta_{ij}\delta_{ab},\qquad
 c_{0k}=\frac{\langle I,H_k\rangle}{\langle I,I\rangle}.
\]

Subtracting the singlet part leaves a traceless adjoint tensor. Its two partial
traces must vanish. For non-singlet scalar representation \(R\), the source
RGBeta generator normalization gives

\[
 \big\|\sum_{A=1}^3 t_H^A\otimes t_R^A\big\|_F^2
 =3T(\boldsymbol{2})T(R),\quad
 T(\boldsymbol2)=\tfrac12,\quad
 T(\boldsymbol3)=2.
\]

The audit chooses the orientation of the first *nonzero native traceless*
portal CG and rescales it to that generator norm. It then produces a **real
conditional matrix** mapping each Matchete coupling to `lambdaH{which}` and a
`lambdaH{which}Adj_MatcheteOriented` parameter. It verifies every component,
plus independent random complex-field potentials, not merely matching one
vacuum field direction. It includes all 16 models and 58 portal coupling slots
(32 singlet and 26 additional adjoint directions).

**Crucial limitation:** the first-native-CG adjoint orientation need not equal
the sign/orientation chosen by native RGBeta's abstract `tGen` generators.
The output labels this sign `UNDETERMINED_PLUS_OR_MINUS`; it is *not* a
physical `lambdaHAdj` identity and must not feed the Wilson background.
Numerical RGBeta generator components or an independently verified common
SU(2) intertwiner are still required for that certification.

## Installation (PowerShell, repository root)

Install the ZIP **into**:

`C:\thesistemp\T3neutrinomassmodelpipeline`

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
$zip = "C:\path\to\t3_workstream_A_native_portal_audit_v5.zip"
Expand-Archive -Path $zip -DestinationPath . -Force

.\.venv\Scripts\python.exe -m unittest discover `
  -s tests -p "test_t3_parallel_a*.py" -v
```

Run both native exporters. Use **fresh** output paths; they refuse overwrites:

```powershell
wolframscript -file scripts\RunT3ParallelAExportPortalMatchete.wl `
  output\parallel_t3\wilson\matchete_portals_native_v5.json

wolframscript -file scripts\RunT3ParallelAExportPortalRGBeta.wl `
  output\parallel_t3\wilson\rgbeta_portals_native_v5.json
```

Then audit:

```powershell
.\.venv\Scripts\python.exe -m scripts.RunT3ParallelAPortalAudit `
  --matchete output\parallel_t3\wilson\matchete_portals_native_v5.json `
  --rgbeta output\parallel_t3\wilson\rgbeta_portals_native_v5.json `
  --project-root . `
  --output output\parallel_t3\wilson\portal_dictionary_conditional_v5.json
```

Expected **only if the live exporters succeed**: `58/58` portal slots,
`16/16` models and status
`NativeMatchetePortalPhasesAndFullRank_Verified_RGBetaAdjointOrientationUncertified`.
The Wolfram exporters have not been executed in this container (WolframScript
is unavailable). Existing public five-class tensors were used to construct
isolated test fixtures, NOT presented as a new all-16 native run.

Do not run a physically normalized EFT/Wilson evolution or identify `Adj`
signs from this v5 output. The 50 extra-quartic v3 dictionary and v4 real-
action audits remain separate, valid **conditional** evidence. Baseline scalar
self/cross invariant CG normalization, scalar masses, Yukawas and spinor
GammaCC-to-Weyl remain distinct open work.
