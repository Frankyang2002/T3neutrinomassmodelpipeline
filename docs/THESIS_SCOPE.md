# Scientific scope of this repository

**Reference:** *Calculating Quantum Corrections in Generalised Scotogenic Radiative Neutrino Mass Models*, Frank Hao Ran Yang, UNSW Honours thesis, November 2026 (109-page uploaded copy). This is a description of the work reported in that thesis, not a claim that every term is implemented or independently verified in the present software.

**Repository inspected:** `Frankyang2002/T3neutrinomassmodelpipeline`, read-only `main` commit `91f8d7b0f81417df751632bacb436156a04872de` (15 November 2026). The user's uncommitted local checkout was not inspected.

## Central physics

The thesis concerns **sixteen ordinary T3 one-loop radiative neutrino mass models**, not an unrestricted survey of new particle physics. Their new fields are colour singlets and Z2-odd, with electroweak representations of dimension 1, 2 or 3. The five representation families are

| Family | `(d_phi, d_phi_prime, d_psi)` | Allowed alpha in the thesis |
| --- | --- | --- |
| A | `(1, 3, 2)` | `-4, -2, 0` |
| B | `(2, 2, 1)` | `-3, -1, +1` |
| C | `(2, 2, 3)` | `-3, -1, +1` |
| D | `(3, 1, 2)` | `-2, 0, +2` |
| E | `(3, 3, 2)` | `-4, -2, 0, +2` |

The thesis uses `Q = T3 + Y` and

\[
Y(\phi)=\frac{\alpha}{2},\quad
Y(\phi')=\frac{\alpha+2}{2},\quad
Y(\psi)=\frac{\alpha+1}{2}.
\]

The topology interactions schematically involve `Y1 L psi^c phi`, `Y2 L psi phi'^dagger`, and `lambdaT3 H H phi phi'^dagger`. The symmetry-allowed additional scalar and Yukawa interactions are relevant to Appendix C; in particular, B(-1) and C(-1) have exceptional Majorana-flavour contributions.

## Required effective-field-theory order

1. Start from **UV = SM + all T3 fermions + all T3 scalars**.
2. **Integrate out every heavy fermion**, retaining all heavy scalars and the appropriate induced EFT operators.
3. Evolve the intermediate scalar-containing EFT (renormalisable RGEs and the Wilson coefficients relevant at the chosen order).
4. **Integrate out the complete scalar sector**, then use SMEFT Weinberg running to the low scale.
5. Extract `m_nu = -(v^2/2) C5`, Takagi masses, mass-squared splittings and PMNS parameters.

This is the *scientific EFT prescription*. The current repository contains both sequential-theory building blocks and conditional one-scale numerical routes. An output from a conditional route must not be relabelled as a complete sequential calculation merely because it produces the same final observable names.

## What the thesis actually investigates

| Thesis section | Subject | Required repository capability |
| --- | --- | --- |
| §3.1–3.2; Appendix C | Gauge representations, neutral-compatible catalogue, allowed UV interactions | Model catalogue, CG/group factors, Lagrangian construction, scalar-operator inventory |
| §3.3–3.4 | UV RGE, fermion threshold, intermediate Wilson mixing, scalar threshold and low-energy observables | Matchete/RGBeta, threshold interface, scalar-EFT Wilson evolution, SMEFT, Takagi/PMNS |
| §4.1–4.2 | UV benchmark and 18-real-Yukawa fitting | Configuration, accepted benchmark source, Sobol/least-squares controls and fit targets |
| §5.1; Appendix B tables | Transfer of A(-4) and D(-2) fitted Yukawas to all 16 models; Weinberg matrix proportionality; B(-1)/C(-1) departures | Fitted reference matrices, saved C5 and observables, transfer/summary scripts |
| §5.2; Figs. 5.2–5.4 and Appendix B | UV/scalar-EFT/SMEFT running toggles; Wilson mixing comparison | Frozen-stage comparisons, source data, heatmap generators |
| §5.3; Fig. 5.5 and Appendix B | Compensated `T1/T2/T3` coupling sensitivity with `gamma` from `0.99` to `1.01` | Saved sensitivity points, plot scripts, NuFIT/JUNO comparison inputs |

The thesis benchmark in Table 4.1 specifies `muUV=10^7`, `muF=10^5`, `muS=10^3`, `muLow=10^2 GeV`, `M_F=diag(1,1.1,1.2)*10^5 GeV`, scalar mass-squared inputs `1.30*10^6` and `1.00*10^6 GeV^2`, `lambdaT3=0.10`, `lambdaOther=0.01`, and real 3x3 Yukawas. It sets the mixed scalar bilinear to zero at the UV scale. Its fitting objective uses five NuFIT observables with 18 real Yukawa parameters and a Sobol/least-squares approach (thesis §4.2).

## Important source-to-thesis distinctions

**The uploaded 109-page thesis is the selection authority for this cleanup, but it does not override actual source behaviour.** The GitHub tree includes separate benchmark generations, and some are not the same as the PDF:

- The saved `uniform001` reference and transfer documentation assigns `lambdaT3=0.05` to family E, whereas the thesis Table 4.1 says `lambdaT3=0.10` for each class. Do **not** delete either benchmark as redundant or silently combine their results; reconcile the exact fitted configuration against the figure/table source.
- The `uniform001` running-stage README explicitly describes **conditional one-scale C5 matching** rather than independently transported scalar-EFT LLSS coefficients, whereas thesis §3.3 describes the intended sequential approach. These are different levels of physical validation.
- The older running-stage heatmaps and the new `uniform001` figure-generation route both exist. The new Section 5.2 directory in the inspected GitHub commit contains instructions and TeX captions, **but not the regenerated 16-model heatmap PNGs**. Retain the older source figures and their provenance until the specific PDF figures can be matched to them.
- The current source includes broad experimental hard matching, finite one-loop closure, side projects and migration utilities that are **not necessary to read the thesis**. Their presence does not mean they are safe to delete: a source script may own a single figure, accepted fit or saved coefficient.

## Scope exclusions

The public-facing documentation does not need to reproduce individual developer batches, worker assignments, failed experiments, abandoned roadmaps or the chronological sequence of implementations. Store unreferenced historical notes in an **external, reversible archive**, not in the main reader-facing project tree.

Do **not** drop a scientific module, generated beta function, Wilson tensor, fitted Yukawa matrix, result table, plot generator, image or calculation manifest simply because its filename suggests a prior development phase. Until a complete dependency/provenance audit has linked the exact thesis figure to its exact upstream data, such items are **review**, not **delete**.
