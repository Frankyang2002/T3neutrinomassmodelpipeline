# T3 scalar interactions: Batch 1 relevance and perturbative-order audit

**Scope:** Latest inspected GitHub `main`: `a02e04c6def5cf8522b2f612ccd6bff7d9963b91`. No production files are changed by this batch.

## Objective and perturbative accuracy

We want the scalar interactions needed for a **specified leading one-loop neutrino-mass amplitude**, including an explicit accounting of the possible impact of RG running. This is **not** a demand to remove irrelevant terms from the full symmetry-allowed potential. The UV scalar Lagrangian and the interaction subset needed for one coefficient are distinct objects.

Use the existing unbroken-electroweak, fermion-first LLSS -> LLHH construction. Schematically the fermion threshold produces

\[
C_{LLSS} \sim \frac{y_1 y_2}{M_F},\qquad
C_5^{(1)} \sim \frac{1}{16\pi^2}\,C_{LLSS}\,\lambda_{HHSS}\,F(M_{S_1},M_{S_2},M_F).
\]

Here \(\lambda_{HHSS}\) means any correctly contracted Higgs-bearing scalar vertex that actually closes the two relevant scalar lines. **Field content and loop order only identify candidates:** they do not certify the SU(2) contraction, scalar flavour flow, lepton-number spurion content or physical normalization. The existing topology \(\lambda_{T3}\) is required for the selected traditional T3 contribution, not necessarily the only possible one-loop quartic in all quantum-number assignments.

### Distinct routes into the observable

| Interaction | One-insertion LLSS -> LLHH counting | Defensible omission |
|---|---|---|
| Required \(HH S_1 S_2^\dagger\) quartic | One loop | **Not** if studying the stated topology |
| Additional quartic containing two H/H† and two BSM scalars | Candidate at one loop (subject to tensor/flavour selection) | Only after coefficient/diagram-level proof, or explicit zero boundary |
| Quartic with four BSM scalars and no Higgs | Cannot be the **sole** closing quartic in a one-loop LLSS -> LLHH scalar loop with two external Higgs legs | May omit its *direct single-vertex amplitude* at this order; not necessarily the UV coupling, RG effects or all diagrams |
| Allowed extra quadratic bilinear | Tree-level modification of scalar masses, splittings or mixing | **Not** by loop counting; can be zero by explicit model assumption or controlled small-mixing approximation |
| Ordinary Higgs portals | Candidate Higgs-bearing interactions | Needs contraction-level analysis, not a blanket exclusion |
| Scalar mass terms | Already appear in propagators | Retain |

## Strict fixed-order versus RG-improved interpretation

If an extra quartic obeys \(\lambda_i(\Lambda)=0\) but has a one-loop beta-function source, then

\[
\lambda_i(\mu) \simeq \frac{\beta_i^{(1)}}{16\pi^2}\ln\!\frac{\mu}{\Lambda},\qquad
\delta C_5\big|_{\lambda_i\text{ generated}} \sim
\frac{\beta_i^{(1)}}{(16\pi^2)^2}\ln\!\frac{\mu}{\Lambda}\,\mathcal K_i.
\]

This is **formally higher fixed-loop order** than a single \(C_5\) scalar loop. It is **not identically zero** and need not be numerically small when logarithms, cancellations or unusual hierarchies amplify it. For \(10^5\to10^3\) GeV, \(|\ln(\mu/\Lambda)|=\ln100\approx4.605\), and \(\ln100/(16\pi^2)\approx0.0292\). This alone is not a bound on \(\delta C_5/C_5\).

The saved `output/general_rgbeta_source_proof.json` already exhibits nonzero formal sources for two B/C(-1) quartics when those coefficients vanish. Thus **zero at a boundary scale is not RG closure**. Moreover, a quartic that is **independently nonzero** at the matching scale enters at its first allowed loop order: one must not reclassify it as two loop merely because it was absent from the original code.

Extra quadratic terms are different: \(\mu^2S_iS_j\) modifies the tree-level mass matrix. Its mixing rotation can change the couplings used inside a leading one-loop diagram. A quadratic term radiatively induced from a zero boundary can have a perturbative suppression, but a nonzero input quadratic has no automatic extra loop factor. In near-degenerate scalar spectra even numerically small mixing elements can have appreciable effects.

## Evidence inherited from current checked-in reports

- The 16-model quartic audit reports 50 extra complex quartic **slots across separate models**, not 50 fields/couplings in a single theory.
- The independent quadratic inventory has eight extra model-specific slots.
- The saved symbolic matching report explicitly contains two new quartic couplings in **B(-1)** and two in **C(-1)**; the other new Higgs-bearing couplings are not present in the available extracted C5 texts. Lack of appearance **in that output** is not proof in all mass bases or threshold schemes.
- `output/general_t3_physics_summary.json` reports algebraic recovery of the original coefficient when all additional quartics are set to zero. This is a source-level regression result; it is not a check of nonzero scalar quadratic mixing.
- Existing B/C reduced scalar-threshold tests justify omitting four-heavy-scalar quartics from the **specific single-insertion LLSS -> LLHH one-loop matching call**, but explicitly do not verify full EFT Hermiticity. See `docs/T3_GENERAL_EFT1_LEADING_C5_SCALAR_THRESHOLD.md`.

## Batch 1 implementation and outputs

`scripts/AuditT3ScalarRelevance.py` reads five existing repository JSON reports; it **does not** run new Matchete or RGBeta jobs and it does not assume absent fields/couplings are dynamically irrelevant. It fails if any required input, 16-model audit record, algebraic zero-extra proof, quartic count, quadratic count, C5 file-availability flag or observed coupling mapping is missing or inconsistent. Source SHA256 hashes appear in outputs.

From the project root:

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
python -m pytest -q tests\python\test_t3_scalar_relevance_batch1.py
python -m scripts.AuditT3ScalarRelevance
```

Generated locally under `output/scalar_relevance/`:

- `operator_classification.json`: all model-specific entries with classifications, levels of evidence and scope-specific omission rules.
- `operator_classification.csv`: flattened review table.
- `summary.md`: model counts, status summary, limitations, and source hashes.

The original quartic portal/self families are listed by family (`*`) because this batch does not independently enumerate their per-model real SU(2) invariant multiplicities. The original builder already constructs these; claiming exact invariant counts for them would require a separate report. Extra quartic entries **do** use the saved independent contraction counts and actual Matchete coupling symbols. Extra quadratic audit identifiers are **not** claimed to be production coupling names.

## Batch 1 stop condition / following work

The output identifies:

- **Retain:** known one-loop topology and any extra quartic observed in the saved C5 for nonzero boundary values; propagator mass parameters.
- **Needs Batch 2 proof:** two-Higgs additional quartics *not* observed in saved symbolic C5, and ordinary Higgs portal effects. The missing term may be genuinely absent from the fixed-order amplitude, cancelled or convention/basis-dependent; Batch 1 never chooses without a diagram/matching check.
- **No direct single-vertex contribution to the targeted LLSS -> LLHH graph:** zero-Higgs four-BSM-scalar quartics. Their other roles remain.
- **Needs Batch 3 physical treatment or a stated zero input:** eight extra scalar quadratic terms.
- **Needs Batch 4 quantitative control:** running-induced couplings, benchmark stability and large-log/cancellation sensitivity.

**Do not label this Batch 1 report as a complete one-loop C5 proof.** It is the input inventory and evidence-controlled power-counting classification for that proof.

## Verified manifest count by model

These are operator-slot counts directly calculated from the checked-in 16-model manifest, **not** fresh one-loop diagram computations. "Higgs-bearing" means exactly two Higgs or conjugate Higgs legs, and "four BSM" means no Higgs legs.

| Model | Higgs-bearing extra quartics | Four-BSM-scalar extra quartics | Extra quadratic slots |
|---|---:|---:|---:|
| A(-4) | 0 | 0 | 0 |
| A(-2) | 2 | 3 | 1 |
| A(0) | 2 | 3 | 1 |
| B(-3) | 1 | 0 | 0 |
| B(-1) | 4 | 3 | 1 |
| B(+1) | 1 | 0 | 0 |
| C(-3) | 1 | 0 | 0 |
| C(-1) | 4 | 3 | 1 |
| C(+1) | 1 | 0 | 0 |
| D(-2) | 2 | 3 | 1 |
| D(0) | 2 | 3 | 1 |
| D(+2) | 0 | 0 | 0 |
| E(-4) | 0 | 0 | 0 |
| E(-2) | 2 | 4 | 1 |
| E(0) | 2 | 4 | 1 |
| E(+2) | 0 | 0 | 0 |
| **Total** | **24** | **26** | **8** |

Four of the 24 Higgs-bearing slots are observed in B/C(-1) symbolic C5 artifacts (two per model). The other 20 must not be marked "higher order" solely on the basis of the absence of their names in that report.

## Additional vacuum and EWSB caveats

Two-Higgs quartics can also induce electroweak-symmetry-breaking scalar-mass corrections of order \(\lambda v^2\); whether these are retained in a particular benchmark depends on the \(v/M\) expansion and the scalar spectrum. The 26 four-BSM-scalar operators can constrain boundedness from below and feed scalar beta functions, even though their *direct* insertion into the selected LLSS -> LLHH one-loop scalar closure is excluded by Higgs-leg counting. Batch 1 does not claim a stable scalar vacuum, bounded-from-below potential, or a physically complete broken-phase scalar spectrum.
