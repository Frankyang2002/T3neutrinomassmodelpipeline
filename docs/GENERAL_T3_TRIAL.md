# General T3 scalar interaction trial — batched validation package

**Source snapshot:** `Frankyang2002/T3neutrinomassmodelpipeline` main
`41194fcd78fc7fe3c9378c9f8ee41228ad6473e7` (24 Oct 2026).

**Scope: opt-in experimental extension. This is NOT yet a verified fully
interaction-general Matchete/RGBeta/numerical pipeline.** The working model's
ordinary interactions remain unchanged unless specifically enabled. The ZIP
only edits local files, never GitHub.

## What the user's first run proved, and what it did not

The first run of `RunGeneralScalarAudit.wl` generated all **seven** additional
complex B(-1) quartic candidates. `FrontEndObject::notavail` occurred at every
`CheckLagrangian` call, and every check returned `$Failed`. The previous
script called `UsingFrontEnd` only for package loading, then performed the
validation outside that scope. This was an **audit-runner execution bug**;
zero accepted candidates was not evidence that their SU(2) contractions were
invalid. A second issue was that even the combined baseline check happened
outside that scope, so it also failed. The corrected runner encloses its **full
field/model construction and all validation inside**
`UsingFrontEnd[T3TrialAuditRun[trialModel]]`, matching the pattern of the
normal `RunModel.wl` builder.

Wolfram's frontend must actually be available. If the service frontend cannot
start, the new audit reports `Status: FrontendUnavailable` and returns null
accepted/valid counts, *not* seven invalid operators. The `BaselineLagrangianValid`
check runs before new candidates are tested. The final combined check always
includes every extra quartic, even if individually rejected.

## Implemented (new/revised relative to first ZIP)

- **Full replacement** `Lagrangian/RunGeneralScalarAudit.wl` — corrected
  frontend scope; baseline-first and full combined checks; machine-readable
  status and per-candidate outcomes.
- **Full replacement** `Lagrangian/interactions/GeneralScalarOperators.wl`
  — SU(2) invariant enumeration and symmetrisation, hypercharge and Z2 selection,
  complex Matchete couplings, and **fail-fast** behavior when an invariant or
  coupling definition fails (previously a structure could be dropped silently).
- **Full replacement** `Lagrangian/interactions/T3Topology.wl` — same opt-in
  `T3_EXPERIMENTAL_GENERAL_QUARTICS=1` production hook from first package.
- `scripts/RunGeneralScalarTrial.ps1` — drives an audit of B(-1), B/C(-1),
  or all 16 ordinary model points. Can separately run baseline and extended
  Matchete matching into isolated directories, with logs and JSON status.
- `scripts/ReportGeneralScalarTrial.py` — checks observed audit JSON against
  the 16-model independently derived operator manifest.
- `scripts/InspectGeneralScalarMatching.py` — compares *symbolic textual*
  Weinberg artifacts from baseline and extended runs, reports occurrences
  of `lambdaT3Extra...` couplings. It **does not** prove numerical equality
  or physical completeness.
- `Numerical/experimental/ScalarQuadraticPrototype.py` — isolated
  diagonalisation of the B/C(-1) complex-doublet mass matrix and splitting
  of Y=0 complex singlet/triplet into two real multiplets. This illustrates
  the quadratic-mixing extension **without wiring invalid non-diagonal
  propagators into the present matching**.
- Static operator manifest, Python operator counter, standalone unit tests.

**Predicted scalar-sector additions**, counting independent couplings per
separate supported ordinary model: 50 complex quartic couplings over 16 model
points and 8 complex quadratic terms over the same set. These are not 50
couplings of a single model. B/C(-1) each need seven new complex quartic
contractions relative to the Matchete builder. Counts use the original
complex-scalar field definitions and the gauge plus Z2 symmetries.

## Explicitly not implemented

- The additional **quadratic terms in the live Matchete/RGBeta/numerical model**.
  This needs propagator/basis changes, mass-matrix RGEs, rotated Yukawas,
  thresholds and matching. They are catalogued and prototyped, not enabled.
- **RGBeta beta functions for the complete added quartic basis.** The
  B/C(-1) RGBeta model already has six special quartics in a different basis,
  but that does not mean the new Matchete couplings have been connected to
  them. No numeric state-vector or beta-file changes are included.
- Newly allowed fermion-sector terms at special quantum numbers, or a
  repaired physical shared-scalar Matchete matching treatment.
- Matchete/Mathematica tests in the build environment: Wolfram is not
  available here. All actual Matchete results require your machine.

## Installation (PowerShell, repository root)

If you have local uncommitted changes, save or commit them first. The ZIP
replaces `T3Topology.wl` and the previous trial files. Consider a local-only
branch before expanding the ZIP:

```powershell
# Optional: creates a local experiment branch, does not push to GitHub.
git switch -c experiment/general-t3
Expand-Archive -LiteralPath .\T3_General_Scalar_Batch2.zip -DestinationPath . -Force

python -m pytest -q tests\python\test_general_scalar_trial.py tests\python\test_scalar_quadratic_prototype.py
python tests\python\check_t3_scalar_quartics.py --output output\general_scalar_static_audit.json
```

## Run many models with one command

First run the corrected B/C(-1) Matchete validity checks, in separate
`output/general_t3_trial` subfolders, with full logs:

```powershell
powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralScalarTrial.ps1 -Scope BC -Phase Audit
python scripts\ReportGeneralScalarTrial.py --root output\general_t3_trial
```

If those pass, run **all 16** ordinary supported cases:

```powershell
powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralScalarTrial.ps1 -Scope All -Phase Audit
python scripts\ReportGeneralScalarTrial.py --root output\general_t3_trial --all --require-complete --output output\general_scalar_audit_report.json
```

Some models have zero extra couplings. They still must pass the baseline and
combined checks. A failed result exits nonzero and is retained as structured
JSON plus a log. `-Scope All -Phase Audit` runs each model in its own
`wolframscript` process so an individual failed session does not contaminate
other model registry states.

To test symbolic one-loop matching **separately** for B(-1), first ensure
its audit passed, then:

```powershell
powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralScalarTrial.ps1 -Scope B -Phase Matching
python scripts\InspectGeneralScalarMatching.py --root output\general_t3_trial --output output\general_matching_symbolic_report.json
```

This runs original (flag OFF) and extended (flag ON) matching in different
output folders and restores the previous environment flag afterward. Matching
can be expensive. To audit and then automatically match only passing model
points in one command, use `-Scope BC -Phase Both` or `-Scope All -Phase Both`.
**Caution: extended symbolic matching is not yet consistent with the present
RGBeta and numerical state. Do not feed these experimental matching artifacts
into production numerical results.**

## Quadratic matrix prototype only

```powershell
python -c "from Numerical.experimental.ScalarQuadraticPrototype import mixed_t3_bc_doublet_spectrum; s=mixed_t3_bc_doublet_spectrum(700**2,900**2,(200**2)+0j); print(s.mass_squared,s.mass)"
```

For B/C(-1) the gauge-basis matrix used is

\[
M_S^2 = \begin{pmatrix}m_{S_1}^2&-\mu_{12}^2\\
-(\mu_{12}^2)^*&m_{S_2}^2\end{pmatrix},
\quad \Phi=(\widetilde S_1,S_2)^T.
\]

This is an **unbroken-gauge-phase scalar spectrum prototype**. Quartic Higgs
VEV effects and genuine RG mass mixing are not included.

## How to interpret status

- `Success`: baseline and all individual extra operators and their combined
  Lagrangian passed Matchete `CheckLagrangian`.
- `FrontendUnavailable`: Mathematica could not provide a usable frontend;
  *no gauge-invariance conclusion follows*.
- `BaselineValidationFailed`: old/ordinary Lagrangian failed before any
  extra candidate was tested. Diagnose Matchete installation/baseline first.
- `ExtraQuarticGenerationFailed`: CG/tensor/coupling construction failed.
- `IndividualOperatorValidationFailed`: at least one new operator failed
  `CheckLagrangian`, names/results are in JSON.
- `CombinedLagrangianValidationFailed`: individual checks passed but the
  full combined expression did not.

## Reversion

With changes safely committed on an experiment branch, return to the previous
branch with `git switch <original-branch>`; no GitHub push is needed. If the
experiment files are uncommitted, use `git status --short` to inspect them
before reverting. In particular, `git restore --
Lagrangian/interactions/T3Topology.wl` discards **all** your uncommitted edits
to that file, not only edits from this ZIP.

## Batch 3: Windows PowerShell parser correction (26 October 2026)

The Batch 2 PowerShell orchestration file contained an interpolation error:
`"Skipping matching for $name: ..."` is not valid Windows PowerShell because
`$name:` is parsed as a scoped variable reference. Batch 3 replaces the
entire runner, using `${name}:`. This defect prevented *any* audits or
matching from starting, regardless of the previously passing Python tests.

Batch 3 also clarifies the no-audits-found result of
`ReportGeneralScalarTrial.py` and returns a nonzero exit code for that case.
The physics modules, symbolic operator catalogue, and mass prototype are
unchanged from Batch 2. This is a runner hotfix, not a newly verified
Matchete result.

Install the Batch 3 ZIP **over** Batch 2 from the repository root:

```powershell
Expand-Archive -LiteralPath .\T3_General_Scalar_Batch3_Fix.zip -DestinationPath . -Force
python -m pytest -q tests\python\test_general_scalar_trial.py tests\python\test_scalar_quadratic_prototype.py

# One large run: audit all 16 representations; if a point passes, compare
# baseline and extended symbolic one-loop matching for that point.
powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralScalarTrial.ps1 -Scope All -Phase Both

# Summarise all 16, including failed/unrun points.
python scripts\ReportGeneralScalarTrial.py --root output\general_t3_trial --all --output output\general_scalar_audit_report.json
python scripts\InspectGeneralScalarMatching.py --root output\general_t3_trial --output output\general_matching_symbolic_report.json
```

**Important:** `-Scope All -Phase Both` may initiate expensive one-loop symbolic
matching in up to 16 baseline and 16 extended runs, conditional on the audits
passing. If you prefer the faster physics validation first, choose
`-Scope BC -Phase Both` before running `-Scope All -Phase Both`.
A failed audit skips matching for that model and records its log.

## Batch 4: findings from the first successful 16-model Matchete run (26 October 2026)

User-provided Windows log showed **all 16 baseline and extended symbolic matching jobs exited with code 0** and all attempted extra quartics passed individual and combined Lagrangian checks. This is encouraging but not an independent end-to-end physics verification; tree matching in the sample emitted non-fatal `FrontEndObject::notavail` warnings.

The independent Python operator inventory and the observed Mathematica generator agreed for 12/16 models. The missing four were B/C at `alpha=-3,+1`, with one *valid* quartic invariant each:

- B/C(`-3`): `Hdag^2 S2dag^2`, counts `{0,2,0,0,0,2}`.
- B/C(`+1`): `Hdag^2 S1^2`, counts `{0,2,2,0,0,0}`.

This package now uses the independently enumerated **SU(2)-singlet verified monomial manifest** for the 16 production scalar representation points, and requires Matchete to reproduce the *number of independent contractions* of every listed structure. The original dynamic hypercharge enumeration is retained as an uncertified fallback for unsupported models. Crucially, a successful `CheckLagrangian` on zero generated operators does **not** count as completeness verification anymore. The direct cause of the previous four empty Mathematica lists was not reproduced locally (no Matchete installed here); the new approach avoids relying on those lists for supported cases.

The logs also showed *symbolic* `C5` dependence in B/C(`-1`) on precisely two new couplings:

- `lambdaT3ExtraQ020020I1` (`Hdag^2 S2^2`, up to Hermitian conjugation)
- `lambdaT3ExtraQ020200I1` (`Hdag^2 S1dag^2`, up to Hermitian conjugation)

Other tested points showed no newly named coupling in the **available extracted C5 artifacts**. This is **not** a general proof of zero effects in all EFT threshold configurations or at higher orders. No UV numerical input has yet assigned these new couplings nonzero values.

### One-command local retest of the four incomplete points

Install this full consolidated ZIP at your repository root and run:

```powershell
Expand-Archive -LiteralPath .\T3_General_Scalar_Batch4.zip -DestinationPath . -Force
python -m pytest -q tests\python\test_general_scalar_trial.py tests\python\test_scalar_quadratic_prototype.py
python scripts\GenerateGeneralScalarWLManifest.py --check
powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralScalarTrial.ps1 -Scope Problems -Phase Both
python scripts\ReportGeneralScalarTrial.py --root output\general_t3_trial --all --require-complete --output output\general_scalar_audit_report.json
python scripts\InspectGeneralScalarMatching.py --root output\general_t3_trial --require-complete --output output\general_matching_symbolic_report.json
```

`-Scope Problems` selects only B/C(-3,+1). The other 12 completed audit/matching artifacts in the same output folder are retained. A successful aggregate report should read **16 passed / 0 failed**, if the new contractions pass Mathematica validation. Run the full 16-model audit as a further check if desired:

```powershell
powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralScalarTrial.ps1 -Scope All -Phase Audit
```

The PowerShell runner now independently gates matching on the manifest count, even when the Mathematica output mistakenly says `Status=Success`. The symbolic report also checks the number of extra interactions named in the extended matching logs and records text equality of baseline/extended `c5_coefficient.txt` after whitespace removal. **Text equality is not an algebraic equivalence test.**

### Actual implementation boundary

This batch completes the **Matchete candidate generation plus count enforcement for the scalar quartic sector**, subject to the four new local checks. It does not claim to complete the physical UV theory. No active scalar quadratic-mixing propagators, missing Yukawa sectors, corresponding general RGBeta beta functions, RG-to-Matchete CG normalisation map, or numerical boundary/state/threshold propagation were implemented. These are substantial independent changes, and must not be treated as verified solely because symbolic matching succeeded.

## Batch 5: complete scalar-quartic symbolic matching has passed all 16 models

**User-supplied local run on 26 October 2026:** 16/16 Matchete operator
completeness audits passed; all original and extended one-loop matching runs
exited successfully. The text-level `C5` comparison reports new formal
couplings only in ordinary B/C at alpha=-1:

- `lambdaT3ExtraQ020020I1`: `(Hdag)^2 S2^2`, the same *field structure* as
  RGBeta's `lambdaHHdagS2S2`.
- `lambdaT3ExtraQ020200I1`: `(Hdag)^2 (S1dag)^2`, the same *field structure*
  as RGBeta's `lambdaHHdagS1barS1bar`.

The RGBeta structures are **not proven to have the same tensor normalisation**
(or phase conventions) as the dynamically normalised Matchete CGs. Do not
identify the numerical values of these pairs until this is checked.

The user's output also contains `FrontEndObject::notavail` messages during
some tree-matching steps, but every one-loop matching run completed and
extracted C5. This message alone does not invalidate a successful result;
formal regression tests are added below to verify the matching artifacts.

### Batch 5: ONE analysis command (reuses existing matching; no expensive reruns)

```powershell
powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralT3Analysis.ps1 -RequireProof
```

It runs four analyses: verifies all 16 operator counts, inspects textual C5
outputs, enumerates Matchete/RGBeta field-structure coverage, then uses Wolfram
to test formally that setting all **new** quartic couplings to zero in the
extended coefficient gives exactly the original coefficient.

The zero-extra check must either prove the identity or explicitly say the
symbolic reduction was not proven; a missing or timed-out result is not
considered a pass. If an identity cannot be shown, inspect its residual
rather than asserting it must vanish. This is an output-level perturbative
regression, **not** a proof that the extra couplings remain zero under RG flow.

Outputs relative to project root:

- `output/general_scalar_audit_report.json`
- `output/general_matching_symbolic_report.json`
- `output/general_rgbeta_coverage.json`
- `output/general_t3_trial/c5_algebraic_report.json`
- `output/general_t3_trial/<MODEL>/general_c5_algebraic/delta_c5_inputform.txt`
- `output/general_t3_trial/<MODEL>/general_c5_algebraic/delta_c5_atoms.json`
- `output/general_t3_physics_summary.json`
- `output/general_t3_physics_summary.csv`

The last two files are compact machine-readable and human-readable tables
for your thesis notes. All newly generated files go to `output/`.

### RGBeta coverage result and why numerical work remains separate

The 16-model manifest contains 50 **per-model** additional complex quartic
couplings. Counting gauge-invariant field structures, 12 of them have existing
B/C(-1) RGBeta counterparts (six in each of two models), and 38 do not.
Within each B/C(-1), one of the two independent mixed Higgs portal
contractions is absent in RGBeta.

This count is *structural* only; it does not verify CG tensor normalisation or
beta functions. `scripts/CheckGeneralT3RGBetaCoverage.py` explicitly labels
these distinctions and checks that the six known special RGBeta coupling
names are actually defined in the local source. This is the next needed
input to a fully general RGE state extension. The quadratic scalar mixing
terms are still prototypes, not a working generalized matching formalism.

## Batch 6: formal C5 response, Matchete tensor normalisation, and RG source terms

The complete 16/16 Matchete scalar-quartic audit and 16-model zero-extra C5
proof have passed **on the user's Windows environment**. Batch 6 does not
rerun those expensive matching calculations. It adds the next diagnostics;
**it is not a numerical neutrino-mass implementation**.

From the project root, after extracting the full Batch 6 ZIP:

```powershell
python -m pytest -q tests\python\test_general_scalar_trial.py tests\python\test_scalar_quadratic_prototype.py tests\python\test_general_scalar_c5_proof.py tests\python\test_general_rgbeta_coverage.py tests\python\test_general_t3_stage6.py
powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralT3Stage6.ps1
```

This uses the **existing**
`output/general_t3_trial/c5_algebraic_report.json` and matching artifacts,
and the **existing** RGBeta JSON exports
`output/T3_B_alpha_m1/data/uv_rgbeta_rge.json` and
`output/T3_C_alpha_m1/data/uv_rgbeta_rge.json`.

Outputs:

- `output/general_c5_response.json`: validates the independent, linear
  coupling derivatives in the formal B/C(-1) \(\Delta C_5\); labels whether
  their coefficient expressions include poles and logarithms. A failure is
  never interpreted as a proven identity.
- `output/general_rgbeta_source_proof.json`: substitutes zero for all six
  special B/C(-1) RGBeta quartics and records the formal one-loop beta-function
  source for the two C5-active quartics. Demonstrates whether setting all
  these couplings to zero is RG-invariant.
- `output/general_t3_trial/matchete_active_quartic_tensors.json`: computes
  Matchete SU(2) tensors for the two C5-active quartics and compares them to
  a **defined** delta/epsilon symmetrised reference tensor. It reports the
  proportionality factor and an exact component-by-component residual.

**Normalisation limitation:** This last step determines only a Matchete-to-
canonical-index factor. In RGBeta the same operator is written using
`delS2` and RGBeta's own invariant conventions, and `AddQuartic` applies no
implicit combinatorial normalisation. A separate RGBeta-to-canonical check
and Lagrangian sign/factor check are required before identifying numerical
values of `lambdaT3ExtraQ020020I1` with `lambdaHHdagS2S2`, or
`lambdaT3ExtraQ020200I1` with `lambdaHHdagS1barS1bar`. Do **not** feed an
unverified transformation into the production neutrino-mass pipeline.

If the existing RGBeta output files are absent, run the script with `-SkipRGE`
and inspect the Matchete proof separately. If Matchete fails a canonical tensor
check, its JSON preserves the raw tensor to help diagnose SU(2) index
orientation rather than silently declaring a unit conversion.

The one-loop RGBeta normalisations and the missing 38 structurally independent
quartic couplings outside the existing RGBeta set remain unimplemented.
Quadratic scalar mass mixing remains catalogued but is not enabled in matching
or RG evolution. Shared-scalar matching is also a separate physical model.

## Stage 9 tensor exporter regression fix (26 October 2026)

See `docs/GENERAL_T3_STAGE9_STATUS.md`. `LoadModel["SM"]` must run after
`ResetAll[]` in the standalone tensor exporter before BSM representations
are defined. For a targeted rerun only:

```powershell
powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralT3Stage6.ps1 -TensorOnly
```

This does not claim that RGBeta's scalar coupling convention has yet been mapped.

## Stage 10 (experimental numerical pilots)

See `docs/GENERAL_T3_STAGE10_STATUS.md` for the two independent numerical probes and their strict limitations. Run `scripts/RunGeneralT3Stage10.ps1` after the successful Stage 8/9 formal-C5, RGBeta and tensor diagnostics. Stage 10 adds no new production numerics or RGBeta-to-Matchete coupling identification.

## Stage 11 — finite scalar mass-ratio logarithms

See `docs/GENERAL_T3_STAGE11_STATUS.md`. The experimental Matchete hard-C5
pilot now retains finite physical mass-ratio logarithms after applying
`μbar2 = MF^2` (instead of refusing any surviving `Log`). The toy input now
supplies positive `MS1`, `MS2` in GeV. Unrecognized symbolic scales and
unsubtracted UV poles still block evaluation and generate a diagnostic JSON.
The PowerShell runner always attempts all four independent B/C probes before
reporting failure. No production RGE or matching code was changed.

## Stage 12 extra-only mass-matrix/phase diagnostics

After the successful Stage 11 numerical pilots, run:

```powershell
powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralT3Stage12.ps1
```

This reuses the stored two-model Matchete `delta_c5` matrices and reports the additional contribution to `m_nu = -(v^2/2) C5`, the three extra-only Takagi singular values and a 16-point quartic phase scan. No further symbolic matching is performed. Full details and limitations: `docs/GENERAL_T3_STAGE12_STATUS.md`.

## Stage 13: coupled RGBeta UV numerical evolution and SU(2) tensor mapping

See `docs/GENERAL_T3_STAGE13_STATUS.md`. After the Stage 9 tensor audit and
Stage 8 RGBeta-source proof have succeeded, run:

```powershell
python -m pytest -q tests\python\test_general_t3_stage13.py
powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralT3Stage13.ps1
```

This uses real one-loop coupled RGBeta beta functions for B/C(-1) **above the
heavy-fermion thresholds**, including the six complex scalar quartics already
registered in RGBeta. It is an isolated experiment, not a replacement for the
production hierarchical numerical neutrino-mass pipeline. Group tensor ratios
are checked, but no unverified equality of complex Matchete and RGBeta
couplings is asserted.


## Stage 14 (conditional UV-induced C5 sensitivity)

See `docs/GENERAL_T3_STAGE14_STATUS.md`. After Stage 11 and Stage 13
outputs exist, use `python -m pytest -q tests\python\test_general_t3_stage14.py`
then `powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralT3Stage14.ps1`.
This uses hypothetical sign/conjugation maps, **not** a proven Matchete/RGBeta
Lagrangian coupling conversion, and never modifies production files. One-loop
RGE-generated quartics inserted into one-loop C5 are formally of order hbar^2.

## Stage 15 — source-convention map for B/C(-1)

Stage 14's 16 hypothetical sign/conjugation scenarios are narrowed to the
source-derived convention \(\lambda_{\mathrm{Matchete}}=-\lambda_{\mathrm{RGBeta}}\)
for the two C5-active B/C(-1) quartics *in the same named scalar basis*.
This follows the RGBeta MBL quartic sign (Thomsen, EPJC 81 (2021) 408,
Eq. 2.3), the experimental Matchete positive `PlusHc` interactions, and the
previously verified +1 SU(2) tensor factors. The source checker fails closed
on field-word, tensor, or interaction-definition changes; the numerical
calculation must reproduce Stage 14's `negative/negative` case.

```powershell
python -m pytest -q tests\python\test_general_t3_stage15.py
powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralT3Stage15.ps1
```

Results are written under `output/general_t3_stage15/`; no production modules
are modified. **Still not threshold-consistent matching or a complete two-loop
calculation**. See `docs/GENERAL_T3_STAGE15_STATUS.md`.
