# Batch 6 — explicit EFT1 scalar-bubble MSbar counterterm and native matching-scale closure

## Scope

**Repository revision inspected read-only:** `7b827e9f2e05a1ab277cbd0a32f63dfdad6d5bcd` (GitHub `main`, 10 Nov 2026).

The completed Batches 1, 2, 4 and 5 provide the source-resolved heavy-fermion tree LLSS boundary, the hard one-loop triangle, the seven-channel candidate EFT1 anomalous dimension, the 88 native Matchete PL Laurent kernels and the finite EFT1 scalar-loop integral. Batch 6 **reuses these saved sources**; it does not rewrite them or invoke Matchete.

The new genuinely computed object is the **EFT1 scalar-loop MSbar counterterm** for the dimension-five Weinberg operator, in the project's declared `C_ijab` full-sum and real-scalar conventions. It is *not* an unknown UV counterterm guessed to cancel Matchete's native pole.

The EFT sequence is unchanged:

```
UV -> integrate out ALL heavy fermions at one MF = 1e5 GeV
   -> EFT1 with BOTH S1 and S2 propagating
   -> integrate out BOTH scalar multiplets at one MS = 1e3 GeV
   -> SMEFT + Weinberg RG running -> conditional neutrino masses
```

## Calculation

For the existing tree-level EFT1 operator `C_LLSS` and source-resolved quartic `lambda_HHSS`, the zero-momentum scalar loop contains

\[
\mathcal A_{\rm EFT1,soft}^{(1)}
=-\frac{1}{2(16\pi^2)}
 C_{ijab}\lambda_{ABab}
 \left[\frac{1}{\bar\epsilon}+\overline B_0(0;m_a^2,m_b^2;\mu)\right].
\]

The computed **pole coefficient** of its MSbar counterterm is

\[
\delta C_5^{\rm EFT1,CT}
=+\frac{1}{2(16\pi^2)}\mathcal P_5\left[C_{LLSS}\lambda_{HHSS}\right]
\frac{1}{\bar\epsilon}.
\]

`P5` is the existing project's Weinberg tensor projector. We sum all **three heavy fermions**, both heavy scalar representations, all **real heavy scalar components**, and all nonzero **source-named quartic channels**. The heavy scalar eigenbasis is unnecessary in this one ultraviolet divergence because the complete real orthogonal rotation cancels. It remains necessary for the finite B0, which is reused from Batch 5.

The code checks three independently saved/constructed representations of this *same loop*:

1. Direct source-level full `LLSS × HHSS` tensor contraction.
2. The coefficient of the scalar-loop `1/epsilonbar` extracted from the original Batch 5 full massive triangle and B0 calculation.
3. The derivative of the finite B0 result under `mu -> MF/2, MF, 2MF`.

**The native Matchete pole does not get cancelled by setting an unknown counterterm equal to its negative.** The residual native pole is retained as an explicit 3×3 matrix and decomposed by original source quartic label. The seven existing EFT1 ADM terms are converted to *candidate* operator MSbar pole coefficients according to `beta / [2(16π²)]` in the stated d=4−2ε convention; their scalar-pair channel is compared with the explicitly computed scalar bubble. This candidate is **not** a derivation of native UV counterterms or an independent anomalous dimension.

## Matching-scale scenarios and physical limits

For the 3 scale factors `f = 0.5, 1.0, 2.0`, Batch 6 constructs

\[
C_5^{\mathrm{strict}}(M_S; f)
=C_5^{\mathrm{hard}}(fM_F)
+\frac{\ln[M_S/(fM_F)]}{\ln(M_S/M_F)}
\Delta C_5^{\mathrm{EFT1,strict}}(M_S)
+\Delta C_5^{\mathrm{scalar,bubble}}(M_S).
\]

The hard boundary is explicitly selected from three **separate** scenarios:

- `analytic_triangle`: calculated triangle hard local coefficient, slope fixed by its scalar CT. Reproduces Batch 2 at `f=1`.
- `native_formal`: formally pole-dropped matched-EFT PL coefficient, with its actual native logarithmic slope. Reproduces Batch 4 at `f=1`. **Not certified MSbar hard matching.**
- `native_as_UV_minus_soft_hypothesis`: native coefficient minus finite scalar soft, **assuming** the native saved Matchete `MatchedEFT` is a bare UV amplitude. Reproduces Batch 5 at `f=1`. **This assumption is not established; the scenario is not physical.**

All source Yukawas, couplings and scalar masses are **held at their saved nominal MF/MS benchmark points**; varying `f` varies the *subtraction scale*, not the physical scale at which the fermions are integrated out, and does not include consistent UV/background running under a changed threshold. A nonzero derivative is a quantified operator/matching-scale inconsistency, not proof of a missing counterterm formula. The user can skip 144 SMEFT evolutions with `--skip-low-energy-running`; by default all nine branches × 16 models run the existing one-loop SM + Weinberg ODE and export conditional masses and mixing angles.

**Not completed by Batch 6:** native UV counterterms, the UV hard-region decomposition of non-triangle operators, the Matchete `GammaCC` to Weyl sign/phase, physical normalization of all RGBeta quartics, and certified complete one-loop master ADM. Batch 6 does not silently assign values to any of those unknown objects. It should not be described as a completed physical one-loop matching calculation.

## Exact project-relative files

All files in the ZIP are **new** and belong under `C:\thesistemp\T3neutrinomassmodelpipeline`:

```
Numerical/parallel_t3/wilson/msbar_counterterm_batch6.py
scripts/RunT3CountertermBatch6.py
scripts/RunT3CountertermBatch6.ps1
tests/test_t3_counterterm_batch6.py
README_T3_EFT1_MSbar_COUNTERTERM_BATCH6.md
```

No existing source code, numerical engine or output is replaced. Every source-batch input is checked against the saved model and SHA256 source manifests. Every new output includes exact input file SHA256 hashes.

## Windows PowerShell commands

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline

Expand-Archive `
  -LiteralPath "$HOME\Downloads\t3_eft1_msbar_counterterm_batch6_all16.zip" `
  -DestinationPath . -Force

.\.venv\Scripts\python.exe -m pytest -q tests\test_t3_counterterm_batch6.py

.\.venv\Scripts\python.exe scripts\RunT3CountertermBatch6.py

.\.venv\Scripts\python.exe scripts\RunT3CountertermBatch6.py --show
```

For an interrupted run:

```powershell
.\.venv\Scripts\python.exe scripts\RunT3CountertermBatch6.py --resume
```

Fresh all-16 recomputation:

```powershell
.\.venv\Scripts\python.exe scripts\RunT3CountertermBatch6.py --overwrite
```

Counterterm and matching-scale scan without 144 low-energy ODEs:

```powershell
.\.venv\Scripts\python.exe scripts\RunT3CountertermBatch6.py --overwrite --skip-low-energy-running
```

Alternatively use `scripts\RunT3CountertermBatch6.ps1 -Resume`.

## New output locations

```
output/parallel_t3/eft1_msbar_counterterms_all16_batch6_v1/
    summary.json
    model_scenario_comparison.csv
    A_m4/
      counterterms_and_native_sources.npz
      counterterm_and_scale_ledger.json
      matching_scale_scenario_matrices.npz
      mass_and_mixing_scan.csv
    ... sixteen model directories
```

The NPZ includes the **source-quartic × fermion** counterterms, the explicit EFT1 pole and CT cancellation, seven candidate-ADM operator counterterm channels, the raw native pole per Matchete quartic and the unassigned native residual. The scan NPZ stores every `C5(MS)`, and by default the evolved low-energy `C5`. CSV includes masses, splittings and mixing angles for all nine branches per model.

## Validation

Python syntax compilation and **eight isolated algebraic tests passed** in the packaging environment. Those tests cover the actual loop pole/CT contraction, source aggregation, strict anomalous-dimension bookkeeping, matching-scale derivative, and three distinct scenario boundaries. They are not an all-16 Windows integration test. The full real 16-model run must be executed locally; report any blockers from `summary.json` if encountered.
