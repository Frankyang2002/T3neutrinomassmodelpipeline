# Workstream A v14–v15 — native RGBeta package-scope quartic projector repair

## Source diagnosis (checked against installed-package upstream source)

RGBeta v1.2.0 declares these as **PackageScope** rather than PackageExport:

- `RGBeta`PackageScope`Lam` in `Package/FieldsAndCouplings.m`
- `RGBeta`PackageScope`RefineGroupStructures` and `RGBeta`PackageScope`TStructure` in `Package/GroupsAndIndices.m`

The previous v10 and v13 Wolfram scripts called unqualified `Lam` and
`RefineGroupStructures` from a `Global`` script. Both stayed inert. This is
visible directly in the user's v13 witness as the unevaluated structure
`RefineGroupStructures[(Lam[...] * RGBeta`PackageScope`TStructure[...]) ...]`.
The 50/50 v10 smoke failures and the matching v13 isolated/full re-probe
outcomes therefore do **not** constitute evidence of a quartic normalization
failure. They diagnose a script-to-package context error.

## New implementation

The v14 native exporter is a **new file**, rather than altering the historical
v10 exporter and invalidating the immutable failed v10/v13 evidence. Its live
Wolfram calls are `RGBeta`PackageScope`Lam[...]` and
`RGBeta`PackageScope`RefineGroupStructures[...]`; it checks the package-scoped
symbols and DownValues before proceeding. It fails closed on unresolved
`Lam`, `TStructure` or `RefineGroupStructures` after contraction.

For each native source quartic `j` and native target coupling projector `i`:

\[
R_{ij}=P_i[\operatorname{Lam}_j]\big|_{\lambda_j=1},\quad
P_i[\operatorname{Lam}_j]\big|_{\lambda_j=0}=0,\quad
P_i[\operatorname{Lam}_j]\big|_{\lambda_j=2}=2R_{ij}.
\]

Contracted native `TStructure` expressions use the **installed RGBeta
scalar-contraction rules**, including its real/complex doubling matrix; the
probe does not supply custom numerical generator matrices.

- Three-model smoke: `A_-4`, `B_-1`, `E_-2`; 38 native quartic registrations,
  29 real self-projector diagonals, 21 portal response entries, 9 complex probes.
- All-16: 178 native quartic registrations, 150 real diagonals, 110 portal
  entries, 28 complex probes.
- Real response acceptance requires each self diagonal to equal 1 within
  `2e-8`; every portal matrix must be finite, full rank, and unit diagonal.
- Native complex probes independently set `lambda` and `Bar[lambda]` to
  (0,0), (1,0), (0,1), (2,3); they report the direct and conjugate responses
  as exploratory evidence. These do **not** gate the real projector audit.
- The v14 Python full audit additionally recalculates the existing v9 source
  and CG-provenance chain and reports the conditional Matchete observation
  matrices `R * M` without claiming a physical source coupling dictionary.
- Source hashes include this new v14 exporter, the current local T3 RGBeta
  model source, and `RGBetaNativeAccess.wl`. It also checks the installed
  RGBeta main-package file SHA256.

The PowerShell runner clears the four optional experimental flags in its
**own process only**, restoring their process values in `finally`. The
Wolfram exporter retains an independent flag guard. All outputs refuse
overwrite unless an already-existing derived report is identical.

## Install paths

The ZIP contains eight **complete new files**:

1. `scripts/RunT3ParallelANativeQuarticProjectorsV14.wl`
2. `scripts/RunT3ParallelAQuarticNativeContextAuditV14.py`
3. `scripts/RunT3ParallelANativeComplexChannelsV15.py`
4. `scripts/RunT3ParallelAQuarticNativeContextV14.ps1`
5. `Numerical/parallel_t3/wilson/native_context_quartic_v14.py`
6. `Numerical/parallel_t3/wilson/native_complex_channel_v15.py`
7. `Numerical/parallel_t3/wilson/WORKSTREAM_A_NATIVE_CONTEXT_V14_V15.md`
8. `tests/test_t3_parallel_a_native_context_v14.py`

From a PowerShell terminal at `C:\thesistemp\T3neutrinomassmodelpipeline`:

```powershell
$zip = "$env:USERPROFILE\Downloads\t3_workstream_A_native_context_v14_v15.zip"
Expand-Archive -LiteralPath $zip -DestinationPath . -Force
powershell -ExecutionPolicy Bypass -File scripts\RunT3ParallelAQuarticNativeContextV14.ps1
```

Use `-SmokeOnly` to stop after the three-model native smoke and v15
classification if it passes. When a native JSON already exists, the runner
revalidates it without replacing it. If the smoke or full native witness
remains blocked, the script prints per-coupling failures and exits 7.

## Outputs

All paths relative to repository root:

- `output/parallel_t3/wilson/rgbeta_quartic_projectors_native_v14_smoke.json`
- `output/parallel_t3/wilson/rgbeta_complex_conjugate_native_v15_smoke.json`
- `output/parallel_t3/wilson/rgbeta_quartic_projectors_native_v14.json`
- `output/parallel_t3/wilson/rgbeta_complex_conjugate_native_v15.json`
- `output/parallel_t3/wilson/native_quartic_projectors_conditional_v14.json`

The full v14 JSON is written only if its native real-projector stage passes
and the historical v9 source chain can be revalidated.

## Verified vs pending

**Verified here:** source-level PackageScope definitions (RGBeta upstream),
new Python schema/negative-control tests, native model and response counts,
Wolfram bracket balance, and preservation of historical evidence. **Not
run here:** Wolfram execution of the corrected package-scoped contractions;
requires your installed RGBeta on Windows. No v14/v15 native success is
claimed before that execution.

No fully physical scalar action, native numerical S2 frame, nonportal
Matchete–RGBeta source dictionary, B(-1)/C(-1) mixed mass basis, signed
GammaCC-to-Weyl matching, complete one-loop ADM or physical Wilson evolution
is certified by this batch. The EFT order is fixed:

`UV -> integrate out ALL fermions -> integrate out ALL scalars -> SMEFT`.
