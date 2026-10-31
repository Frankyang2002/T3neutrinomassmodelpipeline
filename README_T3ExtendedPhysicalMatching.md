# T3 general-Yukawa physical matching batch

## Scope and prior inputs

This is a **local-only, non-destructive experimental matching batch** for the six A/D/E models with an independent Yukawa to the conjugate of the neutral scalar. It follows the completed 16/16 full coupled RGBeta UV numerical run.

It **does not modify** the installed `T3RGBetaModel.wl`, original Matchete builder, legacy `T3UVState`, existing matching outputs, fitted benchmark references, or GitHub.

**The one-step Matchete calculation is a new candidate matching computation, not the independently renormalized hierarchical two-threshold final C5 of the existing pipeline.** The additional scalar *quadratic* operators are deliberately disabled in this computation, while extra quartics and Yukawas are enabled. B/C mixed scalar mass matching, absolute Matchete–RGBeta field conventions, counterterm completion, and a 3-heavy-flavour C5 reconstruction remain separate gates.

## New files / exact project paths

- `Lagrangian/RunT3ExtendedC5Native.wl` — fresh Matchete native one-loop `EFTOrder->5`, complete 1-step heavy matching, extracted canonical Weinberg coefficient, algebraic zero-extra sector and difference. Tested syntactically here, **not live** (no Matchete kernel in this environment).
- `scripts/RunT3ExtendedMatchingBatch.py` — validates existing source and 16-model results, performs six model-matched numerical sensitivity tests and native one-loop matching in isolated Wolfram processes, writes per-model logs and a unified status report.
- `Numerical/experimental/T3ExtendedC5KernelProxy.py` — derives the unweighted, symmetric ordinary-T3 Yukawa flavour tensor from actual coupled UV trajectories and compares seeded versus zero initial extra Yukawas. **Not a physical C5 or neutrino mass.**
- `Numerical/experimental/T3ExtendedC5Observables.py` — independently computes neutrino masses, mixing and Takagi residuals from an externally supplied source-labelled, symmetric numerical C5 and charged-lepton Yukawa matrix. It **cannot certify** where the supplied C5 came from.
- `scripts/TestT3ExtendedPhysicalStackSynthetic.py` — integration tests with simulated native results; tests reject misidentified reports and overclaimed physical certification.

## PowerShell

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
$zip = "$HOME\Downloads\T3ExtendedPhysicalMatchingStack_Local.zip"
Expand-Archive -LiteralPath $zip -DestinationPath . -Force

python scripts\RunT3ExtendedMatchingBatch.py --self-test
python Numerical\experimental\T3ExtendedC5KernelProxy.py --self-test
python Numerical\experimental\T3ExtendedC5Observables.py --self-test
python scripts\TestT3ExtendedPhysicalStackSynthetic.py

# Check source and existing 16/16 RGE report without running Wolfram
python scripts\RunT3ExtendedMatchingBatch.py --preflight

# Recommended initial live native pilot for A(0), plus the numerical proxy
python scripts\RunT3ExtendedMatchingBatch.py --pilot

# All six real native matches and the six numerical proxies. Results survive
# a failed model; the command exits 4 if any model is blocked.
python scripts\RunT3ExtendedMatchingBatch.py --all
$LASTEXITCODE
```

If the user only wants the six immediate numerical sensitivity computations (no Matchete `Match`), run:

```powershell
python scripts\RunT3ExtendedMatchingBatch.py --proxy-only
```

If the one-heavy-flavour native experiment is incompatible with an installed Matchete version, run a three-heavy-flavour test using `--pilot --heavy-flavours 3`. **This option is not guaranteed to make matching succeed.** Native matching can be resource-intensive; examine the `.progress.json` and log for the last completed stage.

### Output locations

- `Reports/output/workstream_a/extended_c5_native/preflight.json`
- `Reports/output/workstream_a/extended_c5_native/summary.json`
- `Reports/output/workstream_a/extended_c5_native/native/<MODEL>_h1.json`
- `Reports/output/workstream_a/extended_c5_native/native/<MODEL>_h1.progress.json`
- `Reports/output/workstream_a/extended_c5_native/native/<MODEL>_h1.{full_C5,zero_extra_C5,extra_difference_C5}.inputform.txt`
- `Reports/output/workstream_a/extended_c5_native/logs/<MODEL>_h1.log`
- `Reports/output/workstream_a/extended_c5_native/kernel_proxy/models/<MODEL>.json`

## Physics and interpretation

The existing baseline flavor lift uses

\[
P_{pq}=\sum_r \left(\overline y_{1,pr}\overline y_{2,qr}+\overline y_{2,pr}\overline y_{1,qr}\right).
\]

The new numerical proxy computes the change of this **unweighted factor** due to the fully coupled UV feedback of an independent extra Yukawa. It excludes masses, loop functions, scalar mixing, the matching coefficients and cross-package conventions; **never label its numerical entries as C5 or a neutrino mass**.

The native Matchete script instead calls the repository's existing `RunT3Matching`, then `ExtractWeinbergCoefficient`. Its `zero_extra` coefficient is produced by **algebraic substitution in the completed full matching result**, not by an independently computed baseline. It is sensitive to the actual one-loop tensor content, but is only a simultaneous full-heavy one-step match in the restricted field basis.

The **complete physical extension remains blocked** if any of these has not been established:

1. A common, absolute Matchete/RGBeta physical field and coupling dictionary for `y1`, `y2`, extra Yukawa, `MF`, and scalar interactions.
2. The complete scalar quadratic mass basis, where applicable.
3. Consistent hierarchical threshold matching and counterterm/renormalization conventions, rather than one-step matching.
4. Three-heavy-generation flavour reconstruction and matching to the evolved UV matrices.

No native report or numerical proxy is silently promoted to a certified physical cross-package neutrino prediction.

## Optional standalone numerical C5 check

For a **numerical C5 matrix derived independently**, save a JSON document shaped as follows. This example shows the *schema*, **not physical model values**:

```json
{
  "kind": "t3_external_extended_c5_numeric_v1",
  "model": "A_p0",
  "scale_gev": 100.0,
  "vev_gev": 246.22,
  "ordering": "NO",
  "c5_gev_inverse": {
    "real": [[0,0,0],[0,0,0],[0,0,0]],
    "imag": [[0,0,0],[0,0,0],[0,0,0]]
  },
  "charged_lepton_yukawa": {
    "real": [[0.000003,0,0],[0,0.0006,0],[0,0,0.01]],
    "imag": [[0,0,0],[0,0,0],[0,0,0]]
  },
  "provenance": {
    "description": "Specify the actual source of your independently computed C5"
  }
}
```

Run `python Numerical\experimental\T3ExtendedC5Observables.py --input path\to\c5.json`.

The supplied C5's provenance is **not** authenticated; the output status is `NeutrinoObservablesFromSuppliedC5`, not a claim that the physical extended C5 has been established. The example's zero C5 does not represent a fitted or predicted neutrino mass.

## Tests performed by package author here

- Python syntax/compilation checks.
- Exact algebraic flavor tensor tests and negative controls.
- Numerical observable test in a controlled test harness with compatible physics interfaces.
- Six simulated Wolfram execution/report flows; malformed/mismatched/overclaimed native results rejected.
- Wolfram source brackets/comments/strings balanced; **not a live Wolfram syntax/run test**.
- ZIP extraction, installed-path and test runner checks.

**A real Matchete kernel and your local outputs are required to verify this package physically.**
