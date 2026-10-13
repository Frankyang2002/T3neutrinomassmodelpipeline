# Run from the T3neutrinomassmodelpipeline repository root.
# Consolidates INFO notes while preserving the complete existing main guides.
$ErrorActionPreference = "Stop"
$root = (Get-Location).Path
$required = @(
  "INFO_T3_MULTIPLICITY.md", "INFO_RGE.md",
  "INFO_T3_BENCHMARK_VALIDATION.md", "INFO_LLSS_DIAGNOSTIC.md",
  "INFO_T3_NUMERICAL_RESUME.md", "INFO_T3_NUMERICAL_SCENARIO_FIX.md"
)
foreach ($p in $required) {
  if (-not (Test-Path -LiteralPath (Join-Path $root $p))) {
    throw "Missing $p. Run from repository root and retain original INFO files."
  }
}
$archive = Join-Path $root "docs\archive"
New-Item -ItemType Directory -Path $archive -Force | Out-Null
$multPath = Join-Path $root "INFO_T3_MULTIPLICITY.md"
$rgePath = Join-Path $root "INFO_RGE.md"
$mult = Get-Content -Raw -Encoding UTF8 $multPath
$rge = Get-Content -Raw -Encoding UTF8 $rgePath
$multMarker = "## Benchmark validation and reproduction"
$rgeMarker = "## LLSS fixed-order diagnostic and power counting"
if ($mult.Contains($multMarker) -or $rge.Contains($rgeMarker)) {
  throw "One of the consolidation sections already exists. No files changed."
}
$multExtra = @'

## Benchmark validation and reproduction

The fixed comparison grid in `studies/FullT3Study.py` uses
`Y = 0.005` or `0.5` and `lambdaT3_reference = 0.1` or `1.0`.
The benchmark-only factors are:

```python
LAMBDA_T3_NORMALISATION = {"A": 1.0, "B": 1.0, "C": 1.0, "D": 1.0, "E": 0.5}
```

Hence A–D use `lambdaT3_effective = 0.1` or `1.0`, while E uses
`0.05` or `0.5`. This compensates only the two equal-weight E charge
channels in the unbroken-SU(2) multiplet-mass benchmark; it does not
equalize the full representation-dependent Clebsch–Gordan coefficients.

### Validation commands (PowerShell, repository root)

```powershell
python -m pytest -q tests/test_t3_multiplicity_benchmark.py
python -m pytest -q tests/test_full_study_architecture.py
python -m pytest -q
python pipeline.py --help
python pipeline.py --full
```

The full pipeline invokes symbolic/RGE stages and may take considerable time.
Unit-test success does not establish a completed 64-case end-to-end run.

### Output ownership and interpretation

Generated comparison inputs are stored under
`configs/generated_model_sets/comparison/<scenario>/input/`.
The full-study summary is `output/full/full_numerical_summary.json`.
The study records `comparison_point.lambdaT3_reference`,
`comparison_point.lambdaT3_effective`, `.normalisation_factor` and
`.normalisation_basis`, and summary normalization metadata.

The multiplicity-only initialization belongs to `studies/FullT3Study.py`.
It does not require modifying the physical RGE, Lagrangian, numerical
running, or Matchete group factors. Existing generated results should not
be assumed to reflect revised input settings until checked or regenerated.

Historical note: the original benchmark-validation runbook described the
64-case rerun as pending. Later run/recovery notes record numerical jobs,
but do not by themselves establish that every stored scenario used the
intended benchmark inputs. Verify the generated `input` configurations
and scenario-specific diagnostics before interpreting the results.
'@
$rgeExtra = @'

## LLSS fixed-order diagnostic and power counting

The component Wilson implementation in
`RGE/running/intermediate/ScalarOnlyWilsonFlow.py` transports the
dimension-five scalar-lepton coefficient between the fermion and scalar
thresholds at fixed-order leading-log accuracy:

\[
C_{\rm LLSS}(\mu)
=
C_{\rm LLSS}^{(0)}
+
\frac{1}{16\pi^2}\ln\!\left(\frac{\mu}{M_F}\right)
\beta^{(1)}[C_{\rm LLSS}^{(0)}]
+\mathcal O(\hbar^2).
\]

`Numerical/diagnostics/IntermediateWilsonDiagnostics.py` samples this
same symbolic expression at intermediate scales. It neither defines a
new beta function nor solves an RG-improved Wilson-coefficient ODE.

**Flavor scope:** the component expression is a symbolic
one-generation/group-theory diagnostic. Production neutrino calculations
use full complex \(3\times3\) lepton-heavy Yukawa matrices. The diagnostic
does not substitute an arbitrary scalar, matrix element, trace, or norm
for a Yukawa matrix. A numerical LLSS norm requires a separately
verified full-flavor lift of the component Wilson RGE.

**Power counting:** LLSS self-running is a one-loop correction,
\(\delta C_{\rm LLSS}^{(1)}=\mathcal O(\hbar)\). Feeding that correction
through an additional scalar loop produces
\(\delta C_5\sim\hbar\,\delta C_{\rm LLSS}^{(1)}
=\mathcal O(\hbar^2)\). It is therefore excluded from the
authoritative one-loop final matching result,

\[
C_5(M_S)=C_{5,\rm hard}+\Delta C_{5,\rm direct}.
\]

The numerical plot `intermediate_direct_weinberg_running.png` instead
shows the full-flavor direct one-loop \(LLSS\to O_5\) contribution that
enters the final \(C_5\). Do not label the symbolic LLSS fixed-order
diagnostic as the physical Weinberg coefficient.
'@
# Stage both complete replacement files before changing the repository.
$staging = Join-Path $root ".info_cleanup_staging"
if (Test-Path $staging) { throw "Staging directory exists: $staging" }
New-Item -ItemType Directory -Path $staging | Out-Null
try {
  Set-Content -LiteralPath (Join-Path $staging "INFO_T3_MULTIPLICITY.md") -Encoding UTF8 -Value ($mult.TrimEnd() + "`r`n" + $multExtra)
  Set-Content -LiteralPath (Join-Path $staging "INFO_RGE.md") -Encoding UTF8 -Value ($rge.TrimEnd() + "`r`n" + $rgeExtra)
  # Preserve the two troubleshooting documents and the two absorbed source notes.
  foreach ($p in @(
    "INFO_T3_BENCHMARK_VALIDATION.md", "INFO_LLSS_DIAGNOSTIC.md",
    "INFO_T3_NUMERICAL_RESUME.md", "INFO_T3_NUMERICAL_SCENARIO_FIX.md"
  )) {
    $dest = Join-Path $archive $p
    if (Test-Path $dest) { throw "Archive destination already exists: $dest" }
  }
  Copy-Item -LiteralPath $multPath -Destination (Join-Path $staging "INFO_T3_MULTIPLICITY.original.md")
  Copy-Item -LiteralPath $rgePath -Destination (Join-Path $staging "INFO_RGE.original.md")
  Copy-Item -LiteralPath (Join-Path $staging "INFO_T3_MULTIPLICITY.md") -Destination $multPath -Force
  Copy-Item -LiteralPath (Join-Path $staging "INFO_RGE.md") -Destination $rgePath -Force
  foreach ($p in @(
    "INFO_T3_BENCHMARK_VALIDATION.md", "INFO_LLSS_DIAGNOSTIC.md",
    "INFO_T3_NUMERICAL_RESUME.md", "INFO_T3_NUMERICAL_SCENARIO_FIX.md"
  )) {
    Move-Item -LiteralPath (Join-Path $root $p) -Destination (Join-Path $archive $p)
  }
  Write-Host "Updated INFO_T3_MULTIPLICITY.md and INFO_RGE.md; archived four source notes."
  Write-Host "Review: git diff -- INFO_T3_MULTIPLICITY.md INFO_RGE.md"
  Write-Host "Review moved files: git status --short"
}
finally {
  Remove-Item -LiteralPath $staging -Recurse -Force -ErrorAction SilentlyContinue
}
