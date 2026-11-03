# Track B A/D/E: one nine-model batch

## Installation (PowerShell)

The single ZIP contains **new, uniquely named Track B files** only. It does not modify Matchete source, numerical production, EFT/RGE, GitHub, or any prior result. Extract it into the project root (change the ZIP path as required):

```powershell
Set-Location 'C:\thesistemp\T3neutrinomassmodelpipeline'
$zip = Join-Path $env:USERPROFILE 'Downloads\Track_B_ADE_Nine_Model_Batch.zip'
if (-not (Test-Path -LiteralPath $zip)) { throw "ZIP not found: $zip" }
git rev-parse HEAD
git status --short
Expand-Archive -LiteralPath $zip -DestinationPath (Get-Location).Path -ErrorAction Stop
$py = if (Test-Path '.\.venv\Scripts\python.exe') { '.\.venv\Scripts\python.exe' } else { 'python' }
& $py -c 'import numpy, scipy, sympy, mpmath; print("Track B dependencies OK")'
if ($LASTEXITCODE -ne 0) { throw 'Required Python dependencies: numpy scipy sympy mpmath' }
& $py scripts\RunTrackBADEBatch.py --repo-root . --run-self-tests
if ($LASTEXITCODE -ne 0) { throw 'Track B mathematical self-tests failed' }
```

## All-nine fresh native Matchete extraction (one kernel per model)

The *first* measurement for a physical source certificate is the fresh native restricted action, CG and quadratic witness. There is **no mutable registry parallelism**:

```powershell
$stamp = Get-Date -Format 'yyyyMMdd_HHmmss'
$native = "output\general_observable\physical_normalisation\parallel_B_ADE_native_$stamp"
$models = @('A_-4','A_-2','D_-2','D_+0','D_+2','E_-4','E_-2','E_+0','E_+2')
New-Item -ItemType Directory -Path $native -ErrorAction Stop | Out-Null
foreach ($model in $models) {
    $slug = $model.Replace('-', 'm').Replace('+', 'p')
    $witness = Join-Path $native "$slug.json"
    wolframscript -file scripts\RunTrackBADENativeAction.wl $model $witness
    if ($LASTEXITCODE -ne 0) {
        Write-Warning "$model native action BLOCKED. Inspect terminal output or $witness; remaining models continue."
    }
}
& $py scripts\RunTrackBADEBatch.py --repo-root . --native-dir $native
if ($LASTEXITCODE -ne 0) { Write-Warning 'One or more model gates blocked; inspect latest nine-model ledger' }
```

Each native process sets `T3_EXPERIMENTAL_GENERAL_FERMIONS=0`, `T3_EXPERIMENTAL_GENERAL_QUADRATICS=0`, `T3_EXPERIMENTAL_GENERAL_QUARTICS=0`, `T3_EXPERIMENTAL_FLAVOR_MASSES=1`, and refuses to overwrite its input/output. It checks `CheckLagrangian`, original+H.c. Yukawa/projector/GammaCC placement, portal/mixing coefficients, canonical scalar kinetics, indexed Dirac `MF`, and exports the actual CG component arrays plus six source SHA256 values. If the indexed flavour option fails, that model is marked `BLOCKED` or conditional as appropriate, rather than given a false diagonal-mass certificate.

## Archived-CG fallback, full nine-model numerical batch

If installed Matchete fails or `wolframscript` is unavailable, the following command still performs a **conditional** nine-model batch from previously exported *actual native* CG tensors and the source fitted-reference inputs. It does not substitute any provisional matrices for native source factors:

```powershell
& $py scripts\RunTrackBADEBatch.py --repo-root .
if ($LASTEXITCODE -ne 0) { Write-Warning 'Inspect blocked model rows in the newest output folder' }
```

For a strict fresh-native-only gate (and no fallback):

```powershell
& $py scripts\RunTrackBADEBatch.py --repo-root . --native-dir $native --native-only
```

Every execution creates a fresh timestamped directory:

`output/general_observable/physical_normalisation/parallel_B_ADE_<YYYYMMDD_HHMMSS_microseconds>/`

It contains exactly nine `*_report.json` files (including explicit `BLOCKED` status where applicable), `batch_summary.json`, `nine_model_ledger.csv`, `BATCH_SUMMARY.md`, and `ROADMAP_RUNTIME_STATUS.md`. The installed roadmap itself is not overwritten; statuses are updated in the new run's roadmap copy. Portal-zero and conditional-identity matrices, charges, CGs, Hessians, eigenvalues, Takagi masses, pole/integral checks, complex negative controls, and comparison-only residuals are included in each unblocked JSON.

The ZIP also contains a **pre-run archived-source numerical snapshot**, with nine complete reports and `ARCHIVED_NATIVE_NINE_MODEL_SNAPSHOT.json`. These were computed in the delivery environment using read-only GitHub repository data; they are **not** claimed as new executed native Matchete witnesses and **not physically certified**. Do not copy them over a later native result. The zero-portal and benchmark-identity assignments are distinct; no A_+0 one-channel result is reused.

## Validation and hard stop

Nine-model numerical evaluation is supported by a source-only Python engine, not the provisional one-loop engine. Physical certification stays `NO` until a model-specific fresh native Matchete action, indexed mass and absolute spinor dictionary, and independently reviewed physical benchmark portal dictionary are accepted. A small frozen comparison residual alone has no certification value. No Track A, Track C, B2 cross-package, EFT, or RGE work is part of this deliverable.
