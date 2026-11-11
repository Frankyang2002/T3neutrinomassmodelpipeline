# Worker B: EFT1 scalar quadratic endpoint audit

This is a **read-only** A(-4)/D(-2) numerical audit using **actual local**
`output/general_observable/all16_v1/<model>/running.json` states. It uses the
installed `Numerical/running/T3ScalarRunningSourceBridge.py` from the previous
Worker B delivery. Nothing is recomputed in UV, EFT1, or SMEFT and historical
numerical results are not edited.

The comparison is two independent constructions of the quadratic spectrum:

- A(-4): `mS1Sq` twice and `mS2Sq` six times.
- D(-2): `mS1Sq` six times and `mS2Sq +/- 2*abs(muT3GeneralS2Sq)`.

The matrix eigenvalues, basis orthogonality, diagonalization, and mass trace
are checked. For D(-2), the native convention is
`V2=m2 Sdag S-[mu (Sdag)^2 + h.c.]` and the real Hessian is
`[[m2-2Re(mu),-2Im(mu)],[-2Im(mu),m2+2Re(mu)]]`.

**What this does not check:** source HHSS, quartic invariant completeness,
Matchete-to-RGBeta normalization, all 16 scalar models, or C5 matching.
The two-model cached quartic and HHSS checks remain separate and conditional.

## Install (PowerShell)

From `C:\thesistemp\T3neutrinomassmodelpipeline`:

```powershell
Expand-Archive -LiteralPath "$HOME\Downloads\T3_Worker_B_Endpoint_Audit.zip" -DestinationPath . -Force
python -X utf8 -m unittest discover -s tests -p "test_t3_worker_b_endpoint_audit.py" -v
```

## Run on local saved endpoints

```powershell
$run = 'workerb_mass_endpoint_20261112_v1'
python -X utf8 -m scripts.RunT3WorkerBEndpointAudit `
  --models A_-4,D_-2 `
  --out "output\audits\$run\report.json"
```

Use a **new run ID** if that report exists. The script refuses to overwrite
existing files. It only reads existing running snapshots. The report includes
model-by-model analytic and numerical spectra (GeV^2), residuals, and tachyon
flags. Status `PassedInternalQuadraticChecks` does not certify physics.

Optional: to compare independently saved mass-squared spectra, supply an NPZ
containing arrays named `A_-4` and `D_-2` (each eight real numbers in GeV^2):

```powershell
python -X utf8 -m scripts.RunT3WorkerBEndpointAudit `
  --reference-npz 'C:\path\to\mass_spectra.npz' `
  --out 'output\audits\workerb_mass_endpoint_historical_v1\report.json'
```

For the complete quartic comparison, an independently evaluated, same-scale
`HHSS (4,4,8,8)` array and the Batch3B cached response are still required.
