# Worker A — scalar-soft source evidence collection

This package adds a read-only source collector, **not a physical matching module**.
The available A(-4) direct Laurent calculation supports formal cancellation of the
hard and scalar-soft pole only in its shared ordered operator channel. It cannot
certify a renormalised EFT1 Wilson coefficient, and does not fix D(-2)'s scalar
mass-basis/operator normalization. Do not add the formal remainder to EFT1.

## Install in repository

- `scripts/WorkerACollectScalarSoftEvidence.py`

## Run (PowerShell, from project root)

```powershell
$runId = Get-Date -Format 'yyyyMMdd_HHmmss'
python scripts\WorkerACollectScalarSoftEvidence.py --repo . --output "output\audits\worker_a_scalar_soft_$runId.zip"
```

Upload the resulting ZIP. It contains matched InputForm text files and SHA256
provenance, if existing on disk. If none is found, the exact scalar-soft
Matchete export will need to be regenerated locally by the existing Wolfram
matching workflow; the lack of an export does not establish a vanishing soft
amplitude.
