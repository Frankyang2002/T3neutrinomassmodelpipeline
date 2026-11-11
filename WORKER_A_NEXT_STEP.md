# Worker A: numerical one-loop matching source audit

## Verified upstream finding

`Numerical/parallel_t3/complete/a4_direct_hard_soft.py` already contains an
A(-4) *formal ordered-operator* direct hard/scalar-soft pole cancellation.
Its formal finite remainder is **not** a matched EFT1 Wilson because the
full operator basis, counterterms, CG and GammaCC-to-Weyl conventions remain
uncertified. It is not valid to paste it into `T3FermionOneLoopBoundary.py`.

For the intended restricted theory keep only Y1 and Y2. D(-2) generated
Matchete hard expressions contain an additional y2Conj channel. The restricted
UV Lagrangian and EFT1 operator inventory must be checked **before** discarding
this channel or claiming the D(-2) pole cancellation.

## Source export from the user's actual local copy

Run at repository root in PowerShell:

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
$runId = Get-Date -Format 'yyyyMMdd_HHmmss'
python scripts\ExportWorkerAMatchingSources.py --repo . `
  --output "output\audits\worker_a_sources_$runId.zip"
```

Upload the generated `output\audits\worker_a_sources_<runId>.zip` for an exact
local code-and-output audit. This is a **read-only** exporter: it does not
alter any historical results or the production coupled solver. The ZIP may
include local numerical benchmark inputs; inspect before sharing if needed.

## Remaining certification steps

1. Confirm the exact restricted-Yukawa UV Matchete source for A(-4) and D(-2).
2. Project all generated relevant dimension-five operators using native
   SU(2), flavour, real-scalar and GammaCC conventions.
3. Subtract the EFT1-soft amplitudes and CT contributions operator by operator.
4. Demonstrate pole cancellation and determine finite boundary tensors at muF.
5. Validate A(-4) numerical benchmark and then D(-2), keeping tree and
   finite one-loop tensors separate.

No new physical one-loop coefficient is claimed in this package.
