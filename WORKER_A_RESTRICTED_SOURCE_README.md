# Worker A — restricted raw UV→EFT1 Matchete source

Project file: `scripts/WorkerARestrictedFermionSource.py`.
Generated evidence: `output/parallel_t3/worker_a_restricted_f_source_v1/{A_-4,D_-2}/{restricted_raw_inputform.txt,restricted_manifest.json}`.

The script reads the full original Matchete one-loop operator ledger from the previous source archive and imposes `y1Conj = y2Conj = 0` **only on independent source couplings**. `Bar[Coupling[y1,...]]` and `Bar[Coupling[y2,...]]` survive: these are complex conjugates of retained Yukawas, not extra independent spurions. The operation applies to every additive term, including non-LNV sectors and tree-mass remnants. Every source term's index/hash is retained in the JSON manifest; input SHA256 is verified.

Numerics: A(-4) 135 of 135 retained, LLSS 48, LLHH 6; D(-2) 359 of 1200 retained, 841 zeroed, LLSS 96, LLHH 6. The projected source remains **raw/unrenormalised**. It is neither an EFT1-soft subtraction nor a physical Wilson coefficient. There is no GammaCC/Weyl normalisation, operator mixing, counterterms or pole cancellation.

Run in PowerShell from the project root:

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
python -m py_compile scripts\WorkerARestrictedFermionSource.py
$runId = Get-Date -Format 'yyyyMMdd_HHmmss'
python scripts\WorkerARestrictedFermionSource.py `
  --sources-zip 'output\audits\worker_a_sources_20261112_003127.zip' `
  --output "output\audits\worker_a_restricted_source_$runId.zip"
```

The next physical calculation is to construct the EFT1 loop amplitude from the **restricted tree EFT1 operators**, match against this restricted UV one-loop amplitude, and add a consistent counterterm prescription. Do not feed the restricted raw expression directly into the production boundary.
