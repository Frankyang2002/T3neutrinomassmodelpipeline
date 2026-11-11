# Worker A — fermion-threshold operator inventory and subtraction gate

Install `scripts/WorkerAFermionSubtractionGate.py` into the project, preserving its path. No existing production files are replaced.

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
$zip = Get-ChildItem "$env:USERPROFILE\Downloads\worker_a_fermion_subtraction_gate*.zip" | Sort-Object LastWriteTime -Descending | Select-Object -First 1
if (-not $zip) { throw "Worker A ZIP not found" }
Expand-Archive -LiteralPath $zip.FullName -DestinationPath . -Force
python -m py_compile scripts\WorkerAFermionSubtractionGate.py
$runId = Get-Date -Format 'yyyyMMdd_HHmmss'
python scripts\WorkerAFermionSubtractionGate.py --sources-zip 'output\audits\worker_a_sources_20261112_003127.zip' --output "output\audits\worker_a_f_source_gate_$runId.json"
```

This script reads *only* the existing extracted first-fermion-threshold Matchete source ZIP, checks the complete source/ledger/Laurent SHA256 chain, enumerates the full operator inventory with per-term SHA256, and marks terms containing prohibited independent conjugate Yukawas. This is *not* an amplitude-level subtraction. It deliberately fails closed when UV-minus-EFT1-soft subtraction, counterterms, normalisation, and associated checks have not been implemented or verified.

Verified on the uploaded source ZIP:

| Model | Raw LLSS | Restricted LLSS | Raw LLHH | Restricted LLHH | Blockers |
| --- | ---: | ---: | ---: | ---: | ---: |
| A(-4) | 48 | 48 | 6 | 6 | 8 |
| D(-2) | 352 | 96 | 12 | 6 | 8 |

`restricted` means source term contains neither `y1Conj` nor `y2Conj`. It does **not** mean the term is renormalised or can be added directly to an EFT1 Wilson boundary.

Important input inconsistencies verified from source metadata: `matchete_loop_function_evaluator_called=false`, `one_loop_hbar_linearity_exact=false`, and `unreconciled_tree_mass_coordinate_identity_verified=false`. A full finite threshold coefficient cannot currently be inferred. Those gates reflect the audit's declared metadata; they do not establish that Matchete was never used elsewhere.

No changes to running, scalar threshold, or coupled solver. Intended output `output/audits/worker_a_f_source_gate_<runId>.json`.
