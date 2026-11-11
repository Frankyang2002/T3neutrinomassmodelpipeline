# Worker B — combined source-HHSS reconciliation batch

Project file: `scripts/RunT3WorkerBBatchReconcile.py`.

Implements in **one local run**:

1. Load hashed Batch3B affine candidate and the actual EFT1 endpoint.
2. Re-evaluate all source HHSS terms using existing Matchete native channels, completed portal and general extra invariant machinery.
3. Handle D(-2) using the existing source scalar mass-basis rotation, back to fixed native real coordinates.
4. Audit all coupling real/imag directions, including independent multi-parameter affine checks.
5. Report largest tensor-component disagreements, model-correct species-pair blocks (A: 2+6, D: 6+2), source/candidate direction overlap Gram matrix, best unconstrained linear coordinate map, and candidate out-of-source-span residuals. This is **diagnostic only**; does not apply fitted transforms.
6. Export a rank-(12,12,12,12) hybrid tensor where source HHSS replaces all six HHSS placements, with every non-HHSS entry preserved exactly.
7. Save per-model NPZ files plus source provenance and a single JSON report.

The full output is **candidate-only** because non-HHSS directions and native/RGBeta normalization remain unverified. No physical C5 or matched scalar Wilson coefficient is produced.

## PowerShell

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive "$HOME\Downloads\T3_Worker_B_Batch_Reconcile.zip" -DestinationPath . -Force
python -X utf8 -m unittest discover -s tests -p "test_t3_worker_b_batch_reconcile.py" -v
python -X utf8 -m scripts.RunT3WorkerBBatchReconcile `
  --models A_-4,D_-2 `
  --response-root "output\audits\batch3b_quartic_fixed_20261111_234059" `
  --state-file "A_-4=output\general_observable\all16_v1\A_m4\running.json" `
  --state-file "D_-2=output\general_observable\all16_v1\D_m2\running.json" `
  --out-dir "output\audits\workerb_batch_reconcile_v1"
```

Use a new `--out-dir` if that directory exists; the program refuses to overwrite.

Outputs: `<out-dir>/report.json`, `<out-dir>/A_-4_reconciliation.npz`, `<out-dir>/D_-2_reconciliation.npz`. NPZ fields `source_hhss`, `source_direction_columns`, `candidate_hhss`, `hybrid_full_tensor`. The tool does not change production solver, saved historical results, or GitHub.

**Dependencies:** packages from earlier Worker B ZIPs, especially `scripts/RunT3WorkerBHHSSDirectionAudit.py`, `scripts/RunT3WorkerBSourceHHSSAudit.py`, `scripts/RunT3WorkerBEndpointAudit.py`, and `Numerical/running/T3CoupledQuarticContext.py`. Prior installations are assumed.
