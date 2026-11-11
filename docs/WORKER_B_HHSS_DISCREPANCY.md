# Worker B: diagnosing the observed HHSS mismatch

The previous `RunT3WorkerBSourceHHSSAudit` measured `1.96043697` (A_-4) and
`1.88489560` (D_-2) for `||Q_candidate-Q_source||/||Q_source||`.
Zero Bose exchange residual does **not** imply compatible native/RGBeta normalization.

This independent, read-only diagnostic uses the same cached quartic responses,
selected running endpoints and already exported source HHSS arrays. It reports:

- Global least-squares coefficient `alpha = <Q_candidate,Q_source>/||Q_source||^2`.
- Shape residual `||Q_candidate-alpha Q_source||/||Q_source||`.
- Sign-reversal diagnostic `||-Q_candidate-Q_source||/||Q_source||`.
- Separate S1S1, S1S2, S2S1, S2S2 and 10 independent Higgs-component blocks.
- Largest componentwise discrepancies, zero-based indices and S1/S2 labels.
- Saved source tensor SHA-256. Neither input is changed.

**Interpretation:** A relative difference near 2 *could* arise from an overall
minus sign; the reported best-scale residual determines whether it actually
explains the tensor. Different block scale coefficients signal nonuniform
basis/normalization disagreement. Do not apply fitted scales or signs in the
production matching; the underlying physical convention is not certified.

Run from `C:\thesistemp\T3neutrinomassmodelpipeline`:

```powershell
python -X utf8 -m unittest discover -s tests -p "test_t3_worker_b_hhss_discrepancy.py" -v
python -X utf8 -m scripts.RunT3WorkerBHHSSDiscrepancy `
  --models A_-4,D_-2 `
  --response-root "output\audits\batch3b_quartic_fixed_20261111_234059" `
  --state-file "A_-4=output\general_observable\all16_v1\A_m4\running.json" `
  --state-file "D_-2=output\general_observable\all16_v1\D_m2\running.json" `
  --source-dir "output\audits\workerb_source_hhss_v1" `
  --out "output\audits\workerb_hhss_discrepancy_v1\report.json"
```

Required existing files: Batch3B cache, both `running.json` files, and
`output/audits/workerb_source_hhss_v1/{MODEL}_source_hhss_native.npz`.

This script **does not** integrate fermions, integrate scalars, run the production
solver, or establish complete physical quartic operator inventories.
