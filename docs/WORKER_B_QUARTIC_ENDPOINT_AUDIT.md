# Worker B: read-only A(-4)/D(-2) quartic endpoint audit

**Files delivered** (complete, additive files):
- `scripts/RunT3WorkerBQuarticAudit.py`
- `tests/test_t3_worker_b_quartic_audit.py`
- `docs/WORKER_B_QUARTIC_ENDPOINT_AUDIT.md`

This uses repository `Numerical.running.T3CoupledQuarticContext.CandidateQuarticContext` to verify the cached Batch3B response provenance (SHA-256 hashes of response tensor, registry and seed) and evaluate the 12^4 fixed-real-component quartic tensor at the **running EFT1@muS** state. Both S1/S2 remain in EFT1; no scalar has been integrated out here. Production solver and saved results are unchanged.

It audits all scalar index compositions H4, H3S, H2S2, HS3, S4, all quartic response-direction names, tensor Bose permutation symmetry, and the symmetry of every affine response column. With explicitly supplied **independent, same-scale and same-basis** source-HHSS NPZ, it records discrepancy from the cached candidate HHSS and confirms six-orbit source replacement and exact non-HHSS preservation. With optional independent direct quartic NPZ, compares all 5 tensor sectors and fails the audit if the direct-vs-cached relative difference is >2e-9. This does **not** independently establish all UV invariants or native/RGBeta normalization.

The source-derived HHSS is not automatically synthesized, because it depends on the source Matchete contraction and basis at muS. This tool will not pretend cached HHSS equals source HHSS. Without optional inputs the report flags `source_hhss_evaluated=false` and `direct_quartic_evaluated=false`.

## PowerShell

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive "$HOME\Downloads\T3_Worker_B_Quartic_Endpoint_Audit.zip" -DestinationPath . -Force
python -X utf8 -m unittest discover -s tests -p "test_t3_worker_b_quartic_audit.py" -v
python -X utf8 -m scripts.RunT3WorkerBQuarticAudit --find-responses --models A_-4,D_-2
```

Copy the returned exact response root into `--response-root`. The command below uses **your verified unambiguous `all16_v1` files** and will write a new report if the Batch3B cache is present:

```powershell
python -X utf8 -m scripts.RunT3WorkerBQuarticAudit `
  --models A_-4,D_-2 `
  --response-root "PATH_FROM_FIND_RESPONSES" `
  --state-file "A_-4=output\general_observable\all16_v1\A_m4\running.json" `
  --state-file "D_-2=output\general_observable\all16_v1\D_m2\running.json" `
  --out "output\audits\workerb_quartic_endpoint_v1\report.json"
```

If independently evaluated tensors are available in NPZ files, add:

```powershell
  --source-hhss "A_-4=path\to\A_hhss.npz" `
  --source-hhss "D_-2=path\to\D_hhss.npz" `
  --direct-quartic "A_-4=path\to\A_direct.npz" `
  --direct-quartic "D_-2=path\to\D_direct.npz"
```

Default keys: `hhss` (4,4,8,8) and `quartic` (12,12,12,12). Change with `--hhss-key`, `--direct-key`. These arrays must use the same unbroken real scalar component basis and the **muS=1000 GeV** running state. A mismatched field-basis comparison is not physically meaningful.

If `--find-responses` finds no candidates: the Batch3B output is absent or elsewhere; use `--response-root` explicitly when the output exists. The script does not generate the source quartic. Existing report files are never overwritten. The `all16_v1` running input may be in a different scalar coupling coordinate convention from the Batch3B reference; missing directions are rejected rather than replaced by zero.

## Verification

Seven synthetic tests executed locally; all passed. No saved Batch3B tensors or same-scale source HHSS files are present in the supplied ZIPs, hence actual endpoint quartic residuals have **not** been measured here.
