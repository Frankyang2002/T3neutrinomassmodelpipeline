# Worker B: HHSS coupling-direction audit

**Project files**
- `scripts/RunT3WorkerBHHSSDirectionAudit.py` (new)
- `tests/test_t3_worker_b_hhss_direction.py` (new)
- `docs/WORKER_B_HHSS_DIRECTION_AUDIT.md` (new)

The previous A(-4) and D(-2) endpoint residuals of ~1.96 and ~1.88, and their best-scale shape residuals ~0.39 and ~0.63, show a global sign is insufficient. This implementation evaluates the derivative of source HHSS with respect to each running `lambda*` real/imaginary coordinate, **including directions not present in Batch3B**, at the same 1 TeV endpoint. It compares those source derivatives with cached Batch3B columns. It also checks a separate multi-coupling probe to reject a non-affine source evaluator.

The source evaluator uses the existing `GeneralT3SequentialLeadingLog.scalar_tensors` with Matchete native quartic channels, portal completion and source extra RGBeta HHSS geometry. D(-2) is transformed from the source's special boundary coordinates back to the fixed native real basis using the existing `build_special_basis` map. This native/RGBeta mapping is *candidate*, not independently certified.

No production module is changed. This code reads the saved RGBeta endpoint, cached Batch3B tensor, source release and quartic registry. It writes one fresh JSON report only.

## PowerShell

Run from `C:\thesistemp\T3neutrinomassmodelpipeline` after expanding the ZIP:

```powershell
python -X utf8 -m unittest discover -s tests -p "test_t3_worker_b_hhss_direction.py" -v

python -X utf8 -m scripts.RunT3WorkerBHHSSDirectionAudit `
  --models A_-4,D_-2 `
  --response-root "output\audits\batch3b_quartic_fixed_20261111_234059" `
  --state-file "A_-4=output\general_observable\all16_v1\A_m4\running.json" `
  --state-file "D_-2=output\general_observable\all16_v1\D_m2\running.json" `
  --out "output\audits\workerb_hhss_direction_v1\report.json"
```

Interpret `candidate_over_source_scale`: approximately -1 with small shape residual indicates a sign error *in that direction*. Values with large shape residual mean a basis/tensor-geometry mismatch beyond a scalar rescaling. `source_only_directions` lists running quartic coordinates missing from the cached Batch3B response inventory. `source_affine_interpolation_relative` measures direct independent source evaluation at a multi-axis perturbation. If this fails, execution stops rather than reporting misleading response derivatives.

The response is **not** an independently certified physical tensor, and no automatic reweighting or matching is performed. The benchmark must be run on the user's local repo: its output data are not present in the supplied ZIPs.
