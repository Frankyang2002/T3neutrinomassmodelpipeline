# Workstream B: B(-1)/C(-1) mixed-bilinear verifier v6

## Scope

Fix the *Python verifier only* after the successful live RGBeta v5 run reported
one `SparseArray`, 295 native rules and 24 mixed-mass coupling witnesses
for each model. Do not change RGBeta, Matchete, the generated RGE or the
native evidence files.

## Verified source-level cause

The original `scripts/VerifyT3BCMixedBilinearMap.py` required:

```python
coefficient == '12*Matrix[muT3GeneralS1S2][]'
```

But the live probe stores `ToString[coupFn[], InputForm]`, i.e. the **evaluated**
zero-index scalar-mass coupling function. In RGBeta v1.2.0:

- `Package/FieldsAndCouplings.m`, in `AddScalarMass`, sets
  `normalization = 24` and multiplies by `1/Sqrt[2]` for each complex scalar,
  producing exactly `24 * (1/Sqrt[2])^2 = 12`.
- `Package/FieldsAndCouplings.m` assigns the scalar-mass coupling using
  `Evaluate[normalization Matrix[coupling] @ ##] &`.
- `Package/GroupsAndIndices.m` explicitly defines `Matrix[m_][] = m`.

Therefore the resulting InputForm is `12*muT3GeneralS1S2` (or with the
explicit `Global\`` context qualifier), not an unevaluated `Matrix` call.

The v6 verifier accepts only these exact evaluated representations. It rejects
any other factor (such as 6 or 24), negative sign, conjugated symbol, other
coupling, or unresolved `Matrix` call. It still rejects missing/blocked live
reports, invalid source SHA256, unexpected SU(2) tensors and any failed
independent action/Hessian/spectrum consistency check.

## Files in ZIP

- `scripts/VerifyT3BCMixedBilinearMap.py`: complete replacement.
- `tests/test_t3_bc_mixed_bilinear_map.py`: complete replacement with native
  evaluation regression tests.
- `docs/WORKSTREAM_B_BC_MIXED_BILINEAR_VERIFIER_FIX_V6.md`: this file.

All files are relative to the local repository root
`C:\thesistemp\T3neutrinomassmodelpipeline`.

## Installation and execution in PowerShell

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "C:\path\to\Workstream_B_BC_MixedBilinear_Verifier_Fix_v6.zip" -DestinationPath . -Force
python -X utf8 -m pytest -q tests\test_t3_bc_mixed_bilinear_map.py tests\test_t3_bc_mixed_bilinear_rgbeta_probe_contract.py
python -X utf8 -m scripts.VerifyT3BCMixedBilinearMap `
  --matchete output\audits\workstream_b_bc_mixed_matchete_v1.json `
  --rgbeta output\audits\workstream_b_bc_mixed_rgbeta_v5.json `
  --output output\audits\workstream_b_bc_mixed_bilinear_map_v6.json
```

The v5 native RGBeta JSON is reused unchanged; **do not rerun Wolfram**.
The verifier refuses to overwrite an existing v6 report.

## Interpretation

A successful result establishes a conditional identical-component
`mu_RGBeta/mu_Matchete` map, with Matchete's **live registered** CG tensor,
RGBeta's live massive `Lam` coupling witnesses, and independently tested
complex-field action, Hessian and spectrum. It **does not** independently
identify Matchete's and RGBeta's SU(2) component/generator bases or certify
full matching. `physical_mapping_certified` correctly remains `false`.

## Local validation

`26 passed` using the two listed test modules with the bundled replacement
files. No installed Wolfram kernel is available in this environment, so this
revision has not been tested against the user's actual local v5 JSON. The
expected evaluated coefficient is confirmed from the RGBeta source.
