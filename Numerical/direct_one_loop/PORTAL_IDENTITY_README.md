# Native Matchete ↔ RGBeta Higgs-portal identity geometry (diagnostic)

## Project files

- `scripts/ResolveT3DirectUVPortalIdentity.py` — new full-tensor geometric audit.
- `tests/python/test_t3_direct_one_loop_portal_identity.py` — six new synthetic regression checks.
- `Numerical/direct_one_loop/PORTAL_IDENTITY_README.md` — this note.

No production file is replaced. No EFT, RGE, matching, threshold, frozen benchmark, or provisional neutrino mass file is edited.

## Why

Previous `scripts/AuditT3DirectUVPortalBasis.py` reports only the **first** native Matchete `lambdaH{1,2}Inv1` projection at the neutral Higgs VEV, not the entire native invariant span. Its conditional `max_abs_difference_from_adapter_gev2` assumes `lambdaH{which}Inv1 == lambdaH{which}` with all other native invariants zero; it is not an observed cross-package physical discrepancy.

RGBeta `RGE/running/rgbeta/T3RGBetaModel.wl` defines basic portals using the invariant

\[T_\mathrm{id}(i,j,a,b)=\delta_{ij}\delta_{ab},\]

ordered `(Bar[H], H, Bar[S], S)`.

Matchete's `Lagrangian/interactions/ScalarPotential.wl` uses `1/2 PlusHc[lambdaHInvk * fields * CG]` with real couplings. For each native CG tensor \(T_k\), this script forms \(P_k=(T_k+T_k^\dagger)/2\), where dagger swaps the two `H` indices and the two `S` indices and complex conjugates. It solves the **real** linear system:

\[\sum_k c_k P_k = T_\mathrm{id}.\]

This tests the entire four-index tensor, not just \(\langle H^0\rangle\), so it is stronger than comparing VEV projector diagonals. An exact solution yields a *candidate geometric* native invariant combination for the RGBeta identity. It **does not** fix the physical coupling dictionary's overall sign/normalization, the Majorana/Dirac spinor phases, or other scalar quadratic operators.

The script requires the original 16-model native portal CG registries and the previous source charge/portal audits; it checks the previously projected diagonals before computing the full-tensor identity. It refuses to overwrite an existing different result.

## Commands (PowerShell, repository root)

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline

python -m unittest discover -s tests/python -p 'test_t3_direct_one_loop*.py' -v

$base = 'output/general_observable/lo_direct_one_loop'
python scripts/ResolveT3DirectUVPortalIdentity.py `
  --project-root . `
  --portal-root "$base/source_portal_cg_all16" `
  --audit-dir "$base/source_audits"

$r = Get-Content "$base/source_audits/native_portal_identity_decomposition.json" -Raw | ConvertFrom-Json
$rows = foreach ($model in $r.models) {
  foreach ($s in $model.scalars) {
    [PSCustomObject]@{
      Model = $model.model
      Scalar = $s.scalar
      FullTensorSpan = $s.full_tensor_identity_in_real_native_span
      Coefficients = ($s.native_coupling_coefficients_for_unit_identity -join ', ')
      Residual = ('{0:E3}' -f $s.max_abs_four_index_residual)
    }
  }
}
$rows | Format-Table -AutoSize
```

Expect 16 models (32 scalar portal sectors). The result can legitimately report `BLOCKED` for some sectors; do not change benchmark scalar masses or declare physical normalization based on the least-squares coefficients alone. If a full-tensor span is blocked, first validate the Matchete representation-component mapping and the native CG Hermitian action.

## Output

`output/general_observable/lo_direct_one_loop/source_audits/native_portal_identity_decomposition.json`.

Contains SHA256 source provenance, 32 full-index residuals, numerical real candidate coefficients, neutral-projector residuals, and explicit uncertified-physics flags.
