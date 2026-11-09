# Workstream A v4 — source-conditioned scalar quartic-action audit

**Status:** `Conditional_RGBetaRealActionAlgebraPassed_NativeScalarBasisUnverified`.
This report distinguishes *algebraic action normalization in a declared scalar frame* from *native RGBeta projectors / physical coupling normalization*. No shared project sources or GitHub repository are changed.

## Source checks (read-only GitHub main SHA `a944e0030b8ef24b23342d99e0d8f33cef3cffd5`)

* `RGE/running/rgbeta/T3RGBetaModel.wl`, `T3RGBetaAddBasicQuartics`: explicitly constructs `lambdaH`, `lambdaS1`, `lambdaS2` with invariant coefficient `1/2`, and `lambdaH1`, `lambdaH2`, `lambda12` with coefficient `1`. It declares all six real. These are the **six simple norm quartics only**; it is **not** the full baseline portal inventory.
* Upstream RGBeta `Package/FieldsAndCouplings.m`: `AddQuartic` sets `normalization = 24` and multiplies by `1/Sqrt[2]` for each of the four complex scalar legs, therefore `6*Matrix[coupling][...]` for the exported 50 complex extra quartics. `CouplingPermutations` computes unique arrangements, `AvgPermutations` weights them by their reciprocal, and `Lam` adds a conjugate for `SelfConjugate -> False` but not for `True`.
* `Lagrangian/interactions/ScalarPotential.wl`, `T3HiggsPortalHermitianPhase`: calculates the Hermitian-conjugate symmetry of the *native* Matchete CG and multiplies an anti-Hermitian portal by `I` **inside** `1/2 PlusHc[...]`. This v4 check **does not** extract that native portal tensor or replace this source.
* Actual uploaded source exports: 50 exact Matchete tensors and 50 conditional reconstructed RGBeta expressions, all 16 models. Record-level provenance SHA256 hashes are in the output.

## Independent real-action algebra

For real Cartesian components `x` and a complex scalar `z_a=(x_{a,R}+i*x_{a,I})/sqrt(2)`, let `P` denote the *unnormalized* projection `P x = x_R + i*x_I`. For each extra complex invariant tensor `T_{ijkl}` with coefficient `c`, construct

```
L_abcd = 6 * Sym4[ c T_ijkl P_(a,i) P_(b,j) P_(c,k) P_(d,l)
                     + conj(c T_ijkl P_(a,i) P_(b,j) P_(c,k) P_(d,l)) ]
V_real = L_abcd x_a x_b x_c x_d / 24
V_complex = c I(z) + conj(c I(z))
```

The factor `6` is **the actual source-exported coupling-function coefficient**, and `Sym4` follows RGBeta's `AvgPermutations` convention, explicitly checked for each `native_unique_arrangements` list. The identity

`6 * (sqrt(2))^4 / 24 = 1`

then reconstructs a **unit prefactor** for the complex invariant in this declared component frame. This is checked numerically with random *complex coefficients*, not just real couplings, and with independent random *real* scalar inputs. The calculation expands the 4-tensor in real coordinates rather than reusing the previous v3 raw polynomial holdout.

For the six real `SelfConjugate -> True` baseline norm vertices, the same calculation omits a second conjugate term and verifies the `1/2` or `1` invariants against their explicit source scalar potentials.

Results on the **uploaded evidence**:

* **50/50** extra complex quartics in 16 models; max relative real-action residual **2.565e-15**.
* **96/96** simple norm-quartic model-coupling combinations (`6 × 16`) from their source-defined polynomials; residual reported in the output.
* **38/38** installed `test_t3_parallel_a*.py` tests passed from the combined prior A workspace and this v4 batch.

## What is still blocked

This deliberately **does not execute the native RGBeta `sDelS`/`$scalarContraction`/`Lam` objects** in Wolfram, so it is not an independent native projector normalization check. It works in the **conditional Matchete-aligned Cartesian triplet frame** used in v3, not in a separately verified native RGBeta frame. The action audit establishes consistency with the upstream `6*Matrix` and `AvgPermutations` conventions **given that coordinate identification**; it does not establish the physical Matchete↔RGBeta coupling map.

The six baseline quartics covered are **only norm-type** couplings. The remaining independent Higgs-adjoint portal operators, `lambda12Cross`, special RG-closed quartics, and especially A(-4) `lambdaH2Inv1` with its anti-Hermitian CG phase are not compared numerically here. Their native Matchete CGs and RGBeta component tensors require additional export/evidence. Signed GammaCC→Weyl and the complete one-loop ADM are **also not certified**. No physical Wilson results should be promoted or recomputed from this audit alone.

## Installation, paths and PowerShell

The ZIP contains **four full files**; install them in the named project-relative paths:

1. `Numerical/parallel_t3/wilson/scalar_action_audit.py`
2. `Numerical/parallel_t3/wilson/SCALAR_ACTION_AUDIT_V4.md`
3. `scripts/RunT3ParallelARealQuarticActionAudit.py`
4. `tests/test_t3_parallel_a_scalar_action_v4.py`

From `C:\thesistemp\T3neutrinomassmodelpipeline`:

```powershell
$zip = "C:\path\to\t3_workstream_A_scalar_action_audit_v4.zip"
Expand-Archive -Path $zip -DestinationPath . -Force
.\.venv\Scripts\python.exe -m unittest discover -s tests -p "test_t3_parallel_a*.py" -v

.\.venv\Scripts\python.exe -m scripts.RunT3ParallelARealQuarticActionAudit `
  --matchete output\parallel_t3\wilson\matchete_native_tensors.json `
  --rgbeta output\parallel_t3\wilson\rgbeta_abstract_components_conditional_v2.json `
  --output output\parallel_t3\wilson\real_quartic_action_conditional_v4.json
```

The script refuses to overwrite an existing result. No B/C-owned files or existing A(-4) source are touched.

**Next input required to remove the remaining baseline/portal block:** native Matchete CG components for all baseline independent portals and native RGBeta `Lam` projector evidence, especially A(-4) anti-Hermitian Higgs-triplet portal. The existing extra-quartic 50-slot export does not contain those coefficients.
