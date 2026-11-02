# T3 direct one-loop — integrated source-convention audit (Batch 1)

**Project:** `C:\thesistemp\T3neutrinomassmodelpipeline`

## Purpose and current status

This is a **read-only, additive audit** for the existing 16 restricted fitted-reference benchmark points. It cross-checks the current numerical Weyl component inputs against the existing native Matchete CG evidence and the frozen source benchmark, then performs algebraic scalar-action and fermion-mass checks. It does **not** alter the provisional numerical inputs, RGE, EFT matching, benchmark JSONs, or GitHub.

**It does not constitute a physically normalized UV certificate.** The source can prove that the quadratic scalar mass matrix is consistent with the stated *candidate* interpretation of `+lambdaT3 HH S1 S2†+h.c.` in the Matchete Lagrangian, but without a native, physically anchored spinor and coupling dictionary one must not set `validated_uv_component_conventions=True`.

This package was authored against source files read from the repository's GitHub default branch. The actual local repository was not available in this session and may differ; SHA provenance and source-content guards reject mismatched inputs instead of guessing.

## Files installed

| Project-relative path | Purpose |
|---|---|
| `scripts/AuditT3DirectUVConventionBatch1.py` | Complete standalone read-only all-16 audit executable |
| `tests/test_t3_direct_uv_convention_batch1.py` | Isolated numerical and negative-control tests |
| `docs/T3_Direct_UV_Convention_Batch1_README.md` | This document |

For Batch 1.1, these three previously installed Batch 1 files are replaced in full; no production physics files are modified.

## Install and test (PowerShell)

```powershell
Set-Location 'C:\thesistemp\T3neutrinomassmodelpipeline'
Expand-Archive -LiteralPath "$env:USERPROFILE\Downloads\T3_Direct_UV_Convention_Batch1.zip" -DestinationPath . -Force
python -m unittest discover -s tests -p 'test_t3_direct_uv_convention_batch1.py' -v
```

Or use the actual ZIP path returned by the download, rather than `Downloads`, if it is stored elsewhere.

## Run using the existing native audits and provisional inputs

```powershell
Set-Location 'C:\thesistemp\T3neutrinomassmodelpipeline'
python scripts\AuditT3DirectUVConventionBatch1.py `
    --project-root . `
    --physical-input 'input\lo_direct_one_loop' `
    --audit-dir 'output\general_observable\lo_direct_one_loop\source_audits' `
    --output-dir 'output\general_observable\lo_direct_one_loop\convention_batch1'
```

The audit intentionally returns **exit code 2**, including when all 16 conditional tensor checks pass: physical certification is still blocked. Read the JSON output rather than treating exit code 2 as an accidental failure.

Result:

`output/general_observable/lo_direct_one_loop/convention_batch1/integrated_convention_audit.json`

If the source audits were saved elsewhere, pass their actual directory in `--audit-dir`; do not copy or fabricate missing reports. If `input/lo_direct_one_loop` is absent, generate it only using the repository's existing provisional builder *after verifying the source audits and without changing frozen fitted benchmarks*.

## Exact scope of checks

1. Source version guards on the native topology, field definitions, SU(2) invariants, portal construction and numerical core; component-input SHA provenance tied to the benchmark and source audit files.
2. Reconstruct native charge-basis CG coefficients with an **independent array contraction implementation**, using the *existing* SU(2) generator-candidate covariance result; do not relabel that generator candidate as a live native export.
3. Charge conservation of every nonzero Yukawa, topology VEV-mixing and Majorana charge-pairing channel.
4. For Dirac fermions in the all-left-handed basis, verify the candidate off-diagonal mass matrix
   \[
   M_{\rm Weyl}=\begin{pmatrix}0&M_F\otimes I_d\\M_F^T\otimes I_d&0\end{pmatrix}.
   \]
   For real SU(2) Majorana representations, verify the candidate \(M_{\rm Weyl}=M_F\otimes G\), where \(G=U^TU\) is the charge-pairing metric, symmetric and unitary.
5. For the restricted scalar input, verify \(B=0\) and the Hermitian \(A\) reconstructed from
   \[
   \mathcal L_{\rm mix}=+\lambda_{T3}HH S_1S_2^\dagger+\mathrm{h.c.},\qquad H^0=v/\sqrt2,
   \]
   with
   \[
   (A_{21})_{ts}= -\lambda_{T3}\frac{v^2}{2}(C_{HH S_1S_2^\dagger})_{00st},\quad
   (A_{12})_{st}=(A_{21})_{ts}^*,
   \]
   and **candidate** universal portal shifts \(-\lambda_{H_i}v^2/2\). Both are tested against the expanded polynomial for two complex scalar configurations. Flipping mixing sign, dropping the VEV factor 1/2, and flipping the portal shift are tested as negative controls where the corresponding couplings are nonzero.
6. Reconstruct the existing **conditional** `y1_HC` and `y2_HC` Weyl tensor convention independently and compare to both source-anchored inputs. This tests code and conjugation consistency under the assumed dictionary, **not the physical value of the spinor phase or an absolute Matchete–RGBeta map**.
7. Require all 16 earlier CG/covariance records and all 32 scalar-portal full-tensor identity representations; incorporate, without upgrading, the earlier A(0) conditional charge-conjugation report if it exists.

## Known blocked physical verification gates

- The native abstract Matchete `GammaCC` identities and the chosen explicit 4×4 Weyl gamma-matrix convention do not by themselves tie the native `CConj` phase to the physical input Yukawa matrices.
- The corresponding `y2` phase and Dirac mass sign must be transported using **one common heavy field definition**. The sign of a field-redefined Yukawa by itself has no physical meaning.
- For self-conjugate heavy fermions (`B_-1`, `C_-1`), the full native mass/charge-conjugation term must be projected to canonical two-component fields to settle the Majorana normalization, not merely checked to be symmetric.
- The full-index portal identity already establishes the existence of the universal structure; it does **not** identify which native couplings equal the frozen `lambdaH1` and `lambdaH2` numerically or fix their relative action signs.
- The extra symmetry-allowed interactions excluded from the restricted fitted inputs cannot be silently treated as part of a fully general UV calculation.

**Next implementation gate:** obtain the native Matchete-to-frozen-benchmark 2-component action dictionary (including `y1`, `y2`, `MF`, `lambdaT3` and the portal coefficients) for one consistent canonical field basis, and verify the same dictionary across the 16 native tensor representations. Until that is supplied or demonstrated, keep all computed masses at `uv_mass_evaluated_uncertified_source_mapping`, not `uv_mass_evaluated_source_certified`.

Batch 2 (final all-16 numerical spectra) is deliberately not launched or simulated by this package. It should read corrected inputs only after the physical dictionary is independently established.

## Batch 1.1 — JSON serialization regression correction (2026-11-02)

`candidate_negative_controls` originally returned NumPy boolean scalars
(`numpy.bool_`) for three detected sign/factor mutations. Python's standard
`json.dumps` cannot encode these scalars, causing a `TypeError` at report
writing after the all-model audit had executed. These three values are now
converted to built-in `bool` at their point of creation. The existing test
suite gains a direct regression for the negative controls and a synthetic
16-model report-serialization test. No numerical formulas or physical
certification gates have been changed.

Run from the repository root after extracting this replacement ZIP:

```powershell
python -m unittest discover -s tests -p 'test_t3_direct_uv_convention_batch1.py' -v
python scripts\AuditT3DirectUVConventionBatch1.py `
    --project-root . `
    --physical-input 'input\lo_direct_one_loop' `
    --audit-dir 'output\general_observable\lo_direct_one_loop\source_audits'
Get-Content 'output\general_observable\lo_direct_one_loop\convention_batch1\integrated_convention_audit.json'
```

The script still returns exit code 2 by design because physical normalization
is not certified. The completed report is written before that exit.

If only 15 provisional inputs were generated, inspect the builder summary:

```powershell
Get-Content 'input\lo_direct_one_loop\build_status.json'
Test-Path 'input\lo_direct_one_loop\B_-3.json'
```

A missing `B_-3.json` is a separate upstream input-build issue and does not
account for the JSON serialization failure.
