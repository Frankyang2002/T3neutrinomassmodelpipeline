# General T3 LLSS normalization audit (v1)

## Why B_-1 and C_-1 must be regenerated

Source inspection establishes a classification inconsistency:

- `Lagrangian/RunMatching.wl` labels native `Proj[-1]` as `PL` and `Proj[1]` as `PR`.
- The original `Numerical/running/GeneralT3MixedWilsonTensor.py` had these labels reversed.
- Actual synchronized B_-1 (20 terms) and C_-1 (24 terms) contain *both* projectors, with `Proj[-1]` acting on `Transp[Field[l]] ... Field[l]` and `Proj[1]` acting on barred lepton fields.

Therefore earlier v3 mixed `PL` projections were projections of the conjugate source sector. The totals `16/16`, `0 blocked` are not sufficient to certify the *holomorphic PL* selection. All 14 ordinary sources use `Proj[-1]` as expected.

The replacement `GeneralT3MixedWilsonTensor.py` follows the native convention and rejects projector/spinor-field mismatches. To preserve output provenance, regenerate the B/C tensor source and polynomial bridge into new versioned directories.

## What is established by the general audit

For each source term:

1. Verify the source contains the native `GammaCC` PL bilinear and exactly the projector consistent with its tag.
2. Use the already validated lepton-flavour dummy map, not the Yukawa matrix's source axis order.
3. Rebuild the symmetric formal polynomial representative and compare every numeric independent coordinate.
4. Verify orbit multiplicities of 1, 2 or 4, including identical lepton or scalar coordinates.
5. Preserve the separation between holomorphic `PL` and Hermitian-conjugate `PR` terms; do not double the PL boundary.

The source polynomial is

\[
P=\sum_{Ip,Jq,ab} T_{Ip,Jq,ab}\,\psi_{Ip}\psi_{Jq}\phi_a\phi_b.
\]

For an orbit with `m` distinct ordered members, the symmetric representation assigns the *monomial-sum* coefficient divided by `m` to each member. This is an algebraic identity of formal commuting scalar variables and *symmetric placeholders for Weyl fermion bilinears*. No additional identical-field `2!` or `4` factor is introduced.

The master Weinberg adapter's `diff(...)/4` tensor reconstructs its reference polynomial directly. The project `RunMatching.wl` independently declares a reported *C5* convention

\[
\mathcal L_5\supset\tfrac12 C_5^{ij}\mathcal O_{ij}+\mathrm{h.c.},\quad
C_5^{ij}=A_5^{ij}+A_5^{ji}.
\]

**Do not apply this C5 reporting factor to generic LLSS tensors.** It does not establish the relative normalization of `NCM[Transp[Field[l]],DiracProduct[GammaCC,Proj[-1]],Field[l]]` against the master RGE's two-component Weyl bilinear.

## Explicit remaining gap

The absolute Matchete `GammaCC` -> master Weyl bilinear **sign and normalization factor** are not established by the polynomial identity, native PL classification, or the master Weinberg template. The general operator/coupling dictionary and resulting three-family real-linear `A(mu), B(mu)` are also not certified. The batch intentionally reports these as unresolved.

## Install/run (PowerShell from repository root)

```powershell
Expand-Archive -LiteralPath "$HOME\Downloads\t3_general_wilson_normalization_v1.zip" -DestinationPath . -Force
python -X utf8 -m unittest tests.test_t3_general_mixed_wilson_tensor tests.test_t3_general_wilson_normalization_audit -v

# Preserve v3 as the historical (incorrect mixed-chirality) output.
Copy-Item -LiteralPath 'output\general_observable\ordered_wilson_tensor_v3' -Destination 'output\general_observable\ordered_wilson_tensor_v4' -Recurse
python -X utf8 -m scripts.RunT3GeneralMixedWilsonTensor --all --reuse-cg --prior-ordered-root 'output\general_observable\ordered_wilson_tensor_v2' --output-root 'output\general_observable\ordered_wilson_tensor_v4'

# Rebuild from the corrected v4 tensor input; do not reuse polynomial v1.
python -X utf8 -m scripts.RunT3GeneralWilsonPolynomialBatch --all --boundary-root 'output\general_observable\ordered_wilson_tensor_v4' --output-root 'output\general_observable\general_wilson_polynomial_batch_v2'

python -X utf8 -m scripts.RunT3GeneralWilsonNormalizationAudit --all --boundary-root 'output\general_observable\ordered_wilson_tensor_v4' --polynomial-root 'output\general_observable\general_wilson_polynomial_batch_v2' --output-root 'output\general_observable\general_wilson_normalization_audit_v1'
```

If the cached native CG provenance is unavailable, omit `--reuse-cg` to rerun the **local** Wolfram/Matchete matching into v4. The runner does not change GitHub.

The native-matching and master evidence sources are identified by SHA-256 in the audit summary. The final status is intentionally `AllRequestedPLAndPolynomialStructureVerified_SpinorNormalizationOpen`, not physical EFT1 running success.
