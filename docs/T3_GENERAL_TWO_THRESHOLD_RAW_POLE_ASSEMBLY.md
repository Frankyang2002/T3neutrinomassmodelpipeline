# T3 B_-1/C_-1: two-threshold raw Weinberg-sector pole cancellation

## Scope and verification

This batch consumes **existing**, locally generated Matchete InputForm outputs from

- `output/general_t3_mass_eigen_matchete/eft1_fermion_first_oneloop/`
- `output/general_t3_mass_eigen_matchete/eft1_hard_pole_audit/`
- `output/general_t3_mass_eigen_matchete/eft1_scalar_threshold_probe/`

It does not invoke Wolfram, rerun any matching, update GitHub, or modify source expressions. The scalar-threshold inputs are explicitly certified as `LeadingC5OneLoopSectorOnly`, with 122 four-heavy-scalar quartic terms **omitted**; the complete scalar potential still fails a separate Hermiticity check. This batch does **not** establish complete matching or the physical neutrino mass matrix.

## Exact observed term structure

| InputForm sector | B_-1 | C_-1 |
|---|---:|---:|
| Fermion threshold: UV poles | 68 | 68 |
| Fermion threshold: matching logs | 68 | 68 |
| Fermion threshold: nonlog finite terms | 68 | 68 |
| Scalar threshold: UV poles | 68 | 68 |
| Scalar threshold: scalar logs, MSX | 50 | 50 |
| Scalar threshold: scalar logs, MSY | 50 | 50 |
| Scalar threshold: nonlog finite terms | 64 | 64 |
| **Opposite UV pole-residue pairs, matched** | **68/68** | **68/68** |

On both actual datasets, the complete signed operator residue of the scalar-threshold `1/ϵ` terms is the negative of the fermion-threshold residue:

\[
\operatorname{Res}_{1/\epsilon}\!\left[\Delta\mathcal L_{LLHH}^{F,\;1\ell} + \Delta\mathcal L_{LLHH}^{S,\;1\ell}\right] = 0.
\]

Comparison retains *all* `Field`, `NCM`, `CG`, Yukawa-flavour indices, mixing spurions, and `MF` factors. It does **not** delete the operator structures to make a spurious coefficient equality.

`B_-1`: all 68 opposite signed operator residues match without changing dummy-index labels.

`C_-1`: 34 PR/conjugate terms match by SU(2) epsilon antisymmetry, while the 34 PL/holomorphic terms additionally require the change of *summed SU(2) index labels* `1 <-> 2, 3 <-> 4`. Every `Index[d$$N, SU2L[fund]]` appears exactly twice in the affected terms, so this is an allowed closed-dummy relabelling, not a change to a physical index. The two epsilon matrices are independently antisymmetric. A single epsilon argument swap therefore contributes a minus sign; two swaps cancel. The audit rejects deviations from these exact forms.

## Outputs

Run the script from the project root:

```powershell
python -m pytest -q tests\python\test_general_t3_two_threshold_pole_assembly.py
python -m scripts.RunGeneralT3TwoThresholdPoleAssembly --all
```

All output files are under

`output/general_t3_mass_eigen_matchete/two_threshold_operator_assembly/`.

For each `B_m1_final` / `C_m1_final`, the runner writes seven complete Mathematica InputForm operator expressions and a JSON report.

- `*_hard_pole_residue_canonical.inputform.txt`: first-stage `1/ϵ` residue, with `ϵ` removed but otherwise retaining the exact LLHH operator.
- `*_scalar_negative_pole_residue_canonical.inputform.txt`: negative second-stage pole residue. **Bytewise identical** to the corresponding canonical first-stage file.
- `*_hard_finite.inputform.txt`: raw first-stage terms without pole or log.
- `*_soft_finite.inputform.txt`: raw scalar-threshold terms without pole or log.
- `*_soft_logs.inputform.txt`: scalar-threshold mass-dependent log terms.
- `*_fixed_common_mu_unprojected.inputform.txt`: 300 terms from both matching stages after exact cancellation of all `1/ϵ` terms; `μbar2` remains symbolic.
- `*_fixed_muF_unprojected.inputform.txt`: 232 terms in the formal common-scale convention `μbar2 = MF^2`: the 68 first-stage hard logs vanish; second-stage logarithms become `Log[MF^2/MSX^2]` or `Log[MF^2/MSY^2]`.

**These are entire one-loop operators, not entries of the Weinberg Wilson matrix**. `hbar` remains symbolic. For a numerical Mathematica substitution, Matchete's loop bookkeeping convention in this pipeline is `hbar -> 1/(16 Pi^2)`.

## Numerical scale context (illustrative only)

The existing frozen-reference scenario uses `MSX = 1000 GeV`, `MSY = 1142.857142857143 GeV`, and reference heavy fermion masses `100`, `110`, and `120 TeV` (different fermion flavours). Therefore

| MF (GeV) | `Log[MF^2/MSX^2]` | `Log[MF^2/MSY^2]` |
|---:|---:|---:|
| 100,000 | 9.210340 | 8.943278 |
| 110,000 | 9.400961 | 9.133898 |
| 120,000 | 9.574983 | 9.307921 |

These are potentially significant logarithms, indicating why intermediate-EFT leading-log running should ultimately be included. This table is *not* a calculation of the neutrino masses. The symbolic `MF` in the files must be interpreted per fermion mass eigenstate when constructing three-generation coefficients.

## What is and is not established

**Verified from both uploaded raw expressions:**

1. Pole cancellation of 68/68 full signed LLHH operators per model, with the stated explicit SU(2) identities.
2. Every raw hard/scalar LLHH term belongs to its corresponding saved one-loop correction (accounting for a leading Mathematica unary-minus representation).
3. Both matching-stage metadata certificates and the previously generated fermion-stage pole-audit SHA256 hashes agree.
4. Outputs contain no `ϵ` UV-regulator symbol, with 300 formal common-scale or 232 formal `μbar = MF` terms.

**Not yet established:**

- An *independent* `LLSS -> C5` EFT1 anomalous dimension or full RG consistency (pole cancellation across two computations is distinct from an RGE beta-function calculation).
- The complete renormalized EFT1, as the full four-heavy-scalar quartic potential was excluded at the scalar threshold.
- SU(2), spinor, flavour and Majorana normalization of the extracted operator relative to a physical Weinberg coefficient `C5[p,q]`.
- Running from each fermion mass eigenvalue through the scalar threshold and below it, or comparison against the independently obtained direct one-scale mass eigenvalues.
- A matrix-level physical sequential neutrino mass prediction.

**Next targeted physics step:** project the pole-free full `LLHH` operator into the single consistent gauge-invariant Weinberg tensor, including chirality orientation and symmetry factors; then implement the EFT1 `LLSS` mixing and run between the two thresholds. Only after those checks is numerical matrix assembly defensible.

## Implementation provenance

New source: `scripts/RunGeneralT3TwoThresholdPoleAssembly.py`.

New tests: `tests/python/test_general_t3_two_threshold_pole_assembly.py`.

All prior files and validated matching executables remain unchanged. Do not overwrite the reference matching exports with the composite outputs. Do not present the common-scale expression as an already scheme-renormalized `C5`.
