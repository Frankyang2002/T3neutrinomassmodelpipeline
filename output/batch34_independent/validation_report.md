# Batch 3–4: independent one-loop \(\psi^2\phi^2\) comparison

**Status:** `BLOCKED_CONFLICTING_ONE_LOOP_FORMULAE`. This is an independently implemented, **synthetic** numerical audit; it does not certify the physical A(-4) RGE.

## Source revision and inputs

- GitHub `main`: `12128de602ac6bf109dbad50dbac8446ed629434` (unchanged at this check).
- Aebischer, Bresciani and Selimović, *Anomalous Dimension of a General Effective Gauge Theory II: Fermionic Sector*, arXiv:2512.16890v1, Eqs. (2.15), (4.85): seven term expression with an **unstarred** crossed \(yyC\) contribution.
- Fonseca, Olgoso and Santiago, *Renormalization of general Effective Field Theories: Renormalization of fermionic operators*, arXiv:2512.15866v1, JHEP07(2026)094, Eqs. (5), (6), (17), (18), (22): independent off-shell calculation; **no identifiable crossed** \(yyC^*\) or \(yyC\) term in Eq. (22).
- Normalizations: `a^(5)_ijab = 4 C_ijab`. The beta is real-linear in `a`, so beta[a]/4 equals beta[C] at the same numerical `C`.
- Tested limit: gauge couplings zero, all other nonrenormalizable classes zero, no scalar/fermion masses or scalar cubic insertions, real fully symmetric quartics, **general complex symmetric Yukawa and Wilson tensors**. Three Weyl labels and four real scalars; no lepton-only or real-coefficient restriction.

## Independent contraction mapping

The following eight separate terms from Eq. (22) (with `a^(5)` identified as `C` *only after* the common 4-factor cancels) reproduce the six non-crossed contributions of the existing numerical evaluator:

\[
\begin{aligned}
\tfrac14 C_{ijbd}\operatorname{Tr}(Y_aY_d^\dagger+Y_dY_a^\dagger) &\to \beta_6,\\
\tfrac14\operatorname{Tr}(C_{ab}Y_d^\dagger)Y_{ijd} &\to \beta_1,\\
\tfrac14[C_{ab}Y_d^\dagger Y_d]_{ij} &\to \beta_7,\\
[C_{ad}Y_d^\dagger Y_b]_{ij}+2[C_{ad}Y_b^\dagger Y_d]_{ij} &\to \beta_3,\\
\tfrac14\lambda_{abde}C_{ijde} &\to \beta_2,\\
\tfrac14\operatorname{Tr}(C_{ab}^\dagger Y_d)Y_{ijd}+\tfrac12[Y_d C_{ab}^\dagger Y_d]_{ij}&\to \beta_5.
\end{aligned}
\]

Here **each line is inside the sum over the four independent exchanges** \(i\leftrightarrow j\), \(a\leftrightarrow b\). The notation \(Y_a\) is the matrix \((Y_a)_{ij}=y_{ija}\). This mapping does not rely on gauge identities or on the Batch34 beta implementation.

## Numeric results (seed 311)

| Comparison with independent Eq. (22) | Maximum absolute tensor difference | Difference / max reference coefficient |
|---|---:|---:|
| Six non-crossed terms | 1.25607397e-15 | 2.96694735e-16 |
| Seven terms: source-v1, four exchanges, \(yyC\) | 0.68957205 | 0.16288244 |
| Seven terms: repository, two exchanges, \(yyC\) | 0.34478602 | — |
| Seven terms: gauge-covariant candidate, \(yyC^*\) | 0.87633293 | 0.20699686 |

The independent reference beta has maximum absolute component 4.23355664 and Frobenius norm 15.52424042, in synthetic arbitrary normalization. The six non-crossed channels **each** agree within `8.1e-16`. The mismatch is entirely the nonzero extra crossed term in these tests.

## Interpretation and decision

1. The unstarred \(yyC\) term fails gauge covariance (separately established) and the independent one-loop reference.
2. Replacing \(C\) by \(C^*\) restores covariance but does **not** establish the loop formula; the corrected candidate also disagrees with the independent Eq. (22).
3. The six-term gaugeless beta reproduces the independent published expression to floating-point accuracy. This strongly suggests the extra crossed term must be revisited, but **does not establish the full gauge-sector RGE** or determine which publication is correct.
4. **No physical engine file is changed**; both current beta modes remain explicitly unverified. A definitive physical replacement requires checking the September 2026 version-of-record equation or independently recalculating the disputed loop topology. Chat A/B physical inputs are not consumed and no matched physical SMEFT \(C_5\) is reported.

## Local commands

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
python -m unittest discover -s tests -p "test_batch34_independent_master.py" -v
# Diagnostic-only: report the expected formula mismatch without a failure exit code.
python scripts\CheckBatch34IndependentMaster.py --diagnostic `
  --output output\batch34_independent_local\comparison.json
# Without --diagnostic, exit code 2 is intentional (physical gate closed).
```
