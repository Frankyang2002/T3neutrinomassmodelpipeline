# Batch 4: inherited F-loop Weinberg term and fixed-order pole audit

## Scope

This batch builds on the **actual saved native scalar-only one-loop hard matching**
from Batch 2 and the **three-nondegenerate-heavy-fermion lift** from Batch 3.
It reads the **archived Matchete-generated first-fermion one-loop EFT1** from
`output/full16_symbolic/T3_{D_alpha_p2,A_alpha_m4}/data/eft1_after_F_one_loop_inputform.txt`.

No Wolfram matching or first-threshold rematching is performed; all inputs are
read-only. **This is a source-audited conditional Weinberg-sector calculation,
not a certified complete physical one-loop SMEFT matching.** In particular,
the legacy fermion one-loop archive has a single unindexed `MF` and the
per-mass heavy-flavour lift assumes the diagonal benchmark mass basis. Its
full general three-mass first-threshold one-loop matching has not been
independently redone. The full EOM and other SMEFT operator classes and
absolute SU(2)/GammaCC-Weyl normalization remain unverified.

## Files installed

- `Numerical/running/T3OrdinaryInheritedFixedOrder.py`
- `scripts/RunT3OrdinaryInheritedFixedOrder.py`
- `tests/python/test_t3_ordinary_inherited_fixed_order.py`
- `README_T3_SEQUENTIAL_BATCH4.md`

All paths are relative to `C:\thesistemp\T3neutrinomassmodelpipeline`.

## PowerShell commands

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$HOME\Downloads\T3_Sequential_Batch4.zip" -DestinationPath . -Force

python -X utf8 -m pytest -q tests/python/test_t3_ordinary_inherited_fixed_order.py
python -X utf8 -m scripts.RunT3OrdinaryInheritedFixedOrder --check-only --require-batch3-parity
python -X utf8 -m scripts.RunT3OrdinaryInheritedFixedOrder --run-id inherited_pole_01 --require-batch3-parity
```

Resume safely:

```powershell
python -X utf8 -m scripts.RunT3OrdinaryInheritedFixedOrder --run-id inherited_pole_01 --require-batch3-parity --resume
```

If your earlier run directories differ, use `--hard-root`, `--batch3-root`,
`--running-root`, and `--archive-root` to point to the actual saved files.

## Physics and bookkeeping

In the holomorphic Weinberg sector, with `hbar = 1/(16 pi^2)` and a
**common** renormalisable coupling convention at `muF`, define

\[
P^{(r)}_{ij}=-\frac{2\,\hbar\,\lambda_{T3}(\mu_F)}{\sqrt3\,M_{F,r}}
  \frac{\overline y_{1,ir}\overline y_{2,jr}
       +\overline y_{1,jr}\overline y_{2,ir}}{2}.
\]

The archived F-first one-loop EFT1 has six independent pure-SM `LLHH`
terms: the same finite, simple-pole and `Log[muBarF^2 / MF^2]` coefficient
in the holomorphic sector, plus conjugates. The holomorphic coefficient is
**opposite** to the independent scalar-hard pole:

\[
C^{\rm hard}_{5,\mathrm{pole}}=\sum_r P^{(r)}/\epsilon,
\qquad
C^{\rm inherited}_{5,\mathrm{pole}}=-\sum_r P^{(r)}/\epsilon.
\]

The source audit checks that each complete finite/pole/log triplet has
identical non-logarithmic operator structure, the expected +2/sqrt(3)
coefficient, and the correct `lambdaT3`, `y1`, `y2`, mass and chirality.
The Matchete hard source is independently hashed using Batch 3's restricted
fingerprint. This cancellation is established **within this source-derived
Weinberg operator sector and common-coupling convention**, rather than as
an all-SMEFT counterterm proof.

The strict one-loop `muF`-coupling fixed-order pieces are

\[
\begin{aligned}
C^{\rm hard,finite}_{ij} &= \sum_r P_{ij}^{(r)}
  \left[1+\frac{m_1^2\ln(\mu_S^2/m_1^2)-m_2^2\ln(\mu_S^2/m_2^2)}
  {m_1^2-m_2^2}\right],\\
C^{\rm EFT1,direct}_{ij} &= -2\sum_r P_{ij}^{(r)}\ln(\mu_S/\mu_F),\\
C^{\rm inherited,finite}_{ij} &= -\sum_r P_{ij}^{(r)}
  \left[1+\ln(\mu_F^2/M_{F,r}^2)\right].
\end{aligned}
\]

They simplify to a **matching-scale independent** expression:

\[
C^{\rm fixed}_{5,ij}=\sum_r P^{(r)}_{ij}
\frac{m_1^2\ln(M_{F,r}^2/m_1^2)-m_2^2\ln(M_{F,r}^2/m_2^2)}
{m_1^2-m_2^2}.
\]

`m1sq`, `m2sq`, and `lambdaT3` are taken at `muF` for the consistent
one-loop expression, with three diagonal positive `MF` eigenvalues and
`y1`, `y2` at `muF` from the existing illustrative running snapshot. The
numerical result is the symmetric 3x3 complex `C5` matrix and
conditional Takagi mass singular values via `mν = v² C5 / 2`.

The code also outputs a **distinct diagnostic** combining Batch 3's
`muS`-evolved scalar-subset with the `muF` inherited term. This combination
uses unequal couplings; **its poles do not cancel exactly**. Do not treat it
as a fully renormalized matching result. Its residual indicates the order of
unaccounted resummed/higher-order terms. Do not add the older Batch 1 LLSS
matrix to either Batch 4 subtotal: it would double-count the same sector.

## Indicative numbers for the uploaded snapshots

Conditional masses (eV, ascending), **fixed couplings at `muF`**:

| Model | m1 | m2 | m3 |
|---|---:|---:|---:|
| D(+2) | 0.019254645 | 0.020893089 | 0.024031308 |
| A(-4) | 0.019138817 | 0.020766923 | 0.023885167 |

These are derived from the uploaded numerical RG snapshots and the archived
source coefficients. The numbers are conditional on the source spinor and
SU(2) normalization, the single-unindexed-MF legacy archive and the
truncated set of Weinberg contributions. The complete physical C5 is **not
independently certified**.

## Outputs

Under `output/general_sequential_pilot/ordinary_inherited_pole/inherited_pole_01/`:

- `settings.json` and `summary.json` (run-level provenance)
- `comparison.csv` (masses and pole residuals)
- `D_p2/result.json` and `A_m4/result.json`:
  - `common_hard_pole_residue_gev_inverse`
  - `common_inherited_pole_residue_gev_inverse`
  - `common_total_pole_residue_gev_inverse`
  - `common_hard_finite_gev_inverse`
  - `common_direct_eft1_log_gev_inverse`
  - `inherited_F_one_loop_finite_gev_inverse`
  - `common_fixed_order_inherited_inclusive_c5_gev_inverse`
  - `common_fixed_order_masses`
  - `RG_improved_mixed_endpoint_uncancelled_pole_residue_gev_inverse`
  - per-heavy-eigenvalue contributions, source SHA256 fingerprints, and flags

Existing Batch 2/3 outputs are never overwritten. A resumed Batch 4 run
requires unchanged SHA256 fingerprints on all seven source files per model.

## Next independent verification

For genuinely complete one-loop EFT matching, rematch the full-general
**first fermion threshold at one loop** with a verified three-eigenvalue
heavy-mass convention, match its inherited full SMEFT sector consistently,
and audit operator normalization and counterterms (not just the six LLHH
source terms) in the same scheme. Retain a separate proof for mixed and
self-bilinear scalar bases.
