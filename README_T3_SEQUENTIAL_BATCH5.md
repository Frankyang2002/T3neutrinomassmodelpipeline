# Batch 5 — remaining six ordinary T3 models

Scope: **B(−3), B(+1), C(−3), C(+1), E(−4), E(+2)**. The existing D(+2) and A(−4) Batch 2–4 outputs are untouched.

## Source-checked distinctions

| Models | Fermion tree Wilson operators | Additional UV scalar quartics | Inherited Weinberg pole prefactor (per \(\hbar\,\lambda_{T3}\,\bar y_1\bar y_2/M_F\)) |
|---|---:|---:|---|
| B(−3), B(+1) | 2 | 1 | \(+2/\sqrt3\) |
| C(−3), C(+1) | 4 | 1 | \(-2/3\) |
| E(−4), E(+2) | 4 | 0 | \(-2\sqrt{2/3}\) |

These counts and model-specific prefactors come from the existing source-generated first-fermion thresholds. The inherited InputForm also has separate finite, \(1/\epsilon\), and logarithmic Weinberg terms in two conjugate sectors. The new scalar hard template has to be **independently generated on your local Wolfram installation**. The numerical runner does **not** extrapolate from D(+2): it parses each new hard expression and algebraically checks its full scalar-mass dependence and **opposite** pole coefficient. If the hard output contains another coupling or a different mass kernel, that model is blocked rather than assigned a guessed mass.

## Install

In PowerShell:

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$HOME\Downloads\T3_Sequential_Batch5.zip" -DestinationPath . -Force
python -X utf8 -m pytest -q tests/python/test_t3_ordinary_six_fixed_order.py
```

The numerical runner requires `numpy` and `sympy`. The test runner additionally requires `pytest`. Batch 3 and 4's existing modules must remain installed (`Numerical/running/T3OrdinaryHardFlavor.py`, `T3OrdinaryInheritedFixedOrder.py`).

## Stage 1 — reuse F tree, perform scalar-only hard Matchete

```powershell
python -X utf8 -m scripts.RunT3GeneralOrdinaryScalarHardBatch5 --preflight-only
python -X utf8 -m scripts.RunT3GeneralOrdinaryScalarHardBatch5 --run-id ordinary6_hard_01
```

This computes one-heavy-mass-eigenvalue hard scalar templates, without rerunning any F threshold. It requires local WolframScript/Matchete. A single blocked model is saved with a model-specific `error.json` and `wolfram_stdout.log` / `wolfram_stderr.log`. Retry with the same run-id only after correcting the actual failure:

```powershell
python -X utf8 -m scripts.RunT3GeneralOrdinaryScalarHardBatch5 --run-id ordinary6_hard_01 --resume
```

Hard outputs reside in `output/general_sequential_pilot/ordinary_scalar_hard_batch5/ordinary6_hard_01/<tag>/`, including `hard_c5_candidate.txt`, `hard_scalar_report.json` and `hard_scalar_input_audit.json`. The source fingerprints protect successful stages from silent overwrite.

For a short first pass through the different representation classes, run `--models B_-3,C_-3,E_-4 --run-id ordinary6_hard_pivots_01`, and give the subsequent numerical runner `--hard-root output/general_sequential_pilot/ordinary_scalar_hard_batch5/ordinary6_hard_pivots_01 --models B_-3,C_-3,E_-4`.

## Stage 2 — guarded three-MF flavour sum + inherited F Weinberg + poles

After Stage 1 succeeds:

```powershell
python -X utf8 -m scripts.RunT3OrdinarySixFixedOrder --check-only --verify-fermion-sources
python -X utf8 -m scripts.RunT3OrdinarySixFixedOrder --run-id ordinary6_fixed_01 --verify-fermion-sources
```

Resume safely:

```powershell
python -X utf8 -m scripts.RunT3OrdinarySixFixedOrder --run-id ordinary6_fixed_01 --verify-fermion-sources --resume
```

Results: `output/general_sequential_pilot/ordinary_six_fixed_order/ordinary6_fixed_01/summary.json`, `comparison.csv`, and one `<tag>/result.json` per successful model.

For each nondegenerate fermion mass \(M_{F,r}\), with \(\hbar=1/(16\pi^2)\), the source-derived **hard** pole residue is

\[
P^{(r)}_{ij}=-a_X\,\hbar\,\lambda_{T3}\frac{\bar y_{1,ir}\bar y_{2,jr}+\bar y_{1,jr}\bar y_{2,ir}}{2M_{F,r}},
\quad a_B=\frac{2}{\sqrt3},\;a_C=-\frac23,\;a_E=-2\sqrt{\frac23}.
\]

The common-coupling restricted Weinberg subtotal is checked numerically against

\[
C_{5,ij}^{\rm restricted}=\sum_r P^{(r)}_{ij}
\frac{m_1^2\ln(M_{F,r}^2/m_1^2)-m_2^2\ln(M_{F,r}^2/m_2^2)}{m_1^2-m_2^2},
\]

with scalar parameters at \(\mu_F\). This combines hard finite, direct EFT1 logarithm, and inherited F finite/log contributions. The common-coupling opposite \(1/\epsilon\) poles and scalar matching-scale variation are checked. The independent old LLSS kernel is used **only for comparison**, never added to the new result.

`--require-reference-parity` is optional: earlier reference kernels are exploratory and are **not** authoritative physical normalizations. They may differ because their scalar coupling basis was provisional. Do not use this switch as a prerequisite for a successful native hard coefficient.

## Limits

The new Wolfram stages **have not been run in this environment**. The 27 offline regression tests include mock native coefficients across all three CG classes, malformed-source rejections, and a complete fake-data runner + resume. No new real benchmark masses are claimed for the six models until you execute the native stages. The numerical flavour sum is still an inferred per-eigenvalue lift of the Matchete template. Full EFT1 operator mixing/EOM, full first-fermion three-mass native one-loop re-matching, independent MS-bar counterterms and absolute SU(2)/spinor normalization remain unverified. Accordingly `physical_c5_certified = false` in the results.

## Project-relative files (full files, no patches)

- `Lagrangian/RunT3GeneralOrdinaryScalarHardBatch5.wl`
- `scripts/RunT3GeneralOrdinaryScalarHardBatch5.py`
- `Numerical/running/T3OrdinarySixFixedOrder.py`
- `scripts/RunT3OrdinarySixFixedOrder.py`
- `tests/python/test_t3_ordinary_six_fixed_order.py`
- `tests/fixtures/t3_ordinary_hard_c5_candidate.txt` (same D(+2) regression fixture as Batch 3)
- `README_T3_SEQUENTIAL_BATCH5.md`

No changes to GitHub, earlier batch files, or existing matching output directories.
