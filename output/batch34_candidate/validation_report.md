# A(-4) Batch 3–4: crossed-Yukawa gauge-covariance diagnostic

**Source revision inspected:** `12128de602ac6bf109dbad50dbac8446ed629434` (GitHub `main`, unchanged at check).

**Scope:** Necessary Ward-covariance and anti-linearity tests on synthetic gauge-invariant tensors; **not** a physical T3 A(-4) Wilson calculation. The physical data from Chat A/B are not present here and have not been certified.

## Source-derived conflict

In the uploaded `2512.16890v1.pdf`, Eq. (4.85), printed page 24 (PDF page 26), the crossed-Yukawa contribution contains

`sum_{sigma({i,j} x {a,b})} y_jla y_ikc C_lkbc` (no conjugation), while the conjugate-coefficient contribution directly afterwards contains `C*`.

The repository symbolic `RGE/general/WilsonTensorRGE.py` and numerical `GeneralEFT1C5Mixing.py` implement two *simultaneous* exchanges rather than the four independent exchanges on the crossed term in the v1 paper. The separate-pair-symmetric Batch34 engine uses the four exchanges; **both direct-C variants violate gauge covariance** on a gauge-invariant SM-like toy. No old source was changed.

## Ward-covariance candidate (not derived from loop counterterms)

`crossed = sum_{sigma({i,j} x {a,b})} y_jla y_ikc conjugate(C_lkbc)`.

For the tested gauge charges, this is a manifestly consistent contraction of Yukawas with a conjugate Wilson tensor. It changes the crossed term from complex-linear to anti-linear; in a candidate `A C + B C*` map, crossed-Yukawa and conjugate-coefficient Yukawa would both contribute to `B`. This is a **hypothesis pending independent source/derivation confirmation**, not a certified correction to Eq. (4.85).

## Numerical results

| Test | Source-v1 four-swap direct C | C* candidate |
|---|---:|---:|
| SM-like (L,eC,H): relative total-beta Ward residual | 0.7048458149779733 | 1.940015680834073e-16 |
| SM-like: crossed-term Frobenius norm | 0.1008414597276338 | 3.2685934294986725e-18 |
| U(1) nonzero crossed test: crossed-term Frobenius norm | 3.195479074111431e-18 | 0.07747580267412524 |
| U(1) nonzero crossed test: total-beta relative Ward residual | (not used for source validation) | 9.26209903771952e-18 |

All input Wilson, Yukawa, and quartic tensors used here separately pass gauge Ward identities to `<= 5.56e-17`. The `U(1)` test has charges `(0,0,-1)` for the Weyl fermions, a complex scalar of charge `+1`, and a gauge-invariant initial `C[2,2] ~ X X`; the candidate produces a **nonzero** induced `C[0,0] ~ X*X`, so it is not only a trivial zero test.

## Tests executed in this runtime

- 5/5 candidate-regression tests passed.
- 3/3 previously supplied Ward-gate tests passed (expected failure is *detected*, not fixed in the original engine).
- 14/16 existing seven-term infrastructure tests passed; the other two were skipped because the upstream GitHub repository files are not mounted in this runtime.

**Status: candidate Ward tests pass. Physical EFT1 RGE not certified.** The existing seven-term solver is unmodified. Do not attach physical Chat A/B inputs to the candidate until a loop-level derivation or authoritative corrected equation resolves the source conflict.

## Run locally (PowerShell)

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
python -m unittest discover -s tests -p "test_batch34_crossed_candidate.py" -v
python scripts\CheckBatch34CrossedCandidate.py --output output\batch34_candidate\ward_candidate_local.json
```

Outputs are created exclusively. Choose a new filename if the output already exists.
