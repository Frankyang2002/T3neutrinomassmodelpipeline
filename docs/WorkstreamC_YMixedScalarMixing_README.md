# Workstream C — Actual B/C Y-threshold mixed-scalar support

## Source-level result verified in this response

The repository's **actual current native source expressions**, inspected read-only, contain for B(-1) and C(-1):

| Source inventory | B(-1) | C(-1) |
|---|---:|---:|
| EFT1 `LLSS` additive terms | 20 | 24 |
| EFT1 `HHSS` additive terms | 42 | 42 |
| Native Y/LLHH poles | 50 | 50 |
| Native Y/LLHH logs | 50 | 50 |
| Native Y/LLHH standalone constants | 32 | 32 |
| Pole monomials supported by `LLYY` × `HHYY` | **18** | **18** |
| Pole monomials supported by `LLXY` × `HHXY` | **32** | **32** |
| Pole monomials unsupported by these coupling products | **0** | **0** |

The coupling-support test retains `Bar[Coupling]` complex conjugation, explicit coupling powers, `Proj[1]` versus `Proj[-1]` chiralities, complex-field contraction orientations, and external Higgs conjugation. It expands the two-entry `lambdaT3` mixed-scalar product rather than treating it as one coupling monomial.

**Crucial scope:** This checks *support*, NOT exact coefficient equality, actual Feynman diagrams, SU(2) contractions, cancellation of 1/epsilon, or a renormalized Wilson coefficient. EOM and CG normalization are not separately evaluated in this test. The provenance and expected Git blob hashes are encoded in the audit script, to help identify a different local source version.

The original `ComputeT3BCYQuarticMixing.experimental.py` accepted **YY only**. The updated complete replacement also accepts mixed **XY**, with new optional input key `light_x_scalar_indices`. Its real-scalar formula becomes

```
(16*pi**2)*beta(C_ijab)|quartic,Y-related
   = sum_(c,d in X union Y; c in Y or d in Y) lambda_abcd*C_ijcd
```

with **ordered** real component pairs. A real XY pair therefore appears with multiplicity 2. The pure XX contribution is excluded because it belongs to the low-energy EFT with X still active. This restriction is a proposed piece of the EFT1/EFT2 anomalous-dimension difference, not a derivation of the full threshold counterterm; scheme-matching and other contributions must still be evaluated. No X/Y tensor map from native Matchete CG conventions is implemented.

## Files and target paths

- `scripts/AuditT3BCYMixedScalarPoleSupport.experimental.py`: NEW, reads local archived native expressions already in the repository. Reports matched coupling signatures and source-term coordinates for all 50 B/C poles. No modification of inputs.
- `scripts/ComputeT3BCYQuarticMixing.experimental.py`: FULL REPLACEMENT of the previous experimental toy calculator; backwards-compatible with the old YY-only v1 JSON; adds optional `light_x_scalar_indices` and XY support.
- `docs/T3BCYQuarticMixing_XY_TOY_NOT_BC.json`: small XY exact-arithmetic regression fixture; explicitly **not** the B/C tensor input.
- `docs/WorkstreamC_YMixedScalarMixing_README.md`: these instructions.

No production matching, model or RGE module is modified.

## Reproduce locally in PowerShell

From `C:\thesistemp\T3neutrinomassmodelpipeline` after extracting the ZIP into the repository root:

```powershell
python scripts\AuditT3BCYMixedScalarPoleSupport.experimental.py --self-test

python scripts\AuditT3BCYMixedScalarPoleSupport.experimental.py `
  --root output\general_t3_mass_eigen_matchete `
  --output output\audits\t3_bc_y_mixed_scalar_support_v1\summary.json

python scripts\ComputeT3BCYQuarticMixing.experimental.py --self-test

python scripts\ComputeT3BCYQuarticMixing.experimental.py `
  --input docs\T3BCYQuarticMixing_XY_TOY_NOT_BC.json `
  --output output\audits\t3_bc_y_mixed_scalar_support_v1\XY_TOY_result.json
```

If the source tree is exactly the inspected Git snapshot, you can append `--strict-known-snapshot` to the second command to enforce all six blob hashes. Do not force this option when the repository was checked out with line-ending conversion or when you intend to inspect another recorded snapshot. All input files and SHA-256 hashes are recorded in the normal output regardless.

## Expected terminal evidence

```
B_-1: 50/50 native pole monomials supported; YY=18, XY=32, XX=0; counterterms=NOT CALCULATED
C_-1: 50/50 native pole monomials supported; YY=18, XY=32, XX=0; counterterms=NOT CALCULATED
PASS: source-level support audit
```

The XY synthetic calculator must give

```
(16*pi**2)*beta_Y_related = cXY*lamXY/2
```

This validates the new XY branch only. It does **not** map B/C data.

## Actual next calculation

The next *physical* work is the exact native-component map of both `LLYY/LLXY` and `HHYY/HHXY` sectors, including their registered SU(2) CG arrays and `GammaCC` bilinear normalization, into the same real-scalar/Weyl normalization. Then calculate the **full Y-related quartic beta tensor** and compare the relevant operator coefficient to the hard-pole residue *without* subtracting the EFT2 loop a second time. Keep the UV/IR assignment, other field/operator counterterms, finite matching and physical `C5_Y` separate until independently derived.
