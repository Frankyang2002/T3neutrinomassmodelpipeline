# Workstream C — F/Y/X sequential-matching batch v14

## Scope and boundary

This is a *large, local-only evidence and callable-coefficient batch*, not a
claim that complete renormalized sequential EFT matching has been calculated.
It operates on the current B(-1)/C(-1) Matchete artifacts without changing the
UV Lagrangian, RGBeta conventions, threshold run or numerical benchmark.

Previous gates already passed on the user's machine:

- v11.2 `ExactPlusEquality`: native PL equals Warsaw `Q_llHH` with positive sign.
- v12: independently restricted Y scalar-bubble MSbar CT, light-X subtraction,
  PL/PR conjugacy, and conditional normalized matrix.
- v13.1: original fermion-first 68/68/68 F-hard LLHH terms; pole/log equality,
  PL/PR conjugacy and formal coefficient tensor.

This release makes **five substantive linked checks** in one command:

1. Regenerate/check the v12 **Y restricted** and v13 **F formal** ledgers.
2. Verify saved **Y/LLXX** source and SHA256, exact opaque-factor pole–log
   cancellation and constant-sector inventory. B: 240/240/186; C: 198/198/186.
3. Verify saved **X/LLHH** source and SHA256, exact opaque-factor pole–log
   coefficient equality and 18/18/0 inventory in both models.
4. Compare the saved Y tree versus full EFT2 continuation hashes and explicitly
   record that the original X Matchete match was made **on the Y tree EFT2**.
   Thus zero X standalone constant is *not* zero total X physical matching.
5. Optionally evaluate *separate* 3×3 ordered Warsaw/normalized C5 matrices
   for F-formal and Y-restricted sectors from explicitly supplied couplings.
   Their physical sum is intentionally **not produced**.

## Files and paths

Install both complete files under the existing local root:

- `scripts/BuildT3BCSequentialMatchingBatch.experimental.py`
- `scripts/RunT3BCWorkstreamCBatch.ps1`
- `docs/WorkstreamC_SequentialMatchingBatch_v14_README.md`

Outputs (script does not overwrite earlier native evidence):

- `output/audits/t3_bc_y_conditional_finite_interface_v1/summary.json` (v12)
- `output/audits/t3_bc_f_hard_boundary_v1/summary.json` (v13)
- **`output/audits/t3_bc_sequential_batch_v14/summary.json`** (new)

## PowerShell

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$HOME\Downloads\T3_WorkstreamC_SequentialMatchingBatch_v14.zip" -DestinationPath . -Force
& .\scripts\RunT3BCWorkstreamCBatch.ps1
```

If local policy blocks PS1 execution, use the direct Python alternative:

```powershell
python scripts\BuildT3BCSequentialMatchingBatch.experimental.py --self-test
python scripts\BuildT3BCSequentialMatchingBatch.experimental.py --output output\audits\t3_bc_sequential_batch_v14\summary.json
```

The direct command expects v12 and v13 reports already generated (as they are
in your project). `--require-physical` deliberately fails.

## Optional numerical input for both models

Supply an explicit JSON file at a user-chosen project path, e.g.
`input/workstream_c_bc_couplings.json`, with top-level keys `B_-1`, `C_-1`.
Each value is an object with fields

```json
{
  "model": "B_-1",
  "MF_GeV": 2000.0,
  "y1": [[0.0,0.0],[0.0,0.0],[0.0,0.0]],
  "y2": [[0.0,0.0],[0.0,0.0],[0.0,0.0]],
  "quartics": {
    "lambdaT3": [0.0,0.0],
    "lambdaT3ExtraQ020020I1": [0.0,0.0],
    "lambdaT3ExtraQ020200I1": [0.0,0.0]
  }
}
```

The numbers above only demonstrate the input schema; **they are not model
benchmark inputs**. Duplicate the object under `C_-1`, changing `model` and
using that benchmark's actual couplings. The script reads each model's saved
`unitary_U` from `*_final_parameters.json` and rejects invalid inputs. Use:

```powershell
& .\scripts\RunT3BCWorkstreamCBatch.ps1 -NumericManifest input\workstream_c_bc_couplings.json
```

For a single model instead:

```powershell
python scripts\BuildT3BCSequentialMatchingBatch.experimental.py `
  --numeric-input input\one_model.json `
  --output output\audits\t3_bc_sequential_batch_v14\numeric.json
```

Numeric results are in GeV^-1, with common `hbar=1/(16*pi^2)`, `MF` explicit,
complex pairs `[real,imag]`, and original Yukawa orientation:

`C5 = A_PL + A_PL.T` (not `A_PL + A_PL.conj().T`).

Neither separate tensor is a certified full physical C5. **Do not add their
formal values to predict neutrino masses** without the missing matching.

## Exact outstanding physics work

The next *substantive matching calculation* is a proper high-minus-low
renormalized 1PI amplitude subtraction at the F and Y boundaries:

```
Delta C5^ren = finite{Gamma_high^1loop - Gamma_low^1loop +
                        delta_CT + finite_field_EOM_parameter_shifts}
```

The original fermion-first EFT1 and independently reconstructed canonical
scalar EFT have not been proven algebraically equivalent in the full HHXY
quartic sector. The X calculation used Y *tree* EFT2 as input and omitted Y
one-loop inheritance. Threshold-scale RG evolution and scheme-consistent
fixed-order bookkeeping are also outstanding. The raw X coefficient has no
standalone constant at muX=MX; this is not evidence its physical finite
correction vanishes.

**Certification flags remain false and `physical_C5` remains null by design.**
