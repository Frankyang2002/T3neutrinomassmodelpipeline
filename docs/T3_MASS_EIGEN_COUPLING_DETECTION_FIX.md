# Matchete mass-eigenstate source coupling-presence fix

## Evidence from the local Windows run

The `B_-1` and `C_-1` live scripts reached `original-interaction-assembly-started`
after Matchete generated all seven extra quartic candidates. Both failed with
`AdditionalQuarticCouplingsAbsentInOriginalSource`, *before* eigenfield
registration or any heavy-mass-basis validation.

The previous dynamic check

```wl
With[{coupling = symbol}, !FreeQ[expr, HoldPattern[coupling[___]]]]
```

was evaluated independently in a Wolfram kernel and falsely returned `False`
for a present applied coupling head (`aa[]`). Thus the reported failure was a
false-negative detector, not evidence of a missing physical interaction.

## Implementation

Only `Lagrangian/RunGeneralT3MassEigenMatchete.wl` has physics-adjacent changes.
The helper now obtains all **applied expression heads** via `Cases[..., HoldPattern[head_Symbol[___]] :> HoldComplete[head], ...]` and checks whether the **exact symbolic coupling head** occurs. It distinguishes different symbol contexts, unlike a comparison of bare printed names. An early in-kernel positive/negative control checks this mechanism on each run. When a coupling really is absent, the script still aborts, now listing the missing symbols.

No numerical boundary, RG equation, scalar field transformation, Yukawa, quartic
expression, CG registration or Matchete pass criterion was altered. The existing
SHA-256 input provenance checks remain unchanged.

## Exact project paths

- `Lagrangian/RunGeneralT3MassEigenMatchete.wl` — full replacement
- `tests/python/test_general_t3_mass_eigen_matchete.py` — full replacement
- `docs/T3_MASS_EIGEN_COUPLING_DETECTION_FIX.md` — new note

## Install and run (PowerShell)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath 'C:\path\to\T3_MassEigen_Coupling_Head_Fix.zip' -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_mass_eigen_matchete.py tests\python\test_general_t3_mixed_scalar_mass_basis.py
python -m scripts.RunGeneralT3MassEigenMatchete --model B_-1
python -m scripts.RunGeneralT3MassEigenMatchete --model C_-1
python -m scripts.RunGeneralT3MassEigenMatchete --all
```

Offline tests check source contracts and Python computations; a standalone
Wolfram kernel reproduced both the old false-negative and the replacement's
correct match. Full live Matchete execution on the user's installation remains
necessary. No physical threshold matching or Matchete/RGBeta normalization is
certified by this change.
