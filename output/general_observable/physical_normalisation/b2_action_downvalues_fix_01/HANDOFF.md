# B2 batch-4 native diagnostic: `DownValues` held-symbol fix

## Observed Windows execution

The native Wolfram probe loaded RGBeta 1.2.0, then emitted:

```
DownValues::sym: Argument ToExpression[First[$B4GroupInv]] at position 1 is expected to be a symbol or the name of an existing symbol.
BLOCKED: native evaluation failed
```

The following Python validator correctly blocked because the native JSON did not exist. The earlier 64 Python tests passed.

## Root cause and change

`DownValues` has a holding attribute. `DownValues[ToExpression[First[$B4GroupInv]]]` therefore receives the *unevaluated* `ToExpression[...]` rather than the corresponding package symbol. The replacement helper uses:

```wolfram
B4GroupInvDownValuesString[names_List] := If[Length[names]===1,
 With[{symbol=ToExpression[First[names]]}, B4String[DownValues[symbol]]],
 "unresolved"];
```

`With` lexically substitutes the resolved symbol *before* `DownValues` evaluates. If symbol discovery is ambiguous, it returns `"unresolved"` instead of fabricating DownValues.

Only `scripts/RunPhysicalB2RGBetaInternalActionProbe.wl` was changed. No underlying registry test, normalization expectation, coupling, SU(2) tensor, model source, or production file was changed. The new static Python test guards the regression, but **does not** constitute a live Wolfram test.

## Paths

- **Replacement:** `scripts/RunPhysicalB2RGBetaInternalActionProbe.wl`
- **New test:** `tests/test_physical_b2_rgbeta_internal_action_probe_source.py`
- **This handoff:** `output/general_observable/physical_normalisation/b2_action_downvalues_fix_01/HANDOFF.md`

## Validation

- Local Python: 68/68 tests passing (64 existing + 4 new).
- Native Wolfram: not executable in this environment; user must rerun the installed probe.
- Physical Matchete–RGBeta action-normalization certification: **0/16**.

## Windows run (from repository root)

```powershell
$zip = Join-Path $env:USERPROFILE 'Downloads\T3_B2_Internal_Action_DownValues_Fix.zip'
Expand-Archive -LiteralPath $zip -DestinationPath '.' -Force
.\.venv\Scripts\python.exe -m pytest tests\test_physical_b2_rgbeta_internal_action_probe_source.py tests\test_physical_b2_rgbeta_internal_action.py tests\test_physical_b2_rgbeta_ward_witness.py tests\test_physical_b2_gauge_coupling.py -q
$env:T3_EXPERIMENTAL_GENERAL_FERMIONS='0'
$env:T3_EXPERIMENTAL_GENERAL_QUADRATICS='0'
$env:T3_EXPERIMENTAL_GENERAL_QUARTICS='0'
$env:T3_EXPERIMENTAL_FLAVOR_MASSES='0'
wolframscript -file scripts\RunPhysicalB2RGBetaInternalActionProbe.wl 'output\general_observable\physical_normalisation\b2_rgbeta_action_native_02\internal_action_registry.json'
.\.venv\Scripts\python.exe scripts\RunPhysicalB2RGBetaInternalAction.py --native 'output\general_observable\physical_normalisation\b2_rgbeta_action_native_02\internal_action_registry.json' --output 'output\general_observable\physical_normalisation\b2_action_local_02'
```

Inspect any newly reported diagnostic or blocked registry row rather than assuming the native measurement succeeded.
