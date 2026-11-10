# Authoritative T3 thesis calculation and results — 11 November 2026

## Current decision (supersedes the previous two-model-pilot designation)

**Authoritative numerical RESULTS for Chapters 4–5:** the **original sixteen-model thesis qualitative representation-running study**, including its saved full-running and stage-frozen counterfactuals.

**Canonical saved result input:**
`Reports/output/thesis_qualitative/representation_running/csv/frozen_stage_observable_changes.csv`

**Canonical saved mass-squared-difference heatmap matrix:**
`Reports/output/thesis_qualitative/representation_running/csv/oscillation_heatmap_masssq_table.csv`

**Canonical plot generator (verified source):**
`scripts/PlotT3OscillationHeatmaps.py`

**Original upstream numerical producer:** **NOT YET IDENTIFIED from inspected source.** The plotting script only consumes the saved comparison CSV and does not perform EFT physics. Do not substitute a later numerical driver or assert an unverified original driver name. The original results are authoritative **for this thesis write-up**, not automatically certified physical predictions.

**Previous choice withdrawn for thesis results:** `scripts/RunT3ThesisTwoModelPilot.py` is a different A(-4)/D(-2) strict-order numerical study. It remains in the repository, but its results, approximations and benchmarks must not be mixed into the original sixteen-model comparison. The later `scripts/RunT3GeneralisedRunningComparisons.py` and its ~10^5% EFT1-off outputs are also **not** the source of the selected heatmap.

## Physics sequence represented by the thesis calculation

`UV (SM + heavy fermions + two scalar multiplets) -> integrate out ALL heavy fermions at M_F -> intermediate EFT containing BOTH scalar multiplets -> integrate out ALL heavy scalars at M_S -> SMEFT containing C5 -> low-energy neutrino masses and PMNS observables`.

This is the intended EFT hierarchy; the original numerical driver's exact renormalisable interaction inventory, matching expressions, and anomalous-dimension basis are **not independently source-verified yet**.

## Exact verified stage-switch meanings in saved CSV

The saved CSV contains `on_case` and `off_case` flags:

- `1111`: baseline full-running result.
- `0111` (`UV` switch): UV running disabled.
- `1011` (`EFT1` switch): EFT1 renormalisable-background running disabled; Wilson running is a separate switch, **not** disabled by this column.
- `1101` (`Wilson_conditional` switch): a *conditional Wilson log reconstruction*, not an independently evolved stage-off result and not shown in the standard UV/EFT1/SMEFT heatmap.
- `1110` (`SMEFT` switch): the original calculation's SMEFT running disabled, retaining the high-scale result after the scalar threshold. **Do not assert this means only C5 was frozen while SM parameters continued running** without verifying the original numerical generator.

Signed change for an observable O and selected switch X:

`delta_X(O) [%] = 100 * (O_full - O_X_off) / O_X_off`.

The stored `full_value` and `frozen_value` are numerically evaluated observables (not merely percentages). To use the no-SMEFT-running prediction, read `frozen_value` from the `SMEFT` rows directly. Do not convert a -12.3% quoted delta to a +12.3% rise; the rise relative to the full result is approximately +14%.

## Verified example: A(-4), from original CSV

| Observable | Full (`1111`) | SMEFT off (`1110`) |
|---|---:|---:|
| `delta_m21_sq_ev2` | 7.537009533779758e-05 | 8.594622876037567e-05 |
| `delta_m3l_sq_ev2` | 0.0025109952138040316 | 0.0028633419780036253 |
| `sin2_theta12` | 0.3088001623632854 | 0.30880030985871515 |
| `sin2_theta13` | 0.022479995082291888 | 0.022479990849423166 |
| `sin2_theta23` | 0.4699929478272284 | 0.46999379480939657 |

The SMEFT-on vs SMEFT-off mass-squared percentages in the original full-running comparison are approximately -12.3%, nearly uniform across representations; changes to mixing angles are much smaller. This is an observation in the saved data, conditional on their matching and running conventions.

## Scope and scientific status

- **Verified:** the 16 model names, selected observable values, stage-frozen flags, numerical percentages, and plot script, directly from GitHub source/results.
- **Supported by architecture:** sequential fermion-first then scalar EFT framework and low-energy observables.
- **Not verified:** name and complete call chain of the original numerical generator; exact interaction truncation; complete one-loop finite thresholds; gauge/Ward and Weinberg-projection validation; full general-Wilson mixing closure.
- **Therefore:** describe this as an *exploratory, conditional running-stage sensitivity study in a specified numerical EFT approximation*, not as a complete general one-loop prediction or proof of physical matching accuracy.

## Thesis writing policy

1. Chapters 4–5: all displayed benchmark numerical values and plots must come from the **original** `Reports/output/thesis_qualitative/representation_running/` dataset and its documented inputs, unless explicitly labelled as a separate study.
2. Chapter 3: general symbolic Lagrangian / sequential matching may be explained theoretically but must not be asserted to be the exact source of every numerical coefficient in the older dataset until producer provenance is established.
3. The original heatmap's `EFT1` column refers to **renormalisable-background running**; Wilson running was analysed separately.
4. Report `1110` as **SMEFT-off**, not necessarily selectively **Weinberg-only-off**.
5. Do not use the newer two-model pilot or the newer generalised stage-freeze data to overwrite the selected figures.

## PowerShell: extract selected original SMEFT-off predictions

From repository root:

```powershell
$source = 'Reports/output/thesis_qualitative/representation_running/csv/frozen_stage_observable_changes.csv'
$destination = 'Reports/output/thesis_qualitative/weinberg_running_off.csv'
Import-Csv $source |
    Where-Object { $_.switch -eq 'SMEFT' -and $_.observable -in @('delta_m21_sq_ev2','delta_m3l_sq_ev2','sin2_theta12','sin2_theta13','sin2_theta23') } |
    Select-Object model,observable,full_value,@{Name='smeft_running_off';Expression={$_.frozen_value}},relative_change_percent |
    Export-Csv $destination -NoTypeInformation
```

This extracts the existing values and **does not rerun** or change physics. Output: `Reports/output/thesis_qualitative/weinberg_running_off.csv`.

## Source provenance

Read-only inspection of `Frankyang2002/T3neutrinomassmodelpipeline`, current main as checked on 11 November 2026. The original saved dataset, heatmap CSV, and plotter were verified directly; the original numerical stage generator remains untraced. No GitHub files were changed, and no fresh numerical EFT solve was performed.
