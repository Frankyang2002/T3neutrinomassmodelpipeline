# All-running-stage OFF heatmaps

This script NEVER manufactures EFT1 combined-freeze values. It reuses verified original 16-model UV/SMEFT stage-off data but requires an actual re-evaluation with **both** EFT1 scalar-background and **all** EFT1 Wilson beta functions set to zero.

## Intended files
- `scripts/PlotT3AllStageFrozenComparisons.py`
- `docs/ALL_STAGE_FREEZE_GRAPH_INSTRUCTIONS.md`

## Preliminary graph (EFT1 unavailable, labelled N/A)
```powershell
python -X utf8 scripts/PlotT3AllStageFrozenComparisons.py --allow-incomplete
```
Output: `Reports/output/thesis_qualitative/representation_running_v2_all_frozen/`. Do not put this preliminary figure in the thesis.

## Final graph after rerun
Prepare a CSV with columns `model,observable,full_value,eft1_all_frozen_value` with 80 rows (16 models x 5 observables). Values must be calculated using the ORIGINAL pipeline and its exact fitted UV inputs, original matching, same non-frozen stages, and with all intermediate-EFT running zeroed, including dimension-five operator mixing and any Weinberg tensor component.

```powershell
python -X utf8 scripts/PlotT3AllStageFrozenComparisons.py --combined-eft1 Reports/output/thesis_qualitative/eft1_all_frozen_16models.csv
```

No archive files or original datasets are modified. This script plots data, it does not numerically solve the EFT.
