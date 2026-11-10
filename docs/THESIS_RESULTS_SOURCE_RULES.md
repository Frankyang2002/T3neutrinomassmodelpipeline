# Thesis numerical-result source rules (original 16-model dataset)

**Decision for this submission:** Chapters 4 and 5 use the original 16-model qualitative representation-running study, not the newer strict-order two-model pilot or the later generalised stage-freeze study.

## Canonical data

- Detailed numerical source: `Reports/output/thesis_qualitative/representation_running/csv/frozen_stage_observable_changes.csv`
- Source heatmap matrices: `Reports/output/thesis_qualitative/representation_running/csv/oscillation_heatmap_masssq_table.csv` and `oscillation_heatmap_mixing_table.csv` in the same directory
- Plotting source: `scripts/PlotT3OscillationHeatmaps.py`
- Generated locked hash manifest: `Reports/output/thesis_qualitative/representation_running/THESIS_DATASET_LOCK.json` (created locally using the audit script)

For the thesis, NEVER copy numbers or figures from `output/parallel_t3/thesis_generalised_running_v1`, `thesis_original_style_figures_v2`, `thesis_original_style_figures_v3`, or `thesis_A4_Dm2_pilot_v1` into the old 16-model dataset. Those calculations have different conventions and are not comparable without a separate methodology section.

### Stage definitions from CSV

| Saved label | Baseline | Off case | Interpretation |
|---|---|---|---|
| `UV` | 1111 | 0111 | UV running off |
| `EFT1` | 1111 | 1011 | EFT1 renormalisable background off; Wilson separately controlled |
| `Wilson_conditional` | 1111 | 1101 | Conditional Wilson-log reconstruction, not a standard independently evolved counterfactual |
| `SMEFT` | 1111 | 1110 | Original SMEFT running stage off; not demonstrated to be Weinberg-only running off |

Each row uses `100 * (full_value - frozen_value) / frozen_value`.

## Required review before final thesis submission

1. Run `python -X utf8 scripts/LockT3ThesisDataset.py --write-lock` in the repository root. This checks all 512 detailed rows, exact stage switches, percentages, heatmap matrices and a checkpoint. Then use `--check-lock` to detect accidental substitutions.
2. Replot from the **canonical CSV** only, or keep the existing figures if already correctly linked. Use: `python -X utf8 scripts/PlotT3OscillationHeatmaps.py --input Reports/output/thesis_qualitative/representation_running/csv/frozen_stage_observable_changes.csv --plot-dir Reports/output/thesis_qualitative/representation_running/plots --csv-dir Reports/output/thesis_qualitative/representation_running/csv`.
3. Inspect thesis `.tex` figure paths: they must point to `Reports/output/thesis_qualitative/representation_running/plots/...` or copied images verified against that folder. Do not relabel a newer plot as original.
4. The **upstream numerical producer remains unidentified**. The audit verifies data consistency, not its full one-loop physics correctness or complete matching. Record this honestly.

## Non-destructive archiving

Run `.\scripts\ArchiveNonThesisT3Outputs.ps1` for a preview. Only with `-Apply` will it MOVE five specific, newer **output** directories under `output/archive/non_authoritative_thesis_2026-11-11/`, if present. It leaves all original thesis outputs and source files untouched. Do not move `Numerical/`, `RGE/`, `Lagrangian/`, `scripts/`, or shared `output/parallel_t3/background`/`rough_inputs` data. The source of the old workflow has not yet been fully traced.
