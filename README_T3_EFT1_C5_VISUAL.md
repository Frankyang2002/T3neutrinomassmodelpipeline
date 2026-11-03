# EFT1 C5 yellow-region visual proxy in the ORIGINAL T3 RGE dashboard

## Goal and limitations

This package supplies a numerical **display-only** C5 curve in the original coloured dashboard (`Numerical/plotting/templates/T3OptimisedReferenceRGE.html`). No physical EFT1 Wilson anomalous-dimension calculation, Matchete/Weyl conversion, or certification is implied. The 14 supported nonmixed general models have pre-existing UV/EFT1 dense samples and saved C5 boundaries; B(-1) and C(-1) lack comparable dense UV/EFT1 trajectories.

For each model, let `K_S = numerical_details.c5_at_scalar_threshold` from the **actual saved benchmark** `centre.json`. For each existing sampled EFT1 scale `mu`, draw

```
K_visual(mu) = K_S * ln(M_F/mu)/ln(M_F/M_S)
M_F = 1e5 GeV; M_S = 1e3 GeV
```

The zero at `M_F` is a CHOSEN DISPLAY ENDPOINT, **not** a matching calculation. The endpoint at `M_S` equals the actual saved scalar-threshold matched C5, so all 19 independent `C5` scalar display channels (norm + 6 absolute values + 6 real + 6 imaginary entries) join continuously to the **actual numerical SMEFT solver** at `M_S`. The yellow curves are therefore a logarithmic interpolation; they are **not** a physical EFT1 RGE solution, and may not be quoted as an EFT1 neutrino-mass prediction.

The optional previous externally certified source-only `Delta C5` calculation is preserved and remains separately guarded. The new display proxy requires no certified boundary and does not change any fitting, matching, or RGE equations.

## Install (PowerShell)

In the local checkout (nothing is pushed to GitHub):

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive "$HOME\Downloads\T3_EFT1_C5_Show_In_Original_Dashboard.zip" -DestinationPath . -Force
python -X utf8 -m pytest -q tests/python/test_t3_original_coloured_rge.py
python -X utf8 -m scripts.BuildT3OriginalColouredRGE `
  --run-id shared_A2_general16_gamma08_12_v1 `
  --points 41 `
  --eft1-display-proxy `
  --clean-obsolete
Start-Process "output\studies\general_scalar_extended\shared_a2_general16_reports\shared_A2_general16_gamma08_12_v1\interactive_RGE.html"
```

On dashboard load, the first group reads `C5 — yellow DISPLAY PROXY / green SMEFT solver`. Select `||C5||F` (or an individual matrix element) and `All scales` or `Intermediate EFT` to inspect the yellow section. `SM + Weinberg` zooms into the genuine green numerical SMEFT samples. The original model-line colours, EFT shading, selectors, zoom and y-scale controls are unchanged.

## Files (project-relative locations)

- `scripts/BuildT3OriginalColouredRGE.py` — **full replacement** of the prior extension, now supporting `--eft1-display-proxy`.
- `tests/python/test_t3_original_coloured_rge.py` — **full replacement** tests covering the new visualization and original functionality.
- `README_T3_EFT1_C5_VISUAL.md` — this note.

## Generated files (not in ZIP)

All are within:

`output/studies/general_scalar_extended/shared_a2_general16_reports/shared_A2_general16_gamma08_12_v1/`

- `interactive_RGE.html` — **same original dashboard filename and template**, with yellow C5 visual proxy, green SMEFT C5.
- `<model-tag>/eft1_c5_visual_proxy.json` — 41 numerical rows per displayed model, all matrix components, formula and unphysical flag.
- `<model-tag>/smeft_trajectories.json` and `<model-tag>/smeft_observables.csv` — original genuine SMEFT trajectory results.
- `eft1_c5_readiness.json` — per-model status `VisualProxyPlotted_NoCertifiedBoundary` if visual interpolation generated while no physical certified boundary exists.
- `original_coloured_rge_status.json` — includes `eft1_c5_visual_proxy_model_count`, count of completed SMEFT curves, and provenance flags.

`--clean-obsolete` removes only superseded additive dashboard **scripts**, tests, README files, and `interactive_RGE_with_SMEFT.html`; it does **not** delete existing centre data, full matching/running results, or sensitivity reports.

## How to remove the illustrative interpolation

Re-run without `--eft1-display-proxy` to regenerate the original dashboard with physical SMEFT-only C5, as before. Previously saved `eft1_c5_visual_proxy.json` files may remain on disk as explicitly labelled illustrations, but are not used by the newly generated dashboard.

## Verification

The package offline tests exercise the original-template adapter, numerical log interpolation, exact 3x3 complex-matrix endpoint, threshold continuity checks, and status/HTML export. Local numerical A(-2)–E(+2) execution is required to confirm the complete study in your environment.
