# Batch 17 — All-16 same-endpoint conditional backgrounds

**Project paths**

- `Numerical/experimental/Batch17All16Background.py`
- `scripts/RunBatch17RebuildAll16.py`
- `tests/test_batch17_all16_background.py`

This is a **new** all-16 runner; it does not overwrite the successful
`RunBatch17RebuildAminus4.py` or your original `output/batch17/rebuilt_A_m4`.

From the repository root in PowerShell:

```powershell
python -m unittest discover -s tests -p test_batch17_all16_background.py -v
python scripts\RunBatch17RebuildAll16.py --models A_-4 --out-dir output\batch17\all16_smoke_A_m4
python scripts\RunBatch17RebuildAll16.py --out-dir output\batch17\rebuilt_all16_v1
```

`--background-only` restricts the run to scalar mass/quartic background construction.
For retries, choose a **new** `--out-dir`, because outputs are never silently
replaced. Output contains `<model>_conditional_c5.json`, `<slug>.json`,
`<slug>.npz`, `summary.json`, and per-model `*_error.log` when blocked.

## Physics limitations

- Uses the exact new Batch14 evolved EFT1 coupling values at muS, with the
  source-hash-checked Matchete quartic geometry. The nontrivial RGBeta↔Matchete
  quartic coefficients are **candidate** crosswalks.
- For A_-4 and the other seven no-extra-bilinear models, the original and
  F-threshold scalar bases coincide. For six self-bilinears and B_-1/C_-1,
  the code uses the **existing exploratory** `special_mass_v6` field map.
  It checks its eigenvalue spectrum against Batch16, but does **not** certify
  its CG basis or eigenvectors.
- Uses the existing provisional historical Matchete fermion-tree source;
  this is **not** re-derived from your new Batch13 UV benchmark flavours.
- The downstream `run_conditional` includes the LLSS-induced EFT1 leading log,
  one simultaneous two-scalar finite bubble, and SMEFT evolution. It does
  **not** include all seven coupled EFT1 dimension-five RGE terms, hard F-loop
  matching or the missing source-to-Weyl signed normalisation.
- Outputs are **candidate masses per unit unknown N**, not physical predictions.
  Individual models may be BLOCKED by missing couplings, basis Ward residual,
  or mismatched source inventories. Never relax these checks to force 16/16.

EFT order: **UV → integrate out all F → EFT1(S1,S2) → integrate out both
scalars simultaneously → SMEFT**.
