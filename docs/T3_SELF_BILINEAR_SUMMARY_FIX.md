# Fix: six self-bilinear results completed, 16-model merge failed

The earlier `neutrino_v1` outputs use the explicit field
`full_general_scalar_physical_c5_certified: false`. The six new self-bilinear
outputs use `full_general_scalar_c5_certified: false`. The original merger
required only the latter and rejected the first unchanged source, `A_-4`.

The replacement uses the **correct family-specific fields**, checks the source
formats, and preserves false-certification and missing-contribution guards.
It also adds `--combine-only` to avoid reevaluating the six already-completed
C5 contributions or comparing an old matching-runner hash against this updated
summary code.

Project paths:

- `scripts/RunT3SelfBilinearNeutrino.py` **full replacement**
- `tests/test_t3_self_bilinear_summary_fix.py` new regression tests
- `docs/T3_SELF_BILINEAR_SUMMARY_FIX.md` this note

Install by extracting the ZIP at the project root, allowing replacement of the
single script. Do not change the root README or existing output JSON files.

```powershell
python -X utf8 -m pytest -q tests/test_t3_self_bilinear_summary_fix.py

python -m scripts.RunT3SelfBilinearNeutrino `
    --all `
    --run-id continuous_v1 `
    --combine-only
```

Expected output:

`output/studies/general_scalar_extended/continuous_full_scalar/continuous_v1/self_bilinear_v1/combined_16_summary.csv`

This corrects **only summary aggregation**. It does not calculate the absent
anomalous-propagator or independent-extra-quartic C5 diagrams, and it does not
certify full sequential Wilson matching.
