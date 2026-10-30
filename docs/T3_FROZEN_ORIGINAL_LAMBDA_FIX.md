# Frozen neutrino plots without scalar-dashboard substitutions

The original `BuildFittedModelEnvelopes.HTML` renderer is used without changing its plotting JavaScript, NuFIT/JUNO input values, zoom/axis controls or colours.

The adapter finds frozen model results using the existing frozen-fit discovery, and explicitly writes a 16×2 coverage CSV for missing scans / unsupported cases.

The adapter intentionally **does not** generate a fake frozen neutrino RGE: true C5(mu), mixing and mass scale-dependent histories require computation from the model's matching/running physics. An original historical RGE HTML may be copied separately and is labelled historical, not frozen.

Run from project root:

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
python -X utf8 -m pytest -q tests\test_frozen_original_lambda.py
python -u -X utf8 -m scripts.GenerateT3FrozenSensitivityScans --run-id sensitivity_v2 --run-ids frozen_massnorm_a0_v1,frozen_massnorm_remaining_v1 --source-campaign overnight_six_v1_bulk48
python -X utf8 -m scripts.BuildT3FrozenLambdaOriginal --run-ids frozen_massnorm_a0_v1,frozen_massnorm_remaining_v1 --sensitivity-csv output/studies/general_scalar_extended/frozen_display/sensitivity_v2/frozen_gamma_points.csv --output-dir frozen_original_lambda_v1
Start-Process .\Reports\output\display_frozen_general_scalar\frozen_original_lambda_v1\index.html
```

To include the existing historical Lambda HTML reports as separate links, append:
`--historical-envelopes <path to envelopes.html> --historical-rge <path to original RGE HTML>`.

The two frozen coupling cases are differentiated by `[small]` and `[large]` in the sensitivity legend. The original historical and shared A(-2) studies are separate; this adapter does not manufacture their new frozen evaluations.
