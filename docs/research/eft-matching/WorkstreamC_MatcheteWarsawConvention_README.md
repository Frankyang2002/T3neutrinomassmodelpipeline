# Workstream C v9.1 — Missing `SMEFT_Warsaw` model fallback

The v9 run reached `Needs["Matchete`"]`, but `LoadModel["SMEFT_Warsaw"]` failed.
This is a missing *optional* model file, not a failed B/C threshold or a
failed normalization identity. Matchete supports loading model files by path:
https://matchete.gitlab.io/quickstart.html . The separate model database is
https://gitlab.com/matchete/model-database .

## Default: advance without the optional model

Extract this archive at the repository root. Run:

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$HOME\Downloads\T3_WorkstreamC_MatcheteWarsawConvention_v9_1.zip" -DestinationPath . -Force
python scripts\InspectT3BCWeinbergNormalization.experimental.py --self-test
python scripts\InspectT3BCWeinbergNormalization.experimental.py `
  --output output\audits\t3_bc_y_warsaw_convention_v1\summary.json
```

The Python script reads the existing *actual native Matchete B/C outputs*:

- `output/general_t3_mass_eigen_matchete/eft1_scalar_Y_only_canonical_probe/B_m1_final_y_only_canonical_match_raw_scalar_llhh.inputform.txt`
- `output/general_t3_mass_eigen_matchete/eft1_scalar_Y_only_canonical_probe/C_m1_final_y_only_canonical_match_raw_scalar_llhh.inputform.txt`

It verifies all 132 terms per model and checks that every pole and constant
has two H, two lepton fields, two SU(2) epsilon CGs, one `GammaCC` and one
loop/mass factor, preserving the PR/H.c. orientations. It uses the repository's
existing `Numerical.experimental.T3BCMatcheteFiniteTerms.split_inputform`.

**This proves a structure inventory compatible with the published Warsaw
operator; it does *not* independently fix the native `GammaCC` to Warsaw
operator coefficient, an overall 1/2, or a neutrino mass sign.**

If the failed v9 `installed_warsaw_dim5.json` exists, you may pass it as
`--reference`. The v9.1 Python script recognizes the `WarsawModelLoadFailed`
status and switches explicitly to the native-source fallback.

## Optional: genuinely inspect an official Warsaw model

Once the official `SMEFT_Warsaw.m` is available locally, run the standalone
Wolfram exporter with an **actual model-file path**:

```powershell
$dir = 'output\audits\t3_bc_y_warsaw_convention_v1'
# Replace the path below with the actual official model-file location.
$warsawModel = 'C:\path\to\SMEFT_Warsaw.m'
Test-Path $warsawModel
wolframscript -file Lagrangian\ProbeT3BCWeinbergNormalization.experimental.wl `
  "$dir\installed_warsaw_dim5.json" "$warsawModel"
if ($LASTEXITCODE -ne 0) { throw 'Warsaw model reference export failed' }
python scripts\InspectT3BCWeinbergNormalization.experimental.py `
  --reference "$dir\installed_warsaw_dim5.json" --output "$dir\summary.json"
```

This optional comparison is still structural, rather than a rigorous
independent derivation of the Majorana spinor sign. Do not install a guessed
model file or reinterpret the reference operator as already normalized C5.

## Deliberately uncertified

- `absolute_GammaCC_to_Warsaw_factor_certified = false`
- `Warsaw_C5_to_neutrino_mass_sign_certified = false`
- `all_EFT1_EFT2_counterterms_certified = false`
- `physical_C5 = null`

No changes to production files or the remote GitHub repository.
