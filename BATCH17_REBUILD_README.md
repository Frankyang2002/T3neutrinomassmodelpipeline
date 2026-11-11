# Batch 17 correction: rebuild A(-4) scalar matching background

Project paths:

- `Numerical/experimental/Batch17RebuildBackground.py`
- `scripts/RunBatch17RebuildAminus4.py`
- `tests/test_batch17_rebuild_background.py`
- `BATCH17_REBUILD_README.md`

Existing Batch 13–17 files are **not** modified.

## Run (PowerShell, repository root)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -Path "$HOME\Downloads\T3_Batch17_Rebuild_Aminus4_Background.zip" -DestinationPath . -Force
python -m unittest discover -s tests -p test_batch17_rebuild_background.py -v
python scripts\RunBatch17RebuildAminus4.py --background-only
python scripts\RunBatch17RebuildAminus4.py --out-dir output\batch17\rebuilt_A_m4_with_C5
```

`--background-only` creates `output/batch17/rebuilt_A_m4/A_m4.json` and `.npz`.
The second command uses a **distinct output directory** (the builder refuses to overwrite). It generates the same-endpoint HHSS input and then runs the original Batch17 joint scalar LLSS+SMEFT conditional calculation, writing `A_-4_conditional_c5.json` and `summary.json`.

### Inputs

By default it reads the precise existing source release located in the user's search output:

`output/general_observable/lo_fermion_boundary/releases/source_021dcc1e8fe8ba0f/`

and uses:

- `output/batch14/eft1/A_-4_EFT1_running.json`
- `output/batch15/scalar_threshold/A_-4_scalar_threshold_handoff.json`
- `output/batch16/scalar_mass_audit/A_-4_mass_hessian.json`
- `output/general_observable/lo_joint_scalar_matching/source_quartic_v4/native_quartic_channels/A_m4.json`
- `output/full16_symbolic/T3_A_alpha_m4/data/eft1_after_F_scalar_quartic_seed.json`
- `RGE/general/generated/general_rgbeta_quartic_registry.json`

Override a missing source path with `--quartic-cache`, `--quartic-seed`, `--source-release`, etc.; the program never substitutes values for absent inputs. If the v4 cache is missing, use a source-matched cache or regenerate the native Matchete HHSS channels separately before running this builder.

### Scope and limitations

The A(-4) model has no additional quadratic mass operator; the fixed original/source real scalar basis is identical. Its mass Hessian is `diag(mS1Sq I_2, mS2Sq I_6)` and must agree *entrywise* with Batch16 and *spectrally* with the conditional matcher. The HHSS tensor uses `GeneralT3SequentialLeadingLog.scalar_tensors` with exact existing Matchete native tensors, existing portal completion, and candidate registry-based additional HHSS reference invariants evaluated at the **new Batch14 muS couplings**.

This output is a numerical *conditional* LLSS-induced leading log and joint finite scalar bubble followed by SMEFT evolution. **Not physical C5.** The Matchete GammaCC-to-Weyl signed factor, the native RGBeta/Matchete quartic coupling normalisation, historical Matchete F-tree source flavour values versus the new illustrative UV benchmark, other EFT1 Wilson operators, and independent hard F-loop matching remain uncertified. The builder does not claim to repair these.

The remaining 15 models are **not** enabled here. In particular, the eight special self-bilinear/mixed models require certified component-basis maps and cannot safely be generalized by copying A(-4) matrices.
