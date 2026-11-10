# T3 Batch 9: source-based A(-4)/D(-2) joint scalar handoff

## Exact scope

This **new batch does physics**, rather than rerunning the older all-16 validation packages. It connects the actual saved **general UV/EFT1 numerical trajectories** from thesis Batch1 to the existing **Matchete source-exact F-tree LLSS** component tensor, **complete source-named HHSS scalar-quartic parser**, **simultaneous two-multiplet scalar B0 matching**, and **scalar-pair leading logarithm**. In the D(-2) branch it uses the existing *named* `T3Scalar2Minus`/`T3Scalar2Plus` field map and a complex self-bilinear mass Hessian; it never substitutes an unrotated `NewScalar2` for these split real fields.

It does **not** complete physical full-general one-loop matching. The hard one-loop F-boundary including UV/EFT1 subtraction, the full seven-channel EFT1 Wilson RGE, remaining scalar-soft channels, and signed GammaCC-to-Weyl conversion are not computed by this batch. Its matrix result is a **conditional LLSS one-loop contribution per unknown signed spinor factor**, not a total Weinberg coefficient. **Do not use it for physical masses or fits.**

### Files and installation paths

| File in ZIP (project-relative) | Purpose |
|---|---|
| `Numerical/parallel_t3/complete/two_model_joint_scalar_batch9.py` | New integrated A/D scalar one-loop response driven by *source data* at real threshold scales |
| `scripts/RunT3FullGeneralJointScalarBatch9.py` | Local runner; neither GitHub nor existing physics results modified |
| `tests/test_t3_full_general_joint_scalar_batch9.py` | Minimal implementation checks |
| `README_T3_FULL_GENERAL_JOINT_SCALAR_BATCH9.md` | This explanation |

### Install and run — PowerShell

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$HOME\Downloads\t3_full_general_joint_scalar_batch9.zip" -DestinationPath . -Force
.\.venv\Scripts\python.exe -m pytest -q tests\test_t3_full_general_joint_scalar_batch9.py
.\.venv\Scripts\python.exe -X utf8 -m scripts.RunT3FullGeneralJointScalarBatch9 --all --preflight
.\.venv\Scripts\python.exe -X utf8 -m scripts.RunT3FullGeneralJointScalarBatch9 --all
```

The run reads the repository's saved `output/parallel_t3/full_general_batch1/thesis_full_general_batch1_v1/` and Batch2 F-source outputs, together with `output/full16_symbolic/T3_{A,D}_alpha_m{4,2}/data/eft1_after_F_scalar_quartic_seed.json`. It checks that the archived F-tree Matchete operators and CG registry are byte-for-byte identical to the full-general thesis Batch2 F-tree seed before reusing them; the A and D source files matched in the inspected GitHub revision. It reuses rather than recreates Matchete source tensors. It does not require re-running WolframScript for its pre-existing inputs.

Output folder: `output/parallel_t3/thesis_full_general_joint_scalar_batch9_v1/` with `A_-4/result.json`, `D_-2/result.json`, model-specific native scalar background tensor NPZ/JSON, and `summary.json`. **The runner refuses to overwrite an existing output folder**; use e.g. `--output output/parallel_t3/thesis_full_general_joint_scalar_batch9_v2` to rerun.

### Status and limits

The code reuses source modules checked on 10 Nov 2026 on the public `main` branch. It was syntax-checked in the build environment, but could **not** be executed against the user's Windows project outputs or Matchete environment here. If the source name/CG map differs between the saved full-general thesis run and the archived original-field symbolic source, the runner reports the exact failure instead of inventing a scalar coefficient. The run must be performed locally before declaring successful numerical outputs.

Remaining physical phases (estimated, not already implemented in this package): **(1)** fully renormalized F-hard one-loop, **(2)** full EFT1 Wilson RG and scalar mass-basis transport, **(3)** complete simultaneous scalar matching with low-energy C5 and the two 18-real-parameter Yukawa fits. These should be treated as large batches, not isolated audits.
