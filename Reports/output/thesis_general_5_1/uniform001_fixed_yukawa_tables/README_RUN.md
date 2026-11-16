# Uniform-0.01 general-T3: fixed D(-2) and A(-4) Yukawa tables

## Physics specification

This is the **fixed-Yukawa cross-representation transfer** shown in the user's example tables C.1–C.4, recalculated under the new **general-Lagrangian uniform-0.01 benchmark**. It is **not** the same as the separately refitted 16-model spectrum.

- At the UV scale `1e7 GeV`, set all symmetry-allowed **BSM dimensionless scalar quartics apart from `lambdaT3`** to **real +0.01**; set symmetry-allowed `y1Conj`, `y2Conj` to **zero** in the UV initial state. Do **not** remove these allowed operators from the general Lagrangian or claim their RG values remain zero.
- Retain saved representation-specific `lambdaT3`: **0.1** for A/B/C/D, **0.05** for E. Keep existing SM Higgs quartic (`lambdaH=0.25`), other SM inputs, masses, thresholds, and mixed-scalar dimensionful bilinear prescriptions unchanged.
- Load new *individually fitted* D(-2) and A(-4) `Y1`, `Y2` matrices from **`uniform001_all16_reference_v1/reference/`**. These are 18 real Yukawa entries for each reference.
- For every target model **keep that source Yukawa pair fixed**; do not refit, reseed, or rescale the transferred Yukawas. Run the target's full-general UV/EFT1 RGE, existing conditional one-scale Matchete C5, and low-energy evolution with the existing local evaluators.
- At `100 GeV` report neutrino masses, `delta_m21_sq_ev2`, `delta_m3l_sq_ev2` (NO: `delta_m31_sq_ev2`), the three `sin^2(theta_ij)`, the five-observable Gaussian chi-squared, and the complete `C5` matrix. Large chi-squared values are expected/allowed for unrefitted transferred Yukawas; numerical `Success` is not fit acceptance.

### Weinberg matrix comparison

For each fixed reference `r = D(-2)` or `A(-4)` and target representation `a`, use the **low-scale Weinberg coefficient matrices** in the same numerical flavour basis:

\[
k_{a|r} = \frac{\mathrm{Tr}[(C_5^{(r)})^\dagger C_5^{(a)}]}{\mathrm{Tr}[(C_5^{(r)})^\dagger C_5^{(r)}]},\qquad
\epsilon_{a|r} = \frac{\|C_5^{(a)}-k_{a|r} C_5^{(r)}\|_F}{\|C_5^{(a)}\|_F}.
\]

This least-squares coefficient is **complex** in general. Table C.4 reports both real and imaginary parts, rather than silently dropping the phase. The row `a=r` has `k=1`, `epsilon=0` by construction. This is a comparison of unrotated saved `C5` matrices at `100 GeV`, not a statistical likelihood or a certified gauge-invariant equality between different representation bases.

### Important limits

This uses the project's **existing conditional one-scale C5 matching**, not a certified strictly sequential `UV -> integrate out ALL fermions -> integrate out ALL scalars -> SMEFT` matching calculation. B(-1)/C(-1) retain the existing mixed-scalar matching/sign convention and absolute +0.01 extra-quartic UV boundary.

## Install

Place `t3_uniform001_fixed_yukawa_source_fix.zip` at the project root and extract:

```powershell
Expand-Archive -LiteralPath .\t3_uniform001_fixed_yukawa_source_fix.zip -DestinationPath . -Force
python -X utf8 -m unittest discover -s tests/python -p "test_t3_uniform001_fixed_yukawa_tables.py" -v
```

This correction supplies complete replacements for the standalone report script and its tests, plus this updated README. It does NOT modify or overwrite any RGE, matching, symbolic or fit source files:

- `scripts/BuildT3Uniform001FixedYukawaTables.py` (full executable source)
- `tests/python/test_t3_uniform001_fixed_yukawa_tables.py`
- `Reports/output/thesis_general_5_1/uniform001_fixed_yukawa_tables/README_RUN.md`

## Commands

Read-only preflight validates the two new D(-2), A(-4) source fit configs, their source fingerprints and required target beta/matching inputs. No integration:

```powershell
python -X utf8 -m scripts.BuildT3Uniform001FixedYukawaTables --preflight
```

Generate the saved *individually fitted* reference table (no RGE):

```powershell
python -X utf8 -m scripts.BuildT3Uniform001FixedYukawaTables --reference-table-only
```

Cheap pilot: compute **four reference/model pairs** (two source diagonals reused, two cross-transfers evaluated):

```powershell
python -X utf8 -m scripts.BuildT3Uniform001FixedYukawaTables `
  --models A_-4,D_-2 `
  --run-id uniform001_Dm2_Am4_fixed_v1
```

Complete the remaining 28 cross-representation transfers, without rerunning successful saved points:

```powershell
python -X utf8 -m scripts.BuildT3Uniform001FixedYukawaTables `
  --models all `
  --run-id uniform001_Dm2_Am4_fixed_v1 `
  --resume --frozen-symbolic
```

If a run was interrupted, repeat that `--resume` invocation. For one-shot all-model execution, omit `--models A_-4,D_-2` from the first command and use `--models all`.

Render saved 32 outputs without rerunning the physics:

```powershell
python -X utf8 -m scripts.BuildT3Uniform001FixedYukawaTables `
  --run-id uniform001_Dm2_Am4_fixed_v1 `
  --render-only
```


## Repair for `source fingerprint differs (or missing provenance)`

The error comes from `scripts/RunT3GeneralSymbolicStudy.py:validate_cached`,
which compares the canonical UV/EFT1 symbolic export's saved SHA256 fingerprints
against **all current symbolic source files**, not from the fixed-Yukawa
transfer or its 18 Yukawa parameters. Current GitHub source itself differs
from the cached symbolic provenance for at least five core tracked files:
`T3ModelCatalog.wl`, `ScalarPotential.wl`, `SU2Invariants.wl`,
`T3RGBetaRunnerCommon.wl` and `general_rgbeta_quartic_registry.json`.
It is **not** safe to silently claim that these cached betas were regenerated
from the edited sources.

The standalone `--frozen-symbolic` switch is an **explicit opt-in** to
reproduce the physics of the already-accepted uniform-0.01 fits under their
**original cached beta/C5 snapshot**. It does not edit global symbolic
modules or regenerate betas. Before any nonmixed transfer it:

1. Calls the existing full-general symbolic cache structural validator against
   the cached export's *own* recorded snapshot hash dictionary (not the live
   dictionary). Structural/interaction/beta-inventory validation remains on.
2. Validates that the target model's own accepted uniform-0.01 reference fit
   hashes the same beta files, Matchete coefficient, matching status, original
   source files and merged source manifest currently in the checkout.
3. Evaluates the target's newly refitted 0.01 baseline independently with
   these exact cached exports. All five oscillation observables and the entire
   low-scale 3x3 complex C5 must reproduce the saved reference to within
   a relative tolerance of 1e-5 (observable absolute floor 1e-10).
4. Only then evaluates the fixed D(-2)/A(-4) transferred Yukawas with all
   additional BSM scalar quartics still at +0.01 and conjugate Yukawas at zero
   on the UV boundary.

This is a **pinned symbolic snapshot reproducibility calculation**, not
certification against the current Lagrangian/RGBeta sources. If any hash,
structural or baseline check fails, stop and investigate; do not change
provenance JSON to force acceptance. A genuinely fresh calculation from
current sources requires regenerating and validating the canonical symbolic
exports and then refitting references if the physics changes.

Resume your existing dataset (successful points are skipped, Failed points
recomputed):

```powershell
python -X utf8 -m scripts.BuildT3Uniform001FixedYukawaTables `
  --models all `
  --run-id uniform001_Dm2_Am4_fixed_v1 `
  --resume `
  --frozen-symbolic
```

If it completes, `table_C1_Dm2`, `table_C2_Am4`, `table_C3_PMNS`,
`table_C4_Weinberg_proportionality` and `table_C5_Weinberg_matrix_entries.csv`
are rendered automatically. To re-render completed data without physics:

```powershell
python -X utf8 -m scripts.BuildT3Uniform001FixedYukawaTables `
  --run-id uniform001_Dm2_Am4_fixed_v1 `
  --render-only
```

### Current validation

The repaired script passed 13 isolated unit tests. Full numerical transfers
were not run in the packaging environment, so saved failures may expose
additional physical or local-path problems after the initial fingerprint
check is resolved. No new numerical values are claimed until the user runs it.

## Outputs

Numerical checkpoints (32 points; `D_m2`/`A_m4` are the fixed source references):

`output/studies/thesis_qualitative/uniform_scalar_001/fixed_yukawa_transfer/uniform001_Dm2_Am4_fixed_v1/<target-tag>/<reference-tag>/result.json`

Final tables (CSV, Markdown, LaTeX):

`Reports/output/thesis_general_5_1/uniform001_fixed_yukawa_tables/uniform001_Dm2_Am4_fixed_v1/`

- `table_C1_Dm2` — 16 fixed D(-2)-Yukawa mass/chi-squared/C5 residual rows.
- `table_C2_Am4` — 16 fixed A(-4)-Yukawa mass/chi-squared/C5 residual rows.
- `table_C3_PMNS` — 32 mixing-angle rows.
- `table_C4_Weinberg_proportionality` — 16 rows, both `k` and `epsilon` for both references.
- `table_C5_Weinberg_matrix_entries.csv` — all 3x3 low-scale C5 entries (288 values per real/imag component; 32 x 9 rows).
- `table_R0_individually_refitted` — **separate** 16-model reference-fit table from existing saved results.
- `provenance.json` — result-file hashes and definitions.

The table generator **refuses** to produce the fixed-transfer tables if any source/model result is missing or failed. No results have been fabricated or inferred from the separately fitted models.

## Repair for `Missing full symbolic expression: interaction_inputform_file`

The canonical `lagrangian.json` records the three Matchete InputForm file
locations using **absolute Windows paths from the machine that generated the
symbolic export** (for example, `C:\thesistemp\...`). After moving the
repository, those exact old paths are absent even when all three InputForm
files are present beside the cached `lagrangian.json`.

The repaired standalone driver addresses **only this path relocation** when
`--frozen-symbolic` is explicitly requested:

1. Confirm that each archived pointer's file name is the expected
   `lagrangian.json.interaction.inputform.txt`,
   `lagrangian.json.uv.inputform.txt`, or
   `lagrangian.json.bsm.inputform.txt`.
2. Require each corresponding file in the target's saved canonical export
   directory to be present and nonempty. Verify the recorded SHA256 hashes of
   the original `lagrangian.json`, interaction and UV InputForm files,
   `uv_rgbeta.json`, and `eft1_rgbeta.json` against `status.json`.
3. Copy the saved symbolic metadata into a **temporary directory** and rewrite
   only those three file pointers inside the temporary copy. Call the
   existing `RunT3GeneralSymbolicStudy.validate_cached()` on that copy using
   the originally recorded source fingerprint. The original cached files,
   hash records and canonical validator remain unmodified.
4. Continue the previously implemented per-model accepted-reference source
   hash check and full numerical reproduction of its five observables and
   low-scale Weinberg C5 matrix before any fixed-Yukawa transfer.

This is a **pinned-cache reproducibility run**, not a certification that the
symbolic files reflect the latest edited Lagrangian/RGBeta sources. It refuses
missing, changed or unrecognised cached data instead of silently accepting it.

Install the new ZIP from the project root and resume the **same** run ID:

```powershell
Expand-Archive -LiteralPath .\t3_uniform001_fixed_yukawa_inputform_fix.zip -DestinationPath . -Force
python -X utf8 -m unittest discover -s tests/python -p "test_t3_uniform001_fixed_yukawa_tables.py" -v
python -X utf8 -m scripts.BuildT3Uniform001FixedYukawaTables `
  --models all `
  --run-id uniform001_Dm2_Am4_fixed_v1 `
  --resume `
  --frozen-symbolic
```

Successful points (B(-1)/C(-1) transfers and the two source diagonals) are
skipped; previously failed points are retried. If a hash or baseline physics
check now fails, stop and inspect it rather than skipping those checks. Full
cross-representation physics calculations were not executed in this packaging
environment. 17 isolated tests passed.
