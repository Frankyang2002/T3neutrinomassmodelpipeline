# Original simultaneous-EFT1 LLSS ablation — thesis-consistent source

## Installation

Copy full files into the exact relative paths:

- `scripts/RunT3OriginalSimultaneousLLSSHeatmap.py`
- `tests/test_t3_original_simultaneous_llss.py`
- `README_ORIGINAL_SIMULTANEOUS_LLSS.md`

This is an **additive** script. It does not edit production physics files, old source data, or GitHub.

## Physics definition

The original historical 16-model study is the source of both `1111` and `1001`.
In the original bit convention `U E W S`, `1001` turns off *both* renormalisable EFT1 and EFT1 Wilson running. UV running, simultaneous scalar matching, and SMEFT running remain on.

The new LLSS-off state reproduces the original simultaneous (`U=1,E=1,W=1,S=1`) geometry, apart from removing the **LLSS-to-LLHH** part of the quartic `scalar_pair` Wilson beta at *every* RGE RHS call. The original seven-term beta function, renormalisable RGBeta, mixed Cartesian scalar treatment, threshold matching, and SMEFT code are otherwise reused.

For each of the five observables:

`Running = 100*(1111 - 1001)/1001`

`Mixing off = 100*(LLSS_OFF - 1001)/1001`

The program refuses a missing baseline and refuses to overwrite existing output. It reuses original `cases.json` rather than the newer two-benchmark `thesis_two_benchmarks_all16_stagefreeze_v2` data.

## Run

From project root, PowerShell:

```powershell
.\.venv\Scripts\python.exe -m unittest tests.test_t3_original_simultaneous_llss

# Smoke test the original A(-4) historical source first
.\.venv\Scripts\python.exe scripts\RunT3OriginalSimultaneousLLSSHeatmap.py `
  --models A_-4 `
  --output output\studies\general_scalar_extended\reports\llss_original_smoke_A4

# Then all 16; this reruns original source-matched LLSS-OFF evolution only
.\.venv\Scripts\python.exe scripts\RunT3OriginalSimultaneousLLSSHeatmap.py `
  --output output\studies\general_scalar_extended\reports\llss_original_all16 `
  --font-size 16 `
  --cell-fontsize 12 `
  --top-tick-fontsize 18 `
  --bottom-tick-fontsize 15 `
  --y-tick-fontsize 15 `
  --figure-width 23
```

Defaults require these existing original studies:

- `output/studies/general_scalar_extended/experimental_simultaneous_eft1_seven_term/simultaneous_A_m4_v2`
- `output/studies/general_scalar_extended/experimental_simultaneous_eft1_general_nonmixed/simultaneous_nonmixed13_v1`
- `output/studies/general_scalar_extended/experimental_simultaneous_eft1_mixed_BC/simultaneous_mixed_BC_v1`
- source fit `output/.../thesis16_merged_v1`, as configured in original code

Output (new directory): `all16_three_way.csv`, `provenance.json`, `plots/all16_heatmap_percent.png`, `plots/all16_heatmap_percent.pdf`, and scratch proof-of-recomputation under `raw/`.

Replot after a complete calculation without rerunning physics:

```powershell
.\.venv\Scripts\python.exe scripts\RunT3OriginalSimultaneousLLSSHeatmap.py `
  --plot-only output\studies\general_scalar_extended\reports\llss_original_all16 `
  --output output\studies\general_scalar_extended\reports\llss_original_large_font `
  --font-size 18 --cell-fontsize 14 --figure-width 25
```

## Caveats

The original study is exploratory. The native Matchete-to-Weyl spinor-normalisation and full physical MSbar threshold matching are not independently certified. A successful source-contraction parity check does not certify them. In mixed cases `B_-1`/`C_-1`, this code uses the original fixed Cartesian mixed-scalar solver and refuses a zero LLSS source diagnostic; it does not claim independent certification of the mixed basis.

**Verification scope:** synthetic unit tests can verify subtraction algebra and plotting, not an entire 16-model Wolfram/native pipeline. Run the smoke test and inspect logs before treating the new numerical LLSS-off columns as established.

## Recovery from `source fingerprint differs (or missing provenance)`

The canonical symbolic loader compares current source SHA256 hashes with the
archived `complete_scalar_symbolic/<canonical_run_id>/<model>/provenance.json`.
The runner now makes an explicit **in-memory source-only exception** for
`RGE/running/rgbeta/RunGeneralT3QuadraticRGBeta.wl`, using the archived
fingerprint. All other symbolic files must match exactly. The original
`validate_cached` function still checks Lagrangian, extra quartic/Yukawa
inventories, UV/EFT1 beta exports, and stage statuses.

An audit is written under each group, for example:
`<output>/raw/a4/source_fingerprint_reuse_audit.json`.
If any other source has changed, the run stops with the actual file names;
there is no general-purpose provenance bypass.

For a new A(-4) smoke test, run from the project root:

```powershell
.\.venv\Scripts\python.exe scripts\RunT3OriginalSimultaneousLLSSHeatmap.py `
  --models A_-4 `
  --output output\studies\general_scalar_extended\reports\llss_original_smoke_A4_sourcefix_v2
```

After that succeeds, use a new `--output` directory for all 16 models.
Existing failed `llss_original_smoke_A4` and `llss_original_all16`
directories are not reused or overwritten.
