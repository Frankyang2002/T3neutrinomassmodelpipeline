# Workstream C v10.1 — GammaCC source locator (v10 loader correction)

This replaces the v10 Wolfram source probe, which failed at a bare `Needs["Matchete`"]` call even though previous project scripts loaded Matchete successfully inside `UsingFrontEnd`.

## Change

- First search the **installed package source** using `FindFile["Matchete`"]`, local `PacletFind`, and Matchete-specific `$Path` directories. No package initialization is needed for this path.
- Only if passive discovery fails, retry using the project's previously successful `UsingFrontEnd[Needs["Matchete`"]]` loader and repeat source discovery.
- Permit a second, optional Matchete package source root argument if the automatic search fails.
- Scan only Matchete source roots for `GammaCC`, `CConj`, `Transp`, `DiracProduct`, and `Proj`, with bounded excerpts.
- Preserve the existing JSON contract `t3-bc-installed-gammacc-source-probe-v1`, so the **existing v10 Python script** remains unchanged.
- This is source inspection only: no model is loaded, no matching or counterterm is performed, and no GammaCC sign is inferred.

## Commands (PowerShell, project root)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$HOME\Downloads\T3_WorkstreamC_GammaCCCalibration_v10_1.zip" -DestinationPath . -Force
python scripts\CalibrateT3BCMajoranaNormalization.experimental.py --self-test
$dir = "output\audits\t3_bc_y_gammacc_convention_v1"
New-Item -ItemType Directory -Force $dir | Out-Null
wolframscript -file Lagrangian\ProbeT3BCInstalledGammaCC.experimental.wl "$dir\installed_source.json"
if ($LASTEXITCODE -ne 0) {
    Get-Content "$dir\installed_source.json"
    throw "Installed GammaCC source locator failed"
}
python scripts\CalibrateT3BCMajoranaNormalization.experimental.py --require-installed-probe --output "$dir\summary.json"
```

If automatic location fails, you can provide an explicit installed Matchete source directory as a second argument:

```powershell
wolframscript -file Lagrangian\ProbeT3BCInstalledGammaCC.experimental.wl `
  "$dir\installed_source.json" "C:\PATH\TO\Matchete"
```

Use the actual Matchete package directory, **not** the entire Wolfram `Applications` directory.

The package can print a source-location failure JSON even if the optional front-end loader fails. This keeps the diagnostic explicit rather than reporting a successful physical normalization.

## Interpretation

A successful result means only that source snippets containing `GammaCC` were located. The source excerpts need to be examined before claiming any exact identification between Matchete's `NCM[... GammaCC ...]` and the Warsaw Weinberg operator. The v9.1 native source skeleton and the v10 Grassmann/VEV factors remain valid, but the physical C5 factor is **NOT CERTIFIED**.

## Files

- `Lagrangian/ProbeT3BCInstalledGammaCC.experimental.wl` — **full replacement** of v10 source probe.
- `docs/WorkstreamC_GammaCCCalibration_README.md` — **full replacement** of v10 documentation.

The Python script from v10 is unchanged. No production file, Matchete package installation, or remote repository is modified.
