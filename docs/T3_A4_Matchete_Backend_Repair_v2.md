# A(-4) launcher repair v2

The previous `wolframscript -code` invocation lost the double quotes around
`FindFile["Matchete`"]` when called from Windows PowerShell. Mathematica
received the invalid expression `FindFile[Matchete`]`, so the launcher stopped
before creating any run output.

This full replacement writes Wolfram code to a temporary `.wl` file and uses
`wolframscript -file`. Both the version probe and the Matchete upgrade
avoid passing quoted Wolfram expressions through a Windows command line.
If an upgrade is necessary, the launcher downloads the **official upstream
master archive** and validates its version **before** backing up and replacing
the local installed package. The old package is backed up beneath
`output/local_matchete_backups/`.

Project paths:
- `scripts/InvokeT3A4ScalarMatcheteRepair.ps1` (complete replacement)
- `docs/T3_A4_Matchete_Backend_Repair_v2.md` (new)

Run from project root:

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$HOME\Downloads\T3_A4_Matchete_Backend_Repair_v2.zip" -DestinationPath . -Force
powershell -NoProfile -ExecutionPolicy Bypass -File scripts\InvokeT3A4ScalarMatcheteRepair.ps1 -InstallMatchete
```

Expected run path:
`output/parallel_t3/full_general_final_one_batch/thesis_a4_m061_scalar_v2/A_-4/`

Inspect `native_two_threshold_status.json` (if generated), and the matching
log. A nonzero status does not by itself mean the repair script failed:
the T3 wrapper deliberately refuses to declare a completed physical C5
until the F-hard Wilson boundary, EFT1 Wilson RG and finite matching exist.

This fix does not change any Matchete CG tensor, T3 Lagrangian, or physics source.
It has been statically reviewed here, but **not run on Windows with Wolfram**.
