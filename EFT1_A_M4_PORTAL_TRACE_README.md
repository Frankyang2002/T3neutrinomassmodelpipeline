# A(-4) full-general portal-retention trace

**Project files:** `Lagrangian/TraceAminus4PortalRetention.wl` (new, additive); this README at repository root (new, additive).

**Purpose:** Identify whether `lambdaH2Inv1` disappears during the *assembled general UV Lagrangian*, Matchete's `Match`, `GreensSimplify`, `EOMSimplify`, effective-coupling replacement, or the later compact-tree/export step. The script uses the original checked `T3ModelCatalog`, `LagrangianBuilder`, `RunMatching`, and Matchete routines, in a fresh Wolfram kernel. It retains all four general-operator switches and requires only `F` to be heavy. It does **not** modify Lagrangian, threshold, or RGE production files.

From the repository root in PowerShell:

```powershell
Set-Location 'C:\thesistemp\T3neutrinomassmodelpipeline'
$env:T3_EXPERIMENTAL_GENERAL_FERMIONS  = '1'
$env:T3_EXPERIMENTAL_GENERAL_QUARTICS  = '1'
$env:T3_EXPERIMENTAL_GENERAL_QUADRATICS = '1'
$env:T3_EXPERIMENTAL_FLAVOR_MASSES    = '1'

wolframscript -file Lagrangian\TraceAminus4PortalRetention.wl `
  output\eft1_general_lagrangian_trace\A_m4_fresh_01

Get-Content 'output\eft1_general_lagrangian_trace\A_m4_fresh_01\portal_retention_trace.json' -Raw
```

If that output folder exists, change `A_m4_fresh_01` to a new run ID; the diagnostic refuses to overwrite it.

**Outputs (all under `output/eft1_general_lagrangian_trace/A_m4_fresh_01/`):** `portal_retention_trace.json`, `LUV.inputform.txt`, each of `RawEFT`, `GreenEFT`, `EOMEFT`, `LoopEFT`, `MatchedEFT`, `CompactExplicitEFT`, and `CanonicalTransitionEFT` as `.inputform.txt`.

**Interpretation:**

- If `both_portals_present_in_assembled_LUV` is false, the problem is in assembly/source registration, *before* `Match`.
- If both appear in LUV but one is absent from `CompactExplicitEFT`, compare stage text exports. Intermediate expressions may hide native coupling names in Matchete effective-coupling structures; a direct string absence is not a proof of tensor loss.
- If the fresh compact EFT contains both but the old `eft1_after_F_tree_inputform.txt` does not, the old saved export is stale or its production/export path needs correction.
- `source_hashes_vs_saved_run` compares archived SHA-256 hashes with the **actual local** source to detect stale saved data.

Do **not** invent a term, set a missing coupling to zero, substitute a restricted Lagrangian, or regard a successful run as certified physical `C5`. Once the failing stage is identified, repair that stage and rerun the *complete general F-only matching*, then recheck the matched tree and renormalizable coupling dictionary. The remaining GammaCC-to-Weyl factor and scalar-basis mapping remain separate certification tasks.
