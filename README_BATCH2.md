# T3 Batch 2 — source-evidence collector (read-only)

Install `scripts/CollectT3FermionMatchingEvidence.py` in the project at
`C:\thesistemp\T3neutrinomassmodelpipeline\scripts\CollectT3FermionMatchingEvidence.py`.

This script copies selected Matchete source ledgers / Laurent audit metadata from
local outputs into a ZIP. It deliberately does not construct the finite physical Wilson
boundary, run the solver, refit Yukawas, or touch GitHub. The hard matching is still open.

Run from PowerShell (full commands in conversation). Check the manifest for missing
required files or omissions due to the 12 MiB per-file and 65 MiB collection caps.
Source files are hashed and stored at their original project-relative paths.
