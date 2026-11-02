# Install and run isolated B1 first batch

## Intended project paths (all NEW)

1. `Numerical/lo_certification/physical_normalisation/b1_matchete_weyl.py`
2. `Numerical/lo_certification/physical_normalisation/B1_MATCHETE_WEYL_DERIVATION.md`
3. `scripts/RunPhysicalB1MatcheteAction.wl`
4. `scripts/RunPhysicalB1WeylAudit.py`
5. `tests/test_physical_b1_matchete_weyl.py`
6. `output/general_observable/physical_normalisation/b1_01/matchete_weyl_rules.json`
7. `output/general_observable/physical_normalisation/b1_01/all16_component_check.json` (**pending, not fabricated**)
8. `output/general_observable/physical_normalisation/b1_01/HANDOFF.md`

No shared source, benchmark JSON, mass engine or CG audits are changed. GitHub remains read-only.

## PowerShell, from repository root

```powershell
Set-Location 'C:\thesistemp\T3neutrinomassmodelpipeline'
# Record actual local source revision before installing. Do not use GitHub commit as local truth.
git rev-parse HEAD
git status --short
# Unzip this file to a NEW tree; do not overwrite any existing B1 work.
Expand-Archive -LiteralPath 'C:\PATH\TO\B1_Matchete_Weyl_FirstBatch.zip' -DestinationPath '.'

# Exact synthetic tests; expected 13 tests passing on Python dependencies numpy/sympy.
.\.venv\Scripts\python.exe -m unittest discover -s tests -p 'test_physical_b1_matchete_weyl.py' -v
# Short chosen-convention check; does not need Wolfram.
.\.venv\Scripts\python.exe scripts\RunPhysicalB1WeylAudit.py --self-test

# Native probes: separate fresh wolframscript processes with FOUR experimental flags forced 0.
# These outputs do NOT claim a concrete numerical GammaCC phase.
$native='output\general_observable\physical_normalisation\b1_02\native'
New-Item -ItemType Directory -Force $native | Out-Null
wolframscript -file scripts\RunPhysicalB1MatcheteAction.wl A_p0 "$native\A_p0.json"
wolframscript -file scripts\RunPhysicalB1MatcheteAction.wl B_m1 "$native\B_m1.json"
wolframscript -file scripts\RunPhysicalB1MatcheteAction.wl C_m1 "$native\C_m1.json"

# Execute independent 16-model reconstruction and frozen-candidate comparison.
.\.venv\Scripts\python.exe scripts\RunPhysicalB1WeylAudit.py --native-dir "$native"
# Output: output\general_observable\physical_normalisation\b1_02\all16_component_check.json
```

The native driver outputs `NativeActionTermsIsolated_PhysicalDictionaryPending` **only** when the action is accepted and original/H.c. terms were distinguished. Native failure exits nonzero and yields `Blocked` or fails before writing. Matchete is not available in the generation environment, so these native commands are **untested**; inspect each `mass_terms_inputform` and both `original_terms_inputform` / `hc_terms_inputform` before any physical inference.

The Python all-16 audit exits 0 only when all required native CG support files, benchmarks and provisional inputs exist; their model/provenance and source Git blob guards match; and the fixed-basis tensors reproduce the provisional builder up to a global neutrino-field sign. The successful status is **`All16ConditionalCandidateWeylTensorsReconstructed_PhysicalDictionaryPending`**, never a physical source certificate. Exit 2 writes a `Blocked` report explaining missing inputs or source drift. Exit 3 refuses to overwrite an existing result; choose a new `--output-dir` when rerunning. If local sources are newer than pinned main, supply the exact drifted source file(s) to revise the B1 audit rather than disabling guards.

To run without native files, omit `--native-dir`, and the program explicitly reports `native_action_observed=false`; all-16 conditional results may still be computed. Do not confuse this with native verification. Copy both `b1_02/all16_component_check.json` and the three `b1_02/native/*.json` outputs back to B1 for review. No EFT tasks should be run as part of this work.
