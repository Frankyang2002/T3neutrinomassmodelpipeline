# Filename migration: historical Stage 23 report provenance

The research filename migration changes the Matchete tensor export script from
`scripts/ExportGeneralT3Stage22Matchete.wl` to
`scripts/ExportGeneralT3ScalarTensorMatchete.wl`. The original saved report
`output/general_t3_stage22/stage23_quartic_prefactor_candidates.json` still
records the old pathname in `source_verified.matchete_source` as historical
provenance. The current `scripts/QuarticPrefactorAudit.py` correctly describes
its renamed source instead.

`NativeQuarticActionAudit.py` now permits precisely this **single string
metadata difference** while continuing to:

- recalculate the entire Stage 23 report from Stage 22 inputs;
- validate their byte-level SHA-256 hashes via `stage23_run`;
- require exact equality for every other saved Stage 23 field;
- reject missing or tampered physics arrays, tensor blocks and status flags;
- perform the original native RGBeta assembly verification.

No saved JSON or physical source is rewritten. A failed hash check or any other
change to the candidate report still causes a fail-closed error. This is a
provenance-path compatibility exception, **not** a new Matchete/RGBeta coupling
certification.

## Local verification

```powershell
python -X utf8 -m pytest -q tests/python/test_renamed_quartic_action_provenance.py --tb=line
python -X utf8 -m pytest -q tests/python/test_general_t3_native_quartic_action.py tests/python/test_general_t3_quartic_action_identity.py tests/python/test_general_t3_triplet_intertwiner.py tests/python/test_general_t3_triplet_native_structure.py --tb=line
python -X utf8 -m pytest -q --tb=line
```
