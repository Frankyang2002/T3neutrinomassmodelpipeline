# B2 Git-root provenance guard — isolated correction

Reference GitHub main: `0b210f18f6a2cfb10fff03059e2ec4a73d2195a4` (reported user local HEAD also matches).

## Observed local result, supplied by user

- B2 local audit `b2_local_02`: `CONDITIONAL_NATIVE_CG_SNAPSHOT_LOADED`, `models=16`, `physical certificates=0`.
- User local pytest before this correction: 14 passed, one failed at `current_git_head_measured` in `test_untrusted_value_and_drift_detected`.
- The original assumption that an arbitrary `tmp_path` cannot inherit a Git HEAD is unsafe: Git searches ancestors and may be influenced by environment repository settings.

## Corrective scope

- `Numerical/lo_certification/physical_normalisation/b2_gauge_coupling.py`: source_inventory only trusts `git rev-parse HEAD` if `git rev-parse --show-toplevel` resolves to the exact requested project root. Preserve path-specific Git clean-filter CRLF acceptance, genuine drift rejection, and all physical gate flags.
- `scripts/RunPhysicalB2GaugeDictionary.py`: does not overwrite the validated Git head provenance using an independent unanchored Git command.
- `tests/test_physical_b2_gauge_coupling.py`: add explicit nested-directory-in-ancestor-repository negative control and check root validation in existing fixtures.

## Evidence and limitations

- Synthetic pytest tests in this session: 16/16 passed. Python syntax compile: passed.
- Installed RGBeta/Matchete native execution in this session: not performed.
- Native physically certified coupling ratios: none. `physical_source_certified` remains `false`.
- Safe to integrate: Git root provenance guard and existing conditional CG audit data. No physical action dictionary is certified.

Use a NEW output directory for any rerun. Do not overwrite `b2_local_02`.
