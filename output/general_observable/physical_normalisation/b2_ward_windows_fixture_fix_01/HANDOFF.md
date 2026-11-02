# B2 Ward witness: Windows Git fixture correction

Scope: synthetic pytest fixture only. No production, source provenance, RGBeta registry,
Matchete invariant, or physical coupling calculations have been modified.

## Observed user results

- Existing B2 gauge tests: previously 16/16 passed on Windows.
- B2 extended native RGBeta registry: exported, formal classes 5/5.
- B2 Ward witness: conditional status `LITERAL_RGBETA_DOUBLE_DELTA_WARD_MISMATCH_CONDITIONAL`;
  formal classes 5/5, B/C models 6/6, Matchete comparisons 6/6;
  physical certificates 0/16.
- Ward witness pytest: 36 passed, one Windows-specific synthetic Git fixture failure.

## Verified issue

The old `_mini_git_repo` test fixture used the raw bytes of a newly written
working-tree file (`git_blob(path.read_bytes())`) as the *expected committed*
Git source blob. On Windows, Python text writes can produce CRLF and Git with
`core.autocrlf` converts that content to LF in the index/commit. The expected
fixture SHA then differs from the real HEAD tree and the production-safe
source guard correctly refuses promotion.

## Correction

- The fixture writes explicit bytes, making platform behavior deterministic.
- The fixture sets `core.autocrlf=true` and obtains its expected reference via
  `git rev-parse HEAD:<path>`, after committing.
- It checks that `git hash-object --path=<path> -- <path>` matches the HEAD blob.
- New regression test exercises CRLF raw worktree bytes and verifies the
  canonical Git hash is accepted.

All **38/38** Python tests passed in the local Linux test environment
(22 Ward witness tests + 16 existing B2 tests). The corrected file must be
run on the user's Windows system for local confirmation.

## Files

Complete replacement:

`tests/test_physical_b2_rgbeta_ward_witness.py`

This report:

`output/general_observable/physical_normalisation/b2_ward_windows_fixture_fix_01/HANDOFF.md`

No native Wolfram scripts need rerunning; preserve existing outputs.
Do not infer Matchete–RGBeta physical normalization from these tests.
