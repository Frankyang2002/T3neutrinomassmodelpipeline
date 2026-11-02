# B1 Windows EOL guard fix — handoff

## Source provenance

Reference source commit `0b210f18f6a2cfb10fff03059e2ec4a73d2195a4`.
GitHub `main` at inspection had advanced by one commit to
`5e88158147bfce94e34da534fa0be3c905561e1c`; the original seven
physics source file blobs remained unchanged in that commit.

The user observed for `Lagrangian/model/SU2Invariants.wl`:
`git ls-files --eol` => `i/lf w/crlf`;
raw blob `631781061c700b9b70b3d3f0b7cae051dd43e16a`;
Git-filtered blob `dbe0cd490204dec363fcf4278f079b58c8ddcefc`;
`git status --short` empty. This is an EOL-only mismatch.

## Change

Only `scripts/RunPhysicalB1WeylAudit.py` is replaced. The source guard
now checks the exact pinned Git blob after `\r\n` -> `\n`, refuses
lone `\r` bytes, and records both raw and normalized digests. It does
not call `git checkout`, alter source files, or weaken the guard for
operator/coefficients changes. `tests/test_physical_b1_line_endings.py`
adds focused positive and negative tests. No other production, B2, A,
RGE, EFT, benchmark, CG audit, or native driver file is modified.

## Outcomes / certification

Independent LF/CRLF synthetic fixtures pass; actual Windows repository
audit was not available in this environment. This change is a
source-verification implementation fix only, not a new physical proof.
`physical_source_certified = false`; no models promoted.

## Required run

From the local repository, extract ZIP with `-Force` for the B1 runner,
run the tests, then execute the audit with the previously generated
`b1_02/native` evidence into a new directory. Do not overwrite `b1_02`.
