# B2 G2 stage 1 — conditional pseudoreal SU(2) intertwiner

**Status: read-only exact mathematical witness; physical certification 0/16.**

Read from the six source-controlled Matchete/RGBeta source blobs, the previous
Matchete all-16 CG audit (`output/general_observable/lo_direct_one_loop/source_audits`),
the *existing* native RGBeta formal invariant report, and the successfully
validated `_05` RGBeta internal registry JSON. No Mathematica rerun or EFT.

The Windows `eps[SU2L[fund],1,2]` and `delS2[SU2L,1,1,2]` probe returned
unevaluated symbolic expressions. RGBeta's exact abstract contractions do
**not** uniquely select a numerical basis. Therefore this batch tests a
**hypothesis**, not a measured RGBeta index orientation.

For B/C, the formal mixing tensor (order `H,H,S1,Bar[S2]`) is

\[
D_{h k s b}=\tfrac12(\delta_{hb}\delta_{ks}+\delta_{kb}\delta_{hs}).
\]

With the usual naive `(fund,fund,fund,anti)` interpretation, the exact squared
Ward residuals are `(3,0,3)`. If the **formal** S1 index is instead treated as
antifundamental (so both its conjugates transform accordingly), the same
untransformed tensor is a singlet under `(fund,fund,anti,anti)` with `(0,0,0)`.
The B `y1` formal epsilon on `(Bar[S1],l)` is then interpreted `(fund,fund)`
and has `(0,0,0)`; under the naive `(anti,fund)` interpretation it has `(2,0,2)`.

Choose a pseudoreal intertwiner
\[ J=\begin{pmatrix}0&1\\-1&0\end{pmatrix},\qquad
J T_a^* J^{-1}=-T_a.\]

Transporting the formal S1 slot to a standard fundamental yields
\[
\sum_r D_{h k r b}J_{rs}=
\tfrac12(\delta_{hb}J_{ks}+\delta_{kb}J_{hs}),
\]

exactly the previous delta–epsilon candidate. Likewise transporting the
formal `Bar[S1]` epsilon gives `J J = -I`, or `-delta`, on `(anti,fund)`.
Both transported structures have three exact zero Ward residuals.
The transpose `J.T=-J` reverses their overall signs without spoiling the Ward
identity; this freedom prevents identification of an absolute physical phase.

The report also re-reads the Matchete source generator diagonaliser and checks
`U3.T U3` as the triplet weight-basis invariant metric. It compares six
recorded native Matchete B/C mixing CG tensors to this *conditional*
intertwiner candidate, recording exact proportionality and source hashes.
The 16/16 Matchete CG covariance records are source evidence; they do **not**
measure RGBeta numeric triplet axes.

**Not proved:** that installed RGBeta assigns these field orientations, a
numeric `delS2` or `eps` basis, Matchete↔RGBeta physical coupling conversion,
fermion spinor phase, action factors, all-five-class numeric basis equivalence,
or a source bug. Do not modify production couplings on these results.

Files added only:

- `Numerical/lo_certification/physical_normalisation/b2_su2_intertwiner.py`
- `Numerical/lo_certification/physical_normalisation/B2_SU2_INTERTWINER_METHOD.md`
- `scripts/RunPhysicalB2SU2Intertwiner.py`
- `tests/test_physical_b2_su2_intertwiner.py`

The script defaults to output
`output/general_observable/physical_normalisation/b2_su2_intertwiner_local_01`.
Refuses an existing directory. Use `--output` for a new suffix.

## Later local Git HEADs (G2 continuation)

The `_05` RGBeta native witness was captured at
`1fb8f4696a5568ec80f5f6dd5c8d6c61780ce40c`. Subsequent Git commits
may add unrelated work without altering the six pinned production source blobs.
The initial Windows command's exact-HEAD condition was too restrictive for
**reusing** that existing witness in a later read-only G2 audit.

The G2 runner now enforces all of the following before accepting archived G1
measurements: the actual current checkout passes `source_provenance(root)`
(HEAD/index/clean-filter and six source contracts), the original G1 capture
HEAD is a real Git ancestor of current HEAD, its *historical commit tree*
contains exactly the same six pinned source blob objects, the three source
hashes recorded in the G1 native JSON equal the pinned hashes, and the
current working-tree source files' Git clean-filter hashes still match.
Finally the entire existing G1 JSON is validated for 30/30 registry ratios.
The original G1 `native_snapshot_guard()` stays strict for same-HEAD capture.

This guards **historical evidence reuse**, not a new native RGBeta execution
at the current HEAD. It does not certify installed RGBeta's component-axis
convention or physical coupling factors. The generated JSON explicitly
records both original capture HEAD and current checkout HEAD.

Do not use the old PowerShell `if ($head -ne "1fb8...") { throw ... }`
pre-check. Use the revised batch ZIP and run the G2 launcher; it performs
the stronger ancestry and pinned-source checks internally. If the repo was
rebased, the current worktree has changed sources, or the original capture
commit is missing, this batch deliberately blocks rather than quietly
reusing an unproven native report.
