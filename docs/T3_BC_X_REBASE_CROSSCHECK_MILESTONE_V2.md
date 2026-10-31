
## Experimental B/C X-heavy continuation exact rebase crosschecks (v2)

- **Previous local evidence:** both Y-only full native matching runs completed; B(-1) light-X EOMEFT and MatchedEFT comparator both passed six full Matchete checks, but B(-1) X-heavy rebase rejected. Missing displayed rebase properties left the failed check unresolved.
- **New code:** full replacement `Lagrangian/RunT3BCXYNativeFollowup.experimental.wl` checks the unchanged Y-tree EOMEFT in the heavy-X registry, separately checks the heavy-free-rebased expression and exports exact free-Lagrangian difference. Matching accepts a rebased expression only with six passed checks, or unchanged EOMEFT only when six checks pass **and** the free-Lagrangian difference is exactly zero. No `MatchedEFT` substitution or operator deletion.
- **Runner:** `scripts/RunT3BCXYNativeFollowup.py` adds both-model `--diagnose-x`, fresh v2 outputs, and strict source hashes. `CollectT3BCXYNativeEvidence.py` targets successful v2 native amplitudes.
- **Offline tests:** 24/24 passed (existing 14 regression tests + 10 new tests). Windows Wolfram/Matchete still **not run** for v2.
- **Physical certification:** **false**. Saved F-first EOMEFT equivalence, Y finite loop inheritance, running threshold coupling substitution, and renormalized X/Y matching all remain missing.
