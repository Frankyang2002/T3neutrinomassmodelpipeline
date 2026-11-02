# Gaussian LLSS boundary label fix

The EFT1 fermion boundary releases use two distinct representations of a model
label:

- Directory name / filesystem slug: `A_m4`, `B_m1`, `E_p2`, etc.
- The `model` property inside `boundary.json`: `A_-4`, `B_-1`, `E_+2`, etc.

The previous `gaussian_amplitude.py` mistakenly compared `model` with the
filesystem slug, causing all 16 real-input runs to block with
`F threshold source release model label disagrees`.

The replacement keeps the filesystem slug to find `boundary.json` but checks
its internal `model` value against the signed model label. Inconsistent or
swapped signed model records still cause a hard error. The threshold-ordering
checks are unchanged.

This correction does not modify any Gaussian/loop tensor calculation, frozen
benchmark, existing EFT release, or source-to-Weyl normalization. The bridge
still operates at the frozen UV input scale (normally `1e7` GeV), while the
EFT1 boundary is at `mu_F=1e5` GeV. No same-scale EFT comparison or physical
UV/EFT equivalence is claimed. Hard one-loop fermion matching remains absent.

Install by extracting the ZIP with `-Force`, then run the regression tests and
the bridge with a fresh output directory (e.g. `gaussian_amplitude_02`). Do not
run in `gaussian_amplitude_01`: the runner is intentionally no-overwrite.

The replacement also corrects the existing
`tests/test_lo_direct_uv_gaussian_amplitude.py` fixture: it now writes
`"model": "B_-1"` inside the `B_m1` release directory. Seven additional
regression tests cover all 16 real release label conventions, bad/swapped
labels, simultaneous threshold ordering, non-comparison at unequal scales,
and an all-16 `batch()` integration with the Gaussian numeric kernel stubbed.
The existing 60 unit tests should be rerun in the local repository.
