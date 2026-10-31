
## 2026-10-31 — Experimental B/C running-Y Wilson tree projection (awaiting local physics run)

- New **experimental**, read-only source:
  `Numerical/experimental/T3BCScalarYThresholdProjection.py` and
  `scripts/RunT3BCYThresholdProjection.py`.
- Inspected public source of the existing coupled Wilson runner, 12-real-field
  tensor engine and complex doublet mass basis; no production source changed.
- The new runner reuses saved `continuous_v1` EFT1 boundary and existing
  hard-F / Stage17 source parity; recomputes its coupled RGBeta background,
  retains the running complex `mu2`, locates its own self-consistent `muY`
  and integrates the coupled Wilson tensor down to that scale.
- Explicit \(12\times12\) orthogonal real map for
  \((S_1,\epsilon S_2^*)=U(X,Y)\), preserving epsilon, complex mixing phases,
  full six physical lepton-flavor pairs and `(L1,L2,eC)` fermion indices.
- Implements \(C^{XY}_{ijABf}=T_{aA}T_{bB}C^{old}_{ijabf}\) and extracts
  all `HH`, `HX`, `XX` blocks in an \(3\times3\times8\times8\times6\)
  **tree-level** EFT2 tensor, not just its Weinberg component.
- 11 standalone NumPy tests passed, including epsilon/phase/ordering negative
  controls and SU(2)×U(1) scalar-covariance checks.
- **Pending user-local run:** `continuous_v1` B/C numerical source data and
  exact Y crossings, Wilson ODE, archived NPZ and source-check outputs.
- Standalone `bc_separated_thresholds_v1` Y-crossing numbers are not imported:
  that study uses a different candidate scalar trajectory. Its latest
  local-only result schema could not be inspected in this environment.
- **Unresolved:** finite one-loop Y matching `Delta C_Y`, EFT2 (H,X) beta
  functions, self-consistent EFT2 X crossing, X matching and neutrino masses.
  Explicit flags remain `False`; no physical certification is claimed.
- Full derivation, status and command sequence:
  `docs/T3_BC_Y_THRESHOLD_PROJECTION_TECHNICAL.md`.
