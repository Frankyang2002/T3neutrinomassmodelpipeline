# Workstream D — Stage 2: running fermion threshold diagnostics

Status: **candidate, no sequential matching, no physical certification**.

## Purpose and unchanged numerical sources

Stage 2 adds a separate read-only numerical audit; it does not edit the Stage 1 scalar adapter or these legacy files:

- `scripts/RunT3FullScalarContinuous.py` — complete UV candidate RGBeta state, initial quartics and quadratic, RHS, beta convention, original DOP853 tolerances;
- `Numerical/experimental/GeneralT3ScalarRGState.py` — pack/unpack with full 3x3 complex `MF`;
- `Numerical/core/StateVector.py` — `MF` packed as flattened real part followed by flattened imaginary part;
- `Numerical/running/FermionThresholdBoundary.py` — existing physical singular-mass definition;
- `scripts/AuditGeneralT3MassCrossings.py` — Stage 1 scalar mass crossings and grouped-F reference state.

A replay is necessary because `continuous_trajectory.json` saves only the *endpoint* fermion masses, not their UV history. Replay uses the **same existing input configuration and Stage 29 assignments**; it never adjusts benchmark inputs. The saved input SHA256 inventory is checked against present local files. The `UV` solution is compared with every saved scalar/quartic/bilinear sample and with the complete saved `MF`, `y1`, `y2` matrices at the grouped-F endpoint. Any disagreement fails closed.

## Numerical definition

The running state contains a complex `3 x 3` mass matrix `MF` with mass dimension 1 and a model-dependent fermion character. At any `t=ln(mu/GeV)` compute its ascending singular values

\[
0<M_0(t)<M_1(t)<M_2(t),\qquad \operatorname{svd}(M_F(t)).
\]

For Majorana `MF`, the code first verifies complex symmetry, so these singular values equal the *Takagi masses*. For vectorlike fermions they are the *bi-unitary SVD masses*. No field rotation or Yukawa matching is carried out. The label `fermion:i` means ascending **mass rank** at the specified scale, not a fixed flavour direction.

A candidate UV crossing solves

\[
M_i(t_i) = \exp(t_i)\ \mathrm{GeV},
\]

using the original full-UV `DOP853` dense interpolant, repeated 3x3 singular decompositions, and `scipy.optimize.brentq`. Brackets are discovered on a tunable logarithmic grid (default 65 points); root locations are **not** computed by linear interpolation between frozen mass eigenvalues. A collision/near-degeneracy between singular masses is rejected because the rank label is ambiguous. An inverse/repeated mass crossing is rejected. No crossing is extrapolated outside the original `[1e5,1e7] GeV` UV integration domain.

### Critical physics limitation

These are roots on the **un-decoupled full-UV equations**. Strict sequential RG running must integrate out the heaviest state at its actual threshold, apply certified matching, and then evolve the reduced-field EFT before determining the next threshold. Thus Stage 2 is a candidate mass-based *schedule diagnostic*, **not** a fully self-consistent sequential series of physical thresholds.

The original Stage 1 scalar crossings remain diagnostics on the saved grouped-fermion EFT1. They are not automatically recomputed after the new fermion roots, and the merger into `candidate_event_order` does not authorize removal of any field.

## Input/output contracts

### Inputs from existing code and other workstreams

| Provider | Contract | Status |
|---|---|---|
| A | Ordinary T3 model identifier; SU(2), hypercharge, fermion/Majorana nature, inventory | Taken from existing representation/manifest; no new inventory invented |
| B | Verified basis/coupling map and complete UV/reduced-EFT beta exports, convention `16*pi^2 dX/dln(mu)` | Existing UV candidate-basis beta export is replayed; other EFT beta maps not inferred |
| C | Matching functions and normalization for each threshold | **Missing / not called** |
| D | Continuous UV state (all real/complex blocks); saved source hashes; `[UV start, grouped-F]` interval; 3x3 `MF` | Input validated/replayed |

All mass scales are GeV; `MF` has dimension 1, scalar bilinear `mu2` dimension 2. No implicit normalization map between RGBeta and Matchete is assumed.

### New JSON outputs

Written to `output/studies/general_scalar_extended/continuous_full_scalar/continuous_v1/workstream_d_fermion_audit_v1/audit_summary.json` by default. **No overwrite** if the file exists; pass an alternate `--output` for a new run.

Top level: `format`, `status`, `physical_matching_certified: false`, `models`.

Each model has:

- `uv_fermion_crossings`: ordered records of `particle`, `kind`, `scale_gev`, `mass_gev`, root residual, mass spectrum at root, method, and `candidate_crossing_matching_required_not_applied`;
- `unresolved_fermion_crossings`: mass branches not bracketed in the **saved** UV interval, explicitly including `below_saved_UV_stop_undetermined` if relevant;
- `eft1_scalar_crossings` and `uv_scalar_crossings`: reused Stage 1 candidate diagnostics;
- `candidate_event_order`: decreasing found crossing scales, **not** an executed EFT sequence;
- `full_uv_mass_diagnostic`: replay numerical checks, `uv_replay_nfev`, `original_continuous_source_hashes_verified`, and flags that no matching/decoupling was performed.

`all_fermion_crossings_bracketed` is false when any crossing lies outside the domain; this is not silently treated as successful threshold matching. Missing physical matching is **not** numerically replaced by zero.

## Installation and execution in PowerShell

Unzip the Stage 2 package at the **existing** root of the local checkout:

```powershell
$repo = 'C:\thesistemp\T3neutrinomassmodelpipeline'
$zip = "$env:USERPROFILE\Downloads\T3_WorkstreamD_RunningFermion_Stage2.zip"
Expand-Archive -LiteralPath $zip -DestinationPath $repo -Force
Set-Location $repo

python -m unittest discover -s tests -p 'test_general_t3_running_fermion_thresholds.py' -v
python -m unittest discover -s tests -p 'test_general_t3_spectrum_threshold_adapter.py' -v
python -m unittest discover -s tests -p 'test_general_t3_bc_conditional_regression.py' -v
python -m scripts.AuditGeneralT3FermionCrossings --all
```

If Stage 2 has already written its default result, do not delete the original; choose e.g. `--output 'output/studies/general_scalar_extended/continuous_full_scalar/continuous_v1/workstream_d_fermion_audit_v2/audit_summary.json'`.

## Scope stop

**Done in Stage 2:** running complex MF extraction, physical singular-mass definition, UV replay/source parity checks, candidate F crossings and 16-model schedule diagnostics.

**Not done:** fermion EFTs between individual thresholds, certified matching, revised scalar RG trajectory after removing separate fermions, EW-broken/pole scalar masses, finite Y/X or LLHH matching, physically certified low-energy neutrino masses. Existing B/C conditional outputs are neither overwritten nor retuned.
