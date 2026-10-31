# Workstream D — Stage 3 universal EFT handoff contract

**Status:** Executable numerical scaffolding and candidate-data readiness audit. **Not a physical sequential T3 computation.** No symbolic source, RGBeta export, fitted UV input, or historical B/C runner has been changed.

## Project paths

- `Numerical/experimental/GeneralT3EFTHandoffController.py`: generic stage and handoff execution engine.
- `scripts/AuditGeneralT3EFTHandoffReadiness.py`: read-only candidate crossing inventory and dependency report.
- `tests/test_general_t3_eft_handoff_controller.py`: synthetic controller tests and candidate-audit contract tests.
- `docs/WorkstreamD_Stage3EFTHandoffContract.md`: this contract.

## Design

A `StageContract` declares the current EFT stage, ordered active particle labels and all numerical coordinate/Wilson blocks. Each `BlockSpec` declares its *exact* shape and real/complex nature. A symmetric block checks transpose symmetry `C == C.T`, not Hermiticity. `C5` specifically requires a complex symmetric `(3,3)` flavour array. All Wilson blocks must specify `operator_basis`, `normalization` and `mass_dimension`; these are **declared metadata**, not independently verified physical normalizations. Scalars can use shape `()`, real or complex; scales are explicitly in GeV.

An `EFTState` holds `model`, `stage`, `mu_gev`, ordered `active_particles`, `parameters` and `wilson`. It is a numerical envelope. When coupling to the existing `T3UVState`, `GeneralT3ScalarRGState` or B/C coefficient tensors, a source-aware adapter must convert between the original representation and this envelope; **that adapter has not yet been implemented or certified**.

The controller has **no built-in RG or matching callbacks**. Callbacks are registered explicitly:

- `controller.add_stage(contract)` registers a complete numerical stage.
- `controller.add_runner(stage, Callback(...))` advances only that stage to the requested endpoint.
- `controller.add_matcher(from_stage, particle, next_stage, Callback(...))` handles a boundary that removes **exactly one** active particle. It may change all block names/dimensions, but the new state must exactly satisfy the registered destination contract.
- `controller.run(initial, events, final_scale_gev=...)` checks scales/events before invoking anything, runs to a threshold, and either applies the declared matching callback or returns a `Blocked` result with the *pre-matching* state. There is no implicit tree projection, C5 matching, change of beta function, or field removal.

The returned `HandoffResult` contains the last validated state, operation history (including callback names/sources), the next event when blocked, and `physical_matching_certified=False`. Synthetic tests that successfully execute a chain do **not** certify its physics. The callback `status` field may say `candidate` or `externally_verified` for provenance; neither makes the controller itself a certification authority.

### Crossing roots versus matched thresholds

Stage 2 finds candidate roots using the **undecoupled full UV RGEs** for all fermions and the **grouped-F EFT1 background** for scalar branches. For an actual descending sequence, once F2 is removed, both the next mass root and all couplings must in principle be recomputed in the correct F2-removed EFT. For this reason, the Stage 3 audit importer only produces a `CandidateInventory`; it never constructs an EFT stage assignment on its own.

`candidate_inventory_from_stage2(row)` recomputes descending order from the Stage 2 raw crossing records, checks the cached order, checks expected scalar branch labels for `none`, `mixed`, `self_s1`, `self_s2`, and reports missing branches. `CandidateInventory.bind_stages(...)` requires an explicit next-stage mapping for every event, and rejects incomplete inventories. **Do not use these candidate roots as certified thresholds in an actual matching calculation.**

### Local verified results (uploaded continuous_v1 audit, 2026-10-31)

- 16/16 candidate model rows pass Stage 3 event-inventory validation.
- 47/48 heavy-fermion candidate roots found; B(-1) `fermion:0` remains below the saved UV endpoint, not extrapolated.
- 37/38 scalar candidate roots found; D(-2) `scalar:S2_minus` is below the saved EFT1 interval or otherwise not bracketed.
- Total **84** found candidate mass roots; **14/16** models have complete candidate inventories. This does *not* imply 14 physically sequential trajectories are available.
- 0 matching and 0 stage-running callbacks registered in the readiness audit.

## PowerShell usage (from project root)

```powershell
python -m unittest discover -s tests -p "test_general_t3_eft_handoff_controller.py" -v
python -m scripts.AuditGeneralT3EFTHandoffReadiness --all
```

Default input:
`output/studies/general_scalar_extended/continuous_full_scalar/continuous_v1/workstream_d_fermion_audit_v2/audit_summary.json`

Default new output (refuses overwrite):
`output/studies/general_scalar_extended/continuous_full_scalar/continuous_v1/workstream_d_handoff_readiness_v1/summary.json`

If rerunning, choose a **new** `--output` path; existing reports are not overwritten. `--model B_-1` selects one model. `--input` accepts an explicit saved Stage 2 audit.

## Dependencies remaining

1. Workstream B: verifiable RGEs and state projections for EFTs with one or two heavy fermion eigenstates removed, plus later scalar EFTs; Stage 3 does **not** fabricate these.
2. Workstream C: finite physical threshold Wilson matching and an operator-normalization certificate, including B/C Y/X corrections; no matching formula is invented here.
3. Workstream D Stage 4: adapt certified numerical states and coefficient tensors to this envelope, then validate regression/parity against B/C without changing production runners; add final all-16 capability checks.

The historical conditional B(-1)/C(-1) neutrino mass values remain unchanged and are not newly certified by this stage.
