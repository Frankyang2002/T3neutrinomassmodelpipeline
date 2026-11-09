# Workstream A v25–v27: B/C native mass-eigenstate concordance and EFT-order gate

This batch is **read-only**. It consumes existing Windows native Matchete checks and
prior Workstream A v21/v22 scalar mass transport. It does not call `Match[]`,
reintegrate fermions, run Wilson coefficients, or modify the original Lagrangian.

## Source contracts inspected (GitHub read-only, 2026-11-10)

- `Numerical/experimental/GeneralT3MixedScalarMassBasis.py`: two-doublet
  convention `Phi=(S1,epsilon conjugate(S2))=U(X,Y)`; complex mass matrix
  and full epsilon index-level transport.
- `scripts/RunGeneralT3MassEigenMatchete.py` and
  `Lagrangian/RunGeneralT3MassEigenMatchete.wl`: an *existing* native Matchete
  check with fresh `NewScalarX`, `NewScalarY` fields and diagonal masses.
  Expected status: `MatcheteMassEigenHeavyBasisPassedPhysicalMatchingNotPerformed`.
- `scripts/RunGeneralT3MassEigenOneLoop.py` and
  `Lagrangian/RunGeneralT3MassEigenOneLoop.wl`: previously executed one-step
  *all-heavy* Weinberg matching probe (dimension-five, one-loop). This is
  **not** the physical fermion-first EFT1 threshold order.
- `Numerical/parallel_t3/wilson/bc_mass_eigen_native_bridge_v21.py`,
  `Numerical/parallel_t3/wilson/bc_mass_bridge_gate_v22.py`: the most recent
  native CG and mass rotation audit, executed successfully on Windows.

### v25: native eigenfield ↔ v21/v22 source agreement

Recomputes v22 from its saved source hash chain and compares against the saved
result. For B_-1 and C_-1 independently, validates:

1. SHA256 of the original mixed mass-basis JSON against the Matchete parameter
   file, accepting Windows or POSIX provenance separators but no unsafe paths;
2. identical source and native `A`, `U` and positive ordered squared masses,
   `U†U=I` and `U†AU=diag(mX²,mY²)`;
3. Matchete's **actual previously executed** HeavyMassBasis, 18 interaction
   inventory, exact `NewScalarX/Y` heavy field registry, and 7+1 complex CG
   source identifiers;
4. native epsilon tensor `{ {0,1},{-1,0} }`, independent portal and self
   CG ranks, seven additional quartics and both native Yukawa phase factors;
5. previously saved source numerical unitary matrix is the same exact `U`
   inserted into the native Matchete spurions (`T3MEU11`...`T3MEU22`).

All native **physical matching** and RGBeta↔Matchete physical-coupling-map
flags must remain `false`.

### v26: previous Matchete one-step probe, not SMEFT matching

Cross-compares the saved `one_loop_probe` against the separately executed
native HeavyMassBasis certificate: all 18 registered interactions, 7 extra CGs
and the U spurions must agree. Validates the **118 ordered holomorphic and
118 conjugate** Weinberg terms per model and records exact SHA256 and sizes of
both large Mathematica `InputForm` artifacts. The ordered `A5` is *not*
silently reinterpreted as canonical symmetric `C5`.

These are prior **simultaneous all-heavy** matching diagnostics. The v26
report explicitly rejects claims of `physical_matching_performed`,
`production_sequential_matching_performed`, or `numerical_c5_computed`.

### v27: threshold-scale order, with no synthetic running

The saved B/C candidate EFT1 states are at 1000 GeV, **below** both scalar
mass eigenvalues. For each model, the output records
`ln(mY/mX)`, `ln(mX/1000 GeV)`, `ln(mY/1000 GeV)`, and the explicit missing
fermion-first EFT1 hard boundary. A 1000-GeV EFT1 endpoint cannot simply be
used as a scalar matching boundary at either eigenmass. A future consistent
run must first compute an ALL-fermion EFT1 boundary, evolve **to the chosen
scalar matching scale** with a verified one-loop ADM, integrate out all scalar
eigenstates, then run SMEFT to the low scale.

The EFT order remains:

> UV → integrate out **ALL fermions** → integrate out **ALL scalars** → SMEFT.

The previous one-step all-heavy Matchete probe remains comparative evidence
only and cannot be inserted as a scalar-threshold Wilson coefficient.

## Files and commands

The delivered ZIP contains **7 new files** (project-relative paths):

- `Numerical/parallel_t3/wilson/bc_native_mass_eigen_reconcile_v25.py`
- `Numerical/parallel_t3/wilson/bc_one_scale_probe_audit_v26.py`
- `Numerical/parallel_t3/wilson/bc_threshold_order_ledger_v27.py`
- `Numerical/parallel_t3/wilson/BC_NATIVE_MATCHING_CONCORDANCE_V25_V27.md`
- `scripts/RunT3ParallelABCMatcheteThresholdBridgeV25V27.py`
- `scripts/RunT3ParallelABCMatcheteThresholdBridgeV25V27.ps1`
- `tests/test_t3_parallel_a_bc_native_matchete_v25_v27.py`

Install the ZIP by expanding it at project root, then:

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
powershell -ExecutionPolicy Bypass -File scripts\RunT3ParallelABCMatcheteThresholdBridgeV25V27.ps1
```

Outputs written **once**, with source-checked replay on reruns:

- `output/parallel_t3/wilson/bc_native_matchete_mass_eigen_concordance_v25.json`
- `output/parallel_t3/wilson/bc_existing_one_scale_matchete_probe_v26.json`
- `output/parallel_t3/wilson/bc_fermion_first_scalar_threshold_ledger_v27.json`

Use the PowerShell `-VerifyOnly` flag to require and independently revalidate
all three. Refuse overwriting altered evidence; delete *only the new result*
when you intentionally want to recompute after changing any source.

## Verification classification

- **Inspected source contracts:** live GitHub repository (read-only).
- **Executed in this package:** synthetic numerical tests of audit entrypoints,
  strict schema and source-hash guards; no Wolfram available in the container.
- **Pending local execution:** your saved native Matchete and v21/v22 evidence
  integration. The test suite includes an auto-detected live integration test.
- **Not certified:** physical Matchete/RGBeta parameter mapping, signed
  GammaCC→Weyl, fermion hard boundary, physical one-loop EFT1 ADM, scalar hard
  threshold, or SMEFT neutrino mass. No physics result is produced.
