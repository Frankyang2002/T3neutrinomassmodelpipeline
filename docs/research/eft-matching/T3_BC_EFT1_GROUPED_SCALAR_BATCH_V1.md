# B/C grouped scalar EFT1 batch v1

## Scope

One intermediate EFT only: UV (F,S1,S2) -> EFT1 (S1,S2) -> SM+C5. This
package uses the pre-existing historical `continuous_v1` EFT1 coupled Wilson
solution, with a **grouped scalar threshold** at its saved EFT1 endpoint.
It intentionally does NOT create a Y-only/X-only intermediate EFT or run
`RunT3BCEFT2Sequential`.

## Important provenance limitation

The historical `continuous_v1` EFT1 trajectory is **not demonstrated to be
identical** to the later full-five-mass B/C running trajectory audited by
Batch 1. This batch checks Batch 1's status but does not numerically splice
its boundary into historical running. Never describe output as a complete
same-trajectory UV->EFT1->SM calculation.

The existing `RunT3BCCoupledWilsonSequential` implementation includes the
one-generation L/eC seven-term Wilson closure, six physical lepton flavour
pairs, conditional additional-quartic mapping and existing grouped hard
boundary. Its scalar-endpoint C5 projection is a *candidate*; independent
finite scalar threshold corrections are missing. Existing output mass
spectra may include a historical SM continuation, which is NOT certified.

## Run (from repository root)

```powershell
& .\scripts\RunT3BCEFT1GroupedScalarBatch.ps1
```

Prerequisites: Batch 1 summary at
`output/studies/t3_bc_fermion_threshold_batch_v1/summary.json`;
existing `output/.../continuous_v1` trajectory and coupled Wilson outputs
from the historical run. If the coupled outputs do not exist, first ensure
that existing direct and selected Wilson prerequisite reports exist, then
run:

```powershell
& .\scripts\RunT3BCEFT1GroupedScalarBatch.ps1 --run-coupled
```

Use `--resume` to permit updating this batch's summary. Source data are
never overwritten by this package. The historical coupled runner may create
its own outputs only when explicitly requested using `--run-coupled`.

## Outputs

`output/studies/t3_bc_eft1_grouped_scalar_batch_v1/summary.json`

The file reports candidate C5 matrices, threshold scales, candidate
neutrino masses, source parities, solver counts and flags specifying all
unresolved physics gates. `physical_C5_certified=false` always.
