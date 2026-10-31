# B/C batch 3: same-UV SM-only reference and historical EFT1 comparison

## Correct EFT chain

`UV(F,S1,S2) -> EFT1(S1,S2) -> SM + C5`. There is **one** intermediate EFT.
Both scalar eigenstates are integrated at a grouped threshold for this restricted study.

## What is calculated

- Revalidates the *new five-mass* UV state and archived 480-term one-scale source parity.
- Runs the **actual installed SM+Weinberg numerical RGE** twice with `rtol=1e-8,1e-9`, using the **same 1e5 GeV five-mass UV reference**.
- Uses `continuous_v1` historical EFT1 coupled Wilson report only for **side-by-side comparison**. Does not insert its coefficients into the new UV history.
- Emits comparable numerical spectra and a solver convergence residual with SHA256 source provenance.

## What is NOT calculated

The *correct* same-five-mass sequential EFT1 Wilson integration requires a UV->EFT1 hard LLSS operator boundary and the actual five-mass EFT1 scalar history. This batch does not construct that new Wilson history, compute a renormalized finite scalar threshold, or certify the inter-package coupling map. SM-only continuation from 1e5 to 100 GeV neglects dynamical scalars over much of the interval; it is a **proxy**, not an EFT calculation.

**Stop gate:** If same-trajectory EFT1 running is absent, the report states `same_trajectory_EFT1_wilson_running_completed=false`. Do not use the historical-versus-new mass difference as a physical RG correction. Physical C5 is always false.

## Output

`output/studies/t3_bc_final_benchmark_consolidation_v1/summary.json`

To rerun with new output ID:

`python -u -X utf8 -m scripts.RunT3BCFinalBenchmarkConsolidation --output-id t3_bc_final_benchmark_consolidation_v2`

No files in `RGE/`, `Lagrangian/`, or production matching modules are changed.
