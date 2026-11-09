# Workstream A v18–v20: complex scalar-action source bridge (conditional)

## Purpose and audited source contracts

The EFT1 RGBeta baseline registers 178 renormalisable quartic coupling slots across the 16 T3 models: 150 real and 28 complex. The complex inventory comprises 16 instances of `lambdaT3` and, for B(-1) and C(-1), six special quartics each. The **optional** general-quartic Matchete basis in B(-1)/C(-1) has **seven** extra complex directions per model, of which the baseline RGBeta's six special couplings span a rank-six subspace. The seventh direction is deliberately **not** misreported as a baseline RGBeta parameter.

This batch only reads the existing source and evidence, writes new immutable diagnostic JSON, and does not call `Match[]`, evolve Wilson coefficients, modify GitHub, or alter B/C files. Its complex coupling comparison is conditional on the declared Matchete-aligned SU(2) frame. It does not independently obtain numerical S2 component axes from native RGBeta.

## New files (project-relative paths)

- `scripts/RunT3ParallelAExportComplexMatcheteV18.wl` — native Matchete topology-CG export in `TOPOLOGY` mode (16 models) and optional general-quartic B(-1)/C(-1) CG export in `BCEXTRA` mode (14 CGs). Each export records the current nine source-file SHA256 values and refuses overwrites. B/C mixed-mass matching is **never** attempted.
- `Numerical/parallel_t3/wilson/complex_scalar_bridge_v19.py` — reads the source CGs; verifies all 16 topology symmetry and action normalizations against the independent real-Lam construction; checks that current B/C CG components agree with saved 50-slot native evidence; independently reruns the old conditional tensor fit and extracts the baseline six-of-seven rank/singular-nullspace check; validates v14/v15 installed RGBeta conjugate-channel evidence.
- `scripts/RunT3ParallelAComplexBridgeV19.py` — v19 CLI, new immutable output or independently verified `--resume`.
- `Numerical/parallel_t3/wilson/full_scalar_census_v20.py` — joins source-verified v17 real census, v19 complex bridge and native RGBeta v14 registration evidence; checks the 178 = 150 + 28 census without claiming a physical coupling map.
- `scripts/RunT3ParallelAFullScalarCensusV20.py` — v20 CLI and guarded resume.
- `scripts/RunT3ParallelAComplexScalarBatchV18V20.ps1` — runs regression tests, two native exports, v19 and v20; clears and restores optional environment flags in the diagnostic PowerShell process only.
- `tests/test_t3_parallel_a_complex_scalar_v18_v20.py` — standalone synthetic and negative-control unit tests.
- `Numerical/parallel_t3/wilson/COMPLEX_SCALAR_BRIDGE_V18_V20.md` — this guide.

## Install and execute (Windows PowerShell)

From the repository root:

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
$zip = "$env:USERPROFILE\Downloads\t3_workstream_A_complex_scalar_v18_v20.zip"
Expand-Archive -LiteralPath $zip -DestinationPath . -Force
powershell -ExecutionPolicy Bypass -File scripts\RunT3ParallelAComplexScalarBatchV18V20.ps1
if ($LASTEXITCODE -ne 0) { throw "v18-v20 batch failed" }
```

To revalidate without running Matchete again after the native v18 JSON has already been written:

```powershell
powershell -ExecutionPolicy Bypass -File scripts\RunT3ParallelAComplexScalarBatchV18V20.ps1 -VerifyOnly
if ($LASTEXITCODE -ne 0) { throw "v18-v20 verification failed" }
```

**Outputs** under `output/parallel_t3/wilson/`:

- `matchete_topology_native_v18.json` (16 topology CGs)
- `matchete_bc_extra_native_v18.json` (14 B/C extra CGs)
- `complex_scalar_bridge_conditional_v19.json` (28 baseline complex quartic slots and six-of-seven B/C source-conditioned bridge)
- `full_scalar_census_conditional_v20.json` (all 178 EFT1 quartic registration slots)

Both v18 output modes refuse overwrite. The v19/v20 CLIs refuse overwrite unless `--resume`; resume recomputes or independently checks the evidence and refuses altered results. The PowerShell runner uses `--resume` for both stages.

## Interpretation and limitations

A successful v20 run **certifies registration counts and conditional source-action algebra only**, not a complete physical Matchete–RGBeta coupling dictionary. Native RGBeta numerical triplet component orientation remains unexposed; native complex `Lam` source and conjugate response are observed symbolically but not equal to a full Matchete component action match. B(-1)/C(-1) heavy-scalar mixed-mass diagonalisation is absent; full signed GammaCC-to-Weyl conversion and the physical one-loop ADM remain unresolved. Physical Wilson evolution is prohibited by the written flags.

The mandatory EFT ordering is **UV → integrate out ALL fermions → integrate out ALL scalars → SMEFT**. Never count the scalar threshold before the complete fermionic threshold or add a leading-log LLSS→C5 result to an already evolved C5.

## Validation provenance

The native exporter was designed against current public repository source (`Lagrangian/model/LagrangianBuilder.wl`, `Lagrangian/interactions/GeneralScalarOperators.wl`, `Lagrangian/interactions/T3Topology.wl`) and the upstream RGBeta `Package/FieldsAndCouplings.m` conventions. Standalone synthetic Python unit tests pass in the authoring container. **Actual native Matchete export and all16 live integration require your Windows environment**; no live Wolfram kernel is installed in the authoring container. Do not treat synthetic tests as native results.
