# Workstream C — v20: source-level one-loop LLHH topology completeness

## Scope

The source-derived B/C(-1) field content has S1,Y=-1/2; S2,Y=+1/2; F,Y=0.
The new optional conjugate-scalar Yukawa family is absent for these models,
because neither S1 nor S2 has zero hypercharge. The original lepton Yukawas
source S1 (`y1`) and S2dag (`y2`). A neutral heavy Majorana F connects them.

At leading O(y_i y_j lambda), in the unbroken EW phase and with diagonal free
scalar propagators and zero additional quadratic mass insertions, an LLHH
1PI one-loop graph needs exactly two original Yukawa vertices and one
four-scalar HHSS vertex. The resulting graph has I=3, V=3, L=I-V+1=1.
Both Hermitian orientations of the scalar quartic must be tested.

For a renormalizable one-loop four-point graph, the vertex valence identity is
`sum_v (valence(v)-2) = 4`. Two required BSM Yukawa cubics contribute 2,
leaving either **one quartic** or **two additional trivalent vertices**.
The latter class, potentially involving SM Yukawa/gauge vertices, is explicitly
**not fully classified** by this batch. Thus three diagrams prove completeness
only within the single-quartic branch.

| Scalar pair at quartic end | Yukawa source pair | UV quartic |
|---|---|---|
| S1dag S1dag | y1 y1 | lambdaT3ExtraQ020200I1 |
| S1dag S2 | y1 y2 | lambdaT3 (h.c. orientation) |
| S2 S2 | y2 y2 | lambdaT3ExtraQ020020I1 |

The B/C operator manifest contains 7 additional quartics and 1 allowed
additional bilinear per model. The bilinear `mu2_extra_000101_I1` has species
`S1dag S2dag + h.c.` and can affect tree-level propagators. A nonzero value
requires re-diagonalization, not loop-counting suppression. The zero-mixed-mass
boundary is an assumption of the earlier scalar relevance benchmark, not a
general constraint imposed by the UV symmetries.

## What the batch does

- Imports and **reuses**, without modifying, the project's prior exhaustive
  16-model bosonic operator inventory / source hash gates and Wick audit.
- Verifies current UV topology, hypercharge and operator enumeration contracts.
- Independently enumerates all direct one-quartic triangles under the stated
  assumptions, rejecting missing/extra channels, wrong orientations and
  unauthorized source or operator-family changes.
- Scans the **entire saved Matchete one-loop EFT correction** at each of the
  F, Y and X stages for LLHH GammaCC candidates, then compares quartic and
  Yukawa-channel coverage with each saved raw sector. Counts must be 204/132/36
  in each model.
- Requires v19's restricted MS-bar hard-soft subtraction and original Matchete
  evidence to be present, and records file hashes in one JSON.

No existing source/matching/RGE file is edited; the only new output is a JSON
under `output/audits/t3_bc_one_loop_topology_completeness_v20/`.

## Explicit limits

This is **exhaustiveness within a vertex truncation**, not a complete
renormalized all-coupling UV+SM diagram calculation. In particular, possible
SM gauge/fermion insertions, nonzero off-diagonal scalar mass insertions,
operator-basis transformations, gauge/EOM effects, and consistent EFT running
are not certified by counting Wick contractions. Source selection using
`GammaCC` is not an independent Fierz proof. `physical_C5` is always null.

The older SU(2) component equivalence test establishes equivalence of the
renormalizable B/C quartics at the component level; it does *not* make the
Matchete symbolic EOMEFT algebraic equality pass.

## PowerShell

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$HOME\Downloads\T3_WorkstreamC_OneLoopTopologyCompleteness_v20.zip" -DestinationPath . -Force
& .\scripts\RunT3BCOneLoopTopologies.ps1
```

Or, to run the audit alone:

```powershell
python scripts\ClassifyT3BCOneLoopTopologies.experimental.py `
  --output output\audits\t3_bc_one_loop_topology_completeness_v20\summary.json
```

If the script fails, its error identifies the first unverified operator family,
source signature or Matchete sector. Do not work around its provenance gates.
