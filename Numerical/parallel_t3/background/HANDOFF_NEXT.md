# Workstream B -> A/C, saved numerical handoff gate (v1)

**Scope:** Checks completed Workstream B numerical data without rerunning RGBeta,
Matchete, Mathematica, any renormalisable ODE or any EFT Wilson ODE. Does NOT
modify GitHub or files owned by A/C. The only persistent output is
`output/parallel_t3/background/<run_id>/handoff_index.json`.

## Install paths (new files only)

- `Numerical/parallel_t3/background/handoff.py`
- `Numerical/parallel_t3/background/HANDOFF_NEXT.md`
- `scripts/RunT3ParallelBHandoff.py`
- `tests/test_t3_parallel_b_handoff.py`

Existing Workstream B modules from the successful 16/16 batch are reused and
must remain installed. The new ZIP is additive and contains complete files.

## Exact PowerShell commands

From `C:\thesistemp\T3neutrinomassmodelpipeline`:

```powershell
# Extract the ZIP into this repository root first (replace <downloaded-zip>):
Expand-Archive -LiteralPath "<downloaded-zip>" -DestinationPath . -Force

# All native-dependency-free handoff tests (test fixtures are clearly synthetic):
.\.venv\Scripts\python.exe -m pytest -q tests\test_t3_parallel_b_handoff.py

# Read-only validation of the actual saved A(-4) first:
.\.venv\Scripts\python.exe -X utf8 scripts\RunT3ParallelBHandoff.py `
  --model A_-4 --run-id B_all16_qualitative_v1 --no-write

# Full saved all-16 handoff audit; creates ONE NEW B-owned JSON index:
.\.venv\Scripts\python.exe -X utf8 scripts\RunT3ParallelBHandoff.py `
  --all --run-id B_all16_qualitative_v1
```

**Expected from previous user terminal evidence:** 16 successful input
backgrounds and 16 formal source boundaries exist, but the actual handoff gate
result is not known until executed locally. It checks stronger cross-file
invariants than the original B runner's final status summary.

## Numerical contract

For each model (`A_m4`, `B_m1`, etc), consumes without altering:

```
output/parallel_t3/background/B_all16_qualitative_v1/<tag>/
    background_manifest.json
    background.npz
    fermion_boundary_manifest.json
    fermion_boundary_formal.npz
```

The `background.npz` contains descending positive scales and all named
`complex128` UV/EFT1 beta coordinates with leading sample dimension. In
particular, every common renormalisable parameter must be continuous across
`MF` under B's **tree-level renormalisable matching assumption**. The
exact inventory satisfies `UV - {MF,y1,y2,y1Conj,y2Conj} = EFT1`, and
`EFT1` includes all scalar quadratic and quartic keys, SM Yukawa matrices
and gauge couplings. The manifest and source hashes must match the saved
numerical input contract.

The source Matchete `C_LLSS_formal` is
`complex128[3,3,2,2,n_real,n_real]` in the order
`(flavour_p,flavour_q,lepton_SU2_i,lepton_SU2_j,scalar_a,scalar_b)`, with
`n_real=2*(dS1+dS2)=8` (A/B/C/D) or `12` (E).
`C_LLSS_by_source_formal` has one extra leading Matchete source-term axis;
its sum must reconstruct the total. The scalar pair and combined lepton
field pair exchange symmetries are verified. Heavy-mass eigenvalues remain
three strictly positive finite numbers at the **single common** `MF`.

The gate verifies the same-run **UV endpoint SHA** across both files,
source Matchete seed SHA, per-artifact SHA256 and all scale/axis labels.
Post-run *current* B source changes are recorded as warnings, without
pretending they retroactively changed immutable saved numerical outputs.
A modified or corrupted saved NPZ/seed itself fails closed.

## Generated handoff index

`output/parallel_t3/background/B_all16_qualitative_v1/handoff_index.json`

- `schema`: `t3-parallel-B-handoff-index-v1`
- `passed_count`, `requested_count`: number of **provisional B data handoffs** passing
- `model_results[model].status`: `Verified_ProvisionalBInputsForAandC` or a concrete blocker
- `model_results[model].files`: four input paths and SHA256s
- `model_results[model].beta_coordinate_inventory`: exact UV/EFT1 parameter names
- `model_results[model].formal_fermion_boundary`: tensor shape, norm, labels, operators, heavy masses
- `model_results[model].workstream_A`: native background/formal source available; canonical Wilson evolution **not** done
- `model_results[model].workstream_C`: scalar coordinates available; independent component tensors **not** certified
- `source_drift_warnings`: post-run B source path differences, if any

A partially blocked index is still written, with per-model blocker strings.
The command returns nonzero when any selected model is blocked. It does not
regenerate or overwrite any saved physics.

## Read-only consumer example

```python
from pathlib import Path
from Numerical.parallel_t3.background.api import Background, load_formal_boundary
from Numerical.parallel_t3.background.contracts import read_json

root=Path('output/parallel_t3/background/B_all16_qualitative_v1')
index=read_json(root/'handoff_index.json')
assert index['passed_count']==16
folder=root/'A_m4'
bg=Background.load(folder)
f_arrays,f_meta=load_formal_boundary(folder)
native_eft1=bg.native_eft1_at(1e4) # native RGBeta MSbar at mu=1e4 GeV
native_scalar_ms=bg.native_scalar_ms() # scalar coords only, no physical tensor
formal_llss=f_arrays['C_LLSS_formal'] # Matchete commuting-spinor source convention
```

### Explicit unresolved physics ownership

1. **Workstream A** must derive the signed GammaCC-to-two-component-Weyl
   conversion, the native RGBeta/Matchete quartic/Yukawa basis map needed for
   the Wilson ADM, and the full coupled dimension-five system with all
   generated operators required for closure, including `C5`. This gate does
   not assign a numerical signed spinor factor or evolve Wilson coefficients.
2. **Workstream C** must construct the complete scalar component quadratic
   and quartic tensors/eigenbases at `MS`, including mixed B(-1)/C(-1) and
   six neutral self-bilinear cases; integrate **both scalars simultaneously**.
3. An absent hard Weinberg coefficient at `MF` is **unknown**, never a
   physical `C5(MF)=0`. The gate permits a nonzero *independently documented*
   `C5_hard_external` candidate while marking full one-loop hard matching
   incomplete. Any missing other hard dim-five terms remain a documented gap.
4. The output is a **candidate native renormalisable and formal tree F-boundary
   handoff**, not a physical neutrino mass prediction. No scalar matching or
   low-energy Weinberg result is produced here.

A successful 16/16 handoff means you can feed the existing saved trajectories
and formal tensors to A/C without recomputing B, **subject to the unresolved
normalisation, ADM and scalar-basis issues above**.
