# T3 extra-Yukawa hierarchical Matchete matching pilot

## Scope

This workstream takes the six previously verified full one-step native Matchete
coefficients and attempts a **genuine two-threshold symbolic matching** in the
existing Matchete implementation:

```text
native UV Lagrangian with F heavy, S1/S2 light
  -> integrate out F (EFTOrder=5; tree + one loop)
  -> intermediate EFT1 with both scalars dynamical
  -> reclassify S1/S2 as heavy in a fresh Wolfram kernel
  -> integrate out S1/S2 (EFTOrder=5; tree + one loop)
  -> attempt native Weinberg coefficient extraction
```

No production files, RGBeta beta functions or existing numerical results are
changed. The six models are `A_p0`, `A_m2`, `D_p0`, `D_m2`, `E_p0`, `E_m2`.
The full native representation uses **three heavy flavours**, matching the
existing Matchete model-builder and successful one-step calculations.

For every model/mode the runner launches two **independent** Matchete kernels:
`seeded` uses the additional Yukawa as defined by the native UV action;
`zero` sets the `y1Conj`/`y2Conj` coupling to zero in the **UV Lagrangian
before matching**. A third kernel imports the two output EFTs and saves
raw symbolic differences. This contrasts native matching actions: it is **not**
the same control as an initially-zero extra Yukawa under the full coupled RGEs.

### What this does not do

* No EFT1 Wilson-coefficient RGE evolution between `mu_F` and `mu_S`.
* No renormalized, pole-subtracted final hierarchical C5, even if the
  native second-threshold extractor reports a symbolic coefficient.
* No Matchete-to-RGBeta absolute coupling/phase dictionary.
* No universal proof of dummy-index equivalence between independently
  matched Wolfram kernels. Exported differences are **formal**; a failed
  `PossibleZeroQ` is NOT a proof that a Wilson operator is physically nonzero.
* No full mixed scalar quadratic basis; the already-established `GENERAL_QUADRATICS=0`
  exclusion remains in force.
* No low-energy neutrino mass prediction from the hierarchical branch.

The project already has a successful conditional sudden-decoupling numerical
calculation; **this package does not overwrite or silently promote that result**.

## Project paths (all complete new files)

| Path | Role |
|---|---|
| `Lagrangian/RunT3ExtendedHierarchyNative.wl` | Matchete F-only / two-threshold source matching, per model and variant |
| `Lagrangian/CompareT3ExtendedHierarchy.wl` | Cross-run WXF expressions and formal UV-seeded-minus-zero comparisons |
| `scripts/RunT3ExtendedHierarchyBatch.py` | Full six-model orchestration, evidence gate, logging and reporting |
| `scripts/TestT3ExtendedHierarchySynthetic.py` | Synthetic orchestration, source-contract and negative controls |
| `README_T3HierarchicalExtraYukawa.md` | This file |

All outputs go under:

`Reports/output/workstream_a/extended_c5_hierarchy/`

## PowerShell (repo root)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline

$zip = "$HOME\Downloads\T3HierarchicalExtraYukawa_Local.zip"
Expand-Archive -LiteralPath $zip -DestinationPath . -Force

python scripts\RunT3ExtendedHierarchyBatch.py --self-test
python scripts\TestT3ExtendedHierarchySynthetic.py
python scripts\RunT3ExtendedHierarchyBatch.py --preflight

# First substantial live stage: integrate out F, with seeded/zero action controls.
python scripts\RunT3ExtendedHierarchyBatch.py --pilot --stage1

# Only if the F-only native pilot succeeded, test the second threshold.
if ($LASTEXITCODE -eq 0) {
    python scripts\RunT3ExtendedHierarchyBatch.py --pilot --hierarchy
}

# When the pilot succeeds, run the entire 6-model family.
# These are two separate large Wolfram batches; choose one or both.
python scripts\RunT3ExtendedHierarchyBatch.py --all --stage1
python scripts\RunT3ExtendedHierarchyBatch.py --all --hierarchy

$LASTEXITCODE
```

If there is a blocker, look at the actual log:

```
Reports/output/workstream_a/extended_c5_hierarchy/stage1/A_p0/seeded/stdout.log
Reports/output/workstream_a/extended_c5_hierarchy/hierarchy/A_p0/seeded/stdout.log
Reports/output/workstream_a/extended_c5_hierarchy/hierarchy/A_p0/seeded/debug/fresh_kernel_stage_2/stdout.log
```

### Main reports

* `stage1/summary.json`: six-model F-only EFT1 report, seeded/zero comparison.
* `hierarchy/summary.json`: attempted full two-threshold matching report.
* `<mode>/<model>/seeded/native.json`: native result and flags.
* `<mode>/<model>/zero/native.json`: UV-zeroed native result and flags.
* `<mode>/<model>/difference.json`: formal WXF comparisons.
* `<mode>/<model>/stage1_tree_seeded_minus_zero.wxf`: matched tree EFT1 difference.
* `<mode>/<model>/stage1_one_loop_seeded_minus_zero.wxf`: matched one-loop EFT1 difference.
* `<mode>/<model>/stage2_c5_difference.inputform.txt`: ordered final Weinberg
  coefficient difference **only if** the full two-threshold extraction completes.
* `source_gate.json`: SHA-256 source/provenance gate.

Any `NativeHierarchyBlocked` state is a genuine failure or unsupported path,
not a zero coefficient. The Python script returns nonzero when any requested
model or comparator is blocked; completed stage-one artifacts are retained
when a later threshold fails.

### Physical meaning and follow-up

This batch tests whether the already implemented general-fermion Matchete UV
Lagrangian can generate the extra-dependent intermediate EFT and propagate it
through the **existing** sequential matching engine. It preserves the original
project's fixed-order tree/new-one-loop/inherited-one-loop distinction.

The next theoretical/numerical step, after obtaining actual EFT1 operators,
is to derive **their RGEs and source-normalized threshold relations** in the
same Matchete field convention, demonstrate pole/log cancellation, and
calibrate the Matchete↔RGBeta coupling dictionary before substituting the
running 3×3 Yukawa matrices. A successful symbolic matching trace alone does
not accomplish that.

## Current validation

Locally verified in the package build environment: Python compilation,
balanced Wolfram source structure, package load ordering, source contracts,
12 synthetic model/mode runs (24 mocked native match calls and 12 mocked
comparisons), and certification/tamper negative controls. **No installed
Matchete/Wolfram kernel was available** in the package build environment;
actual symbolic matching is to be checked on the user's Windows installation.
