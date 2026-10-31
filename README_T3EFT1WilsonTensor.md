# T3 EFT1 heavy-generation Wilson-tensor bridge (local experiment)

## Purpose and scope

The UV model has **three diagonal heavy fermion masses** `MF[r]`, and all three
heavy fermions are integrated out at **one common fermion threshold**.  The
existing Stage-1 Matchete tree EFT contains, for each external lepton pair,

\[
  W^{12}_{pq}=\sum_{r=1}^{3}\frac{y_{1,pr}y_{2,qr}}{M_{F,r}},
\]

and the additional-Yukawa Wilson coefficient

\[
  W^{\rm extra}_{pq}=\begin{cases}
  \sum_r y_{1{\rm Conj},pr} y_{2,qr}/M_{F,r},&\alpha=0,\\
  \sum_r y_{1,pr} y_{2{\rm Conj},qr}/M_{F,r},&\alpha=-2.
  \end{cases}
\]

The real positive, **nondegenerate** `MF[r]` remain in these *definitions*.
The bars appearing in the native stage-1 EFT are placed on the **entire** W
and correspond to complex conjugation of the product for real masses. The
mathematical definition is not a numerical identification with RGBeta.

The first observed failure of `RunThresholdStage.wl` was the triply occurring
`NFlavor` index on `y1[p,r] y2[q,r] / MF[r]`.  We replace only the **tree-level
monomials** matching this precise pattern by native `Coupling[T3EFT1W12,
{Flavor,Flavor},0]` and `Coupling[T3EFT1WExtra,{Flavor,Flavor},0]`. The
scalar threshold is then fed this two-flavour EFT1 tensor; it does not see the
original F-generation index. The first threshold itself remains the existing
real Matchete single-step matching with three nondegenerate heavy generations.

## Project paths in this ZIP

- `Lagrangian/RunT3EFT1WilsonHierarchyNative.wl` — new native F-to-EFT1
  tensor driver and optional scalar matching.
- `Lagrangian/RunT3EFT1WilsonStageLauncher.wl` — an isolated copy of the
  source-verified `RunFreshThresholdStage` function with a different child
  filename; no modification to `RunMatching.wl`.
- `Lagrangian/RunThresholdStageEFT1Wilson.wl` — full copy of the local
  `RunThresholdStage.wl` with only the two Wilson couplings registered in the
  fresh kernel; no modification to the production file.
- `scripts/RunT3EFT1WilsonTensorBatch.py` — six-model seeded/zero dispatcher,
  source/provenance guards, and native reports.
- `scripts/TestT3EFT1WilsonTensor.py` — regression tests, including the
  optional previous A(0) WXF witness.

**The two production source files were verified using the exact copies from
`T3_Ap0_stage2_failure.zip`** (SHA-256 digests checked by preflight). If they
have since changed, the runner blocks rather than silently using an outdated
fresh-kernel copy. No GitHub changes are made.

## PowerShell (repo root)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
$zip = "$HOME\Downloads\T3EFT1WilsonTensor_Local.zip"
Expand-Archive -LiteralPath $zip -DestinationPath . -Force

python scripts\RunT3EFT1WilsonTensorBatch.py --self-test
python scripts\TestT3EFT1WilsonTensor.py
python scripts\RunT3EFT1WilsonTensorBatch.py --preflight

# Isolate the exact 4 original A(0) tree monomials first.
python scripts\RunT3EFT1WilsonTensorBatch.py --pilot --mode prepare

# Next: attempt scalar tree matching of contracted EFT1.
if ($LASTEXITCODE -eq 0) {
    python scripts\RunT3EFT1WilsonTensorBatch.py --pilot --mode tree
}

# Next: scalar one-loop [B], still WITHOUT inherited hard [C] or EFT1 RGEs.
if ($LASTEXITCODE -eq 0) {
    python scripts\RunT3EFT1WilsonTensorBatch.py --pilot --mode soft
}
$LASTEXITCODE
```

**Only after the A(0) scalar pilots pass**, run:

```powershell
python scripts\RunT3EFT1WilsonTensorBatch.py --all --mode prepare
python scripts\RunT3EFT1WilsonTensorBatch.py --all --mode soft
$LASTEXITCODE
```

Outputs in:

`Reports/output/workstream_a/eft1_wilson_tensor/<mode>/<model>/<variant>/`

Each model and variant writes `native.json`, `stdout.log`,
`stage1_tree_original.wxf`, `stage1_tree_wilson.wxf`,
`stage1_hard_one_loop_omitted_from_stage2.wxf`,
`wilson_definitions.json`, and, for successful scalar modes,
`stage2_weinberg_unrenormalized.wxf`. The existing child stage also writes
`debug/fresh_kernel_stage_2/stdout.log`. Combined status is at
`Reports/output/workstream_a/eft1_wilson_tensor/<mode>/summary.json`.

## Physical limitations — do not promote status

- `prepare` checks Wilson contraction, not scalar matching.
- `tree` checks scalar tree matching, not one-loop C5.
- `soft` checks only the **new one-loop scalar matching** of the **tree F-threshold
  EFT1 Wilson operators**. It is one ingredient of the two-threshold result.
- **The stage-1 hard one-loop correction is saved but NOT propagated into
  stage 2.** It contains additional complicated, multiply indexed terms and
  requires separate operator-by-operator Wilson projection.
- EFT1 Wilson RGE evolution and the UV pole/RGE cancellation are absent.
- The fully renormalized hierarchical C5 and the absolute Matchete–RGBeta
  physical coupling dictionary remain **uncertified**.
- This calculation assumes diagonal, independent heavy masses already in the
  heavy mass basis. It does not justify a common mass or perform a new Dirac
  mass-basis rotation.

The first useful live result is whether the scalar tree threshold accepts the
contracted EFT1 tensor. A successful `soft` run would then supply only a
**partial, unrenormalized scalar contribution** to hierarchical C5.
