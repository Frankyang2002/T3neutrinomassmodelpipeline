# General T3: all-16 leading-one-loop LLSS -> Weinberg sequential bridge

## Deliverable and scope

Install the **full files** listed below at their identical relative locations in
`T3neutrinomassmodelpipeline`. No existing repository physics file is replaced.
The implementation was constructed after inspecting the current GitHub `main`
branch. It is a **provisional numerical implementation**, not a certified
complete sequential physical matching calculation. It uses the most-general
symmetry-allowed renormalisable parameter inventories available in the current
all-16 numerical code; it does **not** silently switch to restricted C5 kernels.

The physical integration/decoupling sequence is always:

    UV (1e7 GeV) -> integrate out ALL F (1e5 GeV)
    -> SM + S1 + S2 + dimension-five LLSS operators
    -> integrate out BOTH scalar multiplets (1e3 GeV)
    -> SMEFT SM + C5 -> low scale (1e2 GeV)

The computed intermediate mixing is the **leading-order scalar-pair
LLSS -> Weinberg source**, at fixed overall one-loop order:

    delta_C5_log = log(MS/MF) / (16*pi^2) * beta_C5_16pi2.

At this order, one-loop running of the tree LLSS coefficient itself inside
another one-loop scalar bubble would be two loop. The general seven-term
master Wilson anomalous dimension and conditional implementations are already
present in the repository (`RGE/general/WilsonTensorRGE.py` and
`RGE/running/intermediate/GeneralEFT1C5Mixing.py`). Our bridge reuses the
repo's scalar-pair projector and finite joint-scalar bubble, without asserting
that the other intermediate operator classes are exhaustively matched.

## Implementation batches (all included in this ZIP)

**Batch 1 — General fermion-tree Wilson boundary.**
`Numerical/running/GeneralT3DynamicFermionWilson.py`: for the 14 nonmixed
models, reconstruct tree LLSS from the actual evolved 3-by-3 heavy fermion
Yukawas and each original Matchete term/CG tensor, preserving allowed
conjugate Yukawa source channels. For B(-1)/C(-1), reweight the source-resolved
native mixed-mass reference tensor with the actual symmetrized heavy-flavour
bilinear; record the reference reconstruction residual and explicitly label
this an approximation, not exact native mixed-flavour matching.

**Batch 2 — All general HHSS tensors.**
`Numerical/running/GeneralT3SequentialLeadingLog.py`: construct source-derived
native HHSS Matchete channels, the **full RGBeta baseline portal singlet and
adjoint tensors**, and the exact SU(2)-invariant reference real-component
basis for **all registry-listed additional HHSS couplings**. Reconstruct
all channels at the current running/frozen EFT scalar endpoint and rotate to
mixed/self-bilinear scalar bases where needed. Do not invent zero couplings.

One convention is necessarily provisional: within each extra HHSS field-count
sector, the j-th RGBeta symbol is paired with the j-th independent exact SU(2)
reference invariant; this is not the derived native Matchete-to-RGBeta basis
change. Similarly, the common source GammaCC -> Weyl factor defaults to 1,
and can be varied with `--spinor-factor`. Neither is certified physics.

**Batch 3 — Independent intermediate LLSS leading log and both-scalar hard
bubble.** The tensor kernel reuses `Numerical/joint_scalar_matching/large_batch.py`
`_response`, which contains (i) the fixed-order logarithmic mixing term and
(ii) the simultaneous two-scalar bubble at the same order. Output includes
complex 3-by-3 beta, log, bubble, sum, per-general-HHSS-coupling contributions,
per-fermion-tree-source contributions, Weinberg projection residuals and a
finite matching-scale cancellation check.

**Batch 4 — Running/freeze ablations.**
Two runners are included:

1. `scripts/RunT3GeneralSequentialLeadingLog.py`: all-16 independent
   reference-source 2^3 comparison (EFT1 parameters, intermediate Wilson
   log, SMEFT), no fitted UV source substitution.
2. `scripts/RunT3GeneralSequentialFitAblation.py`: all-16 fitted 2^4
   comparison (UV, EFT1, Wilson log, SMEFT), for a local previously computed
   `RunT3GeneralAll16FromReference.py` campaign. No refitting. Running/frozen
   states include all symmetry-allowed general couplings carried by their
   respective numerical runner.

### Avoiding double counting in the fitted comparison

The existing **full-general one-scale Matchete C5** remains the baseline when
Wilson log is ON. With matched inputs at each U/E configuration, define:

    C5_residual = C5_one_scale - C5_LLSS_log - C5_LLSS_bubble
    C5(W=1)    = C5_residual + C5_LLSS_bubble + C5_LLSS_log
               = C5_one_scale
    C5(W=0)    = C5_residual + C5_LLSS_bubble
               = C5_one_scale - C5_LLSS_log

This is an **anchored diagnostic decomposition**. `C5_residual` is not
independently computed fermion hard matching. The one-scale expression may
contain other logarithms, and independently deriving the hard sequential
matching remains outside this batch. Do **not** add our computed log to the
full one-scale C5 without subtracting it first.

## Installed paths

    Numerical/running/GeneralT3SequentialLeadingLog.py
    Numerical/running/GeneralT3DynamicFermionWilson.py
    scripts/RunT3GeneralSequentialLeadingLog.py
    scripts/RunT3GeneralSequentialFitAblation.py
    tests/python/test_general_sequential_leading_log.py
    README_GENERAL_SEQUENTIAL_EFT1.md

### Local PowerShell installation

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$HOME\Downloads\T3_General_Sequential_EFT1_LeadingLog_BigBatch.zip" -DestinationPath . -Force
python -X utf8 -m pytest -q tests/python/test_general_sequential_leading_log.py
```

### Prerequisite check (all 16, without running a numerical scan)

```powershell
python -X utf8 -m scripts.RunT3GeneralSequentialFitAblation `
  --source-run thesis16_v1 --run-id general_seq_01 --check-only
```

Replace `thesis16_v1` by your actual saved general fitter run ID.
If inputs are missing, the check prints their exact local paths. The relevant
symbolic and numerical files are generated by the **existing** workstreams;
this batch does not attempt a new Wolfram symbolic release, fabricate source
artifacts, or change GitHub.

### All-16 fitted 16-switch study

```powershell
python -X utf8 -m scripts.RunT3GeneralSequentialFitAblation `
  --source-run thesis16_v1 --run-id general_seq_01 --spinor-factor 1
```

For a model subset, add `--models A_-4,A_-2,B_-1`.
Each four-character case is `U E W S` with 1 for running enabled and 0 for
frozen, for example `0000`, `1100`, `1110`, `1111`.

Outputs:

    output/studies/general_scalar_extended/general_sequential_rge_ablation/general_seq_01/
      summary.json
      all16_cases.csv
      <model-tag>/cases.json   (e.g. A_m4/cases.json)
      <model-tag>/error.json   (if a model is blocked)

The CSV includes the neutrino mass singular values, mass-squared differences,
mixing observables, fit chi-square and `||C5_LLSS_log||/||C5_one_scale||`.
Each case JSON additionally contains complete C5, log, bubble, source and
quartic-channel contributions, with model provenance and all assumptions.

### Standalone source-backed LLSS study (does not require fitted campaign)

```powershell
python -X utf8 -m scripts.RunT3GeneralSequentialLeadingLog --check-only
python -X utf8 -m scripts.RunT3GeneralSequentialLeadingLog `
  --with-smeft `
  --output-dir output/general_observable/lo_eft1_mixing/general_sequential_llss/study_01
```

Outputs include `summary.json`, `comparison.csv` and per-model `cases.json`.
`--skip-general-extra-hhss` and `--skip-portal-completion` are optional
**diagnostic truncations only**; defaults keep the general interactions.

### Required existing local inputs

- A 16-model source-resolved Matchete tree release, default
  `output/general_observable/lo_fermion_boundary/releases/source_021dcc1e8fe8ba0f`.
- Native HHSS v4 cached CG channels and corresponding original quartic seeds
  (`output/general_observable/lo_joint_scalar_matching/source_quartic_v4/native_quartic_channels`
  and `output/full16_symbolic`).
- General scalar registry at
  `RGE/general/generated/general_rgbeta_quartic_registry.json`.
- Saved 16-model all-general fits under
  `output/studies/general_scalar_extended/thesis_general_all16/<source-run>`
  and their canonical general beta exports / Matchete coefficient expressions.
- For fitted B/C models, the original
  `output/general_observable/all16_v1` running snapshot and the existing B/C
  general RGBeta and one-scale Matchete source.
- Mathematica/WolframScript for on-demand evaluation of the 14 nonmixed
  full-general one-scale flavour channels.

## Correct interpretation

**Implemented in this package**: the all-general renormalisable parameter
switching, a source-resolved intermediate LLSS scalar-pair leading logarithm,
its per-extra-quartic decomposition, the simultaneous scalar bubble, frozen
Wilson mixing switch, one-scale-anchored final C5 and low-energy neutrino mass
observables, all using the 16-model source catalogue.

**Explicit approximations**: unknown GammaCC/Weyl overall factor fixed by
`--spinor-factor` (default +1), candidate extra invariant-to-RGBeta basis
assignment, conditional B/C mixed flavor lift, possible field-basis and
scalar-threshold matching scheme mismatches. Independent complete hard
matching, full dim-5 operator closure and physical three-generation all-seven
master Wilson RGE are NOT supplied here. Neither physical precision nor
one-loop scheme completeness is claimed.

**Locally verified**: 9 pure algebra/mock-control tests and Python syntax
compilation. The complete 16-model numerical integration cannot be verified
without the repository's actual local symbolic/running artifacts.

No existing source modules are overwritten and nothing is written to GitHub.
