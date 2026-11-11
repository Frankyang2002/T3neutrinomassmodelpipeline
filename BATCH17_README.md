# Batch 17: conditional Wilson matching and SMEFT

This is a conditional calculation, **not a certified physical neutrino mass**.

## Project paths

- `Numerical/experimental/Batch17ConditionalWilson.py`
- `scripts/RunBatch17ConditionalC5All16.py`
- `tests/test_batch17_conditional_wilson.py`

## Dependencies

Previously successful Batches 13–16, and the repository's existing
`Numerical/joint_scalar_matching/fermion_boundary_bridge.py`,
`Numerical/joint_scalar_matching/large_batch.py`,
`Numerical/running/SMWeinbergEvolution.py`, and `physics/NeutrinoMass.py`.

External model-specific tensor inputs are **mandatory** to calculate conditional
C5: a full source-resolved Matchete LLSS release (`summary.json` and
`<slug>/boundary.json`) and a matching explicit HHSS background
`<slug>.json` + its adjacent `.npz`. They are not created by Batch 13–16.

## PowerShell commands (from project root)

```powershell
python -m unittest discover -s tests -p test_batch17_conditional_wilson.py -v
python scripts\RunBatch17ConditionalC5All16.py --model A_-4
```

With the required existing tensor inputs:

```powershell
python scripts\RunBatch17ConditionalC5All16.py `
  --model A_-4 `
  --source-release "C:\path\to\existing\source_release" `
  --background "C:\path\to\matching\Am4.json"
```

For all 16:

```powershell
python scripts\RunBatch17ConditionalC5All16.py `
  --source-release "C:\path\to\existing\source_release" `
  --background-dir "C:\path\to\matching_backgrounds"
```

Output: `output/batch17/conditional_c5/`.

## Exact interpretation

The existing joint matching code computes only the LLSS quartic-induced one-loop
leading logarithm and **one simultaneous finite scalar bubble**, per unit
unknown signed Matchete GammaCC-to-Weyl factor. The new driver checks against
the numerically evolved Batch 16 scalar eigenvalue spectrum, then integrates
the resulting conditional Weinberg operator down to 100 GeV using the existing
SMEFT runner and reports Takagi masses **per unit unknown factor**.

It does not run the full seven-term EFT1 Wilson system in the new Batch 14
trajectory, compute independent hard one-loop heavy-fermion matching, establish
the Matchete/RGBeta HHSS tensor normalisations, or certify source flavour
parameters against new Batch 13 UV initial conditions. Those prerequisites
cannot be inferred from RG running success. The driver does not identify
unverified scalar mass eigenvectors with Matchete source components.

A BLOCKED model is not a failed RGE; it means the required same-scale source
and scalar matching tensors are missing or inconsistent with the current
numerical endpoint.
