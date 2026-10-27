# T3 mass-eigenstate provenance path correction

## Cause

`scripts/RunGeneralT3MixedScalarMassBasis.py` originally serialized provenance
keys via `str(p.relative_to(ROOT))`. On Windows this produced backslash-separated
keys, whereas `scripts/RunGeneralT3MassEigenMatchete.py` looked up fixed
forward-slash keys and falsely reported a changed registry.

## Fix

- Rotation JSON producer now writes portable POSIX-style relative keys.
- Matchete runner normalizes existing Windows and POSIX path keys before
  comparing each SHA-256 against the current file bytes.
- Missing source records, altered digests, traversal paths, malformed values,
  and duplicate aliases still fail closed. No bypass option is introduced.
- Existing rotation results on disk remain readable without regeneration.

## Installation

Unzip from repository root with `Expand-Archive ... -DestinationPath . -Force`.
All files in this archive are complete source files, not patches.

```powershell
python -m pytest -q tests\python\test_general_t3_mass_eigen_matchete.py tests\python\test_general_t3_mixed_scalar_mass_basis.py
python -m scripts.RunGeneralT3MassEigenMatchete --all --prepare-only
python -m scripts.RunGeneralT3MassEigenMatchete --model B_-1
python -m scripts.RunGeneralT3MassEigenMatchete --model C_-1
python -m scripts.RunGeneralT3MassEigenMatchete --all
```

This only fixes provenance; live Matchete heavy-basis validation and Wilson
threshold matching remain unverified until actually executed.
