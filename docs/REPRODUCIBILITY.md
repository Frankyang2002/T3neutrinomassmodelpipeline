# Reproducing a T3 calculation

A numerical result is reproducible only with its **source calculation**, **configuration** and **numerical environment**. This is especially important when symbolic exports are reused by later numerical work.

## Record for each run

Retain the following together with the result:

- Model class, hypercharge parameter, electroweak field representations and relevant symmetries.
- UV couplings, allowed interactions set to zero, heavy mass matrices, flavour conventions and initial scale.
- Threshold order and scales, including the complete fermion-first and scalar-second field content.
- Source coefficient files (or their exact reproducible generation commands) and fingerprints of the relevant Matchete/RGBeta source inputs.
- Included Wilson operators, perturbative order and the matching/renormalisation convention.
- Fit parameters, observables, solver tolerances, seed and each model's accepted/blocked result status.
- Git source commit and Python/Wolfram/package versions, where available.

## Environment record (PowerShell)

Run from the repository root after activating the environment:

```powershell
New-Item -ItemType Directory -Force output\maintenance | Out-Null
python --version | Out-File -Encoding utf8 output\maintenance\python_version.txt
python -m pip freeze | Out-File -Encoding utf8 output\maintenance\python_packages.txt
git rev-parse HEAD | Out-File -Encoding utf8 output\maintenance\source_commit.txt
```

Do not label a rerun using cached symbolic exports as a newly regenerated matching calculation. Check provenance flags in generated JSON and manifests. Reusing cached coefficient files can be scientifically valid, but the source conventions must match the numerical adapters and benchmark being evaluated.

## Verification commands

```powershell
python -X utf8 pipeline.py --help
python -X utf8 full16_symbolic.py --help
.\scripts\RunCoreTests.ps1
```

These commands verify interfaces and a selected software regression suite; they do not execute or certify every symbolic and numerical calculation.

To rebuild presentation artifacts from saved results, follow [Report regeneration](REGENERATE_CURRENT_REPORTS.md). To understand the difference between tests, numerical evaluations and physical matching evidence, see [Capabilities](IMPLEMENTATION_STATUS.md) and [Verification](TEST_STRATEGY.md).
