# Verification and test coverage

The test suites check different levels of evidence. A successful software test does **not** imply that the full symmetry-allowed model, one-loop matching or numerical phenomenology has been physically certified.

## Core regression suite

Run from the repository root in Windows PowerShell:

```powershell
.\scripts\RunCoreTests.ps1
```

This executes a selected group of tests for model identities, charge and threshold conventions, pipeline orchestration, numerical state representations, Weinberg running, coupling inventories and the fermion-first intermediate EFT.

The last locally reported execution on 8 November 2026 completed **189 tests successfully**. This is a fast compatibility check, not a complete physical validation and not a promise that every future checkout has the same result.

## Extended suite

```powershell
.\scripts\RunCoreTests.ps1 -Full
```

The extended suite includes model-specific calculations, older study interfaces, numerical references, generated-artifact integration checks and source-dependent tests. Some tests assume saved Matchete/RGBeta exports or certificates exist at particular paths. A missing artifact should be reported separately from a wrong numerical result.

At the last local extended run reported on 8 November 2026, **4,453 tests passed, 73 failed, 51 were skipped and one had a collection error**. These results belong to that checkout and are not a claim about the current status of every environment. The historical diagnostic record is in `TEST_FAILURE_TRIAGE_2026_11_08.md`.

## What must be checked for the target physics

1. Lorentz, gauge and discrete-symmetry invariance of the UV operator inventory, including independent SU(2) contractions.
2. Correct fermion mass, Yukawa, charge and flavour conventions in the physical mass basis.
3. The ordered EFT transition: UV theory **to the scalar-only intermediate EFT** by integrating out all heavy fermions; then **to SMEFT** by integrating out all heavy scalars.
4. Matching and RGE normalisations of the renormalisable parameters and induced Wilson coefficients, with explicit perturbative order.
5. Correct Weinberg-coefficient orientation and convention, \(m_\nu=-v^2C_5/2\).
6. Numerical stability, threshold continuity, unit conventions and reproducibility of saved inputs and fits.

## Evidence levels

| Check type | What it establishes |
| --- | --- |
| Unit or algebraic test | A stated local numerical or symbolic identity |
| Synthetic interface test | Input validation and control flow for constructed inputs |
| Saved-artifact regression | Compatibility with a particular preserved source/output set |
| Live Wolfram/Matchete/RGBeta run | Behaviour of the installed symbolic tools for a specified model and input |
| Full physical verification | Requires checking theory completeness, native normalisations, matching, running and controlled numerical convergence together |

Tests belonging to superseded studies can be retained in an external, restorable research archive when their results are no longer part of the supported workflow. Never discard a unique covariance, matching or normalisation check based on its filename alone.
