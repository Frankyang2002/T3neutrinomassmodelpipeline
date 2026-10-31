# B(−1)/C(−1): native Matchete Y-only threshold experiment

## Purpose and scope

The conditional numerical sequence through EFT2 and SM+Weinberg evolution has already been executed at:

| Model | muY [GeV] | corrected muX [GeV] | conditional masses [eV], descending |
|---|---:|---:|---|
| B(−1) | 1187.80578221 | 1035.77775471 | 0.0385336229, 0.0156138215, 0.0004584847 |
| C(−1) | 1185.39334545 | 1032.85523509 | 0.0464057419, 0.0459238479, 0.0216846442 |

A subsequent *candidate* scalar-pair one-loop diagnostic produced:

| Model | finite B0(X,Y;muY) | candidate mass spectrum [eV] |
|---|---:|---|
| B(−1) | 0.1317561 | 0.0393113530, 0.0157895632, 0.0005317403 |
| C(−1) | 0.1324234 | 0.0470776247, 0.0467726114, 0.0223335971 |

These outputs are **not physically certified**. The missing step is a complete EFT1 -> EFT2 one-loop threshold matching, separately from the later EFT2 -> SM threshold.

This experiment runs *native Matchete* on the saved F-first canonical EOMEFT with **Y heavy and X light**, attempting a one-loop EFT2 raw operator extraction. It does **not** recompute the numerical RG trajectory and does **not** convert the raw Matchete expressions into a renormalized correction to the running Wilson coefficient.

## Existing source reused

- `Lagrangian/RunGeneralT3EFT1ScalarThresholdProbe.wl`: currently classifies BOTH X and Y as heavy and integrates them out together. Its native source includes CG reconstruction, mass-basis field registration, canonical EOMEFT import, hermiticity checks, and the `RunT3Matching[scalarEFT,5,0/1]` stages.
- `scripts/RunGeneralT3EFT1ScalarThresholdProbe.py`: existing input preparation and source-certificate checks.
- `Lagrangian/RunMatching.wl`: the original Matchete driver.

All source files above remain untouched. The new builder uses strict single-occurrence source-structure guards to create a **complete standalone** `Lagrangian/RunT3BCYOnlyMatching.experimental.wl` from the installed native program. The generated standalone program is a new file, not a patch to the original. If the original script changes structurally, generation aborts. Its source fingerprint is recorded.

## Exact experimental changes to Matchete setup

1. Re-register X as `Mass->{Light,MSX}`, while Y remains `Mass->{Heavy,MSY}`.
2. Retain the source's charge- and SU(2)-covariant native interactions, epsilon tensors and fermion-first EOMEFT operators.
3. Reconstruct heavy-Y versus light-X free-field Lagrangian (the canonical first-stage input classified both scalars light).
4. Require `CheckLagrangian` with `HeavyMassBasis`, hermiticity, canonical normalization and the original gauge checks before matching.
5. Evaluate the existing `RunT3Matching` machinery with `EFTOrder=5`, `LoopOrder=0` and `1`.
6. Assert that the matched EFT retains X but has no Y or integrated F.
7. Preserve complete matched EFT2 tree and one-loop expressions; extract raw `LLHH`, `LLHX`, `LLXX` one-loop sectors for follow-up operator projection.
8. Retain `physical_matching_performed=False`, along with an explicit uncertified driver certificate.

The optional `--scope leading-c5` deliberately removes 122 four-BSM-scalar spectator terms *using the source's existing exact operator partition*. That scope is only suitable for an isolated LLHH test. **Do not use its LLHX or LLXX output as a complete EFT2 Wilson seed**, because Higgs-free `XXYY` quartics can generate such operators when Y is integrated out.

## What is *not* yet done

- **Numerical boundary parity:** the Matchete probe's saved `final` parameter file is *not* identical to `continuous_v1`'s running values at muY. Do not substitute these symbolic expressions into the current numerical result without matching the masses, scalar rotation, all quartics, input Wilson operators, flavor normalization and subtraction scale.
- **Renormalized one-loop matching:** subtract EFT2 low-energy loop amplitudes from EFT1 amplitudes with the same infrared regulator and basis, apply counterterms, and isolate the finite local Wilson coefficients at `muY`.
- **EFT2 operator inventory:** any surviving renormalizable and higher-dimensional operators must be included consistently in RG running. The existing numerical `(3,3,8,8,6)` Wilson closure is conditional, not a complete Green's-basis EFT.
- **X threshold:** a separate verified `EFT2(H,X) -> SMEFT` Matchete matching is still missing, along with dimension-four threshold shifts, field normalization and systematic power counting.
- **Three-family mixing:** the existing Wilson master treats internal SM Yukawas with an effective one-generation closure, even though it preserves six external lepton flavor pairs.

## Installation and execution (PowerShell)

Run from `C:\thesistemp\T3neutrinomassmodelpipeline`:

```powershell
Expand-Archive `
  -LiteralPath "$env:USERPROFILE\Downloads\t3_bc_y_only_native_matchete_probe_v1.zip" `
  -DestinationPath . -Force

python -m unittest discover -s tests -p "test_t3_bc_y_only_matchete_builder.py" -v

python -u -X utf8 -m scripts.RunT3BCYOnlyMatchete `
  --all --scope full --preflight-only

python -u -X utf8 -m scripts.RunT3BCYOnlyMatchete `
  --all --scope full
```

If the **full** EFT fails the native hermiticity or matching gates, do not bypass the failure. To obtain a *restricted LLHH diagnostic*, use a separate output directory:

```powershell
python -u -X utf8 -m scripts.RunT3BCYOnlyMatchete `
  --all --scope leading-c5 `
  --output-dir "output\general_t3_mass_eigen_matchete\eft1_scalar_Y_only_probe_leading_c5"
```

The native experiment may be computationally intensive. The driver does not run physics in the background or modify GitHub. The result directory is `output/general_t3_mass_eigen_matchete/eft1_scalar_Y_only_probe/` for full mode.

A successful native result is **an important intermediate step, not full physical completion**. Only report finished finite C5 matching after the subtraction, Wilson projection and numerical running-parameter handoff have been certified.

## Source and test provenance

The stage was prepared by inspecting the publicly available repository source on 2026-10-31, including its exact `Mass->{Heavy,MSX}`/`Mass->{Heavy,MSY}` setup, `CheckLagrangian` gates, and raw LLHH extraction. The 12 offline tests use a structural fixture, not a full Wolfram kernel. Native Matchete execution and cross-repository source-anchor validation must be performed on the Windows-local repository.
