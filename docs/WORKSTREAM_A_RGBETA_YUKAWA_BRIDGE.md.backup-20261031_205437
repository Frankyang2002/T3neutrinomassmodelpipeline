# Workstream A — experimental native RGBeta extra-Yukawa registration

## Scope

The inspected production `RGE/running/rgbeta/T3RGBetaModel.wl` invokes the
ordinary `y1,y2` Yukawa registrations. The generalized Matchete Workstream A
Lagrangian also has a Yukawa to the opposite orientation of the one neutral
complex BSM scalar for alpha=-2 or 0. The independent 16-model census finds
**six extra matrix types** across the model catalogue: one per model, not six
new couplings within any single model.

This stage **does not edit production RGBeta code**. It creates a fresh RGBeta
model in a separate Wolfram process and registers the missing interaction with
native `AddYukawa` with the exact internal `CheckInvariance` option resolved from `Options[AddYukawa]` in each special model. It reads
RGBeta's actual `$yukawas` association and compares the scalar orientation,
fermion ordering, native invariant, chirality, and flavour-index function to
those of the baseline native Yukawa. The Python runner also requires the
previous live Matchete `uv_closure_stage/native_z2/<model>.json` report to have
accepted the extra Yukawa and its full Lagrangian.

## New native RGBeta terms

| Models | Neutral scalar | New coupling | Native Weyl fields | Baseline |
|---|---|---|---|---|
| A/D/E alpha=0 | S1 | `y1Conj` | `{S1,l,F}` | `y1` with `{Bar[S1],l,F}` |
| A/D/E alpha=-2 | S2 | `y2Conj` | `{Bar[S2],l,Fc}` | `y2` with `{S2,l,Fc}` |

The SU(2) invariant is copied from `T3RGBetaYukawaInvariant[dS,dF]`, the
already registered counterpart for the same dimension pair. A singlet
neutral scalar uses the SU(2) epsilon contraction between the two fermion
doublets; a triplet uses the native `delS2` symmetric-doublet contraction.
The scalar hypercharge is zero, so both orientations have the same U(1)
charge, and the same two Weyl fermion types remain in the new interaction.

## Certification boundary

A successful run verifies native RGBeta registration and *invocation of its
optional generator-level `CheckInvariance` routine without emitted errors*.
That is more than a static inventory, but it does **not** show that the
Matchete and RGBeta Yukawa CG normalizations, scalar conjugation matrices,
Dirac/Majorana convention, or coupling signs are physically equivalent.
It does not evaluate one-loop beta functions. The RGBeta beta *pipeline* still
uses the unmodified production model, where these couplings are missing.

Failing `CheckInvariance` is a diagnostic, not evidence that a gauge-invariant
Yukawa operator cannot exist: the native symbolic barred-triplet convention
may need explicit orientation rather than reusing the baseline invariant.
The result should be investigated from the `.json` and `.log` rather than
silently disabling the native check.

## Run from repository root (PowerShell)

```powershell
python Lagrangian\RunT3RGBetaExtraYukawaBridge.py --self-test
python Lagrangian\RunT3RGBetaExtraYukawaBridge.py --all
```

To rerun an existing audit, first back up
`Reports\output\workstream_a\rgbeta_extra_yukawa_bridge` and pass `--overwrite`.

Output:

- `Reports/output/workstream_a/rgbeta_extra_yukawa_bridge/native_rgbeta/*.json`
- `Reports/output/workstream_a/rgbeta_extra_yukawa_bridge/native_rgbeta/*.log`
- `Reports/output/workstream_a/rgbeta_extra_yukawa_bridge/rgbeta_extra_yukawa_bridge.json`
- `Reports/output/workstream_a/rgbeta_extra_yukawa_bridge/rgbeta_extra_yukawa_bridge.md`

Successful status:
`SixExtraRGBetaYukawasNativeRegisteredPhysicalNormalizationUncertified`.
