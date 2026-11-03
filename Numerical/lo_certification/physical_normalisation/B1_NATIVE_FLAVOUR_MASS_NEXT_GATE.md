# B1 — live diagonal flavour-mass declaration test

## Why this is the next physical gate

The three independently computed A_+0, B_-1 and C_-1 UV mass matrices use
`MF = diag(100,110,120) TeV` in their benchmark loops. The previously live
Matchete spinor action witnesses were generated with
`T3_EXPERIMENTAL_FLAVOR_MASSES=0`, where `FreeLag` has one unindexed `MF` mass
symbol. This leaves a source-to-benchmark assumption for non-degenerate heavy
flavours, even though the actual Python evaluation already resolves their
individual propagators.

`Lagrangian/model/T3Fields.wl` implements a *separate opt-in branch*:

```wl
T3ExperimentalFlavorMassesQ[] :=
  Environment["T3_EXPERIMENTAL_FLAVOR_MASSES"] === "1";
T3MassSpec[label_String, mass_] := Module[{kind},
  kind = If[MemberQ[$T3HeavyFields, label], Heavy, Light];
  If[label === "F" && TrueQ[T3ExperimentalFlavorMassesQ[]],
    {kind, mass, {NFlavor}}, {kind, mass}]];
```

An accepted native `FreeLag` with `Coupling[MF,{Index[r,NFlavor]},0]`
contracted with both fermion fields at **the same** flavour index would
establish that the action permits independent *diagonal* masses. It would
**not** establish a general off-diagonal Majorana or Dirac mass matrix.

## Scope and acceptance

The new standalone Wolfram script deliberately rebuilds only the three
restricted UV actions with flags `general fermions/quadratics/quartics = 0`,
`flavour masses = 1` in a fresh Matchete kernel. Each action must pass its
native `CheckLagrangian`, and `LFree` must have exactly one mass monomial,
`-MF[r] Bar[F_r] F_r` (Dirac A) or `-MF[r]/2 F_r^T GammaCC F_r`
(Majorana B/C). The test checks the **actual native term**: one `MF` flavour
index matching both field flavour indices. Wrong prefactors, unindexed mass
terms, missing/extra flavour indices, and unsuccessful native validation
all produce a *blocked*, nonzero-exit witness.

The result records Git HEAD and hashes of the six used T3 source files.
It never calls Match, RGEs, RGBeta or the EFT pipeline; no previous report
is overwritten.

**Important:** It does not numerically substitute the three specific
benchmark eigenvalues inside Matchete. A successful witness supports the
possibility of the diagonal nondegenerate assignment but does not certify
the chosen benchmark values or the remaining absolute `GammaCC` phase.

## Windows command

After installing `scripts/RunPhysicalB1NativeFlavourMassContract.wl` under
the repository root, run:

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
$stamp = Get-Date -Format 'yyyyMMdd_HHmmss'
$out = "output\general_observable\physical_normalisation\b1_native_flavour_mass_$stamp\native_flavour_mass_action.json"
wolframscript -file scripts\RunPhysicalB1NativeFlavourMassContract.wl $out
if ($LASTEXITCODE -ne 0) { throw 'Live native diagonal flavour-mass gate blocked' }
```

If it fails, inspect the JSON `reason` and the native console output. Do not
change `T3Fields.wl` or accept a failed build as proof of mass generality.

## Gate ledger

- G1/G3 vertex/scalar witnesses: previously passed for three representatives.
- Complex chosen-Weyl consistency: previously passed; native absolute `GammaCC`
  4x4 phase not exposed.
- G2/G3 nondegenerate benchmark mass: conditional until this native action
  check and any required additional numerical native mass assignment.
- G4 certified: 0/16. This witness alone cannot change that count.
