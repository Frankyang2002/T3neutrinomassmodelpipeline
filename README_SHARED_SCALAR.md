# Shared-scalar / scotogenic mode

Two-number dimension input now means one physical scalar multiplet:

```powershell
python pipeline.py --dims 2 1
```

Current production scope is `dS=2` with `dF=1` or `dF=3`.  The hypercharge
parameter is fixed to `alpha=-1`, so the physical field assignment is

- `S ~ (2,+1/2)`
- formal matching leg `S1 = i sigma_2 S*`
- formal matching leg `S2 = S`
- `F ~ (1,0)` or `(3,0)`

The formal two-leg Matchete topology is retained for matching because its
normalisation/CG bridge is already regression-tested.  Physical field counting
in RGBeta and the component EFT RGE uses one scalar block only.

## Threshold syntax

Shared-scalar mode uses the physical scalar name `S`:

```powershell
python pipeline.py --dims 2 1 --threshold F --threshold S
```

Internally only, `S` is expanded to the formal Matchete pair `S1,S2` at the
scalar matching threshold.  User-facing stage labels are therefore
`EFT_1_after_F` and `EFT_2_after_S`.

Three-number mode is unchanged:

```powershell
python pipeline.py --dims 2 2 1 --alpha -1 --threshold F --threshold S1 S2
```

and continues to mean two independent scalar multiplets.

## Shared-scalar RGBeta basis

The physical one-scalar RGBeta model uses the Ma/scotogenic notation
`h, mSSq, lambdaH, lambdaS, lambda3, lambda4, lambda5`.  For the singlet
fermion branch the gauge coefficients are regression-tested to be

`(b1,b2,b3) = (7,-3,-7)`

in the project convention `Q=T3+Y` and `g1=gY`.

## Validation performed here

- all generated Python files pass `py_compile`;
- `tests/python/test_shared_scalar_mode.py` passes 3 tests;
- Wolfram files have balanced `[]`, `{}` and `()` delimiters;
- full Wolfram/RGBeta/Matchete execution was not available in this environment,
  so run the commands below locally before using the branch for final results.
