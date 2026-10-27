# Native HHXX topology CG reconstruction — experimental B_-1/C_-1

## Evidence from local Windows logs

On both mixed models, the original 11-interaction baseline is valid, and the diagonal `FreeLag` passes. The isolated transformed Higgs portal groups also pass. The first rejected transformed term is `lambdaT3|CG:T3MixCG` (`DefineCG::notInvariant`). There is no evidence in those logs that the additional seven quartic operators have failed independently.

## Physics / source change

Original topology vertex: `PlusHc[lambdaT3[] H H S1 Bar[S2] CG[T3MixCG,...]]`.

For `B_-1/C_-1`, define `Phi=(S1,epsilon Conjugate[S2])^T = U (X,Y)^T`, with `epsilon={{0,1},{-1,0}}`. The original `Bar[S2]` is `-epsilon (U21 X+U22 Y)` in this component convention.

We now *remove only the original lambdaT3-expanded terms* from the interaction expression before applying the generic field substitution. We re-express them in a **new Matchete-native four-scalar invariant**:

```
O_HH := H[i] H[j] CG[T3MENativeHHMassCG,{i,j,k,m}]
V_topology in new fields proportional to
    -lambdaT3 O_HH [ U11 U21 X[k] X[m]
                    +(U11 U22 + U12 U21) X[k] Y[m]
                    +U12 U22 Y[k] Y[m] ] + h.c.
```

The code computes the exact scalar tensor factor `alphaCG` relating the original native `H H S1 Bar[S2]` tensor contracted with `epsilon` to the new Matchete tensor. It explicitly requires full array equality and last-pair bosonic symmetry. No normalization is guessed. The four U entries remain symbolic Matchete couplings with the earlier separate numerical mapping.

The remaining 10 baseline quartics, 7 additional quartics and 2 Yukawa vertices remain in the pre-existing rotation pathway. Nothing is changed in RGBeta, EFT1 running or the numerical mass eigenvalues.

## Standalone Matchete check performed in this environment

Temporarily installed Matchete in a stateless Wolfram kernel and evaluated its native `InvariantTensors` routines for the relevant SU(2) reps. They returned one original tensor and one new tensor. Contracting the old tensor's fourth index with `epsilon` gave exactly the new tensor (factor `1`), and the new tensor was symmetric in the two scalar slots:

```
{oldInvariantCount, newInvariantCount, factor, tensorIdentity, bosonicSymmetry}
= {1,1,1,True,True}
```

The same environment **cannot run the complete `CheckLagrangian`**: its Matchete model-check procedure returned `$Failed` with `FrontEndObject::notavail`. It also has no access to the local Windows project and may use a different Matchete package version. This is a native invariant-tensor check, not an end-to-end certification.

## Validation and limitations

- The Python runner now requires the Mathematica export `native_topology_epsilon_tensor` to state that the exact tensor identity and scalar-pair symmetry were verified. Failing or missing evidence is rejected.
- The existing strict inventories, checks for eliminated old scalar fields, and actual `CheckLagrangian`/`HeavyMassBasis` requirements are unchanged.
- If a different operator now fails, the failure-isolation path uses the same newly built native topology expression and reports the next specific coupling/CG group.
- Matchete-to-RGBeta normalisation, two-scale threshold matching and physical Wilson coefficients remain **uncertified/not performed**.

## Install in project root (PowerShell)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath 'C:\path\to\T3_Matchete_Native_Topology_CG_Fix.zip' -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_mass_eigen_matchete.py tests\python\test_general_t3_mixed_scalar_mass_basis.py
python -m scripts.RunGeneralT3MassEigenMatchete --model B_-1
python -m scripts.RunGeneralT3MassEigenMatchete --model C_-1
```

Results/logs: `output/general_t3_mass_eigen_matchete/`. This is a diagnostic mass-basis validation, **not matching**.
