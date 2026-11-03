# D_-2 — Track B archived-source independent direct UV result

**Status:** Chosen-Weyl, archived-native-CG calculation. **Not physically certified.** Fresh model-specific native action and absolute spinor dictionary pending.

Source: fitted_reference benchmark `D_-2`, native Matchete `T3Y1CG`, `T3Y2CG`, `T3MixCG`; full native identity-portal decomposition (conditional). Frozen matrix used for comparison ONLY after independent calculation.

## Physical charge channels

Q=+0 (charged channels retained; only one channel per A/D model, two per E).

## Zero-portal source point

Matrix in eV, real part:

```text
+0.00263477346883 -0.0013028603072 +0.0107501266068
-0.0013028603072 +0.0207143026159 +0.0264284156451
+0.0107501266068 +0.0264284156451 +0.0213461410691
```

Takagi masses in eV, descending: 0.048478472112, 0.0109968725715, 0.00721361761321

## Conditional benchmark-identity portal point

Matrix in eV, real part:

```text
+0.00263499874409 -0.00130291862453 +0.0107510219271
-0.00130291862453 +0.0207159247743 +0.026430540123
+0.0107510219271 +0.026430540123 +0.0213478698232
```

Takagi masses in eV, descending: 0.0484823724024, 0.0109977539447, 0.00721417488384

Full-matrix relative frozen comparison: 7.5435e-12
Maximum independently evaluated (75-digit B0 vs adaptive Feynman integral) scalar-mass kernel residual: 1.765e-16

**Eligibility: NO.** Missing fresh nine-model native action/Dirac indexed-mass measurements and an independently anchored installed Matchete spinor dictionary. The benchmark-to-native portal dictionary is conditional.

The reusable runner recalculates both points, tests complex phases and negative controls, and produces new timestamped results when installed in the actual repository.
