# D_+2 — Track B archived-source independent direct UV result

**Status:** Chosen-Weyl, archived-native-CG calculation. **Not physically certified.** Fresh model-specific native action and absolute spinor dictionary pending.

Source: fitted_reference benchmark `D_+2`, native Matchete `T3Y1CG`, `T3Y2CG`, `T3MixCG`; full native identity-portal decomposition (conditional). Frozen matrix used for comparison ONLY after independent calculation.

## Physical charge channels

Q=+2 (charged channels retained; only one channel per A/D model, two per E).

## Zero-portal source point

Matrix in eV, real part:

```text
-0.00163357161663 +0.00233278151724 +0.00720237733022
+0.00233278151724 +0.0176807599484 +0.023960355913
+0.00720237733022 +0.023960355913 +0.0194882599873
```

Takagi masses in eV, descending: 0.0435912883619, 0.00758115334107, 0.00047468670177

## Conditional benchmark-identity portal point

Matrix in eV, real part:

```text
-0.00163369049292 +0.00233300858802 +0.00720297644497
+0.00233300858802 +0.0176821980767 +0.02396228129
+0.00720297644497 +0.02396228129 +0.0194898411969
```

Takagi masses in eV, descending: 0.043594819764, 0.00758173780894, 0.000474733174385

Full-matrix relative frozen comparison: 7.8374e-12
Maximum independently evaluated (75-digit B0 vs adaptive Feynman integral) scalar-mass kernel residual: 1.765e-16

**Eligibility: NO.** Missing fresh nine-model native action/Dirac indexed-mass measurements and an independently anchored installed Matchete spinor dictionary. The benchmark-to-native portal dictionary is conditional.

The reusable runner recalculates both points, tests complex phases and negative controls, and produces new timestamped results when installed in the actual repository.
