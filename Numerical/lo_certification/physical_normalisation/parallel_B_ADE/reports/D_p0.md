# D_+0 — Track B archived-source independent direct UV result

**Status:** Chosen-Weyl, archived-native-CG calculation. **Not physically certified.** Fresh model-specific native action and absolute spinor dictionary pending.

Source: fitted_reference benchmark `D_+0`, native Matchete `T3Y1CG`, `T3Y2CG`, `T3MixCG`; full native identity-portal decomposition (conditional). Frozen matrix used for comparison ONLY after independent calculation.

## Physical charge channels

Q=+1 (charged channels retained; only one channel per A/D model, two per E).

## Zero-portal source point

Matrix in eV, real part:

```text
+0.00263385240637 -0.00129901656816 +0.0107607055667
-0.00129901656816 +0.0207426073553 +0.0264645251308
+0.0107607055667 +0.0264645251308 +0.0213774930141
```

Takagi masses in eV, descending: 0.0485459973313, 0.0110061574466, 0.00721411289109

## Conditional benchmark-identity portal point

Matrix in eV, real part:

```text
+0.00263407760736 -0.00129907456835 +0.0107616018258
-0.00129907456835 +0.020744231849 +0.0264666522259
+0.0107616018258 +0.0264666522259 +0.0213792240634
```

Takagi masses in eV, descending: 0.0485499027029, 0.0110070394968, 0.00721467031372

Full-matrix relative frozen comparison: 7.2867e-12
Maximum independently evaluated (75-digit B0 vs adaptive Feynman integral) scalar-mass kernel residual: 1.765e-16

**Eligibility: NO.** Missing fresh nine-model native action/Dirac indexed-mass measurements and an independently anchored installed Matchete spinor dictionary. The benchmark-to-native portal dictionary is conditional.

The reusable runner recalculates both points, tests complex phases and negative controls, and produces new timestamped results when installed in the actual repository.
