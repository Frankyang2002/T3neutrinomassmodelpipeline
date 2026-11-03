# A_-2 — Track B archived-source independent direct UV result

**Status:** Chosen-Weyl, archived-native-CG calculation. **Not physically certified.** Fresh model-specific native action and absolute spinor dictionary pending.

Source: fitted_reference benchmark `A_-2`, native Matchete `T3Y1CG`, `T3Y2CG`, `T3MixCG`; full native identity-portal decomposition (conditional). Frozen matrix used for comparison ONLY after independent calculation.

## Physical charge channels

Q=-1 (charged channels retained; only one channel per A/D model, two per E).

## Zero-portal source point

Matrix in eV, real part:

```text
+0.00263371275851 -0.00129842604289 +0.0107623285797
-0.00129842604289 +0.0207469558781 +0.0264700322119
+0.0107623285797 +0.0264700322119 +0.0213822546125
```

Takagi masses in eV, descending: 0.0485563028921, 0.0110075756786, 0.0072141960356

## Conditional benchmark-identity portal point

Matrix in eV, real part:

```text
+0.00263393794824 -0.00129848399436 +0.0107632249743
-0.00129848399436 +0.0207485807307 +0.0264721597507
+0.0107632249743 +0.0264721597507 +0.0213839860467
```

Takagi masses in eV, descending: 0.0485602090998, 0.0110084578426, 0.00721475346842

Full-matrix relative frozen comparison: 8.5303e-12
Maximum independently evaluated (75-digit B0 vs adaptive Feynman integral) scalar-mass kernel residual: 1.765e-16

**Eligibility: NO.** Missing fresh nine-model native action/Dirac indexed-mass measurements and an independently anchored installed Matchete spinor dictionary. The benchmark-to-native portal dictionary is conditional.

The reusable runner recalculates both points, tests complex phases and negative controls, and produces new timestamped results when installed in the actual repository.
