# A_-4 — Track B archived-source independent direct UV result

**Status:** Chosen-Weyl, archived-native-CG calculation. **Not physically certified.** Fresh model-specific native action and absolute spinor dictionary pending.

Source: fitted_reference benchmark `A_-4`, native Matchete `T3Y1CG`, `T3Y2CG`, `T3MixCG`; full native identity-portal decomposition (conditional). Frozen matrix used for comparison ONLY after independent calculation.

## Physical charge channels

Q=-2 (charged channels retained; only one channel per A/D model, two per E).

## Zero-portal source point

Matrix in eV, real part:

```text
-0.00163357161663 +0.00233284998816 +0.00720345095976
+0.00233284998816 +0.0176823659023 +0.0239634714247
+0.00720345095976 +0.0239634714247 +0.0194917378235
```

Takagi masses in eV, descending: 0.0435971121485, 0.00758214355407, 0.000474436485301

## Conditional benchmark-identity portal point

Matrix in eV, real part:

```text
-0.00163369049292 +0.0023330770605 +0.00720405016534
+0.0023330770605 +0.0176838041814 +0.0239653970597
+0.00720405016534 +0.0239653970597 +0.0194933193139
```

Takagi masses in eV, descending: 0.0436006440369, 0.0075827281046, 0.000474482929934

Full-matrix relative frozen comparison: 8.8736e-12
Maximum independently evaluated (75-digit B0 vs adaptive Feynman integral) scalar-mass kernel residual: 1.765e-16

**Eligibility: NO.** Missing fresh nine-model native action/Dirac indexed-mass measurements and an independently anchored installed Matchete spinor dictionary. The benchmark-to-native portal dictionary is conditional.

The reusable runner recalculates both points, tests complex phases and negative controls, and produces new timestamped results when installed in the actual repository.
