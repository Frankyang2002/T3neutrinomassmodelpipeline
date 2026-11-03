# E_+2 — Track B archived-source independent direct UV result

**Status:** Chosen-Weyl, archived-native-CG calculation. **Not physically certified.** Fresh model-specific native action and absolute spinor dictionary pending.

Source: fitted_reference benchmark `E_+2`, native Matchete `T3Y1CG`, `T3Y2CG`, `T3MixCG`; full native identity-portal decomposition (conditional). Frozen matrix used for comparison ONLY after independent calculation.

## Physical charge channels

Q=+2, Q=+1 (charged channels retained; only one channel per A/D model, two per E).

## Zero-portal source point

Matrix in eV, real part:

```text
+0.000870567794213 +0.00750283318319 -0.000794344403878
+0.00750283318319 +0.0148232694166 -0.0208591177511
-0.000794344403878 -0.0208591177511 +0.0188203801113
```

Takagi masses in eV, descending: 0.0386303618094, 0.00744185809056, 0.00332571360322

## Conditional benchmark-identity portal point

Matrix in eV, real part:

```text
+0.000870618508132 +0.00750346018778 -0.000794421490819
+0.00750346018778 +0.01482448645 -0.0208608468577
-0.000794421490819 -0.0208608468577 +0.0188219371142
```

Takagi masses in eV, descending: 0.0386335600259, 0.00744248479616, 0.00332596684255

Full-matrix relative frozen comparison: 2.0296e-11
Maximum independently evaluated (75-digit B0 vs adaptive Feynman integral) scalar-mass kernel residual: 3.028e-16

**Eligibility: NO.** Missing fresh nine-model native action/Dirac indexed-mass measurements and an independently anchored installed Matchete spinor dictionary. The benchmark-to-native portal dictionary is conditional.

The reusable runner recalculates both points, tests complex phases and negative controls, and produces new timestamped results when installed in the actual repository.
