# E_-2 — Track B archived-source independent direct UV result

**Status:** Chosen-Weyl, archived-native-CG calculation. **Not physically certified.** Fresh model-specific native action and absolute spinor dictionary pending.

Source: fitted_reference benchmark `E_-2`, native Matchete `T3Y1CG`, `T3Y2CG`, `T3MixCG`; full native identity-portal decomposition (conditional). Frozen matrix used for comparison ONLY after independent calculation.

## Physical charge channels

Q=+0, Q=-1 (charged channels retained; only one channel per A/D model, two per E).

## Zero-portal source point

Matrix in eV, real part:

```text
+0.000385550192658 -6.82682301807e-05 +0.00884940050415
-6.82682301807e-05 +0.0228212871414 +0.0184928516205
+0.00884940050415 +0.0184928516205 +0.0232234436978
```

Takagi masses in eV, descending: 0.042466029049, 0.00867806075181, 0.00471380876896

## Conditional benchmark-identity portal point

Matrix in eV, real part:

```text
+0.000385611824112 -6.83277806316e-05 +0.00885010592883
-6.83277806316e-05 +0.0228232170663 +0.0184943537947
+0.00885010592883 +0.0184943537947 +0.0232252864375
```

Takagi masses in eV, descending: 0.0424694776252, 0.00867882875803, 0.00471419105535

Full-matrix relative frozen comparison: 2.1194e-11
Maximum independently evaluated (75-digit B0 vs adaptive Feynman integral) scalar-mass kernel residual: 3.028e-16

**Eligibility: NO.** Missing fresh nine-model native action/Dirac indexed-mass measurements and an independently anchored installed Matchete spinor dictionary. The benchmark-to-native portal dictionary is conditional.

The reusable runner recalculates both points, tests complex phases and negative controls, and produces new timestamped results when installed in the actual repository.
