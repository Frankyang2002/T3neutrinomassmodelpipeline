# E_+0 — Track B archived-source independent direct UV result

**Status:** Chosen-Weyl, archived-native-CG calculation. **Not physically certified.** Fresh model-specific native action and absolute spinor dictionary pending.

Source: fitted_reference benchmark `E_+0`, native Matchete `T3Y1CG`, `T3Y2CG`, `T3MixCG`; full native identity-portal decomposition (conditional). Frozen matrix used for comparison ONLY after independent calculation.

## Physical charge channels

Q=+1, Q=+0 (charged channels retained; only one channel per A/D model, two per E).

## Zero-portal source point

Matrix in eV, real part:

```text
+0.00042517915895 -3.42291696551e-05 +0.00881832939275
-3.42291696551e-05 +0.0228366700569 +0.018514913904
+0.00881832939275 +0.018514913904 +0.0232546723168
```

Takagi masses in eV, descending: 0.0425119055002, 0.00864782382377, 0.00464320779138

## Conditional benchmark-identity portal point

Matrix in eV, real part:

```text
+0.000425243985633 -3.42856330624e-05 +0.00881903273675
-3.42856330624e-05 +0.0228385997072 +0.0185164170612
+0.00881903273675 +0.0185164170612 +0.023256518334
```

Takagi masses in eV, descending: 0.0425153568616, 0.00864858944624, 0.00464358428104

Full-matrix relative frozen comparison: 2.5835e-11
Maximum independently evaluated (75-digit B0 vs adaptive Feynman integral) scalar-mass kernel residual: 3.028e-16

**Eligibility: NO.** Missing fresh nine-model native action/Dirac indexed-mass measurements and an independently anchored installed Matchete spinor dictionary. The benchmark-to-native portal dictionary is conditional.

The reusable runner recalculates both points, tests complex phases and negative controls, and produces new timestamped results when installed in the actual repository.
