# E_-4 — Track B archived-source independent direct UV result

**Status:** Chosen-Weyl, archived-native-CG calculation. **Not physically certified.** Fresh model-specific native action and absolute spinor dictionary pending.

Source: fitted_reference benchmark `E_-4`, native Matchete `T3Y1CG`, `T3Y2CG`, `T3MixCG`; full native identity-portal decomposition (conditional). Frozen matrix used for comparison ONLY after independent calculation.

## Physical charge channels

Q=-1, Q=-2 (charged channels retained; only one channel per A/D model, two per E).

## Zero-portal source point

Matrix in eV, real part:

```text
+0.000870567794213 +0.00747705598962 -0.000788351090202
+0.00747705598962 +0.014767229335 -0.0207807091863
-0.000788351090202 -0.0207807091863 +0.0187503037372
```

Takagi masses in eV, descending: 0.0384852687548, 0.00741631220069, 0.00331914431232

## Conditional benchmark-identity portal point

Matrix in eV, real part:

```text
+0.000870618508132 +0.00747768098906 -0.000788427618294
+0.00747768098906 +0.0147684413674 -0.0207824317754
-0.000788427618294 -0.0207824317754 +0.0187518549514
```

Takagi masses in eV, descending: 0.0384884547869, 0.00741693705784, 0.0033193970979

Full-matrix relative frozen comparison: 1.5003e-11
Maximum independently evaluated (75-digit B0 vs adaptive Feynman integral) scalar-mass kernel residual: 3.028e-16

**Eligibility: NO.** Missing fresh nine-model native action/Dirac indexed-mass measurements and an independently anchored installed Matchete spinor dictionary. The benchmark-to-native portal dictionary is conditional.

The reusable runner recalculates both points, tests complex phases and negative controls, and produces new timestamped results when installed in the actual repository.
