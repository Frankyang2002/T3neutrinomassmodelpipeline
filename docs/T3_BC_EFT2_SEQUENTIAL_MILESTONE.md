
## Experimental B/C EFT2 sequential full numerical continuation (conditional, 2026-10-31)

Added local-only files `Numerical/experimental/T3BCEFT2HX.py`,
`scripts/RunT3BCEFT2Sequential.py`, and standalone regression tests.

The stage consumes existing exact Y-threshold projection outputs and replays
`continuous_v1` EFT1 scalar RG data with unchanged SHA256 sources. It rotates
and fits the complete H,X quartic tensor to the six-real-parameter Z2 inert
doublet basis, evolves the **Y-decoupled** renormalisable one-loop EFT2 and
project seven-term `psi2phi2` Wilson tensor, locates the corrected running X
threshold, and runs the conditional Weinberg coefficient to 100 GeV.

**Explicit limitations:** finite Y and X scalar matching not derived;
additional-quartic physical convention flags remain false; the initial Higgs
quadratic parameter is a tree-level proxy; external three-family flavour
indices are retained but the existing gauge-index master uses projected SM
Yukawa anomalous dimensions. The output is conditional; no physical neutrino
mass prediction is certified.

Standalone local algebra tests passed (14/14) at package creation. The
end-to-end Windows numerical run and the native Wolfram finite matching
have **not** been performed in the packaging environment. Production sources
and old outputs have not been modified.
