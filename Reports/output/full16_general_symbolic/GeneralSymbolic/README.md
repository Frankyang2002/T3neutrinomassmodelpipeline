# Full general-interaction T3 symbolic study

Validated: **16/16** model(s).

The BSM UV Lagrangian and **UV/EFT1 renormalisable beta functions** come from the general-interaction Matchete/RGBeta workflows, not from the restricted `pipeline.py` matching outputs.

The beta expressions contain their group-dependent coefficients, but no independent Matchete-to-RGBeta physical coupling-map or term-by-term group-factor identity is claimed. `group_invariants.csv` gives the exact SU(2) Casimirs and Dynkin indices as a separate representation check.

**Not computed:** general-interaction threshold Wilson matching, EFT1 higher-dimensional Wilson RGEs, and final Weinberg C5. For B(-1)/C(-1), mixed scalar masses require diagonalisation.

| Model | Extra quartics | Extra bilinears | Extra Yukawa | Status | Report |
|---|---:|---:|---|---|---|
| A_+0 | 5 | 1 | y1Conj | SymbolicRGExportedMatchingUnverified | [Detail](models/A_p0.md) |
| A_-2 | 5 | 1 | y2Conj | SymbolicRGExportedMatchingUnverified | [Detail](models/A_m2.md) |
| A_-4 | 0 | 0 | — | SymbolicRGExportedMatchingUnverified | [Detail](models/A_m4.md) |
| B_+1 | 1 | 0 | — | SymbolicRGExportedMatchingUnverified | [Detail](models/B_p1.md) |
| B_-1 | 7 | 1 | — | SymbolicRGExportedMassBasisBlocked | [Detail](models/B_m1.md) |
| B_-3 | 1 | 0 | — | SymbolicRGExportedMatchingUnverified | [Detail](models/B_m3.md) |
| C_+1 | 1 | 0 | — | SymbolicRGExportedMatchingUnverified | [Detail](models/C_p1.md) |
| C_-1 | 7 | 1 | — | SymbolicRGExportedMassBasisBlocked | [Detail](models/C_m1.md) |
| C_-3 | 1 | 0 | — | SymbolicRGExportedMatchingUnverified | [Detail](models/C_m3.md) |
| D_+0 | 5 | 1 | y1Conj | SymbolicRGExportedMatchingUnverified | [Detail](models/D_p0.md) |
| D_+2 | 0 | 0 | — | SymbolicRGExportedMatchingUnverified | [Detail](models/D_p2.md) |
| D_-2 | 5 | 1 | y2Conj | SymbolicRGExportedMatchingUnverified | [Detail](models/D_m2.md) |
| E_+0 | 6 | 1 | y1Conj | SymbolicRGExportedMatchingUnverified | [Detail](models/E_p0.md) |
| E_+2 | 0 | 0 | — | SymbolicRGExportedMatchingUnverified | [Detail](models/E_p2.md) |
| E_-2 | 6 | 1 | y2Conj | SymbolicRGExportedMatchingUnverified | [Detail](models/E_m2.md) |
| E_-4 | 0 | 0 | — | SymbolicRGExportedMatchingUnverified | [Detail](models/E_m4.md) |

See also [exact group invariants](group_invariants.csv) and [machine-readable status](summary.json).
