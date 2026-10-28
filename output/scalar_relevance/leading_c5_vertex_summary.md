# T3 scalar relevance — Batch 2: Wick-connection audit

- Models checked: **16**
- New Higgs-bearing quartic slots: **24**
- Wick-eligible and present in the saved C5: **4** (two each in B_-1/C_-1)
- Excluded from the specific single-vertex graph: **20** (10 H Hdag, 10 scalar-orientation mismatch)
- Four-BSM/no-H quartics: **26**; no direct single-quartic insertion
- Additional quadratic mass parameters: **8**, not excluded

## Essential physical derivation

Original lepton Yukawas source `S1` and `S2dag`. With diagonal complex scalar propagators, the Higgs-bearing quartic must supply `S1dag` or `S2`. A single quartic feeding the matched holomorphic Weinberg sector must supply `Hdag^2`, using the conjugate orientation of the baseline `H^2 S1 S2dag` coupling. Therefore the only scalar field pairs admitted at the quartic end are `S1dag^2` (`y1y1`), `S1dag S2` (`y1y2`) and `S2^2` (`y2y2`). Both hermitian orientations of every quartic were tested.

## Model breakdown

| Model | New 2-H quartics | Wick-eligible | H Hdag mismatch | Scalar mismatch |
|---|---:|---:|---:|---:|
| A_+0 | 2 | 0 | 1 | 1 |
| A_-2 | 2 | 0 | 1 | 1 |
| A_-4 | 0 | 0 | 0 | 0 |
| B_+1 | 1 | 0 | 0 | 1 |
| B_-1 | 4 | 2 | 2 | 0 |
| B_-3 | 1 | 0 | 0 | 1 |
| C_+1 | 1 | 0 | 0 | 1 |
| C_-1 | 4 | 2 | 2 | 0 |
| C_-3 | 1 | 0 | 0 | 1 |
| D_+0 | 2 | 0 | 1 | 1 |
| D_+2 | 0 | 0 | 0 | 0 |
| D_-2 | 2 | 0 | 1 | 1 |
| E_+0 | 2 | 0 | 1 | 1 |
| E_+2 | 0 | 0 | 0 | 0 |
| E_-2 | 2 | 0 | 1 | 1 |
| E_-4 | 0 | 0 | 0 | 0 |

## Direct new vertices in the saved one-loop C5

- `B_-1` — `lambdaT3ExtraQ020020I1`: `Hdag^2 S2^2` (y2_y2)
- `B_-1` — `lambdaT3ExtraQ020200I1`: `Hdag^2 S1dag^2` (y1_y1)
- `C_-1` — `lambdaT3ExtraQ020020I1`: `Hdag^2 S2^2` (y2_y2)
- `C_-1` — `lambdaT3ExtraQ020200I1`: `Hdag^2 S1dag^2` (y1_y1)

## Limitations

- This is a necessary contraction test in the DIAGONAL scalar basis, NOT a proof of irrelevance with nonzero extra quadratic mixing.
- Fermion contractions and group tensors can make a Wick-eligible vertex vanish; the four observed terms are supported separately by saved symbolic Matchete matching.
- Non-observed 20 are excluded only from the defined direct one-quartic-vertex amplitude, not from the general UV scalar potential.
- Higher-loop, RG, v^2/M^2, vacuum stability, and general two-threshold Wilson matching are not addressed here.

- Local Batch 1 source hash verification: **yes**
