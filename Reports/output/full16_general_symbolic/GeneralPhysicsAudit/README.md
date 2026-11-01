# General 16-model symbolic physics audit

Structurally verified: **16/16**.
Complete 16-model inventory totals gate: **PASS**.

## What has actually been checked

- Current-source SHA-256 provenance and existing Matchete detailed-check statuses
- Each independent additional scalar monomial's U(1) hypercharge and Z2 selection rules
- Quartic multiplicities cross-checked among the manifest, Matchete interaction names and RGBeta registry
- Additional scalar mass and Yukawa operator inventory (including six exceptional Yukawas)
- Required UV and EFT1 beta-function *keys*, including all general quartics and massive operators
- Independent exact SU(2) Casimirs and Dynkin indices compared with the symbolic-study report

## What remains NOT CERTIFIED

1. Physical action-level Matchete↔RGBeta coupling dictionary: tensor basis, phase, complex conjugation, signs, factorials and H.c.
2. Equality of each group-factor/beta coefficient after applying that physical dictionary
3. Scalar mass-basis matching for B(-1), C(-1) and complete extended sequential EFT/C5 matching

The existing tensor-Ward and normalization-bridge tests are partial evidence; this audit never upgrades them to physical equivalence.
Inspect `audit.json` -> `saved_bridge_report_checks` for the status of existing 38-Yukawa witness files.

## Per-model structural gate

| Model | Quartics | Mass terms | Extra Yukawa | UV betas | EFT1 betas | Status |
|---|---:|---:|---|---:|---:|---|
| A_+0 | 5 | 1 | y1Conj | 27 | 23 | VerifiedStructural |
| A_-2 | 5 | 1 | y2Conj | 27 | 23 | VerifiedStructural |
| A_-4 | 0 | 0 | — | 20 | 17 | VerifiedStructural |
| B_+1 | 1 | 0 | — | 22 | 19 | VerifiedStructural |
| B_-1 | 7 | 1 | — | 29 | 26 | VerifiedStructural |
| B_-3 | 1 | 0 | — | 22 | 19 | VerifiedStructural |
| C_+1 | 1 | 0 | — | 22 | 19 | VerifiedStructural |
| C_-1 | 7 | 1 | — | 29 | 26 | VerifiedStructural |
| C_-3 | 1 | 0 | — | 22 | 19 | VerifiedStructural |
| D_+0 | 5 | 1 | y1Conj | 27 | 23 | VerifiedStructural |
| D_+2 | 0 | 0 | — | 20 | 17 | VerifiedStructural |
| D_-2 | 5 | 1 | y2Conj | 27 | 23 | VerifiedStructural |
| E_+0 | 6 | 1 | y1Conj | 32 | 28 | VerifiedStructural |
| E_+2 | 0 | 0 | — | 24 | 21 | VerifiedStructural |
| E_-2 | 6 | 1 | y2Conj | 32 | 28 | VerifiedStructural |
| E_-4 | 0 | 0 | — | 24 | 21 | VerifiedStructural |

See `audit.json`, `quartic_invariants.csv`, and `group_invariants.csv` for machine-readable evidence.
