# Qualitative SU(2) representation and running effects: 16 general T3 models

## Scope and provenance

Read-only analysis of saved data (no matching, RG integration, or refits).
Pinned reference main commit: `6ebe9937a72880c6f180bcd026261779800ac71f`.
Run provenance field: `thesis16_merged_v1`; this is **not** a replacement for the final fit-status manifest.
Final fit-status source: `continue_general16_v1/merged_all16.json` (16 accepted exploratory fits).
14 ordinary model results: `general_seq_05`; B(-1)/C(-1): **newer** `general_seq_bc_normalized_01`.
The independent B/C mass-log audit is `log_bc_normalized_01`; local coupling sensitivity is `zee_general_all16_merged_v1`.

## Switch and EFT semantics

`U/E/W/S`: UV beta (1e7→1e5 GeV); EFT1 scalar beta (1e5→1e3 GeV);
reconstructed intermediate LLSS→Weinberg leading log; SMEFT/SM+Weinberg beta (1e3→1e2 GeV).
All fermions integrated at 1e5 GeV, then both scalar multiplets at 1e3 GeV.
`1111` is fully enabled; freeze one stage using `0111/1011/1101/1110` respectively.
The stage deltas are **conditional ablations at fixed fitted couplings**, not separately refitted models.
They need not add or factorize: turning on U or E modifies downstream matching inputs.
The W-on C5 is the full-general **one-scale anchor**; W-off subtracts the reconstructed LLSS log while retaining the bubble and bookkeeping residual.
W-off is **not independent physical hard matching**; its large suppression must not be claimed as an independently certified running prediction.
The GammaCC-to-Weyl factor is assumed 1; the extra HHSS basis and absolute normalization remain provisional.

## Verified qualitative patterns in the saved cases

- UV enabled: heaviest mass changes **+1.08% to +13.05%**, positive in all 16 models; family E is largest.
- SMEFT enabled: changes **-6.47% to -6.30%**, negative and almost universal.
- Intermediate EFT1: suppresses the heaviest mass for every B and C assignment;
  enhances every E assignment; A and D show mixed or small shifts.
- W-off conditional comparison: the full/W-off heaviest-mass ratio ranges from 5.69 to 12.18.
  This is a **bookkeeping sensitivity** to the LLSS log, not independent hard-matched physics.

Family-dependent EFT1 shifts to $m_3$ (other switches enabled):

| Family | SU(2) dimensions $(F,S_1,S_2)$ | EFT1 effect |
|---|---|---:|
| A | 2, 1, 3 | -2.91% to +0.48% |
| B | 1, 2, 2 | -4.24% to -1.99% |
| C | 3, 2, 2 | -4.98% to -3.54% |
| D | 2, 3, 1 | -3.07% to +0.38% |
| E | 2, 3, 3 | +2.47% to +5.52% |

The integer alpha in model labels parametrizes hypercharge assignments.
Pairing models within each family probes alpha dependence at fixed scalar SU(2) dimensions;
however family differences also include distinct fermion sectors and channels, so this is not a controlled single-representation derivative.

## Interpretation limitations

All 16 accepted UV Yukawa configurations were calibrated against the same oscillation targets.
Therefore nearly identical `1111` mass-squared splittings/mixings are *not* evidence that different representations physically predict the same observables,
nor can best-fit agreement alone establish SU(2) or hypercharge effects.
The meaningful representation comparison here is **the response to frozen stages under the same campaign conventions**.
Some `W=0` predictions can have strongly distorted mass spectra or mixing angles; interpret near-zero references using absolute differences.
B(-1)/C(-1) have channel-normalised conditional inputs and should be compared with the other fourteen only as a *qualified* combined exploratory sample.
The independent B/C audit checks the mass-log coefficient, not the complete finite or hard matching.
The compensated-coupling sensitivity CSV records local saved finite differences; it is **not an additional running stage**.

## Reproduction checks

- Selected case count: 256 (16 models × 16 U/E/W/S combinations).
- Max relative oscillation-observable difference between `1111` and final manifest: 3.7e-15.
- Sensitivity points checked: 528 (all reported successful).
- No saved physics outputs written.
