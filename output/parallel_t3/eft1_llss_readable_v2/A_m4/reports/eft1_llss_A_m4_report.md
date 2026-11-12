# EFT1 LLSS mixing isolation: A(-4) donor

Completed exploratory models: **0/16**.

## Method and percentage definitions

All three cases use the same UV-to-fermion matching, renormalisable EFT1
background, simultaneous scalar-threshold matching, and SMEFT evolution.
Only the dimension-five Wilson evolution is changed: frozen (all Wilson beta
functions zero), mixing off (only the quartic LLSS-to-LLHH source removed),
or running (complete candidate Wilson beta).

$\delta_{\rm Running}=100(O_{\rm Running}-O_{\rm Frozen})/O_{\rm Frozen}$.

$\delta_{\rm MixingOff}=100(O_{\rm MixingOff}-O_{\rm Frozen})/O_{\rm Frozen}$.

Both are signed percentages with a **common frozen denominator**.
The difference between these columns isolates the observable response to
restoring the LLSS-to-Weinberg quartic mixing term.

## Numerical comparison

### Mass-squared differences

| Model | $\Delta m^2_{21}$ Running | $\Delta m^2_{21}$ Mixing off | $\Delta m^2_{3\ell}$ Running | $\Delta m^2_{3\ell}$ Mixing off |
| --- | ---: | ---: | ---: | ---: |
| A(-4) | -- | -- | -- | -- |
| A(-2) | -- | -- | -- | -- |
| A(+0) | -- | -- | -- | -- |
| B(-3) | -- | -- | -- | -- |
| B(-1) | -- | -- | -- | -- |
| B(+1) | -- | -- | -- | -- |
| C(-3) | -- | -- | -- | -- |
| C(-1) | -- | -- | -- | -- |
| C(+1) | -- | -- | -- | -- |
| D(-2) | -- | -- | -- | -- |
| D(+0) | -- | -- | -- | -- |
| D(+2) | -- | -- | -- | -- |
| E(-4) | -- | -- | -- | -- |
| E(-2) | -- | -- | -- | -- |
| E(+0) | -- | -- | -- | -- |
| E(+2) | -- | -- | -- | -- |

### Mixing-angle observables

| Model | $\sin^2\theta_{12}$ Running | $\sin^2\theta_{12}$ Mixing off | $\sin^2\theta_{13}$ Running | $\sin^2\theta_{13}$ Mixing off | $\sin^2\theta_{23}$ Running | $\sin^2\theta_{23}$ Mixing off |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| A(-4) | -- | -- | -- | -- | -- | -- |
| A(-2) | -- | -- | -- | -- | -- | -- |
| A(+0) | -- | -- | -- | -- | -- | -- |
| B(-3) | -- | -- | -- | -- | -- | -- |
| B(-1) | -- | -- | -- | -- | -- | -- |
| B(+1) | -- | -- | -- | -- | -- | -- |
| C(-3) | -- | -- | -- | -- | -- | -- |
| C(-1) | -- | -- | -- | -- | -- | -- |
| C(+1) | -- | -- | -- | -- | -- | -- |
| D(-2) | -- | -- | -- | -- | -- | -- |
| D(+0) | -- | -- | -- | -- | -- | -- |
| D(+2) | -- | -- | -- | -- | -- | -- |
| E(-4) | -- | -- | -- | -- | -- | -- |
| E(-2) | -- | -- | -- | -- | -- | -- |
| E(+0) | -- | -- | -- | -- | -- | -- |
| E(+2) | -- | -- | -- | -- | -- | -- |

## Observable-level diagnostics

| Observable | Median abs(Running) (%) | Median abs(Mixing off) (%) |
| --- | ---: | ---: |
| $\Delta m^2_{21}$ | -- | -- |
| $\Delta m^2_{3\ell}$ | -- | -- |
| $\sin^2\theta_{12}$ | -- | -- |
| $\sin^2\theta_{13}$ | -- | -- |
| $\sin^2\theta_{23}$ | -- | -- |

## Figures

- `plots/all16_heatmap_percent.png` (and `.pdf`): 16 x 10 comparison.
- `plots/all16_masssq_heatmap_percent.png` (and `.pdf`).
- `plots/all16_mixing_heatmap_percent.png` (and `.pdf`).

## Scope and limitations

These are conditional exploratory numerical results. The chosen candidate
Wilson anomalous dimension, source-to-Weyl normalisation, scalar basis
maps and independent hard one-loop matching are not fully certified.
The `Frozen` case freezes the Wilson RGEs only: the renormalisable
EFT1 scalar evolution is active in every comparison. The resulting
plots are not numerically interchangeable with the original EFT1-off
figure that also freezes the scalar background.

Large mass-squared changes can coexist with very small mixing-angle
changes when the flavour matrix is approximately rescaled.

A conclusion that LLSS dominates is justified only when the Running
column is much larger than the Mixing-off column for the relevant
observables; exceptions must be reported rather than discarded.

Blocked or missing models: A_-4, A_-2, A_+0, B_-3, B_-1, B_+1, C_-3, C_-1, C_+1, D_-2, D_+0, D_+2, E_-4, E_-2, E_+0, E_+2.
