# Step 4: charged-lepton mass-basis validation

## Field convention

The Ma runner stores $Y_e$ and $h$ with lepton-doublet flavour on the matrix
columns. Write the singular-value decomposition as

$$
Y_e=U_R D_e V_L^\dagger,
\qquad
D_e=\operatorname{diag}(y_e,y_\mu,y_\tau)>0.
$$

The corresponding field definitions are

$$
e_{R,\mathrm{old}}=U_R e_{R,\mathrm{mass}},
\qquad
L_{\mathrm{old}}=V_L L_{\mathrm{mass}}.
$$

They give

$$
U_R^\dagger Y_eV_L=D_e.
$$

Because the lepton-doublet index is the column of the scotogenic Yukawa
matrix, its mass-basis form is

$$
h_{\mathrm{mass}}=hV_L.
$$

For the symmetric Weinberg coefficient,

$$
L_{\mathrm{old}}^TC_5L_{\mathrm{old}}
=L_{\mathrm{mass}}^T(V_L^TC_5V_L)L_{\mathrm{mass}},
$$

so

$$
C_{5,\mathrm{mass}}=V_L^TC_5V_L.
$$

The transpose is required; replacing it by $V_L^\dagger$ would be incorrect
for a Majorana operator.

## Implementation

`MaFullRunningComparison.py` now:

1. computes the full complex SVD of $Y_e$ at the matching scale;
2. sorts the singular values in ascending electron--muon--tau order;
3. applies a deterministic column-phase convention to $V_L$;
4. verifies $U_R^\dagger Y_eV_L=D_e$ numerically;
5. rotates $C_5$ with $V_L^TC_5V_L$; and
6. records the rotation matrix, mass-basis Yukawas, and relative SVD residual
   in the output JSON.

The low-energy SMEFT runner then receives a diagonal, real, non-negative
charged-lepton Yukawa vector.

## Validation checks

A general complex $Y_e$ with known non-degenerate singular values was used.
The following checks pass:

| Check | Result |
|---|---|
| Recovered $y_e,y_\mu,y_\tau$ in ascending order | Pass |
| $U_R^\dagger Y_eV_L=D_e$ | Relative residual $4.8\times10^{-16}$ |
| Rotated matching equals direct matching with $h\to hV_L$ | Pass |
| $C_5=C_5^T$ after rotation | Pass |
| $V_L^\dagger V_L=1$ | Pass |

For the supplied numerical benchmark:

| Quantity | Frozen input | After Ma UV running |
|---|---:|---:|
| SVD relative residual | $0$ | $6.55\times10^{-16}$ |
| Off-diagonal norm of $V_L$ | $0$ | $8.35\times10^{-6}$ |
| $y_e$ | $2.9000\times10^{-6}$ | $2.81294\times10^{-6}$ |
| $y_\mu$ | $6.1000\times10^{-4}$ | $5.91686\times10^{-4}$ |
| $y_\tau$ | $1.0200\times10^{-2}$ | $9.89368\times10^{-3}$ |

The RGE comparison remains

$$
\Delta_{\rm UV}=2.29545\%,\qquad
\Delta_{\rm EFT}=4.29400\%,\qquad
\Delta_{\rm full}=5.58232\%.
$$

## Phase convention

An SVD permits a simultaneous phase change of corresponding columns of
$U_R$ and $V_L$. This phase is physically unobservable but can make two
numerically equivalent $C_5$ matrices appear different entry by entry. The
implementation fixes each $V_L$ column by making its largest component real
and positive, and applies the same phase to $U_R$. This preserves the positive
diagonal $D_e$ and stabilizes comparisons between nearby RGE trajectories.

## Result and limitation

The charged-lepton basis transformation, flavour-index direction, transpose,
ordering, and phase handling are correct for non-degenerate charged-lepton
singular values. The physical charged-lepton spectrum is strongly
non-degenerate, so the remaining mathematical ambiguity associated with an
exactly degenerate SVD is not relevant to the benchmark.
