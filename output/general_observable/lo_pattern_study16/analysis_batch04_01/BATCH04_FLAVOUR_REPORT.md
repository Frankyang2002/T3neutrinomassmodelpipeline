# Chat A — Batch 4: why the 16 neutrino-mass hierarchies differ

**Source:** saved Batch-1/2/3 conditional UV results; no new Feynman diagrams, UV runs or EFT matching. **Physical source normalization remains 0/16 certified.**

## Analytical separation: flavour geometry vs heavy-state cancellation

The exact source-conditional heavy-state decomposition in Batch 3 is
\[M_\nu=\sum_g H_g,\qquad H_g=\frac{1}{16\pi^2}\sum_q\sum_{k\in g}K_q(m_k^2)\,\mathrm{sym}\big[Z_{s_q}(M_F^\dagger P_k)Y_{t_q}^{T}\big].\]
Here `g` labels a **distinct heavy singular-mass eigenspace**. Degenerate individual eigenvectors cannot be assigned a unique contribution. The complete projector is invariant.
For the complex SVD `Mν=U diag(m1,m2,m3) V†`, take singular vectors `u_i`,`v_i` and set
\[a_{gi}=u_i^\dagger H_g v_i,\qquad m_i=\sum_g a_{gi}.\]
Define `G = Σ_g ||H_g||F`, `A_i = Σ_g |a_gi| / G` (geometric projection fraction) and `C_i = m_i / Σ_g |a_gi|` (coherence fraction, 0 = cancellation, 1 = aligned). Then, exactly (up to roundoff),
\[\boxed{\frac{m_i}{G}=A_i C_i}.\]
This disambiguates a small eigenmass caused by **weak flavour projection** (`A_i` small) from one caused by **phase cancellation between heavy eigenspaces** (`C_i` small). No notion of adding group Takagi masses is used.

## All-16 measured mechanism diagnostics

| Model | fermion | heavy mass groups | m1/m3 | A1 (geometry) | C1 (coherence) | lightest-mode diagnosis |
|---|---|---:|---:|---:|---:|---|
| A_-4 | Dirac | 3 | 0.010882 | 0.1025 | 0.057116 | mixed_or_neither_dominant |
| A_-2 | Dirac | 3 | 0.14857 | 0.16843 | 0.48065 | mode_geometry_dominant |
| A_+0 | Dirac | 3 | 0.14875 | 0.16868 | 0.48007 | mode_geometry_dominant |
| B_-3 | Dirac | 3 | 0.14935 | 0.16888 | 0.48012 | mode_geometry_dominant |
| B_-1 | Majorana | 3 | 0.089659 | 0.15701 | 0.20854 | mixed_or_neither_dominant |
| B_+1 | Dirac | 3 | 0.1501 | 0.15949 | 0.50522 | mode_geometry_dominant |
| C_-3 | Dirac | 3 | 0.48958 | 0.26878 | 1 | mode_geometry_dominant |
| C_-1 | Majorana | 3 | 0.60203 | 0.5691 | 0.28578 | mixed_or_neither_dominant |
| C_+1 | Dirac | 3 | 0.0076484 | 0.3192 | 0.013122 | heavy_eigenspace_cancellation_dominant |
| D_-2 | Dirac | 3 | 0.1488 | 0.1693 | 0.47851 | mode_geometry_dominant |
| D_+0 | Dirac | 3 | 0.1486 | 0.16846 | 0.48058 | mode_geometry_dominant |
| D_+2 | Dirac | 3 | 0.01089 | 0.10255 | 0.057118 | mixed_or_neither_dominant |
| E_-4 | Dirac | 3 | 0.086244 | 0.15726 | 0.214 | mixed_or_neither_dominant |
| E_-2 | Dirac | 3 | 0.111 | 0.038958 | 1 | mode_geometry_dominant |
| E_+0 | Dirac | 3 | 0.10922 | 0.038752 | 1 | mode_geometry_dominant |
| E_+2 | Dirac | 3 | 0.08609 | 0.15766 | 0.21322 | mixed_or_neither_dominant |

### C_-1: Majorana, m1/m3=0.602027
- Frozen ||y1||F × ||y2||F = 1.937133e-07; |lambdaT3| = 0.1.
- The full lightest mass is 0.03749726 eV; the corresponding group-projection sum of absolute magnitudes is 0.1312111 eV.
- A1=0.5691; C1=0.28577815508140025; Both projection geometry and inter-eigenspace interference contribute; neither is separated by more than a factor of two in this diagnostic.
- Exact mode amplitudes and complex phases are listed in `batch04_mode_projections.csv`; deleting a group is a **diagnostic counterfactual**, not a new UV calculation.

### C_+1: Dirac, m1/m3=0.00764841
- Frozen ||y1||F × ||y2||F = 6.14575e-08; |lambdaT3| = 0.1.
- The full lightest mass is 0.0003654878 eV; the corresponding group-projection sum of absolute magnitudes is 0.02785293 eV.
- A1=0.319197; C1=0.013122059990926294; Large projections from distinct heavy-mass eigenspaces partially cancel in the lightest singular direction.
- Exact mode amplitudes and complex phases are listed in `batch04_mode_projections.csv`; deleting a group is a **diagnostic counterfactual**, not a new UV calculation.

### C_-3: Dirac, m1/m3=0.48958
- Frozen ||y1||F × ||y2||F = 6.3062e-08; |lambdaT3| = 0.1.
- The full lightest mass is 0.02694687 eV; the corresponding group-projection sum of absolute magnitudes is 0.02694687 eV.
- A1=0.268784; C1=1.0000000000050753; Its smallest singular direction has weak overlap with the heavy-group matrices; geometric suppression dominates the diagnostic split.
- Exact mode amplitudes and complex phases are listed in `batch04_mode_projections.csv`; deleting a group is a **diagnostic counterfactual**, not a new UV calculation.

## Additional cross-family observations

- **C_-1 versus C_+1 is not a pure Yukawa-controlled comparison:** the frozen fitted Yukawas differ **and the fermion construction differs** (neutral Majorana C_-1 versus Dirac C_+1). The group decomposition identifies numerical mechanisms, but does not isolate their causal origin with all other inputs fixed.
- **A/D matched cases:** compare their mode fractions and group interference, not just their close `m3`; near equality of selected frozen matrices is not an all-parameter SU(2) identity.
- **E:** Batch 3 found two perfectly aligned scalar-charge channels. The heavy-eigenspace coherence here is a *different* quantity; aligned charge channels do not preclude cancellation or weak projections among heavy generations.
- The determinant relation `|det Mν|=m1 m2 m3` provides a complementary flavour-invariant hierarchy measure; a small determinant need not mean a small loop amplitude.

## Source and interpretation limits

- Calculations are on the complex **source-conditional** matrix in eV, not physical-source-certified amplitudes.
- The mode diagnostic threshold compares the two exact suppression factors: “dominant” means one is at least twice smaller; no model-selection significance is implied.
- A near-degenerate *light-neutrino* singular pair makes individual singular-vector projections orientation-sensitive; the summed singular masses and group matrices remain unambiguous.
- The direct B1/B2 spinor, symmetry, portal and CG physical normalization is still not certified.

## Per-model measured interpretation

- **A_-4**: `m1/m3=0.0108825`, `A1=0.102495`, `C1=0.057116`. Both projection geometry and inter-eigenspace interference contribute; neither is separated by more than a factor of two in this diagnostic.
- **A_-2**: `m1/m3=0.148573`, `A1=0.168426`, `C1=0.480645`. Its smallest singular direction has weak overlap with the heavy-group matrices; geometric suppression dominates the diagnostic split.
- **A_+0**: `m1/m3=0.148751`, `A1=0.168679`, `C1=0.480074`. Its smallest singular direction has weak overlap with the heavy-group matrices; geometric suppression dominates the diagnostic split.
- **B_-3**: `m1/m3=0.149346`, `A1=0.168877`, `C1=0.480121`. Its smallest singular direction has weak overlap with the heavy-group matrices; geometric suppression dominates the diagnostic split.
- **B_-1**: `m1/m3=0.0896585`, `A1=0.157008`, `C1=0.208537`. Both projection geometry and inter-eigenspace interference contribute; neither is separated by more than a factor of two in this diagnostic.
- **B_+1**: `m1/m3=0.150097`, `A1=0.159493`, `C1=0.505219`. Its smallest singular direction has weak overlap with the heavy-group matrices; geometric suppression dominates the diagnostic split.
- **C_-3**: `m1/m3=0.48958`, `A1=0.268784`, `C1=1`. Its smallest singular direction has weak overlap with the heavy-group matrices; geometric suppression dominates the diagnostic split.
- **C_-1**: `m1/m3=0.602027`, `A1=0.5691`, `C1=0.285778`. Both projection geometry and inter-eigenspace interference contribute; neither is separated by more than a factor of two in this diagnostic.
- **C_+1**: `m1/m3=0.00764841`, `A1=0.319197`, `C1=0.0131221`. Large projections from distinct heavy-mass eigenspaces partially cancel in the lightest singular direction.
- **D_-2**: `m1/m3=0.1488`, `A1=0.169299`, `C1=0.478505`. Its smallest singular direction has weak overlap with the heavy-group matrices; geometric suppression dominates the diagnostic split.
- **D_+0**: `m1/m3=0.148603`, `A1=0.168458`, `C1=0.480577`. Its smallest singular direction has weak overlap with the heavy-group matrices; geometric suppression dominates the diagnostic split.
- **D_+2**: `m1/m3=0.0108897`, `A1=0.10255`, `C1=0.0571179`. Both projection geometry and inter-eigenspace interference contribute; neither is separated by more than a factor of two in this diagnostic.
- **E_-4**: `m1/m3=0.086244`, `A1=0.157262`, `C1=0.214`. Both projection geometry and inter-eigenspace interference contribute; neither is separated by more than a factor of two in this diagnostic.
- **E_-2**: `m1/m3=0.111002`, `A1=0.0389577`, `C1=1`. Its smallest singular direction has weak overlap with the heavy-group matrices; geometric suppression dominates the diagnostic split.
- **E_+0**: `m1/m3=0.109221`, `A1=0.0387518`, `C1=1`. Its smallest singular direction has weak overlap with the heavy-group matrices; geometric suppression dominates the diagnostic split.
- **E_+2**: `m1/m3=0.0860901`, `A1=0.157661`, `C1=0.213222`. Both projection geometry and inter-eigenspace interference contribute; neither is separated by more than a factor of two in this diagnostic.

The exact `m1/G=A1*C1` identity is the analytic explanation of *how* a particular frozen texture suppresses the lightest mode. It does not, by itself, prove which lagrangian coupling made that texture arise.

