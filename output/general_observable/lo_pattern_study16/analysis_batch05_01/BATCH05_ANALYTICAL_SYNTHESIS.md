# Chat A — Batch 5: analytical origin of the 16 neutrino-mass patterns

**Scope:** frozen, source-conditional T3 Weyl amplitudes. All B1/B2 physical normalization is unresolved. No EFT matching, no refit, no production UV-solver calls.

## Exact component and flavour decomposition

The saved candidate Weyl vertices are factorized into **frozen Yukawas** and model-specific component tensors:

\[Z_{iIs}=\sum_r y_{1,ir}^{*} C^{(1)}_{rIs},\qquad Y_{iIt}=\sum_r y_{2,ir}^{*} C^{(2)}_{rIt}.\]

The complete source-conditional mass is rebuilt as the sum of charge and distinct heavy singular-mass eigenspace terms:

\[M_\nu=\frac{1}{16\pi^2}\sum_{q,g}\sum_{k\in g}K_q(m_k^2)\operatorname{sym}\bigl[Z_{s_q}(M_F^\dagger P_k)Y_{t_q}^{T}\bigr],\]

\[K_q(m^2)=\delta_q\frac{B_0(m^2,\lambda_{q+})-B_0(m^2,\lambda_{q-})}{\lambda_{q+}-\lambda_{q-}}.\]

The scalar formula is exact for the restricted isolated 2×2 S1/S2 charge blocks already verified by Batch 3. It is not a claim about missing B1/B2 phases or symmetry factors.

## Frozen physical input versus controlled *counterfactual* inputs

For each model we leave its SU(2), hypercharge, heavy Weyl matrix, scalar masses and candidate component tensors fixed, but substitute exactly the same 3×3 complex Yukawa matrices from four reference frozen benchmarks: `A_-2`, `C_-1`, `C_+1`, `E_-2`. Each replacement is performed at the **vertex level**, with every conjugation preserved. None is a refitted benchmark. A separate common-|lambda| comparison rescales E scalar mixing from 0.05 to 0.10; this is only a controlled source-conditional comparison.

## Observed mass pattern and complex-mode cancellation

| Model | Fermion | m1/m3 | A1 | C1 | cancelled projected strength | strongest opposing pair |
|---|---|---:|---:|---:|---:|---|
| A_-4 | Dirac | 0.010882 | 0.1025 | 0.057116 | 94.288% | 2–3 (cos=-1) |
| A_-2 | Dirac | 0.14857 | 0.16843 | 0.48065 | 51.935% | 1–3 (cos=-1) |
| A_+0 | Dirac | 0.14875 | 0.16868 | 0.48007 | 51.993% | 1–3 (cos=-1) |
| B_-3 | Dirac | 0.14935 | 0.16888 | 0.48012 | 51.988% | 1–3 (cos=-1) |
| B_-1 | Majorana | 0.089659 | 0.15701 | 0.20854 | 79.146% | 1–2 (cos=-1) |
| B_+1 | Dirac | 0.1501 | 0.15949 | 0.50522 | 49.478% | 1–3 (cos=-1) |
| C_-3 | Dirac | 0.48958 | 0.26878 | 1 | 0.000% | 2–3 (cos=1) |
| C_-1 | Majorana | 0.60203 | 0.5691 | 0.28578 | 71.422% | 1–3 (cos=-1) |
| C_+1 | Dirac | 0.0076484 | 0.3192 | 0.013122 | 98.688% | 2–3 (cos=-1) |
| D_-2 | Dirac | 0.1488 | 0.1693 | 0.47851 | 52.149% | 1–3 (cos=-1) |
| D_+0 | Dirac | 0.1486 | 0.16846 | 0.48058 | 51.942% | 1–3 (cos=-1) |
| D_+2 | Dirac | 0.01089 | 0.10255 | 0.057118 | 94.288% | 2–3 (cos=-1) |
| E_-4 | Dirac | 0.086244 | 0.15726 | 0.214 | 78.600% | 1–2 (cos=-1) |
| E_-2 | Dirac | 0.111 | 0.038958 | 1 | 0.000% | 2–3 (cos=1) |
| E_+0 | Dirac | 0.10922 | 0.038752 | 1 | 0.000% | 2–3 (cos=1) |
| E_+2 | Dirac | 0.08609 | 0.15766 | 0.21322 | 78.678% | 1–2 (cos=-1) |

## Controlled reference textures, common scalar mixing parameter

For this comparison, each reference provides **both** `y1` and `y2`. The ratio m1/m3 is invariant under an overall nonzero matrix scaling but not under a change of texture or Majorana/Dirac construction. All entries use |lambdaT3|=0.10, rescaling E by two.

| Model | own frozen m1/m3 | A_-2 texture | C_-1 texture | C_+1 texture | E_-2 texture |
|---|---:|---:|---:|---:|---:|
| A_-4 | 0.010882 | 0.14857 | 0.60203 | 0.0076484 | 0.111 |
| A_-2 | 0.14857 | 0.14857 | 0.60203 | 0.0076484 | 0.111 |
| A_+0 | 0.14875 | 0.14857 | 0.60203 | 0.0076484 | 0.111 |
| B_-3 | 0.14935 | 0.14857 | 0.60203 | 0.0076484 | 0.111 |
| B_-1 | 0.089659 | 0.14857 | 0.60203 | 0.0076484 | 0.111 |
| B_+1 | 0.1501 | 0.14857 | 0.60203 | 0.0076484 | 0.111 |
| C_-3 | 0.48958 | 0.14857 | 0.60203 | 0.0076484 | 0.111 |
| C_-1 | 0.60203 | 0.14857 | 0.60203 | 0.0076484 | 0.111 |
| C_+1 | 0.0076484 | 0.14857 | 0.60203 | 0.0076484 | 0.111 |
| D_-2 | 0.1488 | 0.14857 | 0.60203 | 0.0076484 | 0.111 |
| D_+0 | 0.1486 | 0.14857 | 0.60203 | 0.0076484 | 0.111 |
| D_+2 | 0.01089 | 0.14857 | 0.60203 | 0.0076484 | 0.111 |
| E_-4 | 0.086244 | 0.14857 | 0.60203 | 0.0076484 | 0.111 |
| E_-2 | 0.111 | 0.14857 | 0.60203 | 0.0076484 | 0.111 |
| E_+0 | 0.10922 | 0.14857 | 0.60203 | 0.0076484 | 0.111 |
| E_+2 | 0.08609 | 0.14857 | 0.60203 | 0.0076484 | 0.111 |

## C triplet family: controlled texture swaps and interpolation

C_-3 and C_+1 use charged Dirac heavy fermions; C_-1 has a neutral Majorana heavy fermion. The following Yukawa changes do **not** alter this fundamental difference.

| Frozen model | substitution | m1/m3 | A1 | C1 | m1 (eV) |
|---|---|---:|---:|---:|---:|
| C_-3 | frozen | 0.48958 | 0.26878 | 1 | 0.0269469 |
| C_-3 | C_-1_both | 0.60203 | 0.5691 | 0.28578 | 0.0374973 |
| C_-3 | C_+1_both | 0.0076484 | 0.3192 | 0.013122 | 0.000365488 |
| C_-3 | C_-1_y1_only | 0.081708 | 0.095102 | 0.51652 | 0.0116341 |
| C_-3 | C_-1_y2_only | 0.06682 | 0.13706 | 0.21277 | 0.00358078 |
| C_-3 | C_+1_y1_only | 0.094929 | 0.052237 | 1 | 0.0061553 |
| C_-3 | C_+1_y2_only | 0.57777 | 0.4662 | 0.29599 | 0.00895059 |
| C_-1 | frozen | 0.60203 | 0.5691 | 0.28578 | 0.0374973 |
| C_-1 | C_-1_both | 0.60203 | 0.5691 | 0.28578 | 0.0374973 |
| C_-1 | C_+1_both | 0.0076484 | 0.3192 | 0.013122 | 0.000365488 |
| C_-1 | C_-1_y1_only | 0.60203 | 0.5691 | 0.28578 | 0.0374973 |
| C_-1 | C_-1_y2_only | 0.60203 | 0.5691 | 0.28578 | 0.0374973 |
| C_-1 | C_+1_y1_only | 0.012319 | 0.24644 | 0.031275 | 0.0010263 |
| C_-1 | C_+1_y2_only | 0.20697 | 0.30958 | 0.45519 | 0.0249774 |
| C_+1 | frozen | 0.0076484 | 0.3192 | 0.013122 | 0.000365488 |
| C_+1 | C_-1_both | 0.60203 | 0.5691 | 0.28578 | 0.0374973 |
| C_+1 | C_+1_both | 0.0076484 | 0.3192 | 0.013122 | 0.000365488 |
| C_+1 | C_-1_y1_only | 0.20697 | 0.30958 | 0.45519 | 0.0249774 |
| C_+1 | C_-1_y2_only | 0.012319 | 0.24644 | 0.031275 | 0.0010263 |
| C_+1 | C_+1_y1_only | 0.0076484 | 0.3192 | 0.013122 | 0.000365488 |
| C_+1 | C_+1_y2_only | 0.0076484 | 0.3192 | 0.013122 | 0.000365488 |

### Measured C-family heavy-eigenspace phasors

Each row is the **complex** contribution projected on the full lightest-neutrino SVD mode. The sum is a real positive `m1` within numerical accuracy; the individual terms need not be positive.

| Frozen model | heavy group | mass [GeV] | Re(a_g1) [eV] | Im(a_g1) [eV] | angle [deg] |
|---|---:|---:|---:|---:|---:|
| C_-3 | 1 | 100000 | 0.021121056 | 0 | 0 |
| C_-3 | 2 | 110000 | 0.0020814536 | 0 | 0 |
| C_-3 | 3 | 120000 | 0.0037443648 | 0 | 0 |
| C_-1 | 1 | 100000 | 0.084354162 | 0 | 0 |
| C_-1 | 2 | 110000 | -0.0093686349 | 0 | 180 |
| C_-1 | 3 | 120000 | -0.03748827 | 0 | 180 |
| C_+1 | 1 | 100000 | 0.0048812181 | 0 | 0 |
| C_+1 | 2 | 110000 | -0.01374372 | 0 | 180 |
| C_+1 | 3 | 120000 | 0.0092279896 | 0 | 0 |

### C_-1 versus C_+1: algebraic suppression factors

\[m_1=G A_1 C_1,\quad G=\sum_g\|H_g\|_F,\quad A_1=\frac{\sum_g|a_{g1}|}{G},\quad C_1=\frac{|\sum_g a_{g1}|}{\sum_g|a_{g1}|}.\]

The frozen ratio m1(C_-1)/m1(C_+1) = 102.5951 factorizes as [2.64223] (overall group strength) × [1.78291] (mode projection) × [21.7785] (heavy-eigenspace coherence).
This is an *exact numerical factorization*, not an independent theory prediction or a causal attribution to Majorana versus Dirac structure.

- C_+1 → C_-1 fitted texture: frozen m1/m3 = 0.007648411; endpoint = 0.602027; smallest **sampled** ratio 0.007648411 at t=0.0.
- C_-1 → C_+1 fitted texture: frozen m1/m3 = 0.602027; endpoint = 0.007648411; smallest **sampled** ratio 0.007648411 at t=1.0.

The CSV `batch05_C_texture_interpolation.csv` evaluates 11 intermediate Yukawa pairs in each direction. Its extrema and the mode projection/coherence diagnostics reveal whether the C_+1 cancellation persists after the C_-1 Yukawa texture is substituted. A small mass at an intermediate point is a **counterfactual** algebraic cancellation and is not a fitted or naturalness statement.

## Pairwise A/D, B, C, E patterns

To avoid mistaking a fit for a group-theory identity, both the original frozen mass-matrix proximity and the controlled-common-texture proximity are reported in `batch05_matched_pairs.csv`. A near-perfect alignment of the former alone does not establish representation equivalence. E has two aligned scalar-charge channels in the saved inputs; this is separate from the three heavy-eigenspace interference phases.

### Controlled matrix proportionality tests

A small residual below means both *complex* 3x3 mass matrices are related by a single fitted complex scalar for the indicated shared source Yukawa texture. It tests a conditional benchmark relation, not an operator-level SU(2) identity.

| Pair | frozen residual | common A_-2 | common C_-1 | common C_+1 | common E_-2 |
|---|---:|---:|---:|---:|---:|
| A_-4 / D_+2 | 2.81e-05 | 5.63e-16 | 1.84e-15 | 1.17e-15 | 1.67e-15 |
| A_-2 / D_+0 | 3.31e-05 | 1.87e-16 | 4.84e-16 | 5.6e-16 | 2.54e-16 |
| A_+0 / D_-2 | 5.4e-05 | 1.88e-16 | 8.65e-16 | 1.76e-16 | 2.61e-16 |
| B_-3 / B_+1 | 0.000838 | 1.37e-16 | 0 | 0 | 0 |
| C_-3 / C_+1 | 0.988 | 1.6e-16 | 3.19e-16 | 1.83e-16 | 5.57e-16 |
| E_-4 / E_+2 | 0.000161 | 2.73e-16 | 3.6e-16 | 2.22e-16 | 3.12e-16 |
| E_-2 / E_+0 | 0.00192 | 6.74e-15 | 1.43e-14 | 1.48e-14 | 4.06e-16 |

## Model-specific explanations

### A_-4: Dirac

- Source masses: (0.0004744829, 0.007582728, 0.04360064) eV. m1/m3=0.010882; candidate charge channels=1.
- Lightest mode: A1=0.1025, C1=0.057116; 94.288% of absolute projected strength cancels.
- Original source vertices reconstructed with relative errors 0.00e+00 and 6.49e-17.
- At common scalar mixing, C_+1 reference texture gives m1/m3=0.0076484; C_-1 reference texture gives m1/m3=0.60203.
- The heavy eigenspace amplitudes (magnitudes and phases) are recorded in `batch05_mode_phasors.csv`; the reference-texture counterfactuals isolate flavour input changes only.

### A_-2: Dirac

- Source masses: (0.007214753, 0.01100846, 0.04856021) eV. m1/m3=0.14857; candidate charge channels=1.
- Lightest mode: A1=0.16843, C1=0.48065; 51.935% of absolute projected strength cancels.
- Original source vertices reconstructed with relative errors 9.15e-17 and 1.43e-17.
- At common scalar mixing, C_+1 reference texture gives m1/m3=0.0076484; C_-1 reference texture gives m1/m3=0.60203.
- The heavy eigenspace amplitudes (magnitudes and phases) are recorded in `batch05_mode_phasors.csv`; the reference-texture counterfactuals isolate flavour input changes only.

### A_+0: Dirac

- Source masses: (0.007214275, 0.01100004, 0.04849895) eV. m1/m3=0.14875; candidate charge channels=1.
- Lightest mode: A1=0.16868, C1=0.48007; 51.993% of absolute projected strength cancels.
- Original source vertices reconstructed with relative errors 4.58e-17 and 1.57e-16.
- At common scalar mixing, C_+1 reference texture gives m1/m3=0.0076484; C_-1 reference texture gives m1/m3=0.60203.
- The heavy eigenspace amplitudes (magnitudes and phases) are recorded in `batch05_mode_phasors.csv`; the reference-texture counterfactuals isolate flavour input changes only.

### B_-3: Dirac

- Source masses: (0.007218231, 0.01098145, 0.0483324) eV. m1/m3=0.14935; candidate charge channels=1.
- Lightest mode: A1=0.16888, C1=0.48012; 51.988% of absolute projected strength cancels.
- Original source vertices reconstructed with relative errors 9.15e-17 and 5.86e-18.
- At common scalar mixing, C_+1 reference texture gives m1/m3=0.0076484; C_-1 reference texture gives m1/m3=0.60203.
- The heavy eigenspace amplitudes (magnitudes and phases) are recorded in `batch05_mode_phasors.csv`; the reference-texture counterfactuals isolate flavour input changes only.

### B_-1: Majorana

- Source masses: (0.004522493, 0.009803383, 0.05044131) eV. m1/m3=0.089659; candidate charge channels=1.
- Lightest mode: A1=0.15701, C1=0.20854; 79.146% of absolute projected strength cancels.
- Original source vertices reconstructed with relative errors 4.91e-17 and 0.00e+00.
- At common scalar mixing, C_+1 reference texture gives m1/m3=0.0076484; C_-1 reference texture gives m1/m3=0.60203.
- The heavy eigenspace amplitudes (magnitudes and phases) are recorded in `batch05_mode_phasors.csv`; the reference-texture counterfactuals isolate flavour input changes only.

### B_+1: Dirac

- Source masses: (0.007215827, 0.0109461, 0.04807432) eV. m1/m3=0.1501; candidate charge channels=1.
- Lightest mode: A1=0.15949, C1=0.50522; 49.478% of absolute projected strength cancels.
- Original source vertices reconstructed with relative errors 9.01e-17 and 4.10e-18.
- At common scalar mixing, C_+1 reference texture gives m1/m3=0.0076484; C_-1 reference texture gives m1/m3=0.60203.
- The heavy eigenspace amplitudes (magnitudes and phases) are recorded in `batch05_mode_phasors.csv`; the reference-texture counterfactuals isolate flavour input changes only.

### C_-3: Dirac

- Source masses: (0.02694687, 0.02819229, 0.05504074) eV. m1/m3=0.48958; candidate charge channels=1.
- Lightest mode: A1=0.26878, C1=1; 0.000% of absolute projected strength cancels.
- Original source vertices reconstructed with relative errors 6.93e-17 and 5.12e-17.
- At common scalar mixing, C_+1 reference texture gives m1/m3=0.0076484; C_-1 reference texture gives m1/m3=0.60203.
- The heavy eigenspace amplitudes (magnitudes and phases) are recorded in `batch05_mode_phasors.csv`; the reference-texture counterfactuals isolate flavour input changes only.

### C_-1: Majorana

- Source masses: (0.03749726, 0.03858887, 0.062285) eV. m1/m3=0.60203; candidate charge channels=1.
- Lightest mode: A1=0.5691, C1=0.28578; 71.422% of absolute projected strength cancels.
- Original source vertices reconstructed with relative errors 7.85e-17 and 1.38e-16.
- At common scalar mixing, C_+1 reference texture gives m1/m3=0.0076484; C_-1 reference texture gives m1/m3=0.60203.
- The heavy eigenspace amplitudes (magnitudes and phases) are recorded in `batch05_mode_phasors.csv`; the reference-texture counterfactuals isolate flavour input changes only.

### C_+1: Dirac

- Source masses: (0.0003654878, 0.008268705, 0.04778611) eV. m1/m3=0.0076484; candidate charge channels=1.
- Lightest mode: A1=0.3192, C1=0.013122; 98.688% of absolute projected strength cancels.
- Original source vertices reconstructed with relative errors 1.63e-16 and 1.08e-16.
- At common scalar mixing, C_+1 reference texture gives m1/m3=0.0076484; C_-1 reference texture gives m1/m3=0.60203.
- The heavy eigenspace amplitudes (magnitudes and phases) are recorded in `batch05_mode_phasors.csv`; the reference-texture counterfactuals isolate flavour input changes only.

### D_-2: Dirac

- Source masses: (0.007214175, 0.01099775, 0.04848237) eV. m1/m3=0.1488; candidate charge channels=1.
- Lightest mode: A1=0.1693, C1=0.47851; 52.149% of absolute projected strength cancels.
- Original source vertices reconstructed with relative errors 9.59e-17 and 5.22e-18.
- At common scalar mixing, C_+1 reference texture gives m1/m3=0.0076484; C_-1 reference texture gives m1/m3=0.60203.
- The heavy eigenspace amplitudes (magnitudes and phases) are recorded in `batch05_mode_phasors.csv`; the reference-texture counterfactuals isolate flavour input changes only.

### D_+0: Dirac

- Source masses: (0.00721467, 0.01100704, 0.0485499) eV. m1/m3=0.1486; candidate charge channels=1.
- Lightest mode: A1=0.16846, C1=0.48058; 51.942% of absolute projected strength cancels.
- Original source vertices reconstructed with relative errors 1.01e-16 and 8.92e-19.
- At common scalar mixing, C_+1 reference texture gives m1/m3=0.0076484; C_-1 reference texture gives m1/m3=0.60203.
- The heavy eigenspace amplitudes (magnitudes and phases) are recorded in `batch05_mode_phasors.csv`; the reference-texture counterfactuals isolate flavour input changes only.

### D_+2: Dirac

- Source masses: (0.0004747332, 0.007581738, 0.04359482) eV. m1/m3=0.01089; candidate charge channels=1.
- Lightest mode: A1=0.10255, C1=0.057118; 94.288% of absolute projected strength cancels.
- Original source vertices reconstructed with relative errors 8.12e-17 and 9.41e-18.
- At common scalar mixing, C_+1 reference texture gives m1/m3=0.0076484; C_-1 reference texture gives m1/m3=0.60203.
- The heavy eigenspace amplitudes (magnitudes and phases) are recorded in `batch05_mode_phasors.csv`; the reference-texture counterfactuals isolate flavour input changes only.

### E_-4: Dirac

- Source masses: (0.003319397, 0.007416937, 0.03848845) eV. m1/m3=0.086244; candidate charge channels=2.
- Lightest mode: A1=0.15726, C1=0.214; 78.600% of absolute projected strength cancels.
- Original source vertices reconstructed with relative errors 6.92e-17 and 1.55e-16.
- At common scalar mixing, C_+1 reference texture gives m1/m3=0.0076484; C_-1 reference texture gives m1/m3=0.60203.
- The heavy eigenspace amplitudes (magnitudes and phases) are recorded in `batch05_mode_phasors.csv`; the reference-texture counterfactuals isolate flavour input changes only.

### E_-2: Dirac

- Source masses: (0.004714191, 0.008678829, 0.04246948) eV. m1/m3=0.111; candidate charge channels=2.
- Lightest mode: A1=0.038958, C1=1; 0.000% of absolute projected strength cancels.
- Original source vertices reconstructed with relative errors 8.90e-17 and 2.28e-16.
- At common scalar mixing, C_+1 reference texture gives m1/m3=0.0076484; C_-1 reference texture gives m1/m3=0.60203.
- The heavy eigenspace amplitudes (magnitudes and phases) are recorded in `batch05_mode_phasors.csv`; the reference-texture counterfactuals isolate flavour input changes only.

### E_+0: Dirac

- Source masses: (0.004643584, 0.008648589, 0.04251536) eV. m1/m3=0.10922; candidate charge channels=2.
- Lightest mode: A1=0.038752, C1=1; 0.000% of absolute projected strength cancels.
- Original source vertices reconstructed with relative errors 1.31e-16 and 7.01e-17.
- At common scalar mixing, C_+1 reference texture gives m1/m3=0.0076484; C_-1 reference texture gives m1/m3=0.60203.
- The heavy eigenspace amplitudes (magnitudes and phases) are recorded in `batch05_mode_phasors.csv`; the reference-texture counterfactuals isolate flavour input changes only.

### E_+2: Dirac

- Source masses: (0.003325967, 0.007442485, 0.03863356) eV. m1/m3=0.08609; candidate charge channels=2.
- Lightest mode: A1=0.15766, C1=0.21322; 78.678% of absolute projected strength cancels.
- Original source vertices reconstructed with relative errors 1.58e-16 and 1.01e-16.
- At common scalar mixing, C_+1 reference texture gives m1/m3=0.0076484; C_-1 reference texture gives m1/m3=0.60203.
- The heavy eigenspace amplitudes (magnitudes and phases) are recorded in `batch05_mode_phasors.csv`; the reference-texture counterfactuals isolate flavour input changes only.

## Interpretation constraints

- A common Yukawa **matrix** is not a universal physical effective coupling across different SU(2) representations. The source-conditional group coefficients remain embedded in the extracted component tensors.
- Identical heavy singular-mass eigenvalues do not imply identical fermion propagator numerators or Majorana/Dirac contractions.
- Source Yukawas used to extract the tensors must be identifiable (full-rank); ambiguous tensors trigger an error rather than being inferred.
- Phasor angles refer to projections onto singular vectors of the *full* Mν. They are not physical CP phases and should not be compared as generation eigenvector phases within degenerate heavy subspaces.
- None of the controlled substitutions is used to certify the source Lagrangian, and no EFT calculation is performed.
