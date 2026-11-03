# Chat A — Batch 5: analytical patterns from fixed source-conditional Weyl tensors

**Purpose:** Explain *why* the frozen 16 T3 neutrino-mass matrices display their observed hierarchy patterns. Previous batches established the one-loop numerical result, the exact scalar kernel and the heavy-eigenspace coherence diagnostic. This batch makes the connection to the *actual frozen Yukawa matrices and component tensors*, using controlled counterfactual inputs.

## Assumptions and source boundary

- Fixed restricted T3 benchmark, renormalization scale `1e7 GeV`, two candidate Yukawa vertices `y1_HC` and `y2_HC`, zero scalar holomorphic B matrix, isolated same-charge S1/S2 2×2 mixing blocks. No extra BSM insertions or EFT.
- Uses the **same candidate Matchete/RGBeta-to-Weyl map** as Batches 1–4. Its physical spinor/CG normalization is not certified. No quantity in this batch changes `physical_normalisation_certified=false`.
- The independently stored Batch-3 scalar-kernel mass and three heavy-mass-group matrices provide numerical conservation laws. The file lineage and SHA256s must match the unchanged original source and saved reports.
- The source builder for the examined revision explicitly populates the vertex arrays using `conjugate(y1[:,r])` and `conjugate(y2[:,r])`. Tensor extraction below is justified by that source convention; it is **not** inferred solely from the fact that a 3×3 matrix can always be inverted. The frozen input and builder source hashes are checked.

## 1. Derive Yukawa-controlled components

In the original left-handed Weyl component basis, the saved source vertices are `Z` for the conjugate-S1 Yukawa and `Y` for the S2 field Yukawa. The source index contractions are

\[
Z_{iIs}=\sum_{r=1}^{3}y_{1,ir}^{*}\,C^{(1)}_{rIs},\qquad
Y_{iIt}=\sum_{r=1}^{3}y_{2,ir}^{*}\,C^{(2)}_{rIt}.
\]

Here `i` labels SM lepton flavour, `r` labels heavy flavour, `I` labels Weyl-fermion components and `s,t` label component scalars. The tensors `C^(a)` absorb component CG coefficients, fermion charge/conjugation routing, and source-specific signs. **They are not physically normalized until B1/B2.**

For each frozen model, solve

\[
C^{(1)}_{\cdot Is}=(y_1^*)^{-1}Z_{\cdot Is},\qquad
C^{(2)}_{\cdot It}=(y_2^*)^{-1}Y_{\cdot It}.
\]

A rank or conditioning check is mandatory: no unique tensor can be identified from a rank-deficient frozen Yukawa matrix. Rebuild both vertices using the original Yukawas and confirm the relative Frobenius residual is below `5e-9`. The production analytical mass must then reconstruct the original exact Batch-3 result within `2e-7` relative error, **including each heavy eigenspace matrix**, not just the three Takagi masses.

The extracted tensors are unchanged during controlled Yukawa substitutions, so all group-theory and Weyl-basis conventions in that model remain fixed.

## 2. Exact scalar and heavy-fermion kernel

For each original electric-charge sector `q`, in its `(S1,S2)` scalar basis,

\[
A_q=\begin{pmatrix} a_q&\delta_q^*\\ \delta_q&b_q\end{pmatrix},
\quad
\lambda_{q\pm}=\frac{a_q+b_q}{2}\pm
\sqrt{\frac{(a_q-b_q)^2}{4}+|\delta_q|^2}.
\]

The matrix-function off-diagonal propagator at zero momentum is

\[
\mathcal K_q(m_k^2)=\delta_q\,
\frac{B_0(0;m_k^2,\lambda_{q+})-B_0(0;m_k^2,\lambda_{q-})}
{\lambda_{q+}-\lambda_{q-}}.
\]

The finite, subtraction-scale-independent divided difference is exact for these isolated 2×2 blocks. With `Q_k` the eigenvector of `M_F M_F†` and `P_k=Q_k Q_k†`, the original restricted source-conditional matrix is

\[
M_\nu=\frac{1}{16\pi^2}\sum_{q,k}\mathcal K_q(m_k^2)
\operatorname{sym}\!\left[Z_{s_q}\,M_F^\dagger P_k\,Y_{t_q}^{T}\right],
\quad\operatorname{sym}[X]=X+X^T.
\]

`k` is aggregated into distinct heavy singular-mass eigenspaces, so the group matrix is invariant under rotations within a degenerate eigenspace. We never interpret individual degenerate heavy eigenvectors as unique physical contributions.

At small mixing, the relative correction to the linear scalar mixing law is `O(|lambdaT3|²)`. The existing Batch-3 result verified `3.6e-8`–`1.9e-7` relative first-order errors across the frozen 16 inputs. Batch 5 reuses the **exact** kernel, not the approximation.

## 3. Exact amplitude-phase cancellation identity

For the full complex SVD `Mν=U diag(m1,m2,m3)V†`, define

\[
a_{gi}=u_i^\dagger H_g v_i,
\qquad m_i=\sum_g a_{gi},\qquad
G=\sum_g\lVert H_g\rVert_F.
\]

The phase-coherent sum is a real positive singular value within numerical accuracy. For each mode,

\[
\frac{m_i}{G}=\underbrace{\frac{\sum_g|a_{gi}|}{G}}_{A_i\ \text{(projection)}}
\times\underbrace{\frac{|\sum_g a_{gi}|}{\sum_g|a_{gi}|}}_{C_i\ \text{(coherence)}}.
\]

And the explicit destructive cross terms obey

\[
\left|\sum_g a_{gi}\right|^2=
\sum_g|a_{gi}|^2+2\sum_{g<h}\operatorname{Re}(a_{gi}a_{hi}^*).
\]

Batch 5 records for every group its mass (GeV), complex amplitude (eV), magnitude, angle, and fraction of the sum of magnitudes. The strongest negative cross term is tabulated; this is a *diagnostic* of cancellation, not a CP-invariant phase claim.

Using the earlier user-run results:

- `C_-3`: `C1≈1`, heavy group projections reinforce each other in the lightest mode.
- `C_-1`: `C1≈0.28578`, partial cancellation.
- `C_+1`: `C1≈0.013122`, the projected magnitudes cancel by approximately **98.69%**.

The ratio `m1(C_-1)/m1(C_+1)≈102.6` decomposes into an overall heavy-group strength factor, a flavour projection factor, and a heavy-state coherence factor. The new report computes these factors directly from the saved complex matrices.

## 4. Controlled Yukawa-versus-representation experiments

The user-supplied original frozen benchmark remains **unchanged** throughout. For each candidate model `m`, replace only its two Weyl Yukawa vertex arrays with factorized arrays built from a *shared* `(y1,y2)` texture from one of

`A_-2`, `C_-1`, `C_+1`, `E_-2`.

Use common scalar mixing `|lambdaT3|=0.10`: A–D already use 0.10; E uses 0.05 frozen and its S1/S2 mixing blocks are multiplied by `2` **only in the counterfactual calculation**. All other scalar/fermion matrices, CG contractions, fermion nature, hypercharges, and source conventions stay as recorded in that model's original input. The full complex matrix is recalculated via the analytical scalar kernel. **Do not** interpret equal Yukawa matrices as physically equivalent normalized couplings across different SU(2) multiplets.

This yields `16×4=64` controlled results. For each model save all three masses, `m1/m3`, `m2/m3`, geometry `A1`, coherence `C1` and heavy-group strength. The report can distinguish:

1. A hierarchy inherited from the chosen Yukawa texture, if multiple model variants converge to a common hierarchy after substitution.
2. A hierarchy traceable to source representation/fermion/scalar structure, if controlled results remain substantially different.
3. Mixed effects when results change under both Yukawa and representation substitution.

These are **conditional causal counterfactuals within the recorded Weyl map**, not a physical model-fit inference.

For A/D, B, C and E reflected pairs, compute their best complex proportionality factor `c=<A,B>/<A,A>` and residual `||B−cA||/max(||B||,||cA||)` for both original matrices and four shared-reference textures. If the controlled complex-matrix residual is also small, that supports an *input-level* transformation relation; an all-parameter SU(2) identity still requires an independent analytic group derivation.

## 5. C-family single-Yukawa swaps and interpolation

For each of `C_-3`, `C_-1`, `C_+1`, evaluate seven source-conditional calculations: the original frozen texture; both Yukawas from C_-1; both from C_+1; C_-1 y1 only; C_-1 y2 only; C_+1 y1 only; C_+1 y2 only.

This is `3×7=21` counterfactual scenarios. It can identify when the cancellation is sensitive to one Yukawa matrix versus their joint alignment. The Dirac/Majorana nature **does not change** in any scenario.

Also linearly interpolate the complete *complex* Yukawa matrices of C_+1 toward C_-1, and vice versa:

\[
y_a(t)=(1-t)y_a^{\rm start}+t y_a^{\rm destination},
\qquad t=0,0.1,\ldots,1,\quad a=1,2.
\]

These `22` points are diagnostics for isolated cancellation surfaces in texture space. They are **not** a renormalization-group trajectory, a global fit, a probabilistic naturalness measure, or a scan over physical observables.

## 6. Results and inference limits

The runner writes a fresh, separate `analysis_batch05_01` directory with:

- `BATCH05_ANALYTICAL_SYNTHESIS.md`: formulae, all-16 mass patterns, C phasor table, controlled hierarchy table, C swaps/interpolation, matched-matrix proportionality diagnostics, and model-by-model interpretations.
- `batch05_synthesis.json`: complete machine-readable analytical results.
- `batch05_mode_phasors.csv`: all 48 original complex light-mode heavy-group contributions.
- `batch05_common_texture_controls.csv`: 64 controlled model/texture calculations.
- `batch05_C_texture_interpolation.csv`: 22 points.
- `batch05_matched_pairs.csv`: original versus controlled complex-matrix comparisons.
- `SOURCE_MANIFEST.json`, `HANDOFF.json`: original SHA256s, conditional-only limitations, and outputs lineage.
- Three optional `png` figures (when matplotlib is installed).

Actual model-specific counterfactual numbers cannot be asserted until the script is executed on the exact local frozen source checkout. The synthetic tests verify the analytic identities and report logic, not physical normalization. All files are created only after input lineage and matrix reconstructions pass, and an existing nonempty target directory is never overwritten.
