# T3 scalar-potential relevance: leading-one-loop benchmark conclusion

**Status:** restricted fixed-order one-loop **benchmark approximation**, not the full gauge- and Z2-invariant scalar theory. This consolidates the local Batch 1–7 audits and the locally reported successful Batch 7 B/C runs (28 October 2026). No new fitting or exact nonzero-extra-coupling matching is claimed.

## 1. Declared approximation

For each ordinary T3 model (16 assignments), keep the scalar mass terms, the baseline Hermitian self/portal quartics, and the original mixed T3 topology coupling. With model-specific gauge invariant contractions implicit,

\[
V_{\mathrm{red}}=V_{\mathrm{SM}}+m_{S_1}^2S_1^\dagger S_1+m_{S_2}^2S_2^\dagger S_2
+V_{\mathrm{baseline,quartics}}
+\big[\lambda_{T3}(HH S_1 S_2^\dagger)_{\mathbf{1}}+\mathrm{h.c.}\big].
\]

At the **specified** UV boundary \(\mu_{\rm UV}=10^7\) GeV set *all additional symmetry-allowed scalar quartic and quadratic coefficients* to zero. This is an input boundary condition, **not a new protecting symmetry**. The original \(\lambda_{T3}\), Yukawa couplings and ordinary Hermitian baseline quartics are retained. Perform leading-dimension-five, fixed-order one-loop matching with diagonal tree-level scalar propagators. An arbitrary nonzero additional quartic or quadratic term is outside this approximation.

## 2. Operator and RG evidence from the six relevance audits

- 16 ordinary models, **50 additional complex quartic slots** (24 Higgs-bearing, 26 no-Higgs), **eight additional quadratic slots**.
- Under the original Yukawa scalar orientations and diagonal propagators, only **four** extra quartic slots can directly supply the two-Higgs one-loop Weinberg topology: `lambdaT3ExtraQ020020I1` and `lambdaT3ExtraQ020200I1`, each in B(-1) and C(-1). Their symbolic appearance in saved Matchete \(C_5\) matching was independently observed. The remaining **46** do not supply that specific single-quartic-vertex amplitude; this is not their irrelevance to the complete potential, RG, EWSB or vacuum.
- With every additional quartic set zero at the UV boundary, the saved **candidate RGBeta** one-loop equations source **six** quartic slots at one iterated RGE order (three each in B(-1)/C(-1), including the four directly \(C_5\)-eligible slots).
- In the *separately audited* candidate UV/EFT1 one-loop beta-source chains, **none of the eight additional quadratic slots is reached** from the completely zero-additional-coupling boundary. This does **not** prove vanishing under genuine two-loop beta functions, finite thresholds, or the complete quantum theory.
- Nonzero independently specified quartics or bilinears can affect leading one-loop \(C_5\) and cannot be neglected on loop counting alone. In particular, additional quadratic terms alter scalar propagators and Yukawa rotations at tree level.

With \(L=\ln(\mu/\mu_{\rm UV})\), the first radiatively generated extra quartic behaves as

\[
\lambda_{\rm extra}(\mu)\simeq \frac{L}{16\pi^2}\,\beta_{\lambda,\mathrm{extra}}^{(1)}\big|_{\lambda_{\rm extra}=0}.
\]

Inserting it into a **one-loop** matching graph yields a contribution formally of order \((16\pi^2)^{-2}L\) relative to tree-level coupling factors. This is formal power counting, not a uniform numerical bound, especially with large logarithms, cancellations or changes of basis.

## 3. Locally reported Batch 7 numerical benchmark

The user ran `python -m scripts.RunT3ZeroExtraBoundaryBenchmark --all --c5-preview`, and reported success for B(-1) and C(-1). The script integrates native candidate RGBeta from \(10^7\to 10^5\to10^3\) GeV and reuses a **tree-level** fermion-threshold projection, with all seven additional quartics and the one additional quadratic coefficient set exactly zero at \(10^7\) GeV for each B/C model.

| Candidate RGBeta coordinate at \(10^3\) GeV | B(-1) | C(-1) |
|---|---:|---:|
| `lambdaHHdagS2S2` | \(-3.8118134961\times 10^{-11}\) | \(+1.3346023870\times 10^{-10}\) |
| `lambdaHHdagS1barS1bar` | \(+1.2753856756\times 10^{-11}\) | \(-4.4651950563\times 10^{-11}\) |
| \(\lvert\mu_{\rm extra}^2\rvert/\min(m_{S_1}^2,m_{S_2}^2)\) | 0 | 0 |
| **UNCERTIFIED** one-scale \(\|\Delta C_5\|_F/\|C_5^{(0)}\|_F\) | \(4.982653\times10^{-10}\) | \(3.357537\times10^{-9}\) |
| Frozen-scenario leading Takagi mass (eV) | 0.04331467 | 0.05282823 |

The last two lines **use a frozen Yukawa/scalar-diagonal one-scale scenario** and an **uncertified identity map** between the native RGBeta quartics and Matchete quartic coefficients. They are useful magnitude checks, **not physical \(C_5\) error bounds or new neutrino-mass predictions**. In particular, no nonzero-mixing scalar mass-basis \(C_5\) was calculated here.

## 4. What may be stated in the project

For a deliberately zero-additional-coupling UV benchmark, the existing reduced scalar potential captures the **direct fixed-order leading-one-loop** Weinberg contribution in the specified diagonal tree-level scalar basis. The specific B/C source quartics are RG-induced beyond this fixed order; their candidate numeric values are very small for the chosen benchmark. The original ordinary Higgs/scalar portal quartics remain part of \(V_{\rm red}\).

**Do not state** that the model's complete renormalizable scalar potential contains no other terms, that omitted terms are protected by a symmetry, or that the neutrino-mass prediction is accurate to \(10^{-9}\) for generic UV parameters. Matchete/RGBeta basis normalizations, physical scalar threshold matching, scalar vacuum stability and genuinely higher-order effects remain unverified.

## 5. Reproducibility and next stage

Previous local outputs:

- `output/scalar_relevance/final_leading_scalar_potential.md` (Batch 6 synthesis)
- `output/scalar_relevance/zero_boundary_benchmark/summary.md` (Batch 7)
- `output/scalar_relevance/zero_boundary_benchmark/B_m1/result.json`
- `output/scalar_relevance/zero_boundary_benchmark/C_m1/result.json`

The next **numerical continuation** reuses the existing model-specific hierarchical fixed-order \(C_5\) and SM+Weinberg running (Stage 36), but evaluates them on the reconstituted **zero-additional-coupling** B/C RG backgrounds. The comparison to previously saved Stage 36 results is a **difference between two UV boundary scenarios**, *not* an isolated direct extra-quartic \(C_5\) correction. See `docs/T3_ZERO_BOUNDARY_LOW_ENERGY_BATCH8.md`.

Read-only reference code snapshot inspected: `Frankyang2002/T3neutrinomassmodelpipeline`, commit `a02e04c6def5cf8522b2f612ccd6bff7d9963b91`. No GitHub changes are required.
