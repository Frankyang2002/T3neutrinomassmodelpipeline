# A(-4): exact scalar Hessian and native RGBeta fresh-index experiment

**Scope:** independent analysis of the scalar-only one-loop quartic mismatch in the 45-Weyl/12-real-scalar EFT1, where **all heavy fermions are already integrated out** and both BSM scalars remain active. This is a read-only diagnostic batch. It does not change any production source or assert that the Matchete-to-RGBeta quartic dictionary is certified.

## Files and project paths

- `Numerical/running/Aminus4QuarticHessianIdentities.py` — exact symbolic SU(2) component proof of all 55 scalar-loop quadratic monomials.
- `scripts/RunAminus4RGBetaFreshIndexProbe.wl` — *native* four-mode RGBeta model build; redefines only the A(-4) RG-closed adjoint quartic registrations **in-memory**, never modifying source.
- `scripts/RunAminus4FreshIndexCausalBatch.py` — orchestrator, source/hash gates, all 55 cross-checks for each mode, logs and JSON report, signed-spinor convention check.
- `tests/test_a_minus4_fresh_index_hessian.py` — nine algebraic and convention-regression tests.
- `docs/A_MINUS4_SU2_HESSIAN_FRESH_INDEX_BATCH.md` — this file.

Nothing in `RGE/running/rgbeta/T3RGBetaModel.wl`, the installed RGBeta package, the actual beta exports, the threshold matching code, or the Wilson ODE is overwritten.

## Physics calculation

Use the declared potential `Numerical/running/Aminus4EFT1Couplings.py`:

\[
\begin{aligned}
V_4={}&\tfrac12\lambda_H(H^\dagger H)^2+\tfrac12\lambda_{S1}|S_1|^4
+\tfrac12\lambda_{S2}(S_2^\dagger S_2)^2
+\tfrac12\lambda_{S2Adj}\sum_A(S_2^\dagger T^A S_2)^2\\
&+\lambda_{H1}(H^\dagger H)|S_1|^2
+\lambda_{H2}(H^\dagger H)(S_2^\dagger S_2)
+\lambda_{H2Adj}\sum_A(H^\dagger t^A H)(S_2^\dagger T^A S_2)\\
&+\lambda_{12}|S_1|^2(S_2^\dagger S_2)
+\bigl[\lambda_{T3}(HH)_3 S_1 S_2^\dagger+\mathrm{h.c.}\bigr].
\end{aligned}
\]

All six complex components of \(H,S_1,S_2\) are \((R+iI)/\sqrt2\), with 12 independent canonically normalized real components. The SU(2) generators are the ordinary Pauli/2 doublet and spin-1 triplet matrices. The topology term uses the normalized symmetric doublet×doublet CG.

The scalar one-loop polynomial is obtained from **the Hessian of this potential**, independently of RGBeta beta expressions or the fourth-derivative projection:

\[
16\pi^2 \beta(V_4)\big|_{\mathrm{scalar}}
=\frac12\operatorname{Tr}\left[\left(\frac{\partial^2 V_4}{\partial\phi_i\partial\phi_j}\right)^2\right].
\]

Every one of the ten-potential-invariant coefficients is fixed by an exact coefficient map in the 12-variable polynomial ring; the full polynomial residual is checked to vanish for all 55 quadratic products. Nothing is fitted to RGBeta. The exact SU(2) identities are also verified, including

\[
\sum_A t^A_{ij}t^A_{kl}=
\tfrac12\delta_{il}\delta_{jk}-\tfrac14\delta_{ij}\delta_{kl},
\qquad
\sum_A(S_2^\dagger T^A S_2)^2=(S_2^\dagger S_2)^2-|2S_{2,1}S_{2,3}-S_{2,2}^2|^2.
\]

### Exact reference coefficients in disputed channels

| Input monomial | Output coefficients in \(16\pi^2\beta\) (Hessian reference) | Archived RGBeta output |
|---|---|---|
| \(\lambda_{S2Adj}^2\) | \(\beta_{S2}:8,\ \beta_{S2Adj}:6\) | \(16,20\) |
| \(\lambda_{S2Adj}\lambda_{H2Adj}\) | \(\beta_{H2Adj}:6\) | \(16\) |
| \(\lambda_{H2Adj}^2\) | \(\beta_H:1,\ \beta_{S2Adj}:1,\ \beta_{H2}:2\) | \(3,3,6\) |
| \((\mathrm{Re}\,\lambda_{T3})^2\) | \(\beta_{H2Adj}:8\) | \(8/3\) |
| \((\mathrm{Im}\,\lambda_{T3})^2\) | \(\beta_{H2Adj}:8\) | \(8/3\) |

For the topology squares, other beta components such as \(\beta_H,\beta_{H1},\beta_{H2},\beta_{12}\) **agree**; see `exact_hessian_reference.json` for the full coefficients.

### Specific native RGBeta hypothesis

Thomsen, *Introducing RGBeta*, Eur. Phys. J. C **81**, 408 (2021), section 4.2 explicitly warns that internal dummy-index names must be unique for each user-defined invariant invocation; it recommends `SetDelayed` + `Module` and notes that quartics sharing field content require simultaneous projection. (https://doi.org/10.1140/epjc/s10052-021-09142-4)

The inspected repository source `RGE/running/rgbeta/T3RGBetaModel.wl` uses literal internal adjoint index `A` for both `lambdaS2Adj` and `lambdaH2Adj`, whereas `lambdaT3` already uses an invariant with no such internal sum. The fixed index is a **candidate explanation**, not a demonstrated cause of incorrect output. The existing global `BetaTerm[Quartic,1]` projection is retained in this experiment.

The Wolfram probe runs four different, independent *fresh kernels*:

- `baseline`: original unchanged model.
- `triplet`: only `lambdaS2Adj` uses a fresh `Module[{adj}, ...]` dummy index.
- `portal`: only `lambdaH2Adj` uses a fresh dummy index.
- `both`: both invariants use fresh dummies.

Each experiment retains the same nine source quartic definitions/fields/normalizations. The code checks the source SHA-256; if the user's local `T3RGBetaModel.wl` differs, it **fails closed**. Each native output records a two-call raw invariant witness in addition to the nine projected beta functions. The Python driver independently re-evaluates all 55 scalar-only probes for each mode and checks the unchanged `baseline` numerically against the archived projection. This is more informative than just creating yet another empirically corrected beta.

Interpretation:

- `both=0/55`, baseline `34/55`: direct evidence that fresh dummy names resolve the discrepancy for the tested scalar sector, subject to direct review of the two-call witness.
- Some but not all discrepancies vanish: inspect `triplet` and `portal` to localize the channels; no single-cause conclusion.
- No improvement: literal `A` is not sufficient to explain the observed polynomial discrepancy; investigate internal projection/generator contraction conventions next.
- `baseline` differs from archived beta: **stop**; native installation/source has changed, so comparing interventions to old results is invalid.

This does not compare the full gauge- and Yukawa-dependent quartic beta, nor the physical Matchete invariant tensor map. Do not install a modified production beta on the strength of this experiment alone.

### Spinor sign/physical neutrino masses

The existing `declared_convention_boundary.json` fixes a **working** \(\rho=+1\) Weyl convention, but not the installed Matchete `GammaCC` sign. The two 36-coordinate fixed-native-Yukawa candidates are opposite tensors. For *any homogeneous linear* Wilson evolution \(C\mapsto\mathcal E(C,C^*)\), \(\mathcal E(-C,-C^*)=-\mathcal E(C,C^*)\); a linear scalar-threshold mass map changes \(m_\nu\mapsto-m_\nu\), which preserves Takagi singular values. This explains why an **overall** sign by itself does not change masses in that restricted linear setup. It does **not** establish the relative phase to an omitted fermion hard one-loop amplitude, nor certify the physical Wilson boundary.

## Local Windows PowerShell

```powershell
Set-Location "C:\thesistemp\T3neutrinomassmodelpipeline"

Expand-Archive `
  -Path "C:\Users\popok\Downloads\Aminus4_SU2_Hessian_FreshIndex_CausalBatch.zip" `
  -DestinationPath "C:\thesistemp\T3neutrinomassmodelpipeline" `
  -Force

python -m unittest discover -s tests `
  -p "test_a_minus4_fresh_index_hessian.py" -v

python scripts\RunAminus4FreshIndexCausalBatch.py `
  --output-root "output\a_minus4_next_batch\fresh_index_causal"
```

**Expected outputs:**

- `output/a_minus4_next_batch/fresh_index_causal/exact_hessian_reference.json` — fully exact 55-monomial scalar-only reference.
- `output/a_minus4_next_batch/fresh_index_causal/native/{baseline,triplet,portal,both}.json` — four freshly computed native quartic beta sets.
- `output/a_minus4_next_batch/fresh_index_causal/native/{baseline,triplet,portal,both}.log` — native stdout/stderr and command logs.
- `output/a_minus4_next_batch/fresh_index_causal/causal_comparison.json` — compact headline results and all nonzero residual monomials for each mode.

The Python driver returns a nonzero exit code and records a blocking reason if any native mode fails; `exact_hessian_reference.json` is still written first. It refuses to overwrite existing completed native reports.

The source code and its **nine Python unit tests were executed** in the authoring environment. The native WolframScript/RGBeta experiment **could not be executed** in that environment; results from its four modes must be established on the user's local installation. No `GammaCC` sign is inferred by this package.
