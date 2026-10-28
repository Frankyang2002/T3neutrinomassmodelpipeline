# T3 scalar relevance — Batch 5: RG generation orders from a zero-extra boundary

## Purpose and fixed scope

This batch follows the 16-model Batch 1 inventory, Batch 2 Wick-vertex audit and Batch 3–4 quadratic audits. It determines **whether a coupling set to zero at an input scale is radiatively generated from the retained interactions, and after how many iterated one-loop RG insertions**. It does not try to certify the full scalar potential or produce an accurate general-parameter neutrino mass.

For each of the 16 supported ordinary T3 models, the script reads two already saved candidate RGBeta exports, `UV` and `EFT1`, under `output/general_t3_quadratic_eft1/`. It validates them against the independent scalar manifest, Stage29 RGBeta quartic registry and saved symbolic matching report. All source files are SHA-256 fingerprinted in its output.

**No live Wolfram calls, new RG integrations or changes to production physics are made.** The symbolic reduction uses the exact rational coefficients in saved InputForm arithmetic and treats Yukawa traces as formal independent objects. It cancels identical monomials, and fails on unsupported functions or coupling-dependent denominators. Each symbolically nonzero source is *generic*; it may vanish at a chosen numerical point through a flavour or parameter cancellation.

## The counting rule

Let \(X_i\) be an omitted coupling, with \(X_i(\Lambda)=0\), and \(B_i\) its one-loop beta function:

\[
\frac{dX_i}{d\ln\mu}=\frac{B_i(X_j,g,y,\lambda)}{16\pi^2}.
\]

Retained couplings, mass parameters and flavour traces have RG-insertion order \(0\). An omitted coupling starts with no finite generation order. For a monomial \(m\) appearing in \(B_i\), assign

\[
\operatorname{ord}(X_i)=1+\min_{m\subset B_i}\sum_{X_j\in m}\operatorname{ord}(X_j),
\]

counting repeated omitted factors with multiplicity and relaxing until no order changes. A monomial independent of all omitted couplings sources \(X_i\) at order one. This is a **formal fixed-order iterative-beta expansion**; it is neither a proof of nonzero physical amplitudes nor a full RG solution.

A quartic sourced at order one and inserted into the usual one-loop matching diagram contributes formally at order \((16\pi^2)^{-2}\), multiplied by the appropriate logarithm. A quadratic term sourced only through an order-one quartic first appears after two one-loop beta insertions, with two logarithms. Nonzero tree-level quadratic boundary terms are **not suppressed** this way and remain a separate benchmark assumption.

## Required project inputs

- `RGE/general/generated/general_rgbeta_quartic_registry.json`
- `docs/general_scalar_trial_operator_manifest.json`
- `output/general_matching_symbolic_report.json`
- `output/general_t3_quadratic_eft1/<model>_uv_quadratic_live_rgbeta.json` for 16 models
- `output/general_t3_quadratic_eft1/<model>_eft1_quadratic_live_rgbeta.json` for 16 models
- `scripts/AuditT3QuadraticRGClosure.py` (from Batch 4, for its validated Wolfram InputForm parser)

The program checks the 50 quartic and 8 quadratic *model-specific* symbols and maps the four previously observed B/C matching quartics **by field structure only**, without claiming RGBeta-to-Matchete Clebsch–Gordan coefficient normalization.

## PowerShell

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "C:\path\to\T3_Scalar_Relevance_Batch5.zip" -DestinationPath . -Force
python -m pytest -q tests\python\test_t3_reduced_scalar_rg_batch5.py
python -m scripts.AuditT3ReducedScalarRG
Get-Content -Encoding UTF8 output\scalar_relevance\reduced_scalar_rg_order_summary.md
```

Generated only under `output/scalar_relevance/`:

- `reduced_scalar_rg_order_audit.json`: full stage-by-stage orders, source monomials, SHA-256 fingerprints and limitations.
- `reduced_scalar_rg_order_audit.csv`: one record per model/stage/extra scalar coupling.
- `reduced_scalar_rg_order_summary.md`: readable model breakdown.

## Important limitations

1. The RGBeta scalar basis remains a **candidate**; Matchete/RGBeta coupling normalization and physical component mapping are not certified.
2. The saved UV and EFT1 beta functions are analysed **separately**, not combined into a physical two-threshold trajectory. The direct quartics generated above the fermion threshold may be inherited by EFT1 even if their local EFT1 beta sources vanish.
3. The analysis counts RG insertions, **not** one-loop Feynman diagrams. It does not prove a coupling is forbidden, a vacuum is stable, or an exact RG-invariant subspace exists.
4. The result applies to boundary conditions where **all** additional scalar couplings start zero. If an extra quartic or quadratic is nonzero at tree level, it may affect the one-loop result without an extra loop factor.
5. Large scale logarithms, threshold matching, Wilson coefficients, renormalization conventions and a quantitative \(\delta C_5\) estimate remain outside this batch.

At the benchmark ratio \(\mu_{\rm upper}/\mu_{\rm lower}=100\), \(\ln(100)/(16\pi^2)\simeq0.02916\) per RG segment; the omitted coupling's beta-function coefficient sets the actual size. This number is **not a universal 2.9% neutrino-mass uncertainty**.
