# T3 scalar-sector relevance: final restricted one-loop benchmark (Batch 6)

## Goal and result

This batch consolidates the five saved local audits into a **conditional, restricted leading-order scalar potential**, without changing the production model, extending C5 matching, or doing a new numerical RGE integration.

For the 16 supported ordinary T3 representations, the additional gauge- and Z2-invariant scalar operator inventory contains 50 complex quartic *model-slots* and eight complex quadratic *model-slots*. Of the 50 quartic slots, the Batch 2 Wick audit rules out 46 as **direct single-vertex insertions** in the specified unbroken, leading-dimension-five, diagonal-propagator one-loop LLSS -> LLHH coefficient. The four remaining quartics are two each in `B_-1` and `C_-1`, where the saved symbolic C5 depends on them if they are nonzero.

Crucially, those four are **not** higher-order operators. They enter a one-loop C5 if their UV couplings are nonzero at tree level. The reduced benchmark fixes all additional quartics and quadratic bilinears to zero at a **named UV input scale**, rather than asserting that gauge invariance removes them.

Under that boundary condition, the candidate one-loop RGBeta UV equations contain six immediately sourced additional quartic slots (three per B/C(-1)); four of them are the C5-active pair. Such generated quartics enter a subsequent one-loop C5 graph at **formal two-loop order with a logarithm**. The candidate beta-source graph did not find a finite sequence of iterated one-loop source insertions that regenerates any of the eight bilinears in the separately examined UV or EFT1 systems. This is narrower than proving an exact RG-invariant scalar subsector, including threshold matching and genuine higher-loop beta functions.

## Potential and approximation definition

With independent SU(2) contractions implied, take

\[
\begin{aligned}
V_{\rm reduced}={}&V_{\rm SM}
+ m_1^2 S_1^\dagger S_1 + m_2^2 S_2^\dagger S_2
+ V_{\rm ordinary,quartic}\\
&+ \left[\lambda_{T3}(HH S_1 S_2^\dagger)_{\mathbf1} + \text{h.c.}\right].
\end{aligned}
\]

`V_ordinary,quartic` consists of the Hermitian self and portal families constructed by the existing `ScalarPotential.wl` generator. It is not being discarded simply because it is not a direct leading C5 insertion: those terms can affect scalar stability, physical masses after electroweak breaking, and beta functions.

The full gauge/Z2-allowed scalar potential permits an additional sum over the catalogued quartic and quadratic operators. Define the restricted boundary by

\[
\lambda_{\mathrm{extra},i}(\Lambda)=0,\qquad
\mu^2_{\mathrm{extra},a}(\Lambda)=0.
\]

This is a *parameter-space restriction*, not a protecting symmetry. Yukawa tensors, fermion masses and the retained `lambdaT3` are not fixed by this audit.

The direct leading fixed-order one-loop matching is performed in the unbroken phase, with diagonal free scalar propagators and one scalar quartic insertion. Operators with additional Higgs VEV insertions, dimension-seven contributions, self-energy insertions and higher-loop corrections lie outside this fixed-order leading-dimension-five statement.

## Physics power counting

Using `\epsilon = (16\pi^2)^{-1}`, a direct extra scalar coupling regenerated from a one-loop beta source can satisfy

\[
\lambda_{\rm extra}(\mu_F)
  \sim \epsilon\,\beta_{\lambda,0}^{(1)}\,
  \ln\!\left(\frac{\mu_F}{\Lambda}\right).
\]

The corresponding direct one-loop C5 insertion has **two formal loop factors**:

\[
\delta C_5 \sim \epsilon^2\,
\ln\!\left(\frac{\mu_F}{\Lambda}\right)
\,\mathcal{F}(y_1,y_2,\lambda_{T3},M_F,m_{S_i}) .
\]

At a 100:1 scale ratio, `|ln(100)|/(16*pi^2) = 0.02916`. This is the logarithmic expansion factor *before* beta coefficients, Yukawas, tensor contractions, normalizations and cancellations. **It is not a 2.9% neutrino-mass error bound**.

A nonzero quadratic bilinear is **tree-level in the scalar propagator**, so it may affect a one-loop neutrino-mass coefficient without an extra loop suppression. The zero-boundary assumption is essential for any of the eight cases admitting one.

## What is source-verified, and what is not

**Supported by the Batch 1–5 local reports:** the 16-model operator inventory; the saved matching dependence for four direct quartics; their Wick compatibility; analytical candidate scalar Hessian spectra; saved candidate RGBeta beta-expression source structure; the absence of a path to extra bilinears in the iterated one-loop source graph when all extra boundary couplings vanish.

**Not verified:** absolute Matchete/RGBeta coupling normalization; full two-threshold Wilson matching with additional quartics and arbitrary quadratic mass mixing; physical C5 error bars; threshold-induced extra coupling generation; genuine two-loop beta-function sources; boundedness from below of the chosen scalar potential; general arbitrary-parameter neutrino masses.

The generated synthesis refuses missing, inconsistent or uncertified input reports, and SHA-256 fingerprints each input. It does not re-check live Matchete or RGBeta; those were the sources used by the earlier audits.

## Deliverables and local use

This ZIP adds exactly three files, preserving these project paths:

- `scripts/FinalizeT3LeadingScalarPotential.py` — read-only synthesis of saved Batch 1–5 reports.
- `tests/python/test_t3_leading_scalar_potential_batch6.py` — offline tests, synthetic fixtures and fail-closed checks.
- `docs/T3_LEADING_SCALAR_POTENTIAL.md` — this theoretical scope statement.

In PowerShell (change the ZIP source to your actual download directory):

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "C:\path\to\T3_Scalar_Relevance_Batch6.zip" -DestinationPath . -Force
python -m pytest -q tests\python\test_t3_leading_scalar_potential_batch6.py
python -m scripts.FinalizeT3LeadingScalarPotential
Get-Content -Encoding UTF8 output\scalar_relevance\final_leading_scalar_potential.md
```

Outputs, *generated on your own checkout*, are:

- `output/scalar_relevance/final_leading_scalar_potential.json`
- `output/scalar_relevance/final_leading_scalar_potential.md`

Neither output is a new independent physical calculation. No production code or GitHub files are changed.

## Honours-level stopping criterion

For the stated leading-dimension-five, fixed-order one-loop neutrino-mass objective, this conditional omission argument is sufficient to define a **restricted benchmark potential**. It is not sufficient to claim a *general scalar potential* or an *arbitrary-mixing, general one-loop neutrino-mass prediction*. Before presenting benchmark physics, specify the exact numerical UV boundary, check perturbative running and vacuum viability, and avoid claiming a numerical C5 uncertainty unless matching and RG improvements are actually performed.
