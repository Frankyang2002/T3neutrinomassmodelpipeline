# Workstream B — B(-1)/C(-1) Higgs mass beta and five-real-mass closure

## Scope and verified source context

Read-only source inspection of the latest repository `main` before implementation:

- `RGE/running/rgbeta/T3RGBetaModel.wl` defines the SM Higgs doublet and two BSM doublets, and the complete B/C scalar quartic interactions, but **does not** register a Higgs quadratic mass parameter.
- `RGE/running/rgbeta/RunGeneralT3QuadraticRGBeta.wl` registers the mixed B/C scalar bilinear and exports the BSM masses, but not the Higgs mass beta.
- The existing B/C UV/EFT1 saved exports contain `mS1Sq`, `mS2Sq`, and `muT3GeneralS1S2` betas. They omit `mHSq` and thus do **not** include its back-reaction on these mass betas.
- The completed Workstream B mixed-bilinear verification supports the **conditional** common-component action map `mu_RGBeta = -mu_Matchete` for both B/C models; the independent global relative field phase remains unfixed.

This addition **does not replace any existing production file**. It builds the current RGBeta model, registers the missing Higgs mass in the potential convention, includes the already existing seven extra quartic operators and mixed scalar mass, then obtains all four complex-parameter mass betas from the installed RGBeta runtime.

## Conventions and equations

The relevant quadratic potential is

\[
V_2=m_H^2 H^\dagger H+m_1^2S_1^\dagger S_1+m_2^2S_2^\dagger S_2
+[\mu_R\,S_{1i}^{*}\epsilon_{ij}S_{2j}^{*}+\mathrm{h.c.}],
\qquad \epsilon_{12}=1.
\]

`mHSq`, `mS1Sq`, `mS2Sq` are real; `muT3GeneralS1S2` is complex. They form **five real running variables**:

\[
x=(m_H^2,m_1^2,m_2^2,\mathrm{Re}\,\mu_R,\mathrm{Im}\,\mu_R)^T,
\qquad 16\pi^2\,dx/d\ln Q=A(g,y,\lambda)\,x+b(g,y,\lambda,M_F).
\]

At one loop, with dimensionless couplings and fermion masses frozen, this is an **affine** mass-only system. UV fermion-mass terms can provide the inhomogeneous source `b`; EFT1 has no dynamical heavy fermion `F` in its renormalisable sector.

For a useful independent limit with `mu=0` and B/C's extra complex quartics zero, the expected Higgs beta is

\[
\beta_{m_H^2}^{(1)}=m_H^2\!\left(6\lambda_H-\frac92g_2^2-\frac32g_Y^2
+6\operatorname{tr}(Y_uY_u^\dagger)+6\operatorname{tr}(Y_dY_d^\dagger)
+2\operatorname{tr}(Y_eY_e^\dagger)\right)
+4\lambda_{H1}m_1^2+4\lambda_{H2}m_2^2.
\]

This follows the **current RGBeta** potential normalization `lambdaH/2*(H†H)^2` and ordinary `lambdaH1,lambdaH2` portal definitions. The adjacent-generator portals `lambdaH1Adj,lambdaH2Adj` are trace-free in this limit. The verifier checks the independent numeric limit, plus the four reciprocal portal derivatives. Other mixed-mass-dependent terms in the full Higgs beta are obtained **from the live RGBeta result**; they are *not* guessed or hard-coded.

## Files to install (under project root)

1. `RGE/running/rgbeta/RunT3BCFullScalarMassRGBeta.wl` — **new**, opt-in live scalar mass RGE runner.
2. `scripts/VerifyT3BCFullScalarMassClosure.py` — **new**, checks metadata, five-real-mass closure, numerical reality, portal coefficients, restricted Higgs beta, and a short frozen-coupling mass-only evolution.
3. `tests/test_t3_bc_full_scalar_mass_closure.py` — **new**, synthetic and negative control tests.
4. `docs/WORKSTREAM_B_BC_HIGGS_MASS_AND_CLOSURE_V1.md` — **new**, this guide.

## Installation and commands (PowerShell)

From `C:\thesistemp\T3neutrinomassmodelpipeline`:

```powershell
Expand-Archive -LiteralPath ".\Workstream_B_BC_Higgs_Mass_Closure_v1.zip" -DestinationPath . -Force
python -X utf8 -m pytest -q tests\test_t3_bc_full_scalar_mass_closure.py

$root = "output\audits\workstream_b_bc_full_mass"
foreach ($cls in @("B", "C")) {
  foreach ($stage in @("UV", "EFT1")) {
    $name = "${cls}_m1_$($stage.ToLower())_full_mass_v1.json"
    wolframscript -file RGE\running\rgbeta\RunT3BCFullScalarMassRGBeta.wl `
      $stage $cls "$root\$name"
    if ($LASTEXITCODE -ne 0) { throw "RGBeta $cls/$stage failed: exit code $LASTEXITCODE" }
  }
}

python -X utf8 -m scripts.VerifyT3BCFullScalarMassClosure `
  --directory output\audits\workstream_b_bc_full_mass `
  --output output\audits\workstream_b_bc_full_mass\closure_v1.json
```

**Do not overwrite existing reports:** both scripts reject existing output paths. If rerunning, use a different filename / new output directory and adjust the verifier's `--directory`. The verifier expects filenames `B_m1_uv_full_mass_v1.json`, `B_m1_eft1_full_mass_v1.json`, `C_m1_uv_full_mass_v1.json`, `C_m1_eft1_full_mass_v1.json` within that directory.

## What a successful run establishes

- The installed RGBeta has actually returned `BetaTerm[mHSq,1]` and the three BSM mass betas with *all* registered quartics and bilinear terms.
- The full five-real-parameter scalar mass beta is affine in the mass parameters at a fixed set of renormalisable couplings.
- The expected SM-plus-two-doublet limit and the scalar mass-portal back-reaction agree numerically with the independent convention check.
- The short mass-only flow solves that affine system exactly **at fixed dimensionless couplings**, via a 6×6 augmented matrix exponential.

**This does not yet establish:** full coupled integration of every dimensionless beta, the Matchete–RGBeta global component-phase identification, mass-eigenstate sequential matching, a three-generation Wilson anomalous dimension, or a neutrino mass prediction. No production RGE or threshold code is changed. These remain separate from the rough mass-sector diagnostic.

## Failure policy

Missing native files, changed registry, missing mass betas, complex diagonal-mass betas, failed portal/SM limits, and non-affine mass flow fail the verifier explicitly. A synthetic pytest pass does **not** mean the installed Wolfram runtime has passed; only the user's live four exports and closure report establish that.
