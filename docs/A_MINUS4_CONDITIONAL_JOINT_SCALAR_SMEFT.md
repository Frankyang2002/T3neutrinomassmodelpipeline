# A(-4) conditional EFT1 -> SMEFT: one common S1/S2 scalar threshold

## Goal and authority

This numerical batch begins with the **completed** full 45-Weyl/12-real
conditional A(-4) Wilson trajectory, itself evolved using the corrected native
RGBeta seven-scale EFT1 background. It then integrates out all eight heavy real
scalar components **in one matching operation at** `mu_S = 1000 GeV`.

Ordering (required):

`UV -> integrate out ALL fermions at mu_F -> EFT1(SM+S1+S2) -> integrate out ALL scalars together at mu_S -> SMEFT`

**Conditional only.** No native Matchete `GammaCC`-to-master-Weyl sign/magnitude
has been measured. The general gauge-dependent six-term Wilson RGE is still a
candidate, the Matchete-to-RGBeta invariant mapping has not been certified,
and hard one-loop fermion matching and other possible one-loop scalar-threshold
operator contributions have not been calculated. This reports the explicit
MSbar **O(C_LLSS lambda_HHSS) scalar bubble subset**, not a certified full
one-loop SMEFT result. Values called neutrino masses below are **mass proxies**
at `mu_S`, without SMEFT running to the electroweak scale.

## Actual files (all additive; complete files)

- `Numerical/joint_scalar_matching/Aminus4ConditionalSMEFT.py`
- `scripts/RunAminus4ConditionalSMEFT.py`
- `tests/test_a_minus4_conditional_smeft.py`
- `docs/A_MINUS4_CONDITIONAL_JOINT_SCALAR_SMEFT.md`

Nothing is written to GitHub or to production source files. The previous
A(-4) Wilson and renormalisable output directories are read-only.

## Existing input paths

- `output/a_minus4_next_batch/conditional_full_wilson_v1/summary.json`
- `output/a_minus4_next_batch/conditional_full_wilson_v1/full_conditional_C_trajectory.npz`
- `output/a_minus4_next_batch/fresh_index_production/running/A_m4_eft1_tensor_basis.json`
- `output/a_minus4_next_batch/fresh_index_production/running/A_m4_eft1_tensors.npz`
- `output/a_minus4_next_batch/fresh_index_production/A_m4_fixed_17beta_compatible.json`
- `output/general_observable/a_minus4_boundary_batch1/declared_convention_boundary.json`

`Aminus4ConditionalFullWilson.load_contract` rechecks the original corrected
source hash, 17-beta contract, mass/quartic basis, sample provenance and declared
45-Weyl boundary. The saved Wilson input hashes and boundary/endpoints are
checked again before matching. No new Wolfram or Wilson ODE is run.

## Numerical convention

At the unbroken common scalar threshold, use the actual corrected native running
mass parameters `mS1Sq(mu_S)` and `mS2Sq(mu_S)` in

\[
\mathcal M_{\mathrm{heavy,real}}^2 =
\operatorname{diag}(m_{S1}^2,m_{S1}^2,
 m_{S2}^2,m_{S2}^2,m_{S2}^2,m_{S2}^2,m_{S2}^2,m_{S2}^2).
\]

There are no *invented* mass-mixing terms or manually chosen scalar masses.
Both squared masses must be **positive**. Unbroken SU(2) and U(1) mass
commutators are tested. Any tachyonic/invalid running mass blocks this benchmark;
there is no fallback synthetic mass.

The real-scalar quartic tensor is taken directly from the previously verified
`lambda4(mu_S)[12,12,12,12]` in the convention
`V4 = lambda_abcd phi_a phi_b phi_c phi_d / 4!`. The full evolved
`C_ijab(mu_S)[45,45,12,12]` is extracted by named, mapped lepton positions,
not by incorrectly assuming that the first six Weyl fields are leptons.

Master Weinberg normalization:

\[
\mathcal L_5 \supset \frac12 C_{5,pq}(L_p\cdot H)(L_q\cdot H)
   + \mathrm{h.c.},\qquad
T_{\nu\nu,H_R^0H_R^0}=\frac14,
\qquad C_{5,pq}=4C_{\nu_p\nu_q,H_R^0H_R^0}.
\]

Use the full SU(2) template, **not only the pivot**, to project the already
transported `LLHH` tensor and the scalar bubble. The existing joint scalar
matching module uses a Weinberg template **twice** the master template; the
batch checks `C5_master = 2*C5_joint` numerically.

For each eigenstate pair `a,b`, the existing
`Numerical.joint_scalar_matching.core.scalar_bubble_matching` computes

\[
\Delta C_{ijAB}^{\rm scalar}(\mu_S)
=-\frac{1}{2(16\pi^2)}\sum_{a,b=1}^{8}
 C_{ijab}(\mu_S)\lambda_{ABab}(\mu_S)
 \overline B_0(0;m_a^2,m_b^2;\mu_S),
\]

where the full ordered sum counts off-diagonal pairs twice. All eight heavy
real scalar eigenstates are integrated out together. The loop is evaluated at
zero external momentum in an unbroken MSbar EFT, not in a broken/mixed scalar
spectrum. The code crosschecks at fixed background

\[
\frac{\partial\Delta C_5}{\partial\ln\mu_S}
=-\frac{\operatorname{Project}_{5}\big[\lambda_{HHSS}C_{LLSS}\big]}
 {16\pi^2}.
\]

The conditional matching subset is

\[
C_5^{\rm conditional}(\mu_S)=
 C_{5,\mathrm{direct\;LLHH}}(\mu_S)
 +\Delta C_5^{\rm finite\;8\;scalar}(\mu_S).
\]

For a **rough mass estimate only**, use `v = 246.22 GeV` by default:

\[
 m_\nu^{\rm proxy}(\mu_S)=-\frac{v^2}{2}
   C_5^{\rm conditional}(\mu_S),\quad
 m_i^{\rm proxy}=10^9\operatorname{svdvals}(m_\nu)[i]\ {\rm eV}.
\]

The three singular values are sorted ascending. Both consistent global signed
branches `rho=+1` and `rho=-1` have opposite `C5` and the **same singular
values**. A generic complex source phase does **not** necessarily factor out
of the anti-linear Wilson RGE.

## Run on Windows PowerShell

```powershell
Set-Location "C:\thesistemp\T3neutrinomassmodelpipeline"

Expand-Archive `
    -Path "C:\Users\popok\Downloads\Aminus4_ConditionalJointScalar_SMEFT_Benchmark.zip" `
    -DestinationPath "C:\thesistemp\T3neutrinomassmodelpipeline" `
    -Force

python -m unittest discover -s tests `
    -p "test_a_minus4_conditional_smeft.py" -v

python scripts\RunAminus4ConditionalSMEFT.py `
    --output-root "output\a_minus4_next_batch\conditional_joint_scalar_v1"
```

Run does not overwrite an existing output directory. A re-run requires a new
output directory. If the endpoint or corrected RGBeta source differs from the
prior verified run, the stage blocks rather than fabricate an alternative.

## Outputs

- `output/a_minus4_next_batch/conditional_joint_scalar_v1/summary.json`
  contains input SHA-256, both running scalar mass squares, matching and
  reference residuals, full 3x3 complex direct/bubble/total `C5`, three mass
  proxies in eV, and explicit certification flags.
- `output/a_minus4_next_batch/conditional_joint_scalar_v1/conditional_smeft_c5_and_mnu.npz`
  contains complete direct/bubble/total 3x3 `C5`, signed branch and complex
  neutrino mass proxy matrix.

A passing result means the **specific** numerical matching subset is operational.
It does **not** certify the UV spinor dictionary, the complete SMEFT hard
matching, or experimentally meaningful neutrino masses.
