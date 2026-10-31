# Workstream C v18 — B/C three-heavy-flavour frozen numerical reference

## Project file locations

- `scripts/EvaluateT3BCThreeFermionTriangle.experimental.py` — full new experimental evaluator
- `scripts/RunT3BCThreeFermionTriangle.ps1` — full new PowerShell runner
- `docs/WorkstreamC_ThreeFermionTriangle_v18_README.md` — these notes

Existing v12–v17 scripts are **not replaced**. No UV model, RGE, production matching, or GitHub files are modified.

## What is new

v17 verified a signed scalar-quartic / LLSS three-propagator master separately for the XX, XY and YY scalar pairs, and verified its strict expansion against the saved F/Y/X Matchete sources. Its optional numerical interface accepts only one heavy fermion (a single 3-vector each for y1/y2 and one MF). The existing B/C frozen reference instead contains 3x3 Yukawa matrices and three diagonal fermion masses.

v18 evaluates all three heavy-flavour columns, **with distinct mass denominators and loop functions**, for each B/C model. It reads:

- `output/general_t3_mass_eigen_matchete/{B_m1,C_m1}_final_parameters.json` for the **current** scalar masses and mixing matrix.
- `Numerical/benchmarks/fitted_reference/T3_dS1_2_dS2_2_dF_{1,3}_alpha_m1.json` for the frozen **UV-base-state** 3x3 Yukawas, 3x3 diagonal fermion masses, reference lambdaT3, and VEV.
- Saved v14/v17 evidence and original F/Y/X Matchete raw coefficients, through the existing v16 and v17 source-certified functions.

The reference Yukawas are defined in `base_state` at **1e7 GeV**, whereas MF is approximately 100/110/120 TeV. v18 **does not run** them to the fermion thresholds. It deliberately sets the *two independent native extra quartics* `lambdaT3ExtraQ020020I1` and `lambdaT3ExtraQ020200I1` to zero as an explicit illustrative scenario; no automatic claim of a Matchete-to-RGBeta quartic mapping is made.

## Maths

For each heavy-fermion flavour `a=1,2,3`, let

`L_a = (conj(y1[:,a]), conj(y2[:,a]))` (a 3x2 flavour matrix) and

`W_XX=-P_X`, `W_XY=-K_Y`, `W_YY=-P_Y+K_Y`,

where the three-component Yukawa-structure polynomials are re-extracted using v16/v17, so all representation-dependent Clebsches are retained. For scalar pair `rs`,

`Phi_rs(MFa^2) = - MFa^2 * int_0^infinity t dt/((t+MFa^2)(t+M_r^2)(t+M_s^2))`.

For each of three comparison schemes, `K_a = W_XX*Phi_XX + W_YY*Phi_YY + W_XY*Phi_XY`, with each `W` carrying the `(y1y1,y1y2,y2y2)` coefficients but not Yukawa matrices. The symmetric matrix reference is:

`C_ref = (1/(16*pi^2)) sum_a (1/MFa) L_a [[2*K_a[0], K_a[1]],[K_a[1],2*K_a[2]]] L_a.T`.

This uses the v11.2 validated **PL/Warsaw positive sign** and `C5=A+A.T` convention; **PR is not added again**. The three comparison schemes are:

1. **strict** — leading 1/MF hard expansion and zero-MX light-X limit at Y.
2. **finite_XY** — same leading 1/MF but retaining `r ln(r)/(1-r)` in the XY loop.
3. **full_mass** — independently evaluated exact three-mass triangle for each XX/XY/YY pair, still with fixed couplings and no RGE.

For *reference only*, `m_ref(eV) = -246.22^2 / 2 * C_ref(GeV^-1) * 1e9`, and singular values of the symmetric complex matrix are reported in ascending order. These are **not measured or physically certified neutrino masses**.

## Run

From project root, after downloading and extracting the ZIP:

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$HOME\Downloads\T3_WorkstreamC_ThreeFermionTriangle_v18.zip" -DestinationPath . -Force
& .\scripts\RunT3BCThreeFermionTriangle.ps1
```

If PowerShell execution policy blocks scripts:

```powershell
python scripts\EvaluateT3BCThreeFermionTriangle.experimental.py --self-test
python scripts\EvaluateT3BCThreeFermionTriangle.experimental.py `
  --output output\audits\t3_bc_three_fermion_triangle_v18\summary.json
```

The sole output is `output/audits/t3_bc_three_fermion_triangle_v18/summary.json`. It contains B and C data, source SHA-256 hashes, three 3x3 matrices (strict, finite_XY and full_mass), illustrative mass singular values and relative norm differences.

## Certification boundaries

- **Actually implemented:** three-heavy-column sum, exact scalar mass dependence in signed triangle, numerical reference matrix and singular values.
- **Verified from source:** F/Y/X polynomial relations through v17 and v16; correct B/C scalar masses; reference frozen Yukawa/mass input format; positive Warsaw PL identity v11.2 (upstream evidence).
- **Inferred/assumed:** native extra quartics set to zero, reference UV Yukawas held fixed at lower scale.
- **Not completed:** full UV vertex sign and finite matching scheme, matched EFT running and finite parameter/EOM/wavefunction shifts. No physical C5 certificate.

Do not paste the numerical reference matrix into production `C5`, neutrino fitting, or RG running as a physical threshold boundary.
