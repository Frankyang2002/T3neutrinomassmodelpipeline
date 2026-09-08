# T3 Weinberg-Coefficient RGE Pipeline

Note that RGBeta only supports fundamental (2D) and adjoint (3D) representations, thus we cant use higher representations for now.

This document describes the renormalisation-group (RGE) part of the T3
neutrino-mass pipeline.

The implemented flow is

$$
C_5(M)
\;\longrightarrow\;
\beta_{C_5}
\;\longrightarrow\;
C_5(\mu)
\;\longrightarrow\;
m_\nu(\mu)
\;\longrightarrow\;
\text{neutrino observables}.
$$

Here, $M$ is the heavy-particle matching scale and $\mu$ is a lower scale.
Below $M$, the heavy T3 fields have been integrated out, so the running is
performed in the Standard Model Effective Field Theory (SMEFT).

For the T3-B, $\alpha=-1$ Ma/scotogenic benchmark, an additional
common-threshold calculation can now run the full UV theory before matching:

$$
\text{Ma model}(\Lambda)
\longrightarrow
\text{Ma model}(M)
\longrightarrow
C_5(M)
\longrightarrow
C_5(\mu).
$$

---

## 0. Quick Summary

- `data/c5_coefficient.txt` is the main input from matching.
- The symbolic stages run automatically after a successful match.
- Numerical running only occurs when `--numerical CONFIG.json` is supplied.
- `MaFullRunningComparison.py` separately compares frozen, UV-only,
  EFT-only, and full UV-to-EFT evolution for the Ma benchmark.
- `rge_report.pdf` is the main human-readable RGE result for each model.
- The later stages do **not** read `c5_beta.txt`, `c5_flavor_matrix.txt`, or
  `c5_loop_kernel.txt`; these are reports, not pipeline dependencies.
- Detailed intermediate expressions are only written with `--debug-reports`.

---

## 1. Starting Point: the Matched Coefficient

The matching stage produces the coefficient of the Weinberg operator,

$$
\mathcal L_{\rm EFT}
\supset
C_5\,\mathcal O_5 + \text{h.c.}
$$

The coefficient is saved as

```text
output/<model>/data/c5_coefficient.txt
```

We use this file for RGE, flavor, and neutrino-mass
stages. The model directory also contains `c5_coefficient.tex` and
`c5_coefficient.pdf` for direct inspection.

---

## 2. One-Generation Symbolic RGE

`MatchedEFTRGE.py` reads the Matchete expression and computes its SMEFT beta
function. The code uses

$$
t=\ln\mu
$$

and the Higgs-potential convention

$$
V(H)=\frac{\lambda_H}{2}(H^\dagger H)^2.
$$

For one generation, the internal benchmark is

$$
\frac{16\pi^2\beta_{C_5}}{C_5}
=
-3g_2^2+2\lambda_H
+6|y_u|^2+6|y_d|^2-|y_e|^2.
$$

The stage fails if its symbolic result does not reduce to this expression.
The useful output is `c5_beta.txt`; the ratio
`c5_beta_over_c5.txt` is a compact diagnostic written only in debug mode.

---

## 3. Full Three-Generation Flavor RGE

The one-generation matched coefficient is factorised as

$$
C_5^{(1)}=F_{\rm loop}\,\bar y_1\bar y_2,
$$

where $F_{\rm loop}$ contains the masses, scalar coupling, and loop
functions. The code then restores lepton and heavy-fermion flavor indices:

$$
(C_5)_{pq}
=
\frac12\sum_r F_r
\left[
y_{1,pr}^*y_{2,qr}^*
+y_{2,pr}^*y_{1,qr}^*
\right].
$$

This construction makes $C_5$ a complex symmetric $3\times3$ matrix. It
currently assumes a diagonal heavy-fermion mass basis and common scalar
masses.

The full one-loop SMEFT equation implemented by
`SMEFTWeinbergFlavorRGE.py` is

$$
16\pi^2\frac{dC_5}{d\ln\mu}
=
\left(2\lambda_H-3g_2^2+2T\right)C_5
-\frac32\left[
Y_eY_e^\dagger C_5
+C_5\left(Y_eY_e^\dagger\right)^T
\right],
$$

with

$$
T=\operatorname{Tr}\!\left(
Y_eY_e^\dagger
+3Y_uY_u^\dagger
+3Y_dY_d^\dagger
\right).
$$

The symbolic beta matrix is saved as `c5_flavor_beta_matrix.txt`.

---

## 4. Neutrino-Mass Matrix

After electroweak symmetry breaking, the current implementation converts the
coefficient using

$$
m_\nu=-v^2C_5.
$$

The symbolic result is saved in `neutrino_mass_matrix.txt`.

This normalization has been checked directly against the matched T3-B
operator convention. The Python symbolic and numerical stages use the same
conversion.

---

## 5. Optional Numerical Running

When a numerical configuration is supplied, the code evolves

- $g_Y$, $g_2$, and $g_3$;
- $\lambda_H$;
- the diagonal Standard Model Yukawa couplings;
- the complex symmetric matrix $C_5$.

The independent variable is $\ln\mu$, and the coupled one-loop equations are
solved with SciPy's `solve_ivp` using the `DOP853` method.

At the low scale, the code constructs $m_\nu$ and applies a Takagi
factorisation,

$$
U^T m_\nu U=\operatorname{diag}(m_1,m_2,m_3),
$$

to obtain neutrino masses, mass-squared splittings, and $|U_{\rm PMNS}|$.
Normal and inverted mass orderings are supported. The optional data-comparison
stage compares these quantities with the included NuFIT reference data.

---

## 6. Code Structure

### Main symbolic pipeline

| File | Purpose |
|---|---|
| `MatchedEFTRGE.py` | Parses `c5_coefficient.txt`, calculates the one-generation beta function, and checks the SMEFT benchmark. |
| `FlavorMatchedC5.py` | Extracts the T3 loop kernel and constructs the symmetric $3\times3$ flavor coefficient. |
| `SMEFTWeinbergFlavorRGE.py` | Implements the full three-generation Weinberg-operator RGE. |
| `FlavorMatchedRGEStage.py` | Connects the flavor matching and matrix RGE stages. |
| `NeutrinoMassStage.py` | Converts the symbolic coefficient matrix into $m_\nu$. |

### Numerical and phenomenology tools

| File | Purpose |
|---|---|
| `NumericalWeinbergRGE.py` | Numerically evolves the SM parameters and $C_5$. |
| `MaUVRGE.py` | Evolves the one-loop Ma-model gauge, Yukawa, Majorana-mass, scalar-quartic, and scalar-mass parameters above the threshold. |
| `MaFullRunningComparison.py` | Matches the Ma model at a common threshold and separates UV-only, EFT-only, and combined effects. |
| `MaPhysicalBenchmarkFit.py` | Fits $h(\Lambda)$ through the complete UV-to-EFT calculation to a selected low-energy neutrino target. |
| `NumericalPipelineStage.py` | Reads the numerical configuration and writes the high- and low-scale matrices. |
| `RGEReportStage.py` | Builds the concise per-model RGE report, equations, numerical matrices, and observable tables. |
| `NeutrinoObservables.py` | Takagi-factorises $m_\nu$ and calculates observables. |
| `NeutrinoDataComparison.py` | Compares the calculated observables with NuFIT data. |
| `T3NeutrinoTarget.py` | Builds a target neutrino-mass or $C_5$ matrix. |
| `T3YukawaFit.py` | Fits T3 Yukawa matrices to a target coefficient. |
| `T3PhysicalYukawaPoint.py` | Creates a fitted numerical point using the matched loop kernel. |
| `T3RGECorrectedYukawaPoint.py` | Fits at the matching scale after accounting for RGE evolution. |

### General tensor infrastructure

| File | Purpose |
|---|---|
| `GeneralWeinbergRGEGenerator.py` | Implements the general $\psi^2\phi^2$ master-equation terms in a real-scalar basis. |
| `GeneralT3WeinbergRGEStage.py` | Assembles exported T3 tensors, constructs the Weinberg tensor, and checks the representation-generic one-loop beta function component by component. |
| `T3RGETensors.py` | Constructs and validates the T3 scalar, fermion, generator, and coupling tensors. |
| `T3YukawaAdapter.py` | Converts exported Wolfram Yukawa invariants into real-component tensors. |
| `WeinbergWilsonAdapter.py` | Embeds $C_5$ into the general Wilson tensor and checks its symmetries. |
| `T3RGETensorExport.wl` | Exports exact T3 tensor data from Wolfram/Matchete. |
| `GeneralRGE.wl` | Provides Wolfram-side foundations only; automatic RGBeta model definition and beta extraction are not yet implemented. |

---

## 7. Running the Pipeline

### Symbolic RGE

```powershell
python pipeline.py --dims 3 5 4 --alpha 0
```

After successful construction and matching, this automatically runs the
one-generation RGE, full-flavor RGE, and symbolic neutrino-mass stages.

### Numerical RGE

```powershell
python pipeline.py --dims 3 5 4 --alpha 0 --numerical path\to\config.json
```

### Include detailed diagnostics

```powershell
python pipeline.py --dims 3 5 4 --alpha 0 --debug-reports
```

### Full Ma UV-to-EFT comparison

Place `MaUVRGE.py` and `MaFullRunningComparison.py` under `RGE/running/`, then
run from the project root:

```powershell
python -m RGE.running.MaFullRunningComparison examples\ma_full_running_example.json --output output\ma_full_running_comparison.json
```

The calculation uses one common matching scale and reports four cases:

| Case | UV running | EFT running |
|---|---:|---:|
| `frozen` | No | No |
| `uv_only` | Yes | No |
| `eft_only` | No | Yes |
| `full` | Yes | Yes |

This makes the separate effects on $C_5$, the neutrino mass matrix, the
heavy masses, $h$, $\lambda_5$, and the other UV parameters explicit.

At the common threshold, the code diagonalizes the full complex charged-lepton
Yukawa matrix as $Y_e=U_RD_eV_L^\dagger$ and applies

$$
h\rightarrow hV_L,
\qquad
C_5\rightarrow V_L^TC_5V_L.
$$

The output JSON records $V_L$, the ordered electron--muon--tau Yukawa
eigenvalues, and the diagonalization residual.

The Ma parameters are connected to the generated T3-B Matchete convention by

$$
\lambda_{T3}=-\lambda_5,
\qquad
y_1=y_2=h^*.
$$

With this bridge, the Matchete equal-scalar loop kernel is the negative of the
paper's $f/M$ function, and the two signs cancel in $C_5$. The completely
degenerate result is $C_5=-\lambda_5h^Th/(32\pi^2M)$.

### Fit a physical Ma benchmark

```powershell
python -m RGE.phenomenology.MaPhysicalBenchmarkFit examples\ma_full_running_example.json --output examples\ma_full_running_fitted.json --summary output\ma_full_running_fitted_fit_summary.json --comparison output\ma_full_running_fitted_comparison.json --ordering NO --m-lightest 0.001
```

The fitter varies only the high-scale scotogenic Yukawa matrix. Every objective
evaluation performs the full Ma UV running, threshold matching, charged-lepton
basis rotation, and SMEFT running. The summary includes the low-energy NuFIT
comparison and a 64-point UV-trajectory check of perturbativity, boundedness
from below, and positive inert-scalar squared masses.

Because the NuFIT central values are used as the fitting target, a very small
diagnostic chi-square demonstrates numerical reconstruction rather than a
model prediction. The lightest mass and CP phases remain benchmark inputs.

### Generate the full-running comparison report

```powershell
python -m RGE.phenomenology.MaRunningReport examples\ma_full_running_fitted.json --output-dir output\ma_running_report
```

This samples the UV and EFT solutions once each and writes machine-readable
CSV trajectories, a four-case comparison table, the complete comparison JSON,
and publication-ready PNG/PDF figures. The report retains the fit-target
provenance and states the common-threshold and small-$\lambda_5$ limitations.

---

## 8. Outputs

Each successful model has the following simplified layout:

```text
output/<model>/
├── c5_coefficient.pdf
├── c5_coefficient.tex
├── rge_report.pdf
├── rge_report.tex
├── data/
│   └── machine-readable RGE results
└── debug/                  # only with --debug-reports
    └── intermediate algebra and logs
```

`rge_report.pdf` is the main result to read. The expanded symbolic matrices
remain available under `data/` without filling the report with very long
expressions.

### Normal symbolic outputs

| File | Meaning |
|---|---|
| `data/c5_coefficient.txt` | Matched Weinberg coefficient; main downstream input. |
| `c5_coefficient.tex` | Standalone LaTeX version of the coefficient. |
| `c5_coefficient.pdf` | Rendered coefficient for inspection. |
| `rge_report.tex` | Editable source for the complete human-readable RGE report. |
| `rge_report.pdf` | Rendered equations, checks, conventions, and numerical tables. |
| `data/c5_beta.txt` | One-generation symbolic beta function. |
| `data/c5_flavor_beta_matrix.txt` | Full symbolic $3\times3$ beta matrix. |
| `data/neutrino_mass_matrix.txt` | Symbolic neutrino-mass matrix. |

Small JSON summary files record stage status and output filenames.

### Debug-only outputs

| File | Meaning |
|---|---|
| `debug/c5_beta_over_c5.txt` | One-generation beta function divided by $C_5$. |
| `debug/c5_flavor_matrix.txt` | Expanded symbolic flavor coefficient matrix. |
| `debug/c5_loop_kernel.txt` | Flavor-independent T3 loop factor. |
| `debug/wolfram_stdout.log` | Captured Wolfram output. |
| `debug/wolfram_stderr.log` | Captured Wolfram errors and warnings. |

These files are useful for checking a derivation, but no later pipeline stage
requires them.

### Additional numerical outputs

| File | Meaning |
|---|---|
| `data/c5_flavor_matrix_numeric.txt` | Numerical coefficient at the matching scale. |
| `data/c5_flavor_matrix_low_scale.txt` | Numerically evolved coefficient. |
| `data/neutrino_mass_matrix_low_scale.txt` | Low-scale neutrino-mass matrix. |
| `data/neutrino_observables.json` | Masses, splittings, mixing matrix, and numerical checks. |

The standalone Ma full-running report additionally writes
`ma_case_comparison.csv`, `ma_uv_trajectory.csv`, `ma_eft_trajectory.csv`,
`ma_running_report.md`, and three figures in both PNG and PDF formats.

---

## 9. Present Limitations and Checks

- The evolution below the matching scale is SMEFT only; threshold splitting
  between non-degenerate heavy particles is not yet implemented.
- The Ma UV-to-EFT benchmark uses the small-$\lambda_5$ matching expression
  and one common matching scale. It is not yet the general T3 UV runner.
- The flavor lift assumes a diagonal heavy-fermion mass basis and common
  scalar masses.
- The numerical SM running currently uses diagonal Yukawa matrices.
- `GeneralRGE.wl` is not yet a complete automatic RGE calculator.

For ordinary use, inspect `c5_coefficient.pdf` and `rge_report.pdf`. Open the
files under `data/` only when the expanded expressions or numerical matrices
are needed. Enable debug reports only when checking the intermediate algebra.

---

## 10. Pipeline Detailed for RGE

1. After our RGE pipeline, we obtain the weinberg coefficients from c5_coefficient.txt created by RunMatching.wl
2. Pipeline.py runs the run_uv_rgbeta_stage() which calls run_rgbeta_t3() from RGBetaT3Running.py. We give the model parameters and get the UV thoery RGEs
3. RGBetaT3Running.py uses runner RunT3RGBeta.wl which creates our UV model parameter beta functions and  is given as a JSON in uv_rgbeta_rge.json
4. After UV RGE, pipeline.py runs EFT RGE with run_matched_eft_rge_stage() which gives c5_coefficient.txt to run_matched_eft_rge() in MatchedEFTRGE.py
5. MatchedEFTRGE.py uses parse_matchete_c5() to convert Matchete coefficient into SymPy for symbolic python expressions
6. MatchedEFTRGE.py constructs one generation SMEFT for checking implementation
7. Matched coefficient $C_5$ becomes general Wilson tensor $C_{ijab}$ for our generic RGE. MatchedEFTRGE.py selects one component and gets the RGE from it.
8. We have RGEModel.py for scalar representation, GaugeGenerators.py for SU(2) and U(1) generator matrix, and RGECommon.py for other indices and objects.
9. MasterWeinbergRGE.py collates every term in the RGE together 
10. AnomalousDimensions.py gives scalar and fermion collinear anomalous dimensions and put into the master-equation. Then we call calculate_master_rge() giving us the full one loop-equation
11. MatchedEFTRGE.py sums all these contributions and divide by the Wilson coefficient compoennt, now we try to confirm if our one generation works $$\frac{16\pi^2\beta_{C_5}}{C_5}=-3g^2_2+2\lambda_H+6|y_u|^2+6|y_d|^2-|y_e|^2$$, and we write this to c5_beta.txt to check
12. Now we use run_flavor_rge_stage() which calls run_flavor_matched_rge() in FlavorMatchedRGEStage.py which tries to get all 3 generation Weinberg Coeffciients
13. With 3 generations, we use FlavorMatchedRGEStage.py to call flavor_match_from_c5_file() from FlavorMatchedC5.py. We separate our flavour dependent and flavour independent part of our Coefficient, splitting our one generation into 
$C_5=F_{loop}y^*_1y^*_2$
We can use $F_{loop}(Masses + ScalarCoupling)$ from our one generation coefficient and make it general for 3 generations
14. Now in FlavorMatchedC5.py we have 3 lepton generations where we make $F_loop$ change with heavy fermion mass, giving us the symmetric flavor matrix
$$(C_5)_{pq}=\frac{1}{2}\sum_{r}F_r[y^*_{1pr}y^*_{2pr}+y^*_{2pr}y^*_{1qr}]$$
15. Then FlavorMatchedRGEStage.py makes our symbolic $3\times3$ SM Yukawa matrices $Y_e,Y_u,Y_d$ and we also put the $C_5$ matrix to beta_weinberg_matrix() in SMEFTWeinbergFlavorRGE.py
16. SMEFTWeinbergFlavor.py gets the whole Weinberg RGE with all 3 flavors and FlavorMatchedRGEStage.py saves this matrix in c5_flavor_beta_matrix.txt. We have all the symbolic RGE for all flavor components
17. pipeline.py calls run_symbolic_neutrino_mass_stage() which uses NeutrinoMassStage which uses our $C_5$  matrix to get the Majorana neutrino mass matrix with $m_\nu=-v^2C_5$
18. If the user uses --numerical and CONFIG.json, the pipeline.py uses run_numerical_rge_stage() which gives our initial conditions to NumericalPipelineStage.py
19. NumericalPipelineStage.py reads our T3 masses, scalar coupling, T3 Yukawa matrices and SM parameters form JSON file and gets symbolic T3 loop kernel for each heavy-fermion mass and gets a numerical complex symmetric matrix $C_5(M)$. This, with our couplings g_{Y,2,3},$\lambda$ and diagonal SM yukawa is put into SMInitialConditions object and put into evolve_weinberg() in NumericalWeinbergRGE.py
20. NumericalWeinbergRGE.py has the actual coupled differential equations which evolves our SM gauge couplings, the diagonal Yukawa couplings and the $C_5$ matrix. We use $t=ln\mu$. Then we pack real SM parameters with the real and imaginary parts of $C_5$ into a real ODE vector where _beta() calculates one-loop derivatives, where we use solve_ivp to integrate from $ln(M)$ to $ln(\mu)$ (Uses Runge-Kutta integration method). This gives us $C_5(M)\to C_5(\mu)$ with a low-scale coefficient. 
21. NumericalPipelineStage.py gives us our low-sclae neutrino-mass matrix $m_\nu(\mu)=-v^2C_5(\mu)$
22. We send our low-scale neutrino mass matrix from pipeline.py to NeutrinoObservables.py and we diagonalise it to obtain our Majorana mass matrix to get $U^Tm_\nu U=diag(m_1,m_2,m_3)$. This gives us our normal or invertex ordering and gives us our squared differences, our masses and PMNS matrix, and this si saved in neutrino_observables.json. THen we get summaries

