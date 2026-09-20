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
- The UV RGBeta stage and the EFT1 threshold-running stage are handled by
  `RGE/running/rgbeta/` and `RGE/running/eft1/`, respectively.
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

`MatchedWeinbergRGE.py` reads the Matchete expression and computes its SMEFT beta
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
`RGE/running/weinberg/WeinbergRunning.py` is

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

After electroweak symmetry breaking, the code uses the same convention bridge
throughout the symbolic and numerical stages. With the $v\simeq174$ GeV
neutral-Higgs convention this is

$$
m_\nu=-v^2C_5.
$$

Equivalently, with the $v\simeq246$ GeV electroweak convention used by the
numerical input, the implementation evaluates $m_\nu=-(v^2/2)C_5$. The two
forms are the same convention expressed with $v_{246}/\sqrt2=v_{174}$.

The symbolic result is saved in `neutrino_mass_matrix.txt`.

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

### Matching and final Weinberg stages

| File | Purpose |
|---|---|
| `RGE/matching/MatchedWeinbergRGE.py` | Parses the matched $C_5$, evaluates the general $\psi^2\phi^2$ one-loop RGE, and checks the one-generation SMEFT benchmark. |
| `RGE/matching/FlavorC5Matching.py` | Lifts the matched coefficient to the symmetric three-generation flavor matrix. |
| `RGE/matching/WeinbergTensorAdapter.py` | Embeds $C_5$ into the general Wilson-tensor representation used by the master RGE. |
| `RGE/running/weinberg/FlavorMatchedWeinbergStage.py` | Runs the symbolic full-flavor beta matrix and constructs the symbolic neutrino-mass matrix. |
| `RGE/running/weinberg/WeinbergRunning.py` | Contains the symbolic three-generation Weinberg RGE and the numerical coupled SM+$C_5$ evolution. |
| `RGE/running/weinberg/FinalWeinbergCoefficient.py` | Combines hard matching and threshold-running pieces into the final Weinberg coefficient. |
| `RGE/running/weinberg/NumericalWeinbergStage.py` | Evaluates the matched coefficient numerically and evolves it to the requested low scale. |
| `RGE/phenomenology/NeutrinoObservables.py` | Takagi-factorises the low-scale Majorana mass matrix and writes neutrino observables. |

### UV and intermediate-EFT running

| File | Purpose |
|---|---|
| `RGE/running/rgbeta/RGBetaT3Running.py` | Python interface to the UV and EFT1 RGBeta calculations. |
| `RGE/running/rgbeta/RunT3RGBeta.wl` | Runs RGBeta for the full UV T3 theory. |
| `RGE/running/rgbeta/RunT3EFT1RGBeta.wl` | Runs RGBeta for the first intermediate EFT after the first threshold. |
| `RGE/running/rgbeta/T3RGBetaModel.wl` | Builds the RGBeta model definition from the exported T3 data. |
| `RGE/running/rgbeta/T3RGBetaRunnerCommon.wl` | Shared Wolfram runner utilities used by the UV and EFT1 RGBeta stages. |
| `RGE/running/eft1/EFT1TensorAdapters.py` | Builds the EFT1 Wilson and quartic tensors from exported Matchete data. |
| `RGE/running/eft1/EFT1WilsonRGE.py` | Evaluates the EFT1 Wilson-coefficient RGE with the general tensor machinery. |
| `RGE/running/eft1/EFT1WilsonFlow.py` | Performs EFT1 Wilson transport and exports the flavor seed used at the next threshold. |
| `RGE/running/eft1/EFT1DirectWeinberg.py` | Builds and exports the direct one-loop $C_{12}\to C_5$ running contribution. |
| `RGE/running/eft1/EFT1ThresholdResume.py` | Resumes the second threshold with the transported EFT1 result. |
| `RGE/running/eft1/MatcheteParsing.py` | Parses Matchete CG registries and matching expressions used by the EFT1 stages. |

### General tensor infrastructure

| File | Purpose |
|---|---|
| `RGE/general/RGEModel.py` | Defines the scalar representation/basis model used by the general RGE machinery. |
| `RGE/general/GaugeGenerators.py` | Builds the real-scalar SU(2) and U(1) generators and gauge sectors. |
| `RGE/general/FermionBasis.py` | Builds the Weyl-fermion basis and its gauge generators. |
| `RGE/general/WilsonTensorRGE.py` | Implements the general $\psi^2\phi^2$ master-equation tensor terms. |
| `RGE/general/AnomalousDimensions.py` | Adds the scalar and fermion collinear anomalous-dimension contributions and assembles the complete master RGE. |

### Group factors and independent checks

The production group-factor code is split into `RGE/group_factors/core/`,
`RGE/group_factors/recoupling/`, and `RGE/group_factors/validation/`.
Independent research checks remain under `tests/`, including
`tests/python/check_t3_scalar_quartics.py`, `tests/wolfram/ProbeT3SU2Factors.wl`,
`tests/reference/MaUVRGE.py`, and the three independent T3 RGE exporters in
`tests/reference/`.

---

## 7. Running the Pipeline

### Ordinary three-field T3 model

```powershell
python pipeline.py --dims 2 2 1 --alpha -1
```

The three dimensions are $d_{S_1}$, $d_{S_2}$, and $d_F$. For the ordinary
T3 pipeline the threshold order is the heavy fermion first, followed by
$S_1/S_2$.

### Shared-scalar mode

```powershell
python pipeline.py --dims 2 1
```

Two dimensions select one physical shared scalar $S$ and the fermion $F$.
Internally the formal topology uses $S_1=i\sigma_2S^*$ and $S_2=S$, with
$\alpha=-1$, while the physical scalar is counted once. The threshold order
is $F$ followed by $S$.

### Numerical RGE

```powershell
python pipeline.py --dims 2 2 1 --alpha -1 --numerical path\to\config.json
```

When `--numerical` is supplied, `NumericalWeinbergStage.py` evaluates the
matched coefficient at the configured matching scale, `WeinbergRunning.py`
evolves the SM parameters and $C_5$, and `NeutrinoObservables.py` calculates
the low-scale neutrino observables.

### Include detailed diagnostics

```powershell
python pipeline.py --dims 2 2 1 --alpha -1 --debug-reports
```

`--debug-reports` writes expanded intermediate algebra and logs in addition
to the normal report products.

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

- RGBeta is currently used only for the SU(2) representations supported by
  the present Wolfram/RGBeta model setup (fundamental and adjoint).
- Sequential threshold running is implemented for the supported T3 threshold
  plans. The intermediate EFT1 calculation should therefore be distinguished
  from the final SMEFT Weinberg running below the last heavy threshold.
- The flavor lift assumes a diagonal heavy-fermion mass basis. The legacy
  scalar-$C_5$ path also assumes common scalar masses; the hierarchical final
  Weinberg JSON path keeps the hard and running pieces separated before the
  flavor matrix is constructed.
- The numerical SM running in `WeinbergRunning.py` uses diagonal SM Yukawa
  eigenvalues and evolves a complex symmetric $C_5$ matrix.
- `tests/reference/MaUVRGE.py` is retained as an independent Ma/scotogenic RGE
  reference implementation; it is not a production pipeline module.
- The independent SU(2), recoupling, RGBeta, and Weinberg regressions under
  `tests/` should be kept separate from production code even when they duplicate
  part of a calculation, because they provide cross-checks rather than wrappers.

For ordinary use, inspect `c5_coefficient.pdf` and `rge_report.pdf`. Open the
files under `data/` only when the expanded expressions or numerical matrices
are needed. Enable debug reports only when checking the intermediate algebra.

---

## 10. Pipeline Detailed for RGE

1. `Lagrangian/RunMatching.wl` produces the matched Weinberg coefficient and
   the pipeline organises it as `data/c5_coefficient.txt` (or the corresponding
   final-Weinberg JSON after hierarchical threshold transport).
2. `pipeline.py` calls `run_uv_rgbeta_stage()`, which uses
   `RGE/running/rgbeta/RGBetaT3Running.py` and `RunT3RGBeta.wl` to obtain the
   UV-theory beta functions.
3. When the threshold plan contains an intermediate EFT after integrating out
   the heavy fermion, the pipeline runs the EFT1 RGBeta/Wilson stages. The
   relevant code is in `RGE/running/rgbeta/` and `RGE/running/eft1/`.
4. `EFT1TensorAdapters.py` converts exported matching data into the Wilson and
   quartic tensors needed by `EFT1WilsonRGE.py`.
5. `EFT1WilsonFlow.py` transports the EFT1 Wilson coefficients between the
   first and second thresholds and exports the flavor seed for the resumed
   threshold calculation.
6. `EFT1DirectWeinberg.py` carries the direct one-loop $C_{12}\to C_5$
   contribution separately. At fixed one-loop order, the heavy LLSS
   self-running insertion is not fed back through the scalar loop because that
   would be an $O(\hbar^2)$ effect.
7. `EFT1ThresholdResume.py` supplies the transported EFT1 result to the second
   threshold calculation. `FinalWeinbergCoefficient.py` then combines the hard
   threshold contribution and running contribution into the final physical
   Weinberg coefficient.
8. `pipeline.py` runs `run_matched_weinberg_rge_stage()`, which calls
   `RGE/matching/MatchedWeinbergRGE.py`. This parses the matched coefficient into
   SymPy, embeds it in the general Wilson tensor through
   `WeinbergTensorAdapter.py`, and evaluates the generic master RGE.
9. The general tensor calculation uses `RGEModel.py`, `GaugeGenerators.py`,
   `FermionBasis.py`, `WilsonTensorRGE.py`, and `AnomalousDimensions.py`.
10. The one-generation result is checked against

$$
\frac{16\pi^2\beta_{C_5}}{C_5}
=-3g_2^2+2\lambda_H+6|y_u|^2+6|y_d|^2-|y_e|^2.
$$

11. `run_symbolic_flavor_weinberg_stage()` calls
    `RGE/running/weinberg/FlavorMatchedWeinbergStage.py`. It obtains the symmetric
    three-generation $C_5$ matrix from `FlavorC5Matching.py`, constructs symbolic
    $Y_e$, $Y_u$, and $Y_d$, and evaluates `beta_weinberg_matrix()` from
    `WeinbergRunning.py`.
12. The symbolic full-flavor beta matrix is written to
    `data/c5_flavor_beta_matrix.txt`.
13. `run_symbolic_neutrino_mass_stage()` calls the neutrino-mass stage now
    contained in `FlavorMatchedWeinbergStage.py`. It converts the matched symmetric
    coefficient to the Majorana mass matrix using the project convention bridge
    described in Section 4 and writes `data/neutrino_mass_matrix.txt`.
14. If `--numerical CONFIG.json` is supplied,
    `RGE/running/weinberg/NumericalWeinbergStage.py` evaluates the matched
    coefficient and builds `SMInitialConditions` for `WeinbergRunning.py`.
15. `evolve_weinberg()` integrates the one-loop SM gauge couplings,
    $\lambda_H$, diagonal SM Yukawa eigenvalues, and the complex symmetric
    $C_5$ matrix using $t=\ln\mu$ and SciPy `solve_ivp` with `DOP853`.
16. The low-scale coefficient and neutrino-mass matrix are written under
    `data/`. `RGE/phenomenology/NeutrinoObservables.py` then performs the Takagi
    factorisation and writes the neutrino masses, mass-squared splittings, and
    mixing information to `data/neutrino_observables.json`.





From bottom up
1. RGEModel.py takes our scalar fields and puts their indices into scalar blocks where we are in the real basis, and we get our scalar models
2. Gauge generators uses these indices and builds gauge generators for these blocks,
3. WilsonTensorRGE.py and AnomalousDimensions.py uses these for beta functions