# T3 Weinberg-Coefficient RGE Pipeline

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

---

## 0. Quick Summary

- `data/c5_coefficient.txt` is the main input from matching.
- The symbolic stages run automatically after a successful match.
- Numerical running only occurs when `--numerical CONFIG.json` is supplied.
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

This file is the source of truth for the RGE, flavor, and neutrino-mass
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

> **Convention check:** some function docstrings still say
> $m_\nu=-(v^2/2)C_5$, while both the symbolic and numerical code actually
> use $m_\nu=-v^2C_5$. Before a final phenomenological interpretation, the
> normalization of the extracted Weinberg operator and this conversion must
> be checked together.

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

---

## 9. Present Limitations and Checks

- The evolution below the matching scale is SMEFT only; threshold splitting
  between non-degenerate heavy particles is not yet implemented.
- The flavor lift assumes a diagonal heavy-fermion mass basis and common
  scalar masses.
- The numerical SM running currently uses diagonal Yukawa matrices.
- The $m_\nu$ normalization must be reconciled with the operator convention,
  as noted in Section 4.
- `GeneralRGE.wl` is not yet a complete automatic RGE calculator.

For ordinary use, inspect `c5_coefficient.pdf` and `rge_report.pdf`. Open the
files under `data/` only when the expanded expressions or numerical matrices
are needed. Enable debug reports only when checking the intermediate algebra.
