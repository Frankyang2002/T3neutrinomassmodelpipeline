# Generalised T3 Neutrino-Mass Pipeline

This repository currently implements the model construction and one-loop EFT matching stages for the T3 radiative neutrino-mass topology.

The implemented pipeline is

$$
\text{T3 representations}
\;\longrightarrow\;
\mathcal L_{\rm UV}
\;\longrightarrow\;
\text{one-loop matching}
\;\longrightarrow\;
\mathcal L_{\rm EFT}^{\rm BSM}
\;\longrightarrow\;
C_5 .
$$

The subsequent RGE evolution and conversion of the Weinberg coefficient into the physical neutrino-mass matrix are outside the scope of this README.

---

## 0. Physics Scope

The T3 topology contains three BSM multiplets:

- two complex scalars, $S_1$ and $S_2$;
- one fermion multiplet, $F$.

The code generalises the usual singlet/doublet/triplet T3 models to arbitrary valid $SU(2)_L$ representation dimensions.

The hypercharge convention used by the code is

$$
Q = T_3 + Y.
$$

For the T3 parameter $\alpha$,

$$
Y(S_1)=\frac{\alpha}{2},
\qquad
Y(S_2)=\frac{\alpha+2}{2},
\qquad
Y(F)=\frac{\alpha+1}{2}.
$$

The hypercharge scalar difference comes from our interaction term $HHSS^\dagger$, where an example would be the two scalars with hypercharge $\frac{1}{2}$ and $-\frac{1}{2}$

This differs by a factor of two from papers that define hypercharge through
$Q=T_3+Y/2$.

The historical T3 classes correspond to

| Class | $d_{S_1}$ | $d_{S_2}$ | $d_F$ |
|---|---:|---:|---:|
| T3-A | 1 | 3 | 2 |
| T3-B | 2 | 2 | 1 |
| T3-C | 2 | 2 | 3 |
| T3-D | 3 | 1 | 2 |
| T3-E | 3 | 3 | 2 |

The general representation constraints are

$$
d_{S_1}=d_F\pm1,
\qquad
d_{S_2}=d_F\pm1,
$$

because each scalar must occur in

$$
2\otimes d_F,
$$

and the scalar pair must contain the triplet required to contract with the symmetric Higgs pair,

$$
S_1\otimes S_2^\dagger \supset \mathbf 3.
$$

---

## 1. UV Lagrangian

For a valid T3 model, the builder defines the BSM fields, their gauge representations, and the interactions required by the topology.

Schematically,

$$
\mathcal L_{\rm UV}
=
\mathcal L_{\rm SM}
+
\mathcal L_{\rm free}^{\rm BSM}
+
\mathcal L_{\rm int}.
$$

The important topology-generating interactions are two Yukawa vertices and the scalar interaction that closes the T3 loop,

$$
\mathcal L_{\rm int}
\supset
y_1\,LFS_1
+
y_2\,LFS_2
+
\lambda_{T3}\,HHS_1S_2^\dagger
+
\text{h.c.},
$$

with the exact conjugations and $SU(2)$ contractions chosen according to the gauge representations.

Additional gauge-invariant scalar terms needed for a complete two-scalar UV theory are also constructed by the builder.

### Gauge invariance

Candidate interactions are checked before being accepted.

For $U(1)_Y$, the total hypercharge must vanish,

$$
\sum_i Y_i = 0.
$$

For $SU(2)_L$, the tensor product of the participating representations must contain a singlet,

$$
R_1\otimes R_2\otimes\cdots \supset \mathbf 1.
$$

The project contains a canonical gauge-invariance implementation used by the higher-level T3 checks. The T3 invariance layer additionally supports the model's $Z_2$ assignment.

### Generic Clebsch-Gordan construction

For arbitrary $SU(2)$ dimensions, the builder uses Matchete invariant tensors and dynamically defined Clebsch-Gordan coefficients.

The older T3-A--E contractions are retained as a legacy/reference path. They provide a regression oracle for the generic construction in the singlet, doublet, and triplet cases.

The builder only reports a valid T3 construction when all three required ingredients are present:

1. Yukawa vertex 1;
2. Yukawa vertex 2;
3. the T3 scalar-mixing interaction.

---

## 2. EFT Matching

Once the UV Lagrangian has been successfully constructed, Matchete is used to integrate out the heavy BSM fields at one loop.

The current default settings are

$$
\text{EFT operator dimension}=5,
\qquad
\text{loop order}=1.
$$

The matching flow is

LUV ->Match ->GreensSimplify ->
EOMSimplify ->EvaluateLoopFunctions ->ReplaceEffectiveCouplings ->
Matched EFT


The same matching procedure is also applied to the Standard Model alone.

The BSM contribution is then obtained by subtraction,

$$
\mathcal L_{\rm EFT}^{\rm BSM}
=
\mathcal L_{\rm EFT}^{\rm full}
-
\mathcal L_{\rm EFT}^{\rm SM}.
$$

The result is canonicalised again after subtraction because algebraically equivalent Matchete expressions can otherwise differ only by dummy-index naming.

---

## 3. Weinberg Operator Extraction

The dimension-five EFT is searched for the Weinberg operator,

$$
\mathcal O_5 \sim LLHH.
$$

In the Matchete expression this is identified through the corresponding charge-conjugated fermion structure.

The extraction code separates the holomorphic operator from its Hermitian conjugate. The coefficient $C_5$ is extracted from the holomorphic $P_L$ sector so that the Hermitian-conjugate contribution is not double counted.

Conceptually,

$$
\mathcal L_{\rm EFT}
\supset
C_5\,\mathcal O_5
+
C_5^\dagger\,\mathcal O_5^\dagger.
$$

Operator-presence detection and detailed coefficient extraction are intentionally separate. Therefore, if the Weinberg structure is present but the coefficient parser cannot isolate it, the pipeline reports an extraction problem rather than incorrectly reporting that the operator is absent.

---

## 4. Code Structure

### `pipeline.py`

Python command-line entry point.

It can run either:

- the classfied T3-A--E benchmark models; or
- an arbitrary valid representation assignment $(d_{S_1},d_{S_2},d_F)$.

It also manages output directories, Wolfram subprocesses, logs, and aggregate summaries.

### `T3ModelCatalog.wl`

Defines the T3 model data.

Responsibilities include:

- historical T3-A--E representation dimensions;
- $\alpha$-dependent hypercharges;
- generic model construction from representation dimensions;
- representation-level validity checks.

### `LagrangianBuilder.wl`

Constructs the Matchete UV theory.

Responsibilities include:

- defining arbitrary $SU(2)$ representations when required;
- defining BSM scalar and fermion fields;
- defining couplings;
- constructing Clebsch-Gordan tensors;
- generating candidate Yukawa and scalar interactions;
- checking candidate validity;
- assembling $\mathcal L_{\rm BSM}$ and $\mathcal L_{\rm UV}$.

### `RunMatching.wl`

Performs the EFT calculation.

Responsibilities include:

- one-loop matching;
- Green-basis simplification;
- equation-of-motion simplification;
- loop-function evaluation;
- effective-coupling replacement;
- SM baseline matching and subtraction;
- Weinberg-operator detection;
- $C_5$ extraction.

### `RunModel.wl`

Wolfram-side orchestration layer.

It connects model construction, UV validation, matching, Weinberg extraction, LaTeX conversion, and output generation.

### `PhysicsLaTeX.wl`

Converts Matchete's internal expressions into more readable physics-oriented LaTeX for exported results.

---

## 5. Running the Pipeline

### Historical T3 smoke test

```bash
python pipeline.py --smoke
```

This runs the five validated benchmark points

```text
T3-B  alpha = -1
T3-C  alpha = -1
T3-A  alpha =  0
T3-D  alpha = -2
T3-E  alpha =  0
```

A successful run should report, for each model,

```text
build=Success
match=Success
T3=True
Weinberg=True
```

and successful $C_5$ extraction.

### Generic representation

For example,

```bash
python pipeline.py --dims 3 5 4 --alpha 0
```

runs

$$
d_{S_1}=3,\qquad
d_{S_2}=5,\qquad
d_F=4,\qquad
\alpha=0.
$$

The Python front end first checks that the requested dimensions satisfy the T3 representation conditions before invoking Matchete.



---

## 6. Output

Each run receives its own directory under

```text
wolfram/output/
```

The pipeline writes Matchete stdout/stderr logs and a machine-readable model summary.

When Weinberg extraction succeeds, the important outputs include

```text
comparison_summary.json
c5_raw.txt
c5_raw.tex
c5_coefficient.txt
c5_coefficient.tex
wolfram_stdout.log
wolfram_stderr.log
```

`c5_raw.*` contains the isolated Weinberg-operator contribution before removing the operator structure.

`c5_coefficient.*` contains the extracted coefficient itself.

For multi-model runs, the Python front end also writes an aggregate comparison file,

```text
wolfram/output/t3_model_comparison.json
```

---

## 7. Regression Testing

The UV-to-$C_5$ stage has a dedicated regression gate:

```bash
python regression.py
```

This performs two checks:

1. reruns the five historical smoke models through the full UV construction and matching pipeline;
2. validates the resulting $C_5$ outputs, including the T3-B scotogenic limit.

Existing matching outputs can be reused with

```bash
python regression.py --no-rematch
```

The refactored pipeline has also been tested with the generic higher-dimensional example

```bash
python pipeline.py --dims 3 5 4 --alpha 0
```

to ensure that the arbitrary-representation path remains functional.

---

## 8. Not Covered Here

This README intentionally stops at the extracted Weinberg coefficient.

The following later stages belong to the next part of the project and are not documented here:

- RGE evolution of the EFT coefficients;
- conversion between real-scalar coefficient conventions and $\kappa$;
- electroweak symmetry breaking;
- construction of the neutrino-mass matrix;
- comparison with measured neutrino masses and mixing parameters.
