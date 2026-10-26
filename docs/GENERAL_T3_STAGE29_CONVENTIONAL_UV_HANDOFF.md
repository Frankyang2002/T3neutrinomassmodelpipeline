# Stage 29 — Convention-labelled UV quartic benchmark handoff

## Scope

Produce a **numerical initial-condition candidate**, **not** RGE evolution, for all
50 scalar-quartic slots (38 new + 12 registered structural couplings) across the
16 ordinary T3 models. The 46 exact conditional signed matrices from Stage 23
are used only after verifying the Stage 26–28 evidence. Nothing is written into
production `Numerical/`, model configuration files, or Wolfram sources.

Current source ownership was inspected on GitHub before creating this script:

- `RGE/general/generated/general_rgbeta_quartic_registry.json`: 38 new vertices,
  12 reused structural slots, 50 total, 16 models.
- `Numerical/core/State.py`: canonical `T3UVState` dataclass.
- `Numerical/core/StateVector.py`: real ODE packing for that state.
- `Numerical/core/RGBetaEvaluator.py`: numerical evaluator with explicit
  UV state symbol environment.
- `docs/PHYSICS.md`: canonical chosen UV scale of `1e7 GeV`,
  `16*pi**2*dX/dln(mu)=beta_X^(1)` and Class-E `lambdaT3` benchmark-only 1/2.

The present numerical ODE does **not** include the 38 newly generated quartics.
Stage 29 inspects local source with Python's AST to record coverage for each
layer without importing or changing the solver.

## Mathematical transformation

For each block the Stage 23 tensor geometry is

\[
R_j=\sum_i M_i A_{ij}.
\]

The conditional Lagrangian sign gives

\[
\lambda_R=-A^{-1}\lambda_M,
\qquad
A\lambda_R+\lambda_M=0.
\]

Stage 29 checks this identity with **exact SymPy rational/algebraic numbers**
for every block before creating floating-point candidate values. It chooses

\[
\lambda_{M,i}(\mu_{UV})=0.01
\]

for **each of the 50 additional scalar-quartic slots** by default. This is a
new explicit *illustrative initial-condition choice*, **not** the historical
production benchmark and **not** a fitted or measured value. For the familiar
B/C `110101` two-dimensional example,

\[
A=\begin{pmatrix}\tfrac12&1\\-\tfrac{\sqrt3}2&0\end{pmatrix},\quad
\lambda_R=\begin{pmatrix}\tfrac{2\sqrt3}3\\-1-\tfrac{\sqrt3}3\end{pmatrix}0.01.
\]

Both the exact expressions and floats are included. The `lambdaT3` topology
coupling is **not part of this 50-slot initialisation**. The Class-E factor 1/2
used for that topology coupling is **not** applied to the additional quartics.

The triplet convention is labelled explicitly as the *Stage 27 Matchete CG
reference basis, provisionally identified with abstract RGBeta `S2`*.
Stage 28 did not independently export a numerical RGBeta component basis or
physical barred-triplet intertwiner, so no physical map is certified.

## Run (PowerShell, repository root)

```powershell
Expand-Archive -LiteralPath .\T3_General_Scalar_Stage29_ConventionUVHandoff.zip `
    -DestinationPath . -Force

python -m pytest -q tests\python\test_general_t3_stage29_convention_uv_handoff.py

python -m scripts.Stage29ConventionUVHandoff --extra-quartic 0.01
```

Expected terminal summary, **conditional on current local evidence passing**:

```
Stage29 convention-labelled numerical UV values: 16/16 models, 50/50 slots, 46/46 blocks
New slots needing production one-loop state/beta support: 38
Largest |lambda_RGBeta| at chosen UV point: <computed number>
Actual RGE integration performed: NO
Physical Matchete/RGBeta coupling mapping: UNCERTIFIED
```

Report written to
`output/general_t3_stage22/stage29_convention_uv_handoff.json`.

The default invocation re-runs Stage 23 and Stage 28 checks, with Stage 28
recursively checking Stages 27 -> 26 -> 25. The diagnostic compares the current
saved Stage 23/28 outputs and hashes all source/evidence inputs to prevent use
of stale mapping data. It must be run with the regenerated local Stage 22–28
JSON evidence, including the successful live Wolfram exports.

For a different **illustrative** coefficient, use for example
`--extra-quartic 0.001`. The value must be a finite decimal in `(0,1]`.

## Output structure

The single JSON file contains 16 model entries, each with:

- `matchete_quartics`: exact and numerical chosen input parameters;
- `rgbeta_quartics`: exact and numerical conditionally mapped values, classified
  as `new` or `reused` in the current general registry;
- `blocks`: the exact signed matrix used and an exact zero residual;
- `production_ode_symbol_coverage`: current local state/vector/evaluator
  coverage of the 38 new and 12 reused slots;
- SHA-256 fingerprints for source and upstream evidence.

`physical_coupling_mapping_certified=false`,
`one_loop_beta_functions_for_new_quartics_evaluated=false`, and
`rge_evolution_performed=false` are deliberate. These are **not** solver inputs
for the current production code until that solver has been explicitly extended.

## Next numerical work

Export one-loop RGBeta beta expressions for the 38 new quartic symbols for each
relevant model and extend the numerical UV state/vector/evaluator to include
those parameters, alongside the 12 reused couplings. Check complete beta-key
coverage and pack/unpack consistency before integrating from `1e7` to `1e5`
GeV. Keep these boundary values labelled convention-dependent until the global
barred-field component identification is independently resolved. Do not silently
set missing beta functions to zero.
