# Batch 3–4: full seven-term Wilson RGE covariance continuation

**Status:** `COVARIANCE_PASS_UNVERIFIED_LOOP_FORMULA`. All results below are **synthetic tests**. No physical A(-4) EFT1 Wilson coefficient, finite scalar matching or SMEFT `C5` is certified.

**Read-only source revision:** GitHub `main` `12128de602ac6bf109dbad50dbac8446ed629434` (compared with earlier audit, unchanged). Original symbolic/finite-log source files were not changed.

## 1. What is confirmed from the source

The uploaded arXiv:2512.16890v1, Eq. (4.85), printed page 24, contains a crossed-Yukawa term with `y_jla y_ikc C_lkbc`, without a star on `C`. Its full pair-exchange symbol requires separate interchange of `i,j` and `a,b`; the previously inspected repository symbolic and fixed-order numerical implementations instead used two simultaneous exchanges. The existing independent report tested four exchanges but **still** failed a gauge-Ward test for the crossed term.

The journal article is **Aebischer, Bresciani & Selimović, JHEP 09 (2026) 159, published 11 September 2026, DOI 10.1007/JHEP09(2026)159**. We could verify publication metadata but not independently inspect an updated equation from the journal PDF. An independent, distinct general-EFT calculation is Fonseca, Olgoso & Santiago, JHEP 07 (2026) 094 / arXiv:2512.15866; a detailed basis/normalisation-level equality of its formula with Eq. (4.85) has **not** been established. It would be incorrect to claim that these references certify the candidate formula.

## 2. Exact unitary-covariance constraint

With `L5 = C_ijab psi_i psi_j phi_a phi_b + h.c.` and `L_Y = -y_ija psi_i psi_j phi_a /2 + h.c.`, under `psi' = U psi` and `phi' = R phi`,

- `y' = (U* ⊗ U* ⊗ R)y`;
- `C' = (U* ⊗ U* ⊗ R ⊗ R)C`;
- `t' = U t U†`, `theta' = R theta R^T`.

The printed crossed contraction `yyC` is not covariant under a generic complex unitary `U`. The candidate contraction **`yy C*`** is covariant, for both independent external pair swaps. As a trivial analytic check, under a global rephasing `U=e^{iα}I`, `y` and `C` pick up `e^{-2iα}`; `yyC` picks up `e^{-6iα}` (incorrect), while `yyC*` picks up `e^{-2iα}` (correct).

The fermion collinear gauge term is another basis-index issue: the operator-index contraction requires the gauge Casimir `C2(t)^T`, which transforms as `U* C2(t)^T U^T`, rather than the untransposed `C2(t)`. The transpose changes nothing in the canonical SM irreducible blocks, where the Casimir is real and diagonal, but matters after complex unitary rotations between blocks.

These are **necessary covariance requirements**. They do not establish a diagrammatic coefficient, counterterm sign, or an authoritative correction to the published master equation.

## 3. Numerical covariance checks

Using one generic complex unitary rotation of four left-Weyl fields, one orthogonal rotation of three real scalars, complex Yukawa/Wilson tensors, a real-symmetric quartic, and a Hermitian nontrivial gauge generator:

| Term | Original relative covariance residual | Candidate residual |
| --- | ---: | ---: |
| Yukawa wavefunction | 4.32e-16 | 4.32e-16 |
| Scalar pair | 4.42e-16 | 4.42e-16 |
| Mixed Yukawa–gauge | 5.05e-16 | 5.05e-16 |
| Crossed Yukawa | **1.2018** | **4.39e-16** |
| Conjugate-coefficient Yukawa | 3.93e-16 | 3.93e-16 |
| Scalar collinear | 3.31e-16 | 3.31e-16 |
| Fermion collinear | **0.2094** | **3.62e-16** |

The previous physical-selection-rule tests also remain passing: SM-like `SU(2)×U(1)` beta Ward residual `1.94e-16`, and a nontrivial `U(1)` crossed mixing beta Ward residual `9.26e-18` with crossed beta norm `0.077476`.

## 4. Matrix-free ODE diagnostics

The opt-in candidate mode now runs via `evolve(..., beta_mode="ward_covariant_candidate")` or `--beta-mode ward_covariant_candidate`. It retains all seven terms; the crossed and conjugate-coefficient Yukawa terms are both anti-linear in `C`. This means `B(mu)` receives two terms in `A(mu)C+B(mu)C*`. Other terms remain complex-linear. The solver never constructs dense matrices `A` and `B`.

- Synthetic gauge-invariant `U(1)` ODE, `MF=1e5 GeV` to `MS=1e3 GeV`: two zero-boundary symmetric coordinates generated, including nonzero `C[0,0]` from the crossed term; maximum Ward residual along the stored trajectory is zero to reported precision.
- Full-size deliberately **non-gauge-invariant synthetic A(-4)-sized** test, `MF=1e5 GeV` to `MS=2e3 GeV`: 45 Weyl components, 12 real scalars, 80,730 independent symmetric Wilson entries; 36 generated initially-zero coordinates; 44 RHS evaluations. `||C_ODE-C_fixed_order||/||C_initial||=0.00115746` for this synthetic benchmark, which is a resummed-log difference rather than one-loop hard matching.
- Outputs under `output/batch34_covariant/A_m4_full45_synthetic/` retain the full independent trajectory, complete endpoint Wilson tensors, all seven beta tensors at both endpoints, full pair-coordinate inventory, generated modes, and source diagnostics.

## 5. Test status, blockers, and safeguards

Local combined tests: 18/18 covariance, ward, and independent SM reference tests passed. Original seven-term tests: 14 passed, 2 skipped because the upstream repository is not mounted into this evaluation container (they can be run in your local full checkout). New six-test unitary/ODE group: 6/6 passed.

**The revised source and CLI expose two explicit variants:**

- `source_v1_unstarred` (default; follows the four-exchange printed formula, **known to fail general gauge/unitary covariance**).
- `ward_covariant_candidate` (opt-in; crossed uses `C*`, fermion Casimir uses a transpose, **mathematically covariant but still unverified at loop level**).

The exact source-formula conflict is unresolved, and both modes explicitly report `beta_physics_certified=false` and `physical_C5_certified=false`. Physical (`authority="certified"`) input contracts now **fail closed** in this continuation, instead of running a coefficient equation with a known unresolved one-loop convention. Chat A's canonical matching normalization and Chat B's complete running general EFT1 inputs have not been independently verified here. The mandatory EFT sequence remains `UV -> integrate out all heavy fermions -> EFT1(SM+S1+S2) -> integrate out all scalars simultaneously -> SMEFT`.

No code owned by Chat A or Chat B, no original Matchete boundary, no original symbolic master, and no GitHub files were edited.

## 6. PowerShell commands

From the project root:

```powershell
python -m unittest discover -s tests -p "test_batch34_unitary_covariant_evolution.py" -v
python -m unittest discover -s tests -p "test_general_eft1_seven_term_flow.py" -v

# Existing deliberately synthetic 6-Weyl fixture (source and running couplings)
python scripts\RunGeneralEFT1SevenTermFlow.py `
    --boundary tests\fixtures\synthetic_A_m4_boundary.json `
    --couplings tests\fixtures\synthetic_A_m4_couplings.json `
    --output-root output\batch34_candidate_local `
    --allow-synthetic --beta-mode ward_covariant_candidate `
    --source-decomposition --convergence-check --points 11

# The full three-generation scaling test, synthetic only
python scripts\RunBatch34Full45Synthetic.py `
    --output-root output\batch34_candidate_full45_local `
    --beta-mode ward_covariant_candidate
```

The target output directories must be new because the runners refuse to overwrite output files.
