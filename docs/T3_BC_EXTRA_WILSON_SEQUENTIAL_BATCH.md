# B/C(-1) selected extra-quartic Wilson source on continuous full-scalar EFT1 running

## Installation

Extract the ZIP at the repository root:

`C:\thesistemp\T3neutrinomassmodelpipeline`

The archive adds **three new files** and changes no existing file:

| ZIP path / project destination | Role |
|---|---|
| `scripts/RunT3BCExtraWilsonSequential.py` | Numerical runner, source contracts and outputs |
| `tests/test_t3_bc_extra_wilson_sequential.py` | Physics/source/runner contract tests |
| `docs/T3_BC_EXTRA_WILSON_SEQUENTIAL_BATCH.md` | This document |

## Exact physics scope

Adds the **selected additional-quartic scalar-pair source** from the already-derived Stage17 and Stage19 B/C(-1) artifacts to the **saved direct-only fermion-first sequential candidate**. No existing UV trajectories or previous matching outputs are overwritten.

At the grouped fermion boundary, the tree-level PL `LLSS` Wilson seed gives the fixed heavy-flavour contractions

\[
F_{i,pq}=\sum_{r=1}^3\frac{\overline{Y_{i,pr}}\,\overline{Y_{i,qr}}}{M_{F,r}},\qquad i=1,2.
\]

For B(-1) and C(-1), respectively, the previously certified Weinberg-subspace projection gives

\[
16\pi^2\frac{d\,\delta K_{pq}}{d\ln\mu}\bigg|_{\rm extra}=
\begin{cases}
2\left[\overline{\lambda_{HH^\dagger S_1^\dagger S_1^\dagger}}F_{1,pq}
+\overline{\lambda_{HH^\dagger S_2 S_2}}F_{2,pq}\right], & B(-1),\\[3pt]
\dfrac{2\sqrt{3}}{3}\left[\overline{\lambda_{HH^\dagger S_1^\dagger S_1^\dagger}}F_{1,pq}
+\overline{\lambda_{HH^\dagger S_2 S_2}}F_{2,pq}\right], & C(-1).
\end{cases}
\]

These factors are **not re-derived by guessing a group factor**. The runner loads:

- `output/general_t3_stage17/{B_m1,C_m1}_extra_eft1_scalar_pair_beta.json` (four independent mixed-Higgs real components);
- `output/general_t3_stage19/{B_m1,C_m1}_selected_extra_weinberg_transport.json` (full 30-component source SU(2) projection and physical Weinberg normalization certificate);
- the model's saved `eft1_after_F_wilson_seed.json`, `eft1_wilson_rge.json` and `final_weinberg_coefficient.json`.

It verifies **symbolically** that all four current Stage17 components match the saved full-tensor Weinberg projection, compares the saved direct normalization against the current authoritative direct coefficient, and reuses Stage18's actual source integrator plus Stage19's projection for the numerical 3×3 matrices.

The running scalar coefficients, including the two extra B/C Higgs quartics, are sampled on a fresh EFT1 integration starting from the **saved** fermion endpoint. Endpoint parity is checked against `continuous_v1`. This repeats **EFT1 only**, not UV running, and retains full backreaction of the candidate scalar quadratic and quartic beta functions on the scalar background.

The resulting update is

\[
\begin{aligned}
\delta K(\mu_S)&=\int_{\ln\mu_F}^{\ln\mu_S}
   d\ln\mu\,\beta_K^{\rm extra}(\mu),\\
K_{\rm candidate}(\mu_S)&=K_{\rm previously\ saved\ direct}(\mu_S)+\delta K(\mu_S),\\
K_{\rm candidate}(100\,\mathrm{GeV})&=\text{existing SM+Weinberg RGE}(K_{\rm candidate}(\mu_S)).
\end{aligned}
\]

The same SM+Weinberg routine is applied independently to the saved direct-only input to cross-check the baseline masses, and the correction at 100 GeV is saved as the difference. All 3×3 matrices are complex and symmetric.

## Scientific limits — **not full sequential matching**

1. **Only one selected scalar-pair source insertion** from the tree-level PL `LLSS` boundary is included. The other LLSS tensor components are not evolved under their coupled homogeneous/mixing RGEs.
2. The heavy-fermion Yukawas and masses in the tree-level `LLSS` seed are frozen at the fermion threshold; their scalar-EFT Wilson self-running is absent.
3. Other new quartic mixing sources, nonzero bilinear insertions and complete one-loop scalar-threshold matching are not included.
4. The B/C mass-eigenstate one-scale results are a **different approximation/scheme** and are not added to this result. The new output builds strictly on `direct_wilson_sequential_v1`.
5. A saved Stage19 SU(2) and normalisation certificate is reused; to repeat its expensive 30-component full-tensor check with current repository sources, pass `--reverify-tensors` on a **fresh output ID**.
6. The broader Matchete–RGBeta general-scalar conventions remain candidate/provisional. This result is **not a complete physical C5 or an oscillation refit**.
7. The Stage29 illustrative B/C extra quartic boundary values are already **nonzero**; the source is therefore a one-loop Wilson insertion, not necessarily an O(hbar²) effect. Avoid applying the earlier zero-UV-extra-coupling loop-order label to this run.

## Dependencies and execution

The previous successful commands have already produced:

- `output/studies/general_scalar_extended/continuous_full_scalar/continuous_v1/{B_m1,C_m1}/continuous_trajectory.json`;
- `output/studies/general_scalar_extended/continuous_full_scalar/continuous_v1/direct_wilson_sequential_v1/{B_m1,C_m1}_direct_wilson.json`.

Run each command **separately** in PowerShell from the project root:

```powershell
python -X utf8 -m pytest -q tests/test_t3_bc_extra_wilson_sequential.py
```

```powershell
python -m scripts.RunT3BCExtraWilsonSequential --model B_-1 --run-id continuous_v1
```

```powershell
python -m scripts.RunT3BCExtraWilsonSequential --all --run-id continuous_v1 --resume
```

For unchanged files, `--resume` uses saved results. To force a fresh independent run without overwrite:

```powershell
python -m scripts.RunT3BCExtraWilsonSequential --all --run-id continuous_v1 --output-id bc_extra_wilson_sequential_v2
```

If your previous direct-source run used a nondefault output ID or sample count, supply `--direct-output-id <name>` and matching `--samples <odd-integer>`.

### Outputs

Under `output/studies/general_scalar_extended/continuous_full_scalar/continuous_v1/bc_extra_wilson_sequential_v1/`:

- `B_m1_extra_wilson.json`, `C_m1_extra_wilson.json`: input hashes, source and projection certificates, complex 3×3 matrix correction, frozen-quartic control, direct and corrected threshold matrices, 100 GeV results, neutrino masses, and honest completeness flags.
- `summary.csv`, `summary.json`: **two** model rows with relative corrections and mass comparisons.

`|delta C5|/|direct C5|` refers to Frobenius norms at 1 TeV. The difference between the new neutrino masses and the saved direct-only masses can be attributed to the selected extra source because they share the same hierarchical baseline and SM evolution. It is not the full general-scalar correction.

## Validation performed on the artifact-generation environment

- Python compilation succeeded.
- **22 isolated unit/contract tests passed**, and **2 source-dependent tests were skipped** because the complete repository is not mounted here. Those source tests will run on the Windows checkout.
- The full 2-model EFT1 source integration and numerical masses were **not run in this environment**; your Windows checkout is required.
- GitHub was read only; nothing was pushed or altered remotely.
