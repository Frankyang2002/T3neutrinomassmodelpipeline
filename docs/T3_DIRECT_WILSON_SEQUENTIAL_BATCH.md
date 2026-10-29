# Batch: direct Weinberg Wilson transport on the continuous general-scalar EFT1 background

## Files and placement

This archive adds **new files** (no existing repository files are replaced):

| ZIP path | Repository path |
|---|---|
| `scripts/RunT3DirectWilsonSequential.py` | `scripts/RunT3DirectWilsonSequential.py` |
| `tests/test_t3_direct_wilson_sequential.py` | `tests/test_t3_direct_wilson_sequential.py` |
| `docs/T3_DIRECT_WILSON_SEQUENTIAL_BATCH.md` | `docs/T3_DIRECT_WILSON_SEQUENTIAL_BATCH.md` |

Install by extracting at the root of `C:\thesistemp\T3neutrinomassmodelpipeline`.

## Purpose and exact physics scope

This stage implements **one strictly source-derived part** of the sequential matching chain:

\[
\begin{aligned}
K(\mu_F)&=\sum_r \left[K_{\mathrm{hard},r}^{\overline{\rm MS}}(M_{F,r})+\beta_{K,r}^{\mathrm{direct}}(\lambda_F)\ln\frac{\mu_F}{M_{F,r}}\right],\\
K(\mu_S)&=K(\mu_F)+\int_{\ln\mu_F}^{\ln\mu_S}dt\,
\beta_{K,\mathrm{direct}}(t),\\
K(\mu_{\rm low})&=\mathrm{SM+Weinberg\,RGE}[K(\mu_S)],\\
M_\nu(\mu_{\rm low})&=-\frac{v^2}{2}K(\mu_{\rm low}).
\end{aligned}
\]

- `K_hard` is parsed from the saved **model-specific hierarchical one-loop MSbar hard matching** JSON at each heavy eigenmass `MF_r`. The finite `ln(mu_F/MF_r)` relocation is included in the grouped F-threshold boundary at `mu_F=1e5 GeV` so that later Wilson running does not lose or double-count the heavy-mass logarithms.
- `beta_direct` is obtained by differentiating the matching JSON's **ordered full-flavour running kernel** with respect to `ln(MS)` and checked against the separately saved `eft1_wilson_rge.json` one-generation `beta_kappa_16pi2`.
- The same model's UV-evolved heavy-fermion Yukawas and masses from the saved `continuous_trajectory.json` are rotated into the physical heavy mass basis and *held fixed in EFT1*. That is appropriate for the source term built from the tree-level heavy-fermion Wilson seed, but this is **not** full Wilson tensor evolution.
- The additional scalar quartics and quadratic coupling fully backreact on the EFT1 scalar RG evolution through the existing `RunT3FullScalarContinuous.coupled_rhs`. `lambdaT3(t)` is sampled and Simpson-integrated against the fixed one-loop direct source kernel.
- A **second EFT1 scalar integration** starting from the saved fermion-threshold state is necessary because the previous trajectory saved full additional quartics but not intermediate `lambdaT3` samples. The **UV stage is not rerun**. The new EFT1 endpoint is compared with the saved old endpoint (masses, topology coupling, all additional quartics, bilinear).
- Frozen-`lambdaT3` hard+log results, **including the separate heavy-mass-to-grouped-boundary logs**, are numerically compared with the original `MatchedKernel.evaluate` for each model. A failed parity check rejects the model.
- SM+Weinberg running to 100 GeV uses the existing numerical pipeline and reports 3 Majorana masses and five oscillation observables.

## Available model sources

The **12 models** with existing saved hierarchical hard boundary *and* direct EFT1 Wilson RGE artifacts are:

- A(+0), A(-2)
- B(+1), B(-1)
- C(+1), C(-1)
- D(+0), D(+2), D(-2)
- E(+0), E(+2), E(-2)

The following **four** have legacy one-scale `c5_coefficient.txt` but no saved hierarchical MSbar hard threshold at the inspected repository paths:

- A(-4), B(-3), C(-3), E(-4)

Those four get explicitly marked `LegacyOnlyNoSavedHierarchicalHardBoundary` and no new sequential neutrino masses are invented. They remain available in `neutrino_v1` and `self_anomalous_massdep_v1` as **different approximation/scheme outputs**. All 16 rows are included in the summary; 12 have new partial sequential calculations.

## Exclusions and scientific status

**NOT a complete physical general-scalar sequential match**:

1. The complete fermion-threshold `LLSS` operator tower is not evolved; only the source that projects directly onto Weinberg is included.
2. Full homogeneous/mixed EFT1 Wilson anomalous dimensions (and possible extra-quartic mixing source terms) are not included.
3. The six nonzero self-bilinear anomalous/extra-HH one-scale contributions are **not** added to these sequential results. Their original mass-dependent one-scale outputs remain separate, because a consistent mixed-order matching prescription is not yet certified.
4. B/C mixed-scalar one-scale extra-quartic corrections are not included in the new sequential values.
5. Finite one-loop scalar threshold and complete scalar mass-eigenstate matching are not derived here.
6. RGBeta–Matchete coupling basis normalisation is not certified for all general scalar interactions.
7. Grouped thresholds `mu_F=100000 GeV`, `mu_S=1000 GeV` are used; the script records hierarchy diagnostics rather than silently claiming a multi-threshold physical approximation.

This is a **candidate source-only sequential Wilson improvement**, intended to produce rough numerics and to establish an actual time-dependent EFT1 Wilson transport before attempting complete general-scalar matching.

## PowerShell commands

From the project root:

```powershell
python -X utf8 -m pytest -q tests/test_t3_direct_wilson_sequential.py

python -m scripts.RunT3DirectWilsonSequential --all --run-id continuous_v1
```

If interrupted, or when the outputs exist and sources have not changed:

```powershell
python -m scripts.RunT3DirectWilsonSequential --all --run-id continuous_v1 --resume
```

To generate a new independent output without overwriting an earlier run:

```powershell
python -m scripts.RunT3DirectWilsonSequential --all --run-id continuous_v1 --output-id direct_wilson_sequential_v2
```

A single model first (recommended for checking API integration):

```powershell
python -m scripts.RunT3DirectWilsonSequential --model B_-1 --run-id continuous_v1
```

Output directory:

`output/studies/general_scalar_extended/continuous_full_scalar/continuous_v1/direct_wilson_sequential_v1/`

Files:

- `<model_tag>_direct_wilson.json`: per-model hard, direct transport, source fingerprints, frozen-kernel parity, low-scale masses and relevant non-certification flags.
- `summary.json`, `summary.csv`: 16 rows; actual partial-sequential versus legacy-only status unambiguously distinguished.

## Validation

In the artifact-generation environment, Python source and test module compile. Unit/contract tests are executable with stubbed external repository modules. Repository-source tests are skipped in the isolated environment and run in the user's complete repository. The 12 real-source numerical integrations and frozen-kernel parity must be verified using the commands above. No GitHub files were modified.
