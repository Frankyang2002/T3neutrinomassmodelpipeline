# T3 six-model self-bilinear anomalous HH matching — candidate numerical batch

## Install

Extract the ZIP into `C:\thesistemp\T3neutrinomassmodelpipeline`, preserving these paths:

- `Numerical/experimental/T3SelfBilinearAnomalousKernel.py` — NEW spectral one-loop and source tensor implementation
- `scripts/RunT3SelfBilinearAnomalousC5.py` — NEW numerical six-model and 16-model summary runner
- `tests/test_t3_self_bilinear_anomalous_c5.py` — NEW source contraction / numerical regression tests
- `docs/T3_SELF_BILINEAR_ANOMALOUS_HH_C5_BATCH.md` — this document

All files are new. No overwrite of the project README or previous physics outputs.

Requires the **already completed** `continuous_v1`, `neutrino_v1`, and
`self_bilinear_v1` runs (including the fixed 16-model summary). It reads the
existing native tensor sources:

- `output/general_t3_stage22/matchete_tensors.json`
- `output/group_factors/T3_A_uv_cg_registry.json`
- `output/group_factors/T3_D_uv_cg_registry.json`
- `output/group_factors/T3_E_uv_cg_registry.json`

No Wolfram, UV/EFT1 RGE rerun, refit or generation of matching artifacts is needed.

## PowerShell, from project root

```powershell
python -X utf8 -m pytest -q tests/test_t3_self_bilinear_anomalous_c5.py

python -m scripts.RunT3SelfBilinearAnomalousC5 `
    --all `
    --run-id continuous_v1
```

On interrupted execution:

```powershell
python -m scripts.RunT3SelfBilinearAnomalousC5 `
    --all `
    --run-id continuous_v1 `
    --resume
```

If all six model files exist and you only need to recombine the 16 models:

```powershell
python -m scripts.RunT3SelfBilinearAnomalousC5 `
    --all `
    --run-id continuous_v1 `
    --combine-only
```

The new output directory is:

`output/studies/general_scalar_extended/continuous_full_scalar/continuous_v1/self_anomalous_v1/`

- `A_p0_anomalous_hh_c5.json` etc: one complete model result each
- `summary.json`: progressive six-model statuses
- `combined_16_summary.csv` and `.json`: six new + ten existing candidate results

All previous numerical reports remain unchanged.

## Implemented one-loop contribution

For the real singlet/triplet self-bilinear quadratic potential

\[
V_2=m_S^2S^\dagger S+\mu_S^2(S^\dagger)^TC(S^\dagger)+\mathrm{h.c.},
\qquad m_\pm^2=m_S^2\pm2|\mu_S^2|,
\]

we use the **saved and already validated** normal-channel calculation

\[
C_{5,N}=\tfrac12(C_5^{\lambda_T3}(m_+^2)+C_5^{\lambda_T3}(m_-^2)).
\]

There is exactly one extra Stage22 quartic in each of these six models with
two conjugate Higgs fields. Its Hermitian conjugate contains the required
\(HH\) fields. The source-matched quartic names are:

- alpha=0, A(+0)/D(+0)/E(+0): `lambdaT3ExtraQ021010I1` / `lambdaT3GeneralQ021010I1`;
- alpha=-2, A(-2)/D(-2)/E(-2): `lambdaT3ExtraQ020101I1` / `lambdaT3GeneralQ020101I1`.

The source Matchete tensor is conjugated and compared to the existing
`T3MixCG`. The scalar anomalous line is contracted with the **native
Matchete real Cartesian metric** \(\delta_{ab}\), with the self-quartic
source tensor used as an independent metric check. The SU(2) tensor ratio is
computed, not assigned. The same class CG source is used for alpha=0 and
alpha=-2; the SU(2) representation is unchanged, but the cross-package
component convention is not certified.

With \(\lambda_M=s\lambda_R\), \(s=-1\) by the explicit default potential
sign hypothesis, define \(U_\pm=C_5^{\lambda_T3=1}(m_\pm^2)\). Then

\[
 C_{5,A}=\frac{\kappa}{2}\,\lambda_M^*\,P_\alpha(\mu_S^2)
  (U_+-U_-),\quad
 P_0(\mu)=\frac{\mu^*}{|\mu|},\quad
 P_{-2}(\mu)=\frac{\mu}{|\mu|}.
\]

\(\kappa\) is the checked native SU(2) CG ratio (expected to equal one
for the saved sources). The result is **exact for the propagator spectral
partial fraction and the assumed source-to-candidate-kernel correspondence**.
It reuses the original matched topology's numerical kernel, rather than
postulating a different loop integral. It vanishes when \(\mu\) or the extra
quartic vanishes. The runner verifies that setting `lambdaT3=1` really
reconstructs the previously saved normal-source \(C_5\), so the method cannot
silently substitute a nonlinear matching kernel.

Finally \(C_{5,N}+C_{5,A}\) is evolved with the existing full-flavour
SM+Weinberg RGE from 1 TeV to 100 GeV, and masses/mixings are recalculated.

## Limitations (important)

1. **The RGBeta-to-Matchete extra coupling normalisation is not certified.**
   `--extra-sign -1` is a declared default hypothesis (as for earlier B/C
   source-oriented explorations), *not* an inferred or fitted value. `--extra-sign 1`
   or `0` is available for independent conditional scenarios; supply a distinct
   `--output-id` for each scenario. These are not new physical fits.
2. The RGBeta quadratic mass basis uses a spherical metric, while source
   Matchete tensors use Cartesian triplet axes. The numerical mass-squared
   eigenvalues agree, but the complete explicit CG/bilinear coefficient map
   across packages is not independently established. The calculation uses
   the native Matchete metric on the source tensor side.
3. The original one-loop matching kernel is a model-specific **existing
   hierarchical/legacy source approximation**. Spectral replacement extends
   that formula; it does **not** produce new, fully renormalised sequential
   fermion-first/scalar-threshold Wilson matching.
4. There are no additional one-vertex HH quartics beyond the one identified
   here among the six Stage22 catalogues, but higher-order effects and
   threshold operator mixing are not certified. Thus the output status
   says `Candidate`, not `PhysicalFullGeneralScalarMatched`.
5. All numbers are conditional on the saved 1%-quadratic and Stage29
   extra-quartic benchmark choices. No re-fitting of the Yukawa sector.

The intended next piece after this numerical batch is full sequential
fermion-first Wilson matching with verified cross-package conventions.
