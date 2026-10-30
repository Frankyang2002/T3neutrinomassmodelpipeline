# Six self-bilinear mass / RGE / one-scale C5 audit

Models: `A_-2`, `A_+0`, `D_-2`, `D_+0`, `E_-2`, `E_+0`. All other T3 models remain unchanged.

**New files**

- `Lagrangian/RunSixSelfBilinearNativeTensor.wl`: extracts each actual Matchete holomorphic mass invariant, including singlet cases; records the tensor and factor `k_native`. No matching operation.
- `scripts/VerifySixT3SelfBilinearBatch.py`: checks the existing 16-model symbolic summaries, UV/EFT1 scalar mass betas, operator names and signs, native tensors, an independent one-loop triangle integral, and both C5 channels.
- `tests/test_six_t3_self_bilinear_batch.py`: pure algebra and certificate-schema regression tests.
- This document.

## Command lines (PowerShell, repository root)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "C:\path\to\T3_six_self_bilinear_batch.zip" -DestinationPath . -Force
python -X utf8 -m pytest -q tests/test_six_t3_self_bilinear_batch.py
python -u -X utf8 -m scripts.VerifySixT3SelfBilinearBatch --model A_-2 --run-id six_self_a2_preflight_v1
python -u -X utf8 -m scripts.VerifySixT3SelfBilinearBatch --all --run-id six_self_all_v1
```

All six native CG probes invoke local Mathematica/Matchete. The earlier successful `all16_v1` symbolic export must exist at `output/studies/general_scalar_extended/complete_scalar_symbolic/all16_v1`. A failed model is reported in the combined summary; a new `--run-id` prevents stale result reuse. The runner does **not** overwrite production Lagrangians, RGBeta source, benchmark configurations or prior certificate files. It may take appreciable time because six Matchete initializations and 36 high-precision loop integrations are performed.

Optional illustrative inputs (independent of fitted couplings):

```powershell
python -u -X utf8 -m scripts.VerifySixT3SelfBilinearBatch --all --run-id six_self_phase2_v1 --mu2-fraction 0.02 --mu2-phase 1.1 --extra-rgbeta-real 0.02 --extra-rgbeta-imag -0.01
```


### Optional: also evaluate the **existing evolved** one-scale candidates

If the already generated `continuous_v1` 16-model trajectories are present locally, the same batch can additionally invoke the repository's existing `scripts.RunT3SelfBilinearMassDependentC5.calculate` for each of the six models. That reuses their evolved UV/EFT1 parameters and low-energy SM+Weinberg stage; it does not duplicate or rederive that code. The separate *illustrative* static calculations above are clearly distinguished in each certificate.

```powershell
python -u -X utf8 -m scripts.VerifySixT3SelfBilinearBatch --all --run-id six_self_with_running_v1 --continuous-run-id continuous_v1
```

The input trajectories are expected under `output/studies/general_scalar_extended/continuous_full_scalar/continuous_v1/<model-tag>/continuous_trajectory.json`. No trajectories are regenerated. The additional candidate outputs are `<model-tag>/<model-tag>_massdependent_c5.json`, and the summary includes the reported relative change in low-energy neutrino mass matrix. They remain **conditional one-scale matching predictions**, not certified sequential threshold EFT results. If any trajectory is absent, the preflight fails before creating the output directory.

Reports: `output/studies/general_scalar_extended/six_self_bilinear_batch/<run-id>/summary.json`, `summary.csv`, and `<model-tag>/certificate.json` plus `native_bilinear.json` / `native_matchete.log`.

## Physics

With `V2 = m² S†S + [mu² (S†)^T C S† + h.c.]`, and `C` the spherical singlet/triplet singlet metric, the real eigenfield masses are

\[
 m_\pm^2=m^2\pm2|\mu^2|.
\]

We obtain the native Matchete mass tensor rather than assuming its normalisation. The source-level conventions are `L_M=+mu_M Q+h.c.` and `V_RGB=+mu_RGB Q+h.c.`; for a verified *common* field basis, `mu_RGB = -k_native mu_M`. The triplet field intertwiner is checked for the quadratic metric, but **a full RGBeta/Matchete gauge-vertex and quartic convention map has not been proven**.

For the six models, the Matchete HH anomalous vertex tensor, contracted with the native real scalar metric, is checked against the original topology Clebsch–Gordan tensor (source-based group-factor check). It must have factor +1 and residual below `2e-12`.

For each of the three heavy eigenmasses and each scalar branch, the saved one-scale Matchete kernel `K(M,m1,m2)` is compared against

\[
 J(M,m_1,m_2)=\int_0^\infty\frac{t\,dt}{(t+M^2)(t+m_1^2)(t+m_2^2)}.
\]

The coefficient ratio `K/[hbar M J]` must be constant at six points. This independently checks mass dependence and detects branch errors, but **the absolute SU(2) factor is taken from the saved source, not independently rederived**, except the previously verified A(-2) identity `-2/sqrt(3)`.

After lifting the unit kernel into the symmetric three-lepton-flavour matrix with nondegenerate heavy masses,

\[
 C_5^N=\tfrac12\lambda_{T3}(K_+ +K_-),\quad
 C_5^A=-\tfrac12\lambda_{\rm extra,RGB}^* \eta_\alpha (K_+-K_-),
\]

where `eta_-2 = mu²/|mu²|`, `eta_0 = conj(mu²)/|mu²|`. The *minus* sign is the **conditional** Lagrangian/potential map; the report includes the opposite sign as a sensitivity diagnostic. `mu²` is chosen illustratively from the benchmark scalar mass, `mu²=0.01 m_self² exp(0.35 i)`, and `lambda_extra=0.02+0.01i` by default. Saved Yukawa matrices and nondegenerate heavy masses are used as source inputs. They are **not refitted**.

## Precisely what is not claimed

- No new Matchete one-loop mass-eigenstate `Match[]` is executed.
- No full RGBeta spherical triplet to Matchete Cartesian triplet **vertex** map is certified; only quadratic metric and original-vs-extra native Matchete tensors are checked.
- The UV/EFT1 RGEs are checked for presence/completeness of beta symbols. The optional trajectory mode **reuses previously evaluated** numerical running but does not integrate or certify it anew, nor match all boundary couplings to a physical component basis.
- No full sequential Wilson threshold matching, new physical neutrino mass fit, or comparison of schemes between models is claimed.
- B/C(-1) mixed-doublet mass-basis matching remains separate and blocked.

A successful run is a substantive six-model **one-scale candidate verification**, not proof that the full general scalar pipeline is physically complete.
