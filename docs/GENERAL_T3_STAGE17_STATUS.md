# Stage 17: additional EFT1 scalar-pair Wilson mixing

The Stage 16 B/C(-1) renormalisable UV and scalar-EFT RGE trajectories, together with the **pre-existing baseline** SM+Weinberg evolution, are already running on Windows.

Stage 17 extracts the **new quartic-dependent part** of the one-loop EFT1 psi2phi2 Wilson anomalous dimension. The existing matched F-threshold tree Wilson seed contains both PL y1*y1 and PL y2*y2 operators. The new real-scalar quartic tensor is constructed from the independent SU(2) tensors proven in Stages 9 and 13:

- `lambdaHHdagS2S2`: `V_extra = lambda2 (Hdag_1*S2_1 + Hdag_2*S2_2)^2 + h.c.`;
- `lambdaHHdagS1barS1bar`: `V_extra = lambda1 (Hdag_1*S1dag_2 - Hdag_2*S1dag_1)^2 + h.c.`;

where `lambda1`, `lambda2` are RGBeta **potential** coefficients and Stage 15 established `lambda_Matchete = -lambda_RGBeta` in the common scalar basis. `V_extra = lambda_abcd phi_a phi_b phi_c phi_d/24`; `lambda_abcd` is computed by exact real-field derivatives, with the source's identical-boson projector already incorporated.

The one-loop master tensor RGE's scalar-pair term, holding existing Wilson coefficients fixed, gives

`(16*pi^2) d C_ijab/dln(mu) |_extra = (1/2)*sum_{swap(a,b)} sum_{c,d} lambda_extra_abcd C_ijcd`.

The new script uses the existing `build_eft1_wilson_tensor` Matchete PL adapter, not an invented analytic Yukawa prefactor. It computes the nonzero real Higgs components of the **new piece only** for one lepton generation and checks that both Yukawa families occur.

**Not completed:** full 3-generation flavor Wilson ODE, inclusion of all already existing anomalous dimensions in a coupled Wilson evolution, scalar-threshold one-loop matching, pole cancellation and the final corrected C5. Do NOT add these beta components directly as the final C5; the stage is an operator-mixing audit. All production files stay unchanged. The Stage 16 legacy low-energy masses remain valid baseline values.

Requires the user's current repository `RGE/running/intermediate/ScalarOnlyTensorAdapters.py`, `output/T3_{B,C}_alpha_m1/data/eft1_after_F_wilson_seed.json`, `eft1_rgbeta_rge.json`, and Stage 15 source-convention certificate.

PowerShell from project root:

```powershell
Expand-Archive -LiteralPath .\T3_General_Scalar_Batch17.zip -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_stage17.py
powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralT3Stage17.ps1
```

Reports: `output/general_t3_stage17/{B_m1,C_m1}_extra_eft1_scalar_pair_beta.json` and `stage17_summary.json`.
