# Stage 13: RGBeta projector mapping + coupled B/C(-1) UV running

This is an **experimental UV-only** extension. The original T3 pipeline is not modified.

## New physics calculation

`Numerical/experimental/GeneralT3CoupledUV.py` calls the *existing production* `Numerical.running.UVRunner.run_uv_segment` with the complete `output/T3_{B,C}_alpha_m1/data/uv_rgbeta_rge.json` payload, evolved between the toy UV scale and a scale **above every initial heavy fermion mass**. The full state includes all six previously supported complex RG-closure quartics. For the two C5-active quartics it reports both (i) the coupled one-loop RGBeta solution and (ii) the first-order frozen-background leading-log estimate over the same interval.

The new code is **not** a separate partial beta-function implementation. It invokes the existing project beta evaluator and validates the complete beta/state key layout. It does not run below a heavy-particle threshold or compute the fully RG-improved neutrino mass.

The example `docs/general_c5_toy_benchmark.json` is not a fitted benchmark. For the UV-only run the default endpoint is **3500 GeV** (initial masses 1500, 2100, 2800 GeV). This differs deliberately from the earlier frozen-source pilot's 1500 GeV endpoint, which lies below the heaviest heavy fermion and was not a proper full-UV-running interval.

## RGBeta symmetric-projector audit

Source-verified against [RGBeta `Package/GroupsAndIndices.m`](https://github.com/aethomsen/RGBeta/blob/master/Package/GroupsAndIndices.m) and the project's `RGE/running/rgbeta/T3RGBetaModel.wl`:

```
P[a,b,c,d] = (delta[a,c] delta[b,d] + delta[a,d] delta[b,c])/2
```

For the all-conjugated interaction, transforming the two last indices by the standard SU(2) pseudo-reality intertwiner `epsilon[1,2]=+1` gives

```
P[a,b,u,v] epsilon[u,c] epsilon[v,d]
 = (epsilon[a,c] epsilon[b,d] + epsilon[a,d] epsilon[b,c])/2.
```

The two tensors match exactly the canonical choices measured in the local Stage 9 Matchete exporter, with factor `+1` **at the SU(2) group-tensor level**. That does **not** prove an equality of the complex coupling coefficients in the physical Lagrangian: the relative field-conjugation/basis orientation, potential-vs-interaction sign, and hermitian-conjugate phase must still be checked at the vertex level. Accordingly `full_lagrangian_coupling_identification_verified=false`.

## Install and run

From project root:

```powershell
Expand-Archive -LiteralPath .\T3_General_Scalar_Batch13.zip -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_stage13.py
powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralT3Stage13.ps1
```

Expected paths:

```
output/general_t3_stage13/doublet_tensor_map.json
output/general_t3_stage13/B_m1_coupled_uv.json
output/general_t3_stage13/C_m1_coupled_uv.json
output/general_t3_stage13/stage13_summary.json
```

If existing `output/T3_{B,C}_alpha_m1/data/uv_rgbeta_rge.json` exports are not available, the runner attempts the corresponding `output/general_t3_trial/{B_m1,C_m1}/baseline/data/uv_rgbeta_rge.json` path and explicitly reports absence if neither exists. Refresh missing exports using the project's normal RGBeta production runner, **not** by fabricating beta functions.

Important: numerical beta-running and matched extra C5 responses are still reported separately. Do not interpret their sum or product as a fully renormalized general T3 neutrino-mass prediction until the complete vertex-level Matchete–RGBeta coupling convention is verified.
