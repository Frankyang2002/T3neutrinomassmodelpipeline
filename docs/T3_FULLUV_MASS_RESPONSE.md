# T3 full-UV scalar mass boundary response

This **diagnoses** the 10 currently supported full-UV-running cases; it does not alter beta functions, matching, or acceptance rules, and does not fit Yukawas.

For historical full-UV Yukawas and quartics, vary only the initial two scalar **mass-squared** values by factors `0.25,1,4`. Evolve the exact repository UV and EFT1 RGEs. Record the UV beta derivatives, masses squared at the fermion threshold and scalar threshold, physical scalar masses, rhoF, and screening outcome. A comparison across factors identifies sensitivity to initial scalar masses versus running-induced contributions; it is not a mathematical decomposition of the nonlinear RG flow.

Run from repository root:

```powershell
python -X utf8 -m pytest -q tests\test_t3_fulluv_mass_response.py
python -u -X utf8 -m scripts.DiagnoseT3FullUVMassResponse --run-id fulluv_mass_response_v1 --source-campaign overnight_six_v1_bulk48 --factors 0.25,1,4 *> .\fulluv_mass_response_v1.log
Get-Content .\fulluv_mass_response_v1.log -Tail 30
```

Results: `output/studies/general_scalar_extended/bulk_48/fulluv_mass_response_v1/mass_response.jsonl`. Resume interrupted runs by adding `--resume`. Single model: `--models A_+0`. Do not treat negative or failed running masses as fit accepted. Scalar TeV mass window remains [0.5,5], rhoF <=0.5.

For a compact independent check of change at fermion threshold, inspect `uv_at_muF_mSsq_GeV2`, and for later running inspect `eft_at_muS_mSsq_GeV2`. The baseline UV beta derivative `beta_mSsq_at_UV_GeV2_per_lnmu` is signed with respect to increasing ln(mu), whereas integration proceeds downward in mu.
