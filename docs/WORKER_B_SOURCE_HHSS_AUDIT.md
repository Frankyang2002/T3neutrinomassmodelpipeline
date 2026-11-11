# Worker B: Direct source HHSS at the EFT1 scalar endpoint

**Scope:** A(-4) and D(-2), `muS=1e3 GeV`, all eight BSM real scalar components retained. Read-only. No UV/fermion evolution and no scalar threshold matching is performed.

The script uses the *current repository source functions*, not fabricated `hhss` arrays:

1. `load_release_manifest` and `load_source` get the source Matchete boundary labels.
2. `RunT3GeneralSequentialFitAblation._load_channels` loads v4 native quartic channels.
3. `GeneralT3SequentialLeadingLog.scalar_tensors` evaluates native channels, completed RGBeta portals and registry-derived extra HHSS invariants at the *given* EFT1 `running.json` scalar threshold.
4. D(-2) requires `build_special_basis`: it returns the quartic in a signed matching boundary basis. The script applies `R Q_boundary R^T` on heavy legs to restore the fixed native 8-real-component candidate basis before comparison.
5. `CandidateQuarticContext.load` verifies the Batch3B response, registry hash, and quartic seed hash.
6. `audit_tensors` reports source versus candidate HHSS relative discrepancy and the hybrid source-HHSS/non-HHSS checks.

**Important:** The D(-2) field-basis alignment is source-convention-derived but not an independently certified native-to-RGBeta identification. Disagreement between source and candidate HHSS is a *measurement*, not an automatic failure; the normalization is already known to be uncertified. This does not establish the full 16-model quartic potential.

Install ZIP from the repository root. Run:

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive "$HOME\Downloads\T3_Worker_B_Source_HHSS_Comparison.zip" -DestinationPath . -Force
python -X utf8 -m unittest discover -s tests -p "test_t3_worker_b_source_hhss_audit.py" -v
python -X utf8 -m scripts.RunT3WorkerBSourceHHSSAudit `
  --models A_-4,D_-2 `
  --response-root "output\audits\batch3b_quartic_fixed_20261111_234059" `
  --state-file "A_-4=output\general_observable\all16_v1\A_m4\running.json" `
  --state-file "D_-2=output\general_observable\all16_v1\D_m2\running.json" `
  --out-dir "output\audits\workerb_source_hhss_v1"
```

This writes a **new** `report.json`, `A_-4_source_hhss_native.npz`, `D_-2_source_hhss_native.npz`. If `--out-dir` already exists, the script refuses to overwrite it. `source_hhss_vs_candidate_relative` is the numerical comparison. The NPZ files can also be supplied to the existing quartic audit via `--source-hhss MODEL=NPZ`.

If a required release, native seed, portal, registry direction, or basis convention is unavailable, the script fails explicitly instead of inserting zero couplings. It does not require re-running the physics solver.
