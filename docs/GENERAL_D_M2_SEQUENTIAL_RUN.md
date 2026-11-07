# General D(-2) local sequential matching — integrated runnable entry point

**Run from the repository root in PowerShell** after extracting the ZIP.

```powershell
python -X utf8 -m unittest discover -s tests/python -p test_general_d_m2_sequential_matching.py
python -X utf8 -m scripts.RunGeneralDMinus2SequentialMatching
```

Default performs the full chain in separate processes:

1. Rebuild unrestricted D(-2) UV and match heavy fermions using the existing first-threshold driver. The original UV registry and first-stage code are unchanged.
2. Recompute one-generation EFT1 beta from the exact real-scalar CG seeds.
3. Regenerate full-flavour direct EFT1 running, which is kept OUT of scalar one-loop `Match[]`.
4. Rebuild exact physical scalar registry and match all heavy scalars using `Lagrangian/RunGeneralDMinus2SequentialScalarThreshold.wl`. Reversibly promote ONLY `MS1`, `T3S2RMass`, `T3S2IMass` from EFT-order 1 to 0 on a temporary full-EFT1 copy; compare the quadratic sector to `FreeLag`, validate `CheckLagrangian`, and run full one-loop `Match[]` plus all four reductions.
5. Check hard plus inherited EFT1 first-threshold O(hbar) terms, and full LLHH operator/chirality inventory.
6. Verify direct-running flavour symmetry modulo bound summation dummy renaming.
7. Assemble the provisional symmetric hard + inherited + running C5 with exactly one `hbar` per term and no double symmetrization of direct running.
8. Validate the extractor against Matchete's Warsaw reference operator.
9. Evaluate the finite-part `NF=1` numerical benchmark with physical real-singlet mass splitting, Takagi singular values and explicit pole reporting.
10. Check exact symbolic cancellation of the simple pole for `NF=1`.

If the first threshold/RGE/full-flavour running are already verified and you only want to rerun the **complete scalar threshold and coefficient assembly**:

```powershell
python -X utf8 -m scripts.RunGeneralDMinus2SequentialMatching --reuse-boundary
```

Increase the per-operation scalar limit (seconds) if your machine is slower:

```powershell
python -X utf8 -m scripts.RunGeneralDMinus2SequentialMatching --scalar-stage-timeout-seconds 180
```

Results: `output/general_sequential_pilot/D_m2_stage1/data/general_d_m2_scalar_threshold/sequential_run/sequential_matching_manifest.json`. Each step also has a `sequential_run/<stage>.log`. Under `profile_tree/`: `profile_full_scalar_one_loop_processed.wxf` (hard scalar-matched EFT), `profile_full_scalar_one_loop_hard_hbar_sector.wxf`, `profile_c5_flavour_candidate_total_provisional.wxf`, `profile_c5_rough_benchmark_result.json`, and `profile_c5_pole_cancellation_result.json`.

The orchestrator compares SHA256 hashes of all four physical first-threshold WXF files before/after scalar processing and fails on modifications; it does not modify the legacy `Lagrangian/RunGeneralDMinus2ScalarThreshold.wl` (which may differ locally from the public snapshot). Only the separately named consolidated scalar driver is used. This avoids overwriting local work.

## Status boundaries

A successful local run reports `SequentialC5BenchmarkReproduced_ProvisionalNormalization`. This demonstrates a reproducible UV -> F EFT1 -> all-scalar one-loop matching and provisional C5/Takagi spectrum for the selected benchmark, **not** an independently certified complete SMEFT operator basis or independently certified model-specific Clebsch–Gordan normalization. The Warsaw reference test does not prove the latter; exact NF=1 pole cancellation is not a full renormalization proof.

Prerequisites include your existing project files (`Lagrangian/RunGeneralDMinus2Stage1.wl`, `RunMatching.wl`, builder and model catalog; Python stage1/RGE/running runners, Matchete and Wolfram). No GitHub actions occur. The ZIP contains full files under their intended project paths, not patches. The final numerical configuration is `benchmarks/general_d_m2_rough_benchmark_finitepart.json` and can be passed via `--config`.
