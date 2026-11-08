# Experimental B/C X-only matching from saved Y EFT2 in a fresh Matchete kernel (v3)

## Verified preceding numerical results (provided by local Windows runs)

| Model | Running muY (GeV) | Running muX (GeV) | Conditional masses descending (eV) |
|---|---:|---:|---|
| B_-1 | 1187.80578221 | 1035.77775471 | 0.0385336229, 0.0156138215, 0.0004584847 |
| C_-1 | 1185.39334545 | 1032.85523509 | 0.0464057419, 0.0459238479, 0.0216846442 |

These are conditional sequential numerical calculations **without finite one-loop scalar matching**. The bubble diagnostic is a candidate sensitivity study, not certified physical matching.

## What the previous failure does and does not demonstrate

The `--diagnose-x` output in v2 verified that all six Matchete validity checks passed for both the *direct heavy-registry* and *heavy-free-rebased* EFT2 expressions, in both models. The light-versus-heavy free Lagrangians differ, so the *rebased* expression is the appropriate selected matching source. The previous `--match-x` stopped at B_-1, with X report status `XFollowupNotStarted` even though the validations passed. Inspection of the v2 source showed that `XFollowupNotStarted` was an initial placeholder and that the X tree/loop matching calls lacked durable intermediate X-report writes. It does **not** show whether the tree call aborted, an internal Matchete stage failed, or the kernel terminated unexpectedly. The precise cause requires the actual matching log.

## What v3 changes

1. Adds `Lagrangian/RunT3BCXFromSavedY.experimental.wl`, a **complete** experimental Wolfram script derived from the v2 X native file. Existing v2 code remains unchanged.
2. Establishes the same independently checked *complete* native F-first-source reconstruction with all scalar quartics, Yukawa-derived structures and exact saved LLSS. This is a consistency and object-definition step only, **not a proof of EOMEFT equivalence**.
3. Instead of running Matchete's native Y-only tree and loop Match in this kernel, loads the **saved tree and one-loop Y-only EOMEFT** from a previous successful v2 `--diagnose-x` run. The Python preflight verifies the original v2 Wolfram source hash, prepared input hash, F-first EOMEFT hash, and native Y/X check reports; unexpected source changes fail closed.
4. Invokes X-only matching on the **saved Y-tree EOMEFT** after reclassifying X as heavy and replacing its free Lagrangian. The Y-loop artifact is retained for provenance but intentionally **not fed into the X matching source**. This matches the earlier structural experiment, but without retained earlier Match state in the current Wolfram kernel.
5. Before each potentially expensive X-only `RunT3Matching` call, exports the X report status `XFollowupTreeMatchingInProgress` or `XFollowupOneLoopMatchingInProgress` and checkpoints `.progress.json`. On a returned Matchete failure, records the actual stage status.
6. Runs B and C independently even if one fails. Each model produces its own `.log`, `.progress.json`, `.provenance.json`, and raw Matchete X report, if available.
7. `CollectT3BCXFreshEvidence` packages all available reports/logs/raw operators with SHA256 hashes, **without requiring a successful X match**. Missing evidence and physical-certification flags are explicitly recorded.

### Important uncertainty and abort interpretation

- If `XFollowupYTreeEOMCheckFailed`, the saved EOMEFT does not validate under the independently reconstructed fresh Matchete registry. New effective couplings/CGs introduced during prior matching may be unavailable; this is an input-object/registry problem and must not be bypassed.
- If `XFollowupTreeMatchingInProgress` persists, the next action is to inspect the actual log and `.progress.json` for the precise Matchete failure, crash, resource exhaustion or timeout.
- `XFollowupOneLoopMatchingInProgress` means tree matching returned and one-loop matching was entered, but did not complete.
- `XFollowupTreeMatchingFailed` / `XFollowupOneLoopMatchingFailed` mean Matchete returned a failure; look at the recorded stage status and log.
- `XFollowupYTreeOnlyXRawMatchingNotCertified` indicates **raw native structural X matching** completed. It does *not* establish full renormalized matching.
- The original v2 log for B_-1 may already contain the decisive diagnostic. Inspect it before running an expensive fresh-kernel continuation.

## Install and commands (PowerShell)

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$env:USERPROFILE\Downloads\t3_bc_x_fresh_kernel_v3.zip" -DestinationPath . -Force

python -m unittest discover -s tests -p "test_t3_bc_x_fresh_kernel_v3.py" -v

# Inspect the existing v2 failure first:
Get-Content "output\general_t3_mass_eigen_matchete\eft1_scalar_xy_native_followup_v2\B_m1_final_xy_native_match_x.log" -Tail 80

# Verify source hashes and previously passed native checks:
python -u -X utf8 -m scripts.RunT3BCXFreshKernel --all --preflight-only

# Fresh-kernel structural check:
python -u -X utf8 -m scripts.RunT3BCXFreshKernel --all --check-x

# If both structural checks pass, perform X tree and loop Match in fresh kernels:
python -u -X utf8 -m scripts.RunT3BCXFreshKernel --all --match-x

# Always usable even when either X match fails:
python -u -X utf8 -m scripts.CollectT3BCXFreshEvidence --all
```

Output directory: `output/general_t3_mass_eigen_matchete/eft1_scalar_xy_fresh_kernel_v3/`

Evidence ZIP: `output/audits/t3_bc_xy_native_followup/xy_fresh_native_evidence_v3.zip`

If the fresh run fails, upload the evidence ZIP. It contains the existing v2 B matching log and both new model-specific logs and source hashes, where available.

### Source-faithfulness and physics limitations

The new script is based on the **verified exact v2 file as installed in the previously supplied ZIP**; the launcher requires its SHA256 to match the source hash recorded during the user's v2 diagnostic. The stored Y EOMEFT is an **independent canonical reconstruction** and has not been demonstrated equivalent to the original saved F-first EOMEFT. No higher-order Y-loop matching is inherited by the X source; no running `continuous_v1` boundary couplings are substituted. The conditional one-generation SM Yukawa approach and additional quartic mapping status are unchanged. No statement of physical C5/neutrino mass certification is made. The v3 scripts have been syntax-compiled and offline-tested, not run with a Windows Matchete installation.
