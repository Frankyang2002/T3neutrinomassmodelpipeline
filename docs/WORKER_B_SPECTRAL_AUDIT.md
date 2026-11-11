# Worker B — orthogonal-basis spectral obstruction audit

Read-only and physics-provisional. Uses already created, hash-bound Batch3B response and Worker B reconciliation. Does not modify production files or historical results.

## Installation

Extract ZIP at `C:\thesistemp\T3neutrinomassmodelpipeline`.

## Run (PowerShell, from repository root)

```powershell
python -X utf8 -m unittest discover -s tests -p "test_t3_worker_b_spectral_audit.py" -v
python -X utf8 -m scripts.RunT3WorkerBSpectralAudit `
  --models A_-4,D_-2 `
  --response-root "output\audits\batch3b_quartic_fixed_20261111_234059" `
  --reconcile-dir "output\audits\workerb_batch_reconcile_v1" `
  --out "output\audits\workerb_spectral_v1\report.json"
```

## Reading the results

- `directions.lambdaT3.real.normalized_spectral_lower_bound`: a nonzero value proves **no real orthogonal Higgs/heavy scalar rotation and overall scale** can identify those **fixed** coupling-direction tensors. Even this lower bound allows an unphysically broad independent O(16) x O(64) matrix rotation, so it is a necessary, not sufficient, condition.
- `lambdaT3_source_phase_scan.real.best_sampled_phase.spectral_distance`: smaller means some real/imag source-direction mixture has a similar SVD fingerprint. This is a finite angular search; do not treat its minimum as proof of a global no-go for arbitrary coupling redefinitions.
- `directions.lambdaH2.real.signflip_rel_to_candidate`: checks the known portal sign result.
- `candidate_species_sectors` and `source_species_sectors`: pure S1S1 and S2S2 block spectra (and mixed block) in the fixed model-specific species ordering. These are invariant only under transformations preserving S1/S2 splitting, not arbitrary scalar mixing.

The operator is the HHSS coefficient Q_ABab reshaped as a 16 x 64 matrix. For real orthogonal H, S changes, M transforms as `(H kron H)^T M (S kron S)`; therefore its singular values cannot change. Normalized singular-value distance is an exact necessary obstruction for fixed direction, though vanishing distance does not prove existence of an admissible field map.

This does **not** establish Matchete-to-RGBeta invariant normalization or physical matching. Source-HHSS hybrid remains the provisional running context. The next action is to trace the individual `lambdaT3` source potential and candidate Batch3B potential if spectral obstruction is nonzero.
