# Worker 1 — Scalar Quartic Tensor Robustness

## Status and source

Reviewed public main repository tree `bab0061547c0d850eef04783c85dfb8b2cab32b9` (tree SHA, not commit SHA), and existing Worker B discrepancy, direction, reconciliation, basis-search and source reports. Their audits use the same running muS=1000 GeV endpoints with the D(-2) special basis mapped back to native components. Results are **conditional** on native/Matchete/RGBeta coordinate identifications.

**Numerical findings:** A(-4) candidate/source coefficient -0.9216565588633524, post-scale residual 0.38800668487045953; D(-2) -0.7764157057450806, post-scale residual 0.6302211134771419. A cross S1S2 scale -0.09763107293781742 with 0.9952226753832595 residual; D S1S2 scale 0.41178128083801047 with 0.9112827095646059 residual and S2S2 -0.011865399420561747 with 0.9999296036704738 residual. Thus the mismatch is **not merely a sign or factor**.

Direction audits: ordinary Higgs portals show candidate/source -1, residual near zero; lambdaT3 real/imag shows coefficient -0.097631 and residual ~0.995223. D extra general HHSS directions Q020101I1 and Q110002I1 show +1 and zero residual. 384 discrete basis transforms per model still leave lambdaT3 best real-scale residual ~0.82231; the best unconstrained 2x2 lambdaT3 coordinate fit residual is ~0.78781. These transformations do not exhaust continuous field rotations.

## Which tensor is in the selected solver?

- `scripts/RunT3SimultaneousEFT1A4.py`: the `HHSSLinearBasis` calls `Numerical.running.GeneralT3SequentialLeadingLog.scalar_tensors` for running HHSS, then uses the explicit seven-term non-HHSS potential. Coupled DOP853 runs renormalisable and Wilson state together.
- `scripts/RunT3SimultaneousEFT1All16.py`: analogous running HHSS with `general_quartic`, whose non-HHSS is provisional.
- `scripts/RunT3SimultaneousEFT1MixedBC.py`: original fixed Cartesian basis for the mixed B/C cases; not part of this A/D investigation.
- `Numerical/running/T3CoupledQuarticContext.py` is a candidate cached Batch3B context, **not** proof that its quartic is the production source. Worker B reconciliation exports a candidate hybrid without connecting it to production.

The physically justified independent replacement has **not** been established. The source and candidate differ in coupling-coordinate geometry; neither is shown to be an independently correct complete physical 12-real-scalar potential. Consequently it would be scientifically invalid to run a fitted sign-flip candidate and present it as an alternative prediction.

## Scope and interpretation

At one scalar endpoint, tensor geometry is demonstrably non-proportional, but flavour may be insensitive if the matrix contributions are aligned with the same LLSS Wilson-flavour tensor. That alignment **has not been verified for the fully coupled, running problem**, and cannot be inferred from a tensor norm. Observed stage-freeze predictions from the original solver are baseline predictions, not discrepancy sensitivities. No current evidence establishes invalidation or robustness of the quantitative neutrino flavour conclusions. Their **qualitative status remains conditional**, rather than disproved.

No full 16-model reconstruction warranted before physically identifying the A/D lambdaT3 map and source quartic normalization.

## Files added (project-relative)

- `scripts/RunT3Worker1ScalarRobustness.py` — existing-audit loader and independent NPZ numeric cross-check, optional gated *saved* paired-IVP comparison
- `tests/test_t3_worker1_scalar_robustness.py` — synthetic numerical tests
- `docs/WORKER_1_SCALAR_ROBUSTNESS.md` — evidence and execution instructions

The tool only writes to a fresh `output/audits/worker1_scalar_robustness_v1/` directory, never edits production, GitHub or historical files. It does **not** run the coupled IVP itself or claim that a physical alternative exists. The optional comparison only reads two already generated `cases.json` files and rejects attempts without `--physical-alternative-certified`; the operator must independently verify identical UV input and scalar conventions.

## PowerShell

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
$zip = "$HOME\Downloads\T3_Worker_1_Scalar_Robustness.zip"
Expand-Archive -LiteralPath $zip -DestinationPath . -Force
python -m pip install numpy scipy
python -X utf8 -m unittest discover -s tests -p 'test_t3_worker1_scalar_robustness.py' -v
python -X utf8 -m scripts.RunT3Worker1ScalarRobustness --project . --out-dir output/audits/worker1_scalar_robustness_v1
Get-Content output\audits\worker1_scalar_robustness_v1\REPORT.md
```

To reproduce independently sourced cached results first, where available in the local project (choose fresh output paths if historical paths already exist):

```powershell
python -X utf8 -m unittest discover -s tests -p 'test_t3_worker_b_batch_reconcile.py' -v
python -X utf8 -m unittest discover -s tests -p 'test_t3_worker_b_basis_search.py' -v
```

An optional paired *certified* solver-results comparison is gated; both inputs must be independently generated under identical UV benchmark, matched convention and evolution flags:

```powershell
python -X utf8 -m scripts.RunT3Worker1ScalarRobustness `
 --baseline-cases 'PATH_TO_BASELINE_CASES.json' `
 --alternative-cases 'PATH_TO_CERTIFIED_ALTERNATIVE_CASES.json' `
 --physical-alternative-certified `
 --out-dir output/audits/worker1_certified_pair_v1
```

**Required local artifacts** for full remeasurement: all five Worker B JSON audit reports under `output/audits/workerb_{hhss_discrepancy,hhss_direction,basis_search,batch_reconcile,source_hhss}_v1/` (exact paths in `AUDITS` in the script), plus `output/audits/workerb_batch_reconcile_v1/{A_-4,D_-2}_reconciliation.npz` for independent tensor cross-check. Existing source arrays are `output/audits/workerb_source_hhss_v1/{A_-4,D_-2}_source_hhss_native.npz`; cached Batch3B arrays `output/audits/batch3b_quartic_fixed_20261111_234059/{A_-4,D_-2}/quartic_response.npz`. To establish final neutrino impact additionally require two **physically certified** and UV-identical model `cases.json` sets, not merely candidate-hybrid simulations.
