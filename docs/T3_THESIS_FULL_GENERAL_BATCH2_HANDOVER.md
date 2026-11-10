# A(-4) / D(-2) — Batch 2 F-threshold source checkpoint

## Scope and run status

Source base: GitHub `main` commit `c054b8f6012485cc1be2351ef4d3a903a94cb8fa` (read-only inspection, 10 November 2026). Batch 1 input is `output/parallel_t3/full_general_batch1/thesis_full_general_batch1_v1`.

**This package does NOT yet complete Batch 2's finite numerical Wilson matching or EFT1 Wilson evolution.** It adds a source-preserving adapter that actually invokes the existing full-general Matchete all-heavy-fermion one-loop matching, retains **every** operator expression, inventories every term without imposing a triangle-only restriction, and checksums the connection to the completed Batch 1 numerical endpoints. Matchete's symbolic couplings are **not** assigned the RG-evolved numerical values without an explicit, validated physical coupling/basis dictionary.

Do not interpret a successful raw one-loop EFT output as a renormalized finite MSbar Wilson boundary, even when the source is described as "full-general". The Match is general in its requested interaction inventory; the numerical EFT matching is not yet connected.

## Files in this ZIP (complete files)

- `scripts/RunT3ThesisGeneralBatch2.py` — source-driven runner, no production source changes.
- `Numerical/parallel_t3/wilson/general_batch2_source.py` — checksummed numerical input, mass diagnostics, operator inventory/source safety gates.
- `Lagrangian/RunT3Batch2OneLoopLedger.wl` — lossless Matchete InputForm expansion, every term with its classification, exact reconstruction flag. Does **not** drop unclassified sectors.
- `tests/test_t3_thesis_general_batch2.py` — seven offline focused tests.
- `docs/T3_THESIS_FULL_GENERAL_BATCH2_HANDOVER.md` — this handover.

All files are NEW in this delivery; no existing working module is replaced and no old result is overwritten. Existing Matchete calls and existing CG exports are reused.

## PowerShell commands

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
$zip = "C:\path\to\T3_FullGeneral_Batch2_2026_11_10.zip"
Expand-Archive -LiteralPath $zip -DestinationPath (Get-Location) -Force

python -X utf8 -m unittest discover -s tests -p test_t3_thesis_general_batch2.py -v
python -X utf8 -m scripts.RunT3ThesisGeneralBatch2 --all --preflight

python -X utf8 -m scripts.RunT3ThesisGeneralBatch2 `
  --all `
  --run-id thesis_full_general_batch2_raw_v1 `
  --accept-raw-source-only
```

`--accept-raw-source-only` returns success **only if** the Matchete raw full-general F one-loop export and complete operator ledger succeeded. It does not certify finite Wilsons. Omit the flag to receive exit code 3 on otherwise successful raw matching, emphasising the remaining physical-stage blocker. A failed Matchete invocation returns exit code 2 and records the failure. Pick a NEW `--run-id` for a rerun; old results are never deleted.

Required local dependencies: existing `wolframscript`/Matchete, NumPy/SciPy environment already used in Batch 1, the existing repository Python package, and the user's completed local Batch 1 folder.

## Expected outputs per model

Root: `output/parallel_t3/full_general_batch2/thesis_full_general_batch2_raw_v1/`.

For each `A_-4` and `D_-2`:

- `numerical_fermion_boundary.json`: **all** complex native values at MF, 3 heavy masses and left/right SVD rotations (diagnostic only), candidate scalar quadratic spectrum, D extra quartics and bilinear, and explicit `matchete_to_rgbeta_dictionary_applied=false`.
- `native_matching/native_status.json`: status from the existing full-general F-threshold matched EFT, all three heavy flavours integrated together.
- `native_matching/data/eft1_after_F_tree_inputform.txt`: full Matchete tree EFT1 expression.
- `native_matching/data/eft1_after_F_one_loop_inputform.txt`: full Matchete bare one-loop EFT1 correction, not just LLHH triangle.
- `native_matching/data/eft1_after_F_compact_inputform.txt`: tree + one-loop EFT1.
- `native_matching/data/*_wilson_seed.json`: exact tree-level lepton-number-violating operator/CG data, including D mass-eigenfield representation.
- `one_loop_operator_ledger.json`: every one-loop term, classified but **not discarded**, plus source checksum and exact reconstruction flag.
- `manifest.json`: success or explicit failure, matching source hashes, input hashes, exact remaining physics.
- `summary.json`: status for both models.

## Relevant verified Batch 1 running

Both A_-4 and D_-2 completed `CandidateGeneralRenormalisableRGIntegrated` according to committed outputs. D_-2 has five extra quartics, `y2Conj` (UV boundary zero), and a **generated** self-bilinear: approximately -4338.3146 GeV² at MF and -9617.6887 GeV² at MS, with negligible imaginary numerical residue. The existing D_-2 symbolic mass splitting uses a **real** self-bilinear; a genuine nonzero complex phase is blocked rather than erased.

## Explicit remaining components needed to finish requested Batch 2

1. **Matchete/native RGBeta dictionary** with all dimensions, heavy flavour chiral orientation, mass/eigenfield conventions and 5 D_-2 quartic normalisations. The source currently marks this `physical_matchete_rgbeta_mapping_certified=false`.
2. **Complete one-loop coefficient extraction** from the entire matched EFT1 operator set, in a declared, nonredundant Wilson basis. A source-level all-term ledger is not this numerical map.
3. **Subtraction prescription**, including UV-hard minus EFT1-soft contributions, local counterterms, the exact dimensional-regulator residue, and finite MSbar threshold coefficients. No pole is simply set to zero and no UV triangle replaces the full match.
4. **Complete EFT1 one-loop operator anomalous dimensions** and evolution to 1e3 GeV, with both heavy scalars active, all D quartics, generated bilinear and relevant higher-dimensional operator mixing. Existing `Numerical/parallel_t3/wilson/adm.py` is a *candidate psi² phi² subset*, not a certified complete EFT1 matrix.
5. **Full mass eigenbasis map for the joint scalar threshold** in Batch 3.

### Numerically executable interface for the next Batch 2 step

```python
from Numerical.parallel_t3.background.general_batch1_api import load_batch1
from Numerical.parallel_t3.wilson.general_batch2_source import native_numbers

b1 = load_batch1("output/parallel_t3/full_general_batch1/thesis_full_general_batch1_v1/D_-2")
values = b1.fermion_matching_input()["native_uv_rgbeta_couplings_at_muF"]
f_boundary = native_numbers(values, "D_-2")
# f_boundary contains ALL UV couplings and the running-generated bilinear;
# this is a parameter boundary, NOT yet a finite Wilson coefficient.
```

No numerical C5, PMNS/chi², or optimised Yukawa is claimed. Do not connect the current raw loop ledger directly to the restricted Sobol fitter.
