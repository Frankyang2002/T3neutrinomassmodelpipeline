# Worker C — D(-2) restricted saved-Matchete finite candidate (diagnostic)

## What this does

Processes the **existing** D(-2) full-EFT1 scalar threshold output, **without calling `Match[]`**. At the saved, assembled Weinberg-expression level it sets the independent `y2Conj` coupling to zero and preserves `y1`, `y2`, all scalar couplings and any Matchete conjugation of the retained couplings. For the **NF=1** test case, it extracts the simple dimensional-regulator pole and checks that the *assembled* restricted candidate has zero pole after exact contraction. Only if the pole cancels does it export the corresponding finite symbolic `C5` candidate and its hard/inherited/running components.

This **is not** a new matching calculation for the restricted UV action; evaluation at zero `y2Conj` is a conditional restriction of a previously matched general source. It does not certify the finite renormalization scheme. It does not compute or subtract the overlap with the existing finite `LLSS*HHSS` scalar bubble. The exported candidate **must not be added** to the Worker C joint result.

The full saved expression has independent `y2` and `y2Conj` sectors. The D(-2) saved reference contains a successful `Match[..., EFTOrder->5,LoopOrder->1]`, 68 one-loop Weinberg candidate terms, and an NF=1 pole-cancellation audit for the *unrestricted* provisional sum. The restricted pole check must be run independently; that unrestricted proof does not suffice.

## Files and destinations

- `Lagrangian/RunWorkerCRestrictedDMinus2Finite.wl`
- `scripts/RunT3WorkerCRestrictedDMinus2.py`
- `tests/test_t3_worker_c_restricted_d_minus2.py`
- `WORKER_C_RESTRICTED_D_MINUS2_README.md`

Existing local scripts/solvers and historical outputs are not modified. Public GitHub was inspected read-only; your Windows local source has not been directly inspected.

## PowerShell commands

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$HOME\Downloads\T3_worker_C_restricted_D_m2_20261112.zip" -DestinationPath . -Force
python -m unittest discover -s tests -p "test_t3_worker_c_restricted_d_minus2.py" -v
python -m scripts.RunT3WorkerCRestrictedDMinus2 --preflight
python -m scripts.RunT3WorkerCRestrictedDMinus2
```

Output is placed in a new `output/audits/worker_c_restricted_D_m2_*` directory (never overwrite). A successful run writes `result.json`, `restricted_total_finite_candidate.wxf`, corresponding text, and component WXF files. A nonzero pole, missing input, or any unsupported regulator structure prevents a finite coefficient being exported.

### Immediate next comparison

After local symbolic extraction succeeds, normalize its operator and flavour conventions against the already-run EFT1 `C5`, evaluate both at the **same** scalar threshold snapshot, and compare against the existing joint bubble. Determine the non-overlap algebraically before introducing any numeric `delta C5` to the production pipeline. The current script does not claim to have done that.

## Local Wolfram import fix (2026-11-12)

The original diagnostic incorrectly passed an Association of file paths to
`AssociationMap`, which expects keys rather than a key-to-path association.
The resulting `AssociationMap::invrlf` errors caused all WXF imports to fail,
followed by a misleading `y2Conj` source mismatch and exit code 8.
The replacement uses `Map` to transform association values for the four WXF
imports, the restricted expressions, pole separation, and source hashes.
The runner now rejects status 0 without an actual `result.json` and finite
WXF candidate, and preserves nonzero Wolfram exit codes.

Eight Python tests passed in the build environment. The Wolfram calculation
requires a Windows-side rerun; no restricted-pole cancellation is yet claimed.

PowerShell from project root after extracting this ZIP:

```powershell
python -m unittest discover -s tests -p "test_t3_worker_c_restricted_d_minus2.py" -v
python -m scripts.RunT3WorkerCRestrictedDMinus2
```

The runner prints the actual report path on success. No physical contribution
is connected to production `C5`.
