# A(-4): CheckLagrangian False and non-UTF-8 Matchete output

Source inspected: repository `Frankyang2002/T3neutrinomassmodelpipeline`, commit `80dcee89d9a7d9565e76013b4e5c78735b6193b6`. The subsequent `thesis_physical_scalar_massrepair_v2` result was supplied by the user; it was not committed to the repository at inspection time.

## Verified results and limitations

- Existing `tests/test_t3_thesis_final_one_batch.py` passed 13 tests, but these are not a native symbolic or numerical physics verification.
- A(-4) passed scalar mass promotion but `CheckLagrangian` returned **False**. The older committed fresh-threshold run also printed `CheckLagrangian::Hermiticity`; whether the new failure has the identical residual is not established without the new native log.
- D(-2) successfully matched three-flavour tree Wilsons across both scalar eigenfields (`baseline=4`, `extra=4`), returning 68 scalar-soft Weinberg terms. This is not an RG-improved renormalised `C5`.
- The old Python runner uses `subprocess.run(text=True)` and fails to decode arbitrary Wolfram stdout under UTF-8; this creates an unrelated `UnicodeDecodeError` and then `TypeError` when `.stdout` becomes `None`.

## New files / local paths

1. `Lagrangian/RunT3AMinus4FullGeneralPhysicalScalar.wl` — **full replacement**. The existing mass promotions remain; if validation fails, it attempts only exact dummy-index relabelling or a proof-checked alpha-equivalent SM-sector rewrite. Neither projected nor symmetrized Lagrangians are used in physical matching. If validation still fails, it records a two-stage diagnostic suppression (`W12` absent; then `W12` and `lambdaT3` absent), only to isolate the problematic operator sectors. That diagnostic does **not** imply those interactions may be dropped.
2. `scripts/RunT3ThesisFinalOneBatchSafe.py` — **new** encoding-safe entrypoint. Runs the existing `RunT3ThesisFinalOneBatch.py` unchanged, with subprocess output captured as binary and decoded with UTF-8 backslash replacement. No thread-based text decoding. Preserves readable Wolfram logs including undecodable byte escapes.
3. `tests/test_t3_a4_canonical_and_safe_capture.py` — **new** static regression checks.
4. `docs/T3_A4_CHECKLAGRANGIAN_AND_ENCODING_v4.md` — **new** instructions.

## Commands (PowerShell)

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$HOME\Downloads\T3_A4_CheckLagrangian_Encoding_v4.zip" -DestinationPath . -Force
python -X utf8 -m unittest discover -s tests -p test_t3_a4_canonical_and_safe_capture.py -v
python -X utf8 -m scripts.RunT3ThesisFinalOneBatchSafe --model A_-4 --run-id thesis_a4_hermiticity_diagnosis_v4 --debug
$root = 'output\parallel_t3\full_general_final_one_batch\thesis_a4_hermiticity_diagnosis_v4\A_-4'
Get-Content "$root\A_joint_scalar_lagrangian_audit.json"
Get-Content "$root\native_two_threshold_status.json"
Get-Content "$root\native_two_threshold.log" -Tail 80
```

The Python process may return exit code 2 and `Blocked_CompletionBlocked` despite successful scalar-soft matching, because the existing driver still demands a complete physical `C5`. Check the **native** `status` in `native_two_threshold_status.json`.

## Interpretation

- `CheckLagrangian_passed=true` and `accepted_canonicalisation=ProvenDummyIndexRelabel` or `ProvenEquivalentSMRepresentation`: physically equivalent tree input passed native validation; a joint scalar match will be attempted.
- `CheckLagrangian_passed=false`, `CheckLagrangian_without_W_DIAGNOSTIC_ONLY=True`: the problematic hermiticity/normalization lies in the `W12` sector or its interference with other couplings. This does NOT authorise zeroing `W12` in matching.
- Both zero-W and zero-W-plus-lambdaT3 fail: the issue is outside those two operators or in interactions between remaining terms; examine native log rather than suppress the checks.
- A new successful scalar-soft match still lacks the one-loop F-hard Wilson boundary, EFT1 Wilson tensor anomalous dimensions, finite hard+soft subtraction and low-energy physical `C5`.

No GitHub edits are performed. The D(-2) physical scalar continuation is untouched.
