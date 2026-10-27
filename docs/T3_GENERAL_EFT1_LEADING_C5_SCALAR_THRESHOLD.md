# Experimental leading-C5-only scalar-threshold mode (`B_-1`, `C_-1`)

## Why it exists

The **complete** saved first-threshold EOMEFT fails Matchete's Hermiticity
check for both models. The uploaded component diagnostics isolate the failure
in the 122 heavy-scalar-only quartic terms. Background (21 terms), the
tree-matched LLSS sector (20 / 24 terms) and all Higgs-containing quartics
(52 terms) separately pass their sector checks. The full Hermiticity problem
has **not** been fixed.

The option `--c5-relevant-only` is a separate, explicit **leading one-loop
LLHH subtheory** used to unblock a raw scalar-threshold contribution. It
retains the exact saved EFT terms in the following inventory:

- 21 background/kinetic/SM terms;
- 20 (`B_-1`) or 24 (`C_-1`) matched LLSS terms;
- all 52 quartics with two external Higgs scalar fields and two X/Y fields.

The 122 quartic terms involving **four heavy-scalar fields and no Higgs**
are excluded from this particular matching call. Such interactions cannot
contribute to the leading one-loop LLSS -> LLHH graph: with one LLSS
insertion, an additional Higgs-bearing vertex is needed to provide two
external H fields, and inserting an extra four-scalar vertex raises the loop
order. This is a power-counting restriction for the selected LLSS -> LLHH
amplitude, **not** an equality to the full EFT. Contributions to other EFT
operators are not claimed to be reproduced by the reduced matching.

No source matching coefficients or UV input files are rewritten or deleted.
The full EOMEFT remains available unchanged, and normal runs continue to
fail on its genuine Hermiticity violation. With this option, the reduced
Lagrangian must independently pass every Matchete `CheckLagrangian` check.
There is **no** automatic fallback to a restricted theory.

## PowerShell

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
python -m pytest -q tests\python\test_general_t3_eft1_c5_relevant_sector.py tests\python\test_general_t3_eft1_scalar_threshold_probe.py
python -m scripts.RunGeneralT3EFT1ScalarThresholdProbe --model B_-1 --c5-relevant-only
python -m scripts.RunGeneralT3EFT1ScalarThresholdProbe --model C_-1 --c5-relevant-only
```

Results are saved under
`output/general_t3_mass_eigen_matchete/eft1_scalar_threshold_probe/`
with `_c5_relevant` in their names so the original full-theory diagnostics are
not overwritten. The summary is `summary_final_c5_relevant.json`.

## Scientific status

**Implemented:** exact partition of the saved tree expression; check that
122 excluded terms, 20/24 LLSS terms, 52 Higgs-bearing quartics and 21
background terms reproduce the full inventory; independent validation of
the reduced subtheory; 0- and 1-loop scalar matching; raw LLHH export.

**Not implemented/verified:** full-EFT Hermiticity, general scalar-EFT
matching, independent hard pole--EFT1-RGE consistency, EFT1 Wilson running,
MSbar renormalization of the scalar-threshold contribution, or a physical
sequential three-flavour C5. The previous one-scale numerical masses are
unaffected.

**Failure handling:** an unexpected operator inventory, failed independent
`CheckLagrangian`, or unsuccessful Matchete match stops execution. The
restricted mode never marks the complete EFT as validated.
