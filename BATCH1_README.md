# T3 Batch 1: executable completion assessment

This package is deliberately **read-only**. It does not replace numerical model modules or modify results: the local repository source is not available in this environment, so edits to those files cannot be validated here. It does include an executable completion gate, using the actual modules revealed in the user's local audit.

## Install
Extract into `C:\thesistemp\T3neutrinomassmodelpipeline` preserving the `scripts/maintenance/` directory.

## Run (PowerShell, from the repository root)
```powershell
python scripts/maintenance/Batch1CompletionGate.py --run-tests --run-self-tests
```

Outputs:
- `output/audit/batch1/batch1_execution_report.json`
- `output/audit/batch1/batch1_execution_summary.md`

Only the new audit output files are written. Test suites may create their own temporary files according to repository test behaviour. The gate runs individual test files, avoiding an unmanageable all-test transcript.

## Source-level conclusions already supported

1. The actual `T3GeneralYukawaCoupledUV.py` enumerates 16 ordinary representations, checks each RGBeta beta key, and builds model-specific fully coupled UV equations. It reads the generated one-loop beta JSON and benchmark configuration.
2. `T3RGAssistedExtendedC5.py` implements a six-model UV-to-extended-C5 numerical bridge. Its scenario accepts explicit, *unverified* RGBeta -> Matchete complex conversion factors; this is an assumption, not a physical convention proof.
3. `T3ExtendedLowEnergyStudy.py` computes charged-lepton-basis neutrino observables from C5 and evolves SM+Weinberg below a single matching scale. This is not a full separated-threshold EFT sequence.
4. `T3BIntegratedCoupledPipeline.py` invokes the existing general pipeline and optionally the B-specific exploratory coupled-leading-log calculation. It is not a universal 16-model full-scalar endpoint.
5. Existing test audit: 3,065 passed, 33 failed, 11 skipped. Physics-priority failures include EFT1 hard pole and SU(2)/quartic convention integration checks. Many other failures concern documentation, reporting, encoding or path assumptions. The supplied audit gives the exact failed test names.

## Batch 2 priority

Do not add another framework. Verify the exact UV -> C5 mapping for one representative model, run the corresponding full numerical benchmark from a single state, compare against a known limit and solver tolerance, then expand the proven path to the supported models. This will require the real local source, input JSON and generated beta/C5 artifacts for direct execution and repair.

## Limitations

A static source gate and selected tests do not certify physics or reproduce the user's benchmark. Do not present this ZIP as a complete numerical repair. No existing production source file has been changed.
