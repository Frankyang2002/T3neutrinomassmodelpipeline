# Simultaneous EFT1 scalar-quartic preflight batch

This batch **does not change physics**. It adds a standalone, read-only inventory audit and a combined PowerShell preflight. It does not claim to implement the complete quartic tensor.

## Project paths

- `scripts/AuditT3SimultaneousQuarticReadiness.py` — reads the 16-model registry and current simultaneous source; emits a JSON readiness report.
- `scripts/RunT3QuarticPreflight.ps1` — runs Stage 22 static checks and the ordinary/mixed solver self-tests.
- `docs/QuarticBatchStatus.md` — this document.

## Local execution

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
powershell -ExecutionPolicy Bypass -File scripts\RunT3QuarticPreflight.ps1
python -X utf8 -m scripts.AuditT3SimultaneousQuarticReadiness --require-ready
```

The second command intentionally exits with code `2` while the complete Wilson quartic input is not demonstrated. The JSON report is written to `output/general_t3_stage22/simultaneous_quartic_readiness.json`.

## Blocking technical issue

- The ordinary simultaneous solver packs the baseline renormalisable state and assembles its non-HHSS quartics through provisional pair contractions.
- The mixed B/C solver packs additional RGBeta quartics, but the Wilson sector calls the same provisional `general_initial_context` tensor constructor.
- The experimental quartic registry contains 50 additional model slots (38 new and 12 reused). The registry expressly labels the new coupling mappings `NotNormalised`.
- The Stage 22 audit expressly states that a static audit is not a tensor normalization proof. The native RGBeta invariant-expression exporter reports `component_tensor_exported=False`; a complete, convention-aligned tensor basis is not established by that file.

**Required before physical replacement:** reliable per-coupling real-scalar quartic tensor response in the same field and potential convention as the Wilson evaluator, and its running RGBeta state at each solver RHS evaluation. Do not silently substitute invariant counts or guessed factors.
