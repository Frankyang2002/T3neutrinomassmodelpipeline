# Direct UV all-16: native charge-support gate correction

## Verified defect

The 2026-11-02 all-16 local run completed 7 provisional evaluations and blocked 9
with `Native CG candidate charge support mismatch`. The old runner compared
`charge_resolved_loop_weights.json` topology charges with the **union** of all
nonzero Yukawa vertex charges. This is not a valid equality: source `weights`
are charges in the simultaneous Y1 x Y2 x T3MixCG contraction, whereas each
Yukawa can also have allowed components that do not participate in that
particular mixed topology. The latter components must not be rejected or
removed from the exact full-UV calculation.

The new check:

- Requires each source topology charge to be represented in BOTH ordinary
  Yukawa vertex charge supports.
- Requires independent `charge_channels.json` and
  `charge_resolved_loop_weights.json` source topology support to agree, when
  both are present.
- Records Yukawa-only charges separately, without removing any vertices.
- Retains the source scalar mixing support and existing channel, fermion,
  scalar, suppressed-term, and cancellation decompositions.
- Keeps source-to-Weyl phases, absolute Yukawa/portal signs and extra general
  interactions **uncertified**; the correction does not establish physical
  certification.

## Full replacement files (owned paths only)

- `Numerical/lo_certification/direct_uv/runner.py` — full replacement.
- `Numerical/lo_certification/direct_uv/CHARGE_SUPPORT_FIX.md` — new note.
- `tests/test_lo_direct_uv_charge_support.py` — new regression tests.

No `Numerical/direct_one_loop` engine or old report is modified.

## Installation (PowerShell at repository root)

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
$zip = "$HOME\Downloads\T3_LO_direct_UV_charge_support_fix.zip"
Expand-Archive -LiteralPath $zip -DestinationPath . -Force

python -m unittest discover -s tests -p "test_lo_direct_uv_*.py" -v
if ($LASTEXITCODE -ne 0) { throw "Direct UV tests failed" }

python scripts/RunLODirectUVCertification.py --project-root . `
  --output-dir output/general_observable/lo_certification_direct_uv/charge_support_fix_01
```

The separate `charge_support_fix_01` output path preserves the previous run's
`all_16_certification_status.json` and seven individual model reports. If that
path already exists, use a **different new subdirectory**; the runner explicitly
refuses to overwrite existing certification outputs.

## Tests and limits

In an isolated fixture containing the original core/comparison sources, all
22 `test_lo_direct_uv_*.py` tests passed (17 previous + 5 new). The new
integration test constructs both source-charge audit report formats and checks
16/16 provisional evaluations and 0 blocked with numerical synthetic inputs.
This is NOT an execution of your actual Windows frozen benchmark data.

The status for genuine source-certified mass matrices remains 0/16 until the
independent four-spinor/Weyl, scalar-potential and coupling dictionary witnesses
and any missing general interaction inputs are established. The exact UV
calculation still does not use an EFT heavy-mass expansion. The EFT threshold
sequence remains UV -> integrate ALL BSM fermions at mu_F -> EFT1 (BOTH scalar
multiplets retained) -> integrate BOTH BSM scalars at mu_S -> SMEFT.
