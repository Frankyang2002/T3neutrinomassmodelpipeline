# B/C native RGBeta `Lam` runtime diagnostic

## Status and purpose

The B(-1)/C(-1) `110101` quartic tensor map and reconstructed potential
normalization passed on the user's local project (see resolution log).
Those checks did **not** directly execute the installed RGBeta `Lam` function.
This package performs that next *diagnostic* runtime call and records the
returned expression structure, not an independent action-level coefficient
certificate. It is intentionally designed as a two-stage, fail-closed gate:

1. Wolfram builds each full RGBeta T3 model, loads its existing extra quartic
   source, inspects the two live native `$quartics` records, then calls the
   **installed RGBeta** `Lam[a,b,c,d]` (with independent symbolic indices).
   The expression's hash, sparse tensor shape, and the first 12 nonzero/default
   rules are saved. No coupling constants or beta functions are changed.
2. Python verifies the live evidence schema, two expected coupling symbols,
   Hermitian-conjugation settings, and the presence of an evaluated sparse
   tensor. It **rejects** incomplete diagnostics rather than inventing data.

**Limit:** A symbolic `Lam` call and its sparse structure do not yet establish
`V_RGBeta = - L_Matchete` by an *independent runtime polynomial contraction*.
The full map and physical neutrino matching remain uncertified.

## Files (install at project root)

This correction ZIP contains four complete files: the Wolfram probe, its tests, this document, and the resolution log. The Python verifier remains unchanged from the earlier probe package.

- `scripts/ProbeT3BCNativeRGBetaLam.wl` — live Wolfram exporter (corrected `PackageScope` resolution).
- `scripts/VerifyT3BCNativeLamProbe.py` — existing evidence gate from the previous ZIP (unchanged; not included in this fix).
- `tests/test_t3_bc_native_rgbeta_lam_probe.py` — schema and Wolfram symbol-resolution regression tests.
- `docs/T3_BC_NATIVE_RGBETA_LAM_PROBE.md` — this document.
- `docs/T3_GENERAL_SCALAR_MATCHING_RESOLUTION_LOG.md` — complete updated log.

## Exact PowerShell commands

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath 'C:\path\to\T3_BC_native_rgbeta_lam_scope_fix.zip' -DestinationPath . -Force
python -X utf8 -m pytest -q tests/test_t3_bc_native_rgbeta_lam_probe.py
wolframscript -file scripts/ProbeT3BCNativeRGBetaLam.wl output/audits/t3_bc_native_lam/native_lam_v2.json
python -u -X utf8 -m scripts.VerifyT3BCNativeLamProbe --input output/audits/t3_bc_native_lam/native_lam_v2.json --output output/audits/t3_bc_native_lam/native_lam_v2_verified.json
```

When rerunning, change `native_lam_v2` to a fresh ID. Neither command will
overwrite old reports. RGBeta and Wolfram must be installed locally. The
package makes **no** changes to GitHub or production model calculations.

If the Wolfram output reports `NativeLamEvaluationRequiresInspection`, please
supply the terminal output **and** the model's `lam_expression_head`,
`lam_expression_excerpt_inputform`, `sparse_array_count`, and
`lam_symbol` from the diagnostic JSON. The runtime representation can then
be inspected without guessing a Mathematica internal data layout.

## Physics follow-up

The next independent test is to recover the active `110101` quartic coefficient
from native `Lam`, including real/complex field projection, all permutations,
and both conjugate terms, and compare its explicit polynomial against the
existing Matchete source `PlusHc` contraction. Thereafter the other B/C
quartics and Yukawa CG vertices still require their own checks. **Passing
this diagnostic alone is not a physical coupling-map certificate.**

## 2026-10-31: `Lam` namespace resolution fix

The initial user-local execution failed before any model build:

```text
Cannot uniquely resolve installed RGBeta Lam symbol: {}
FileNotFoundError: output\audits\t3_bc_native_lam\native_lam_v1.json
```

**Root cause:** The old probe searched `RGBeta`` and `RGBeta`Private`` but
RGBeta declares `PackageScope["Lam"]` in `Package/FieldsAndCouplings.m`.
The documented Wolfram structured-package context for this declaration is
`RGBeta`PackageScope`Lam`. The verifier's file error was downstream of the
failed Wolfram probe. No physics calculation was performed by that run.

**Fix:** The full replacement Wolfram script uses the exact declared context
`RGBeta`PackageScope`Lam`, requires it to exist in `Names`, and requires
nonempty `DownValues` before invoking the installed function. The actual
`Lam` implementation remains untouched. A static Python regression test
checks that the previously broken lookup does not return.

**Verification status:** Python tests can verify the source and JSON gate;
the actual Wolfram runtime check must be rerun on the user's installed RGBeta.
Use the fresh `native_lam_v2.json` output to avoid overwriting previous files.
A successful structural probe is **not** a full potential or matching
certificate.
