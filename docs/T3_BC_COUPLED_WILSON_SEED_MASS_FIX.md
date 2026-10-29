# B/C coupled EFT1 Wilson: symbolic fermion-mass substitution fix

## Verified failure

`independent_tensor()` in `Numerical/experimental/T3BCCoupledWilsonTensor.py` previously substituted `sp.Symbol('MF')` into a Matchete-derived Wilson seed that uses `sp.Symbol('MF', nonzero=True)`. SymPy treats these as different symbols, leaving factors such as `0.125*I/MF` unresolved and blocking B(-1) and C(-1).

## Resolution

The full replacement module now identifies the actual `y1`, `y2`, and `MF` symbol objects from each seed expression before the three Yukawa probes. It rejects unknown symbols or multiple symbols with identical names but incompatible assumptions. The physics RGE, mass-basis flavour contraction, threshold scheme, and numerical integrator are unchanged.

## Installation

Extract the ZIP to the repository root, overwriting only these two files:

- `Numerical/experimental/T3BCCoupledWilsonTensor.py` (full replacement)
- `tests/test_t3_bc_coupled_wilson.py` (full replacement, two regression tests added)

This markdown file is new: `docs/T3_BC_COUPLED_WILSON_SEED_MASS_FIX.md`.

## PowerShell

```powershell
python -X utf8 -m pytest -q tests/test_t3_bc_coupled_parser_api.py tests/test_t3_bc_coupled_wilson.py
python -m scripts.RunT3BCCoupledWilsonSequential --model B_-1 --run-id continuous_v1 --resume
python -m scripts.RunT3BCCoupledWilsonSequential --all --run-id continuous_v1 --resume
```

If the runner blocks on a parity check, retain that error message and inspect the underlying physical normalization; do not bypass the check. Existing outputs are not overwritten.

## Validation

The corrected full module and replacement tests compiled; `tests/test_t3_bc_coupled_wilson.py` passed **24 tests**, with **1 repository-dependent skip** in the isolated environment. The complete B/C numerical run has not yet been executed on the user's checkout.
