# F-only LLHH — exact Laurent pole / finite / scale-log bridge

This batch **does not match the full physical C5**. It organizes the *already
executed* native simultaneous-F threshold outputs, without rerunning Matchete,
changing the LLSS release, integrating out either heavy scalar, or promoting
bare EFT expressions into an MSbar physical Wilson coefficient.

## Inputs (already present in the checkout)

- `output/general_observable/lo_certification_spinor/native_hard_e515d82d118284ba_native01/`
- `output/general_observable/lo_certification_spinor/raw_llhh_e515d82d118284ba_inventory01/`

The validator rehashes the native release summary, every model log/report/raw
Matchete InputForm, recomputes each raw inventory from the original native
InputForm and requires exact equality with every saved `terms.json`. It checks
all 16 native simultaneous heavy-fermion matches at `mu_F = 100000 GeV`.

## Verified algebraic identity and what it is not

All source-level kernel components in the frozen committed inventory can be
matched as exact Matchete InputForm triplets:

`K_r`, `K_r / ϵ`, and `K_r * Log[μbar2/MF^2]`,

including original +/− signs, all 3-flavour indices, the native SU(2) CG,
barred/conjugated couplings and separate PL/PR sectors. No group combines
terms with different flavour tensors or scalar-mixing spurions. A source with
missing, duplicate or non-identical pole, constant or log expression **fails**
closed. A numerator scalar numeric factor is **never** silently stripped.

The complete raw symbolic candidate can therefore be *organized* as

`sum_r K_r * (1 + 1/ϵ + Log[μbar2/MF_r^2])`.

The source expressions can contain non-hbar terms in the original native
**full** one-loop matching export; these are excluded from the LLHH candidate
inventory and remain auditable upstream. `K_r` itself contains `hbar`.

**This is not a true EFT hard-region coefficient:** the pole can contain an
EFT1 soft contribution, subdivergences and renormalization-scheme terms.
Dropping `K_r/ϵ` *formally* does not do MSbar or EFT1 subtraction. In
particular `1 + Log[...]` is not a certified physical short-distance result.
Nor does the native four-component `GammaCC` automatically fix the master Weyl
phase; this remains null. The heavy mass in the log is the flavour-dependent
`MF`, not unconditionally equal to the common matching scale `mu_F`.

## Output (new immutable directory only)

`output/general_observable/lo_certification_spinor/pole_log_e515d82d118284ba_bridge01/`

For each model (`A_m4`, etc.):

- `trios.json`: kernel, exactly three authenticated source indices and hashes,
  pole and logarithmic identities, strict physical-block flags.
- `PL_pole_residue_raw.inputform.txt`, `PR_pole_residue_raw.inputform.txt`
- `PL_finite_constant_raw.inputform.txt`, `PR_finite_constant_raw.inputform.txt`
- `PL_scale_log_raw.inputform.txt`, `PR_scale_log_raw.inputform.txt`
- `PL_formal_pole_drop_NOT_MSbar.inputform.txt`, `PR_formal_pole_drop_NOT_MSbar.inputform.txt`

Batch `summary.json` records model counts and all remaining physics blockers.
These are exact **raw** InputForm expressions, not numeric Wilson coefficients.

## Commands (PowerShell, from repository root)

```powershell
$env:PYTHONPATH = (Get-Location).Path
python -X utf8 -m unittest discover -s tests -p 'test_lo_spinor_*.py' -v
python -X utf8 scripts\BuildLOSpinorPoleLogBridge.py --check-only
python -X utf8 scripts\BuildLOSpinorPoleLogBridge.py --run-id bridge01
```

## Next physics calculations

1. Determine the UV hard and EFT1 soft one-loop contributions in **the same**
   regulator/subtraction scheme, including tree-level-LLSS insertions, scalar
   and SM loops, external-field/coupling counterterms, and pole cancellation.
   Verify no double counting and isolate the true local hard remainder.
2. Fix the native connected Dirac/Majorana `GammaCC` spinor/CG matrix-element
   dictionary to the chosen master Weyl basis. The native abstract identities
   are insufficient for an absolute signed factor.
3. Only after these first-threshold checks, run the EFT1 Wilson system with
   **both** scalar multiplets still present; integrate **both** scalars together
   at the later single `mu_S` and enter SMEFT.
