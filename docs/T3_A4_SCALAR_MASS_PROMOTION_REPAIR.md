# Local repair: A(-4) physical joint scalar threshold

Source inspected: GitHub `Frankyang2002/T3neutrinomassmodelpipeline`, main commit `80dcee89d9a7d9565e76013b4e5c78735b6193b6` (2026-11-10).
The repository was **read only**. These files are designed for local replacement, not a GitHub commit.

## Confirmed failure

The committed output `output/parallel_t3/full_general_final_one_batch/thesis_physical_scalar_continuation_v1/A_-4/debug/fresh_kernel_stage_2/stdout.log` contains:

- `CheckLagrangian::Hermiticity` (may be a separate issue that still requires verification)
- Warnings that `NewScalar1` and `NewScalar2` are defined Heavy but the input EFT still treats them as Light
- An attempted one-loop match which returns, but no authoritative Weinberg coefficient is extracted

The exact saved A(-4) F-boundary EFT1 tree has **one** `Coupling[MS1,{},1]` and **one** `Coupling[MS2,{},1]`, rather than EFT order 0, whereas the second matching kernel defines both scalars Heavy. This repair follows the mass-power-counting approach already implemented successfully for D(-2).

## Files

- `Lagrangian/RunT3FullGeneralTwoThresholdNative.wl` — full replacement; keeps F-stage all-three-flavour contraction but dispatches A scalar matching to a dedicated physical continuation, as D already does.
- `Lagrangian/RunT3AMinus4FullGeneralPhysicalScalar.wl` — new physical scalar continuation; fully rebuilds the Matchete registry and checks both scalar mass label promotions are reversible; checks the Lagrangian; executes one simultaneous S1+S2 one-loop scalar Match; exports intermediate WXF stages and scalar-soft Weinberg coefficient if successful.
- `tests/test_t3_thesis_final_one_batch.py` — full replacement; existing tests plus A-specific contracts.
- `tests/test_t3_scalar_physical_threshold_static.py` — dependency-free new static regressions.
- `docs/T3_A4_SCALAR_MASS_PROMOTION_REPAIR.md` — this document.

The existing D physical continuation is unmodified. The successful prior Batch 1 and Batch 2 outputs are not overwritten. The standalone new A file is invoked by the parent, not run directly.

## Local install from ZIP (PowerShell)

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$HOME\Downloads\T3_A4_PhysicalScalar_HeavyMass_Repair.zip" -DestinationPath . -Force
```

Or extract the ZIP manually into the repository root; preserve relative directories.

## Tests and one run

```powershell
python -X utf8 -m unittest discover -s tests -p test_t3_scalar_physical_threshold_static.py -v
python -X utf8 -m unittest discover -s tests -p test_t3_thesis_final_one_batch.py -v
python -X utf8 -m scripts.RunT3ThesisFinalOneBatch `
  --all `
  --run-id thesis_physical_scalar_massrepair_v2 `
  --debug
```

`--fit` should **not** be used: there is not yet a complete, finite, RG-improved C5.

## Expected results and status

Successful A scalar-soft matching will create:

- `output/parallel_t3/full_general_final_one_batch/thesis_physical_scalar_massrepair_v2/A_-4/A_joint_scalar_mass_power_counting_audit.json`
- `.../A_joint_scalar_lagrangian_audit.json`
- `.../A_EFT1_joint_heavy_scalar_match_input.wxf`
- `.../A_scalar_soft_only_raw_matched.wxf`
- `.../A_scalar_soft_only_GreensSimplify.wxf` (and three later processing stages)
- `.../A_scalar_soft_only_unrenormalized.wxf`

A successful scalar-soft stage will report `A_PhysicalScalarSoftOnly_HardInheritanceAndWilsonRGPending` internally. The outer Python driver intentionally labels it `Blocked_CompletionBlocked` and returns exit code 2, because physical full C5 is still absent. This is **not** evidence that the scalar match failed. Inspect `native_two_threshold_status.json` to distinguish the outcomes.

D remains `D_PhysicalScalarSoftOnly_HardInheritanceAndWilsonRGPending` until the independent Wilson evolution and hard inheritance have been computed and renormalised.

If the A `CheckLagrangian` audit still fails, inspect `A_joint_scalar_lagrangian_audit.json` and the full `native_two_threshold.log`. Do **not** symmetrise the Lagrangian or suppress this validation to produce a numerical coefficient: identify and repair any residual CG/hermitian-conjugation mismatch first.

## Physics and limitations

The required EFT order is preserved:

`UV -> (all three heavy fermion flavours, one threshold at 1e5 GeV) -> EFT1 with S1+S2 -> (both scalar multiplets, one threshold at 1e3 GeV) -> SMEFT`.

The Wilson tensor is an independent EFT1 coupling with exact tree definition `W12[p,q] = sum_{r=1}^3 y1[p,r] y2[q,r] / MF[r]` for diagonal physical heavy-fermion masses. Its own one-loop RG evolution is not yet included.

The remaining requirements for full-general physical neutrino masses in both models are: renormalised finite F hard Wilson boundary, complete EFT1 Wilson anomalous-dimension transport, inherited O(hbar) terms through the scalar threshold at tree order, finite subtraction without double counting, SMEFT Weinberg evolution, and fit of 18 real UV Yukawa entries. They are **not** implemented by this repair.
