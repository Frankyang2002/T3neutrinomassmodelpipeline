# Batch 8 — gated B/C canonical reconstruction for scalar-only matching

## Inputs already established

The saved mixed B(-1), C(-1) fermion-first tree `EOMEFT` and native X/Y mass
basis are available; no first-stage fermion matching is repeated. The local
`mixed_hard_06` evidence has already passed the v7 independent component
checker (`NativeCGComponentConditionsVerifiedNoMatch`, 137 terms in 25 scalar
sectors, 32 native cross-CG entries). This is **not** Matchete's full
Hermiticity check or a scalar-threshold matching result.

The v5/v6 `CheckLagrangian` on the **archived** complete EOMEFT failed only
Hermiticity. The independently constructed native EFT1 was not tested in v6
because the code gated it on syntactic equality, although the subsequent
native tensor-component audit established physical equality for the audited
representation. Do not change the saved `EOMEFT` or replace it with
`(L + HC[L])/2`.

## New full project file paths

- `scripts/RunT3MixedScalarGatedBatch8.py` — **new standalone launcher**, uses
  existing `scripts/RunT3MixedScalarHardBatch6.py` and
  `scripts/RunT3MixedComponentAuditBatch6.py` unmodified.
- `tests/python/test_t3_mixed_batch8_gated_scalar.py` — **new tests**.
- `README_T3_SEQUENTIAL_BATCH6.md` — this **full replacement** document.

No source under `Lagrangian/` is edited. No GitHub files are modified. Prior
outputs `mixed_hard_01` through `mixed_hard_06`, all saved first-stage F
thresholds, and all eight ordinary-model results are left untouched.

## Gate sequence

1. Read the original `mixed_hard_06/component_equivalence_audit.json` and
   `mixed_hard_06/settings.json`. Require the recorded 137-term/25-sector
   certificate, the unmodified source hashes, and the previously certified
   saved and reconstructed EFT1 expression SHA-256 digests.
2. Independently rerun **Python only** the established exact SU(2) solver for
   all 25 sectors and the existing native-array verifier for 32 components.
   No Wolfram/fermion matching occurs at this stage.
3. Build a **temporary**, instrumented copy of the established B/C native
   mixed-doublet builder, stopping before its first fermion `Match` call.
   The original `Lagrangian/RunGeneralT3FermionFirstEFT1Tree.wl` stays intact.
4. In the native kernel, remove F, preserve its real mass parameter, promote
   X and Y to Heavy, and rebase their entire free actions from the first-stage
   Light assignment. Recreate the exact native scalar potential plus all saved
   LLSS operators, without dropping any quartic.
5. Match the *fresh* rebased saved EFT1, the *fresh* reconstructed candidate,
   and the *fresh* two cross-CG arrays against the **three certified v6 witness
   files**. Reject mismatches. The saved EFT1 and reconstruction are allowed
   to differ syntactically from **each other** because their component identity
   was separately proved; the current instances must each match their own
   archived witness.
6. Independently run Matchete `CheckLagrangian` and its detailed output on
   the reconstructed candidate, requiring all five checks:
   `Hermiticity`, `CanonicallyNormalized`, `HeavyMassBasis`, `GaugeInvariant`,
   `AllObjectsDefined`. Any false/missing result blocks `Match`.
7. **Only then** run scalar-only `Match[tree, EFTOrder -> 5, LoopOrder -> 1]`,
   with X,Y as the only heavy fields. The existing Batch 6 hard-only extraction
   writes `hard_c5_candidate.txt` and retains the same provisional labels.

The v8 `hard_scalar_reconstruction_audit.json` distinguishes component
identity from Matchete algebraic equality. It must not claim
`algebraic_equality_proven=True` unless the algebraic test itself passed.

## PowerShell

```powershell
Set-Location 'C:\thesistemp\T3neutrinomassmodelpipeline'

Expand-Archive `
  -LiteralPath "$HOME\Downloads\T3_Batch8_GatedScalarMatching.zip" `
  -DestinationPath . -Force

python -X utf8 -m pytest -q `
  tests/python/test_t3_mixed_batch6.py `
  tests/python/test_t3_mixed_batch6_component_audit.py `
  tests/python/test_t3_mixed_batch8_gated_scalar.py

# Recheck prior component proof against current source hashes; no Wolfram.
python -X utf8 -m scripts.RunT3MixedScalarGatedBatch8 `
  --component-run mixed_hard_06 --preflight-only

# NEW run; does NOT recompute the first fermion threshold.
python -X utf8 -m scripts.RunT3MixedScalarGatedBatch8 `
  --component-run mixed_hard_06 `
  --run-id mixed_hard_08
```

If and **only if** both native hard candidates export successfully:

```powershell
python -X utf8 -m scripts.RunT3MixedFixedOrderBatch6 `
  --hard-run mixed_hard_08 --check-only

python -X utf8 -m scripts.RunT3MixedFixedOrderBatch6 `
  --hard-run mixed_hard_08 `
  --run-id mixed_fixed_08 `
  --conditional-quartic-map
```

## Diagnostic files and output paths

- `output/general_sequential_pilot/mixed_scalar_hard_batch6/mixed_hard_08/component_gate.json`:
  current audited SHA-256 witness paths, source component run, and false
  physical-certification flags.
- `.../mixed_hard_08/B_m1/hard_scalar_reconstruction_audit.json` and the
  analogous `C_m1` file: `component_equivalence_proven`, five `checks`, and
  `independent_matchete_check_passed`.
- `.../mixed_hard_08/{B_m1,C_m1}/native_result.progress.json`:
  checkpoints, including `batch8-component-certificate-replay-returned` and
  `batch8-reconstruction-independent-CheckLagrangian-returned`.
- Conditional numerical outputs, only on hard success:
  `output/general_sequential_pilot/mixed_fixed_order_batch6/mixed_fixed_08/`.

Any native failure is a **BLOCKED** result, not evidence that scalar matching
ran. The component gate can succeed while the independent reconstruction's
Matchete check still fails. No Matchete runtime is installed here; local
Wolfram matching is not yet validated.

## Physics limitations

The required EFT ordering is

`UV → integrate out ALL fermions → EFT1 (SM+X+Y) → integrate out BOTH scalars → SMEFT`.

Passing the component and native Matchete gates authorizes the attempted
second-stage scalar one-loop match, not a completed physical neutrino-mass
calculation. The source has a single symbolic MF; three independent Majorana
masses are a separate numerical lift, not independently matched here.
Inherited one-loop terms, EFT1 counterterms/soft subtraction, and full MS-bar
pole cancellation are **not** certified; neither is the final physical C5.
