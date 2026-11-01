# Full 16-model general-interaction symbolic physics audit

This is a **read-only** consistency audit of saved general T3 UV Lagrangians and
UV/EFT1 renormalisable RGBeta outputs. It does not invoke Mathematica,
Matchete or RGBeta, overwrite existing outputs, or change `pipeline.py`.

## What is verified

- All 16 T3 representation/hypercharge assignments in the existing catalogue.
- Current-source SHA-256 fingerprints and prior Matchete detailed diagnostic
  results (via `RunT3GeneralSymbolicStudy.validate_cached`).
- Hypercharge neutrality and Z2-evenness of **each** additional quartic and
  scalar bilinear monomial.
- The per-monomial number of independent scalar invariants as recorded in the
  independently generated manifest, checked against **both** the exported
  Matchete interaction-name list and the RGBeta quartic registry (including
  the **12 reused B/C couplings**).
- Presence of **50 extra quartic**, **8 extra bilinear**, and **6 extra Yukawa**
  interactions in the combined 16-model source inventories.
- SU(2) triangle conditions and U(1)Y gauge invariance for ordinary and extra
  Yukawas.
- UV and EFT1 beta-function **coverage**: general quartics, scalar masses,
  the 2 original Yukawas, and the 6 extra Yukawas in the appropriate UV models.
  Fermion-Yukawa and fermion-mass beta keys must be absent in EFT1.
- Exact SU(2) quadratic Casimirs `C2(d)=(d²-1)/4` and indices
  `T(d)=d(d²-1)/12`, independently compared with `GeneralSymbolic/group_invariants.csv`.
- Presence and limited status of previous Yukawa Ward/normalization reports,
  if those reports exist on disk.

## Not verified / not yet certified

**The audit does not establish the full Matchete–RGBeta physical coupling
normalisation or group-factor coefficient equality.** That needs action-level
comparison of native complex field components, scalar tensor bases, fermion
bilinears, potential-versus-Lagrangian signs, H.c., factorials and
conjugation. It cannot be inferred from matching registry names, group
Casimirs, matching operator counts or RGBeta internal tensor prefactors.

Mixed scalar masses in B(-1) and C(-1) remain a mass-basis matching limitation.
Full general-theory threshold matching, Wilson RGEs and C5 are outside scope.
The existing restricted pipeline is unchanged.

## Installation

Unzip the supplied archive at the repository root. Install exactly:

- `scripts/AuditFull16GeneralSymbolicPhysics.py`
- `tests/python/test_full16_general_symbolic_physics_audit.py`
- `docs/GENERAL_SYMBOLIC_PHYSICS_AUDIT.md`

## PowerShell (from the repository root)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
$env:PYTHONUTF8 = "1"

python -X utf8 -m scripts.AuditFull16GeneralSymbolicPhysics --self-test
python -X utf8 -m pytest -q tests/python/test_full16_general_symbolic_physics_audit.py

# After your full 16-model symbolic run completes:
python -X utf8 -m scripts.AuditFull16GeneralSymbolicPhysics
```

Audit the completed A(-2) pilot alone (uses the default saved data run ID):

```powershell
python -X utf8 -m scripts.AuditFull16GeneralSymbolicPhysics `
    --model A_-2 --study general_symbolic_pilot
```

If the symbolic sweep used a different `--run-id` and `--study`, pass the same
values to this audit. The audit never starts missing Wolfram jobs; absent,
stale or incomplete inputs are marked `Failed`.

Default audit outputs:

- `Reports/output/full16_general_symbolic/GeneralPhysicsAudit/README.md`
- `Reports/output/full16_general_symbolic/GeneralPhysicsAudit/audit.json`
- `Reports/output/full16_general_symbolic/GeneralPhysicsAudit/audit.csv`
- `Reports/output/full16_general_symbolic/GeneralPhysicsAudit/quartic_invariants.csv`
- `Reports/output/full16_general_symbolic/GeneralPhysicsAudit/group_invariants.csv`

A completed inventory should read `Structural inventory: 16/16` and
`full_16_expected_inventory_passed: true`. The physical mapping certificate
**must remain false** until an independent full component/action-level
normalisation proof is obtained. That is a scientific limitation, not an
execution failure.
