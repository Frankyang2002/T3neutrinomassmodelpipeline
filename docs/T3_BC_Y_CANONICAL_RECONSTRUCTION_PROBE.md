# B/C Y-only Matchete: independent canonical full-EFT reconstruction probe

**Scope:** experimental and non-certifying. No production module is changed.

## Verified failure inputs

The 2026-10-31 user-uploaded B(-1) runtime report says:

- `SavedEFT1TreeHeavyScalarCheckFailed`.
- The saved EOMEFT is non-Hermitian **before** X/Y reclassification.
- All Matchete detailed checks except Hermiticity pass in the final saved-EFT check.
- Background + all LLSS and background + Higgs portals each pass Hermiticity.
- Background + scalar self-quartics fails, including individual coupling groups `lambda12Inv2`, `lambdaS1Inv1`, `lambdaS2Inv1`, `lambdaT3ExtraQ001201I1`.
- The existing independent renormalisable reconstruction has `algebraic_equality_proven=false`, but was **never independently checked** with Matchete because the previous script ran `CheckLagrangian[composite]` only *after* literal/relabeled EOMEFT equality.
- Its residual begins in `lambdaT3ExtraQ020020I1` with different index-contraction representations. This does not establish a physically different potential, but means a replacement is **not certified equal** to the serialized EOMEFT.

This stage fixes the **premature rejection of the independent reconstructed candidate as an experimental diagnostic**, not yet the underlying EOMEFT Hermiticity issue.

## Implementation

`Lagrangian/RunT3BCYOnlyCanonicalReconstruction.experimental.wl` is a full standalone copy derived from the precise user-uploaded `RunT3BCYOnlyMatching.experimental.wl` (original SHA256 `f0364a59b73b4f6cefc5e737f611b167946650367ac8e771a999e11d66dd3ead`). It does not need regenerating with the old Y-only builder. The script:

1. Reads the same first-stage saved canonical EOMEFT; checks the field/operator inventory and its failing Hermiticity without modifying the source file.
2. Builds all native X/Y renormalisable interactions from the original input, retaining **all self-quartics** and Higgs portals, and adds the exact 20 (B) / 24 (C) LLSS terms from the saved tree EOMEFT.
3. Records whether this alternative is **algebraically identical** to saved EOMEFT. A false result is retained as false, never overridden or interpreted as proof of field redefinition equivalence.
4. Checks this alternative independently with ordinary `CheckLagrangian` and the detailed Hermiticity, canonical normalisation, heavy-mass basis, gauge-invariance and field-definition checks. No `(L+L†)/2` fix, coupling deletion or weakened acceptance check is applied.
5. Only with all independent checks passing, permits an **explicit opt-in** experimental candidate to replace the failed in-memory EOMEFT.
6. `--check-only` exits before matching. `--match` can run tree and one-loop native Y-only matching and export raw EFT2 LLHH/LLHX/LLXX sectors, but is deliberately marked `YOnlyIndependentCanonicalRawEFT2NotCertified`.

**Critical:** Even if this alternative passes Matchete and produces raw EFT2, it is **not a certified repair** of the saved EOMEFT until its equivalence to the matching-derived EFT1 is established. There is no finite matching coefficient to inject into the numerical chain yet.

## Install paths

- `Lagrangian/RunT3BCYOnlyCanonicalReconstruction.experimental.wl`
- `scripts/RunT3BCYCanonicalReconstruction.py`
- `tests/test_t3_bc_y_canonical_reconstruction.py`
- `docs/T3_BC_Y_CANONICAL_RECONSTRUCTION_PROBE.md`
- `docs/T3_BC_Y_CANONICAL_RECONSTRUCTION_MILESTONE.md`

Outputs: `output/general_t3_mass_eigen_matchete/eft1_scalar_Y_only_canonical_probe/`.

## PowerShell commands

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$env:USERPROFILE\Downloads\t3_bc_y_canonical_reconstruction_probe_v1.zip" -DestinationPath . -Force
python -m unittest discover -s tests -p "test_t3_bc_y_canonical_reconstruction.py" -v
python -u -X utf8 -m scripts.RunT3BCYCanonicalReconstruction --all --preflight-only
python -u -X utf8 -m scripts.RunT3BCYCanonicalReconstruction --all --check-only
# Only if both models report IndependentCanonicalFullEFTCheckPassedEOMEquivalenceUnproved:
python -u -X utf8 -m scripts.RunT3BCYCanonicalReconstruction --all --match
```

The `--match` command requires a successful prior check-only certificate and identical SHA256 hashes for the WL file, parameter JSON, and first-stage tree source. Each run produces a report, complete Matchete log, and `.provenance.json`. To rerun, supply `--overwrite`; previously saved 2026-10-31 runs are not touched.

## Failure interpretation

- `IndependentCanonicalReassemblyFailedMatcheteChecks`: the native full-EFT reassembly **also** fails Matchete. Inspect `*_reconstruction_audit.json` field `independent_hermiticity`. No matching attempted.
- `IndependentCanonicalFullEFTCheckPassedEOMEquivalenceUnproved`: the alternative full Lagrangian is accepted by Matchete but canonical EOMEFT equality remains unproved. No matching attempted.
- `YOnlyIndependentCanonicalRawEFT2NotCertified`: one-loop raw EFT2 sector exports are present but cannot be promoted to physical Wilson matching.

C(-1) has not yet been run as of this package; its success is not inferred from B(-1).
