# Workstream 4 — native basic scalar portal dictionary (v2)

## Reason for the v2 update

The first dictionary runner stopped with `Native portal CG tensors linearly
dependent` after the native Matchete exporter completed successfully.
The solver wrongly required **all raw native invariant tensors to remain
independent after the real-coupling Hermitian projection**. The real source
operator is `1/2 PlusHc[lambdaHiInvN * Hbar H Sbar S * CG]`; these are not
equivalent rank conditions.

Native sources inspected:

- `Lagrangian/interactions/ScalarPotential.wl`: real `lambdaHiInvN`
  and `1/2 PlusHc` for each native Higgs portal.
- `RGE/running/rgbeta/T3RGBetaModel.wl`: basic `lambdaHi`
  quartic with `delta(Hbar,H) delta(Sbar,S)`.
- `scripts/BuildT3DirectOneLoopInputs.py`: conditional diagonal scalar
  blocks `mSiSq - lambdaHi*v*v/2`.
- `scripts/RunLODirectUVNativePortalTensorExport.wl`: existing live
  registered tensors, unchanged in this batch.

## Exact operator problem

For raw registered native tensor `T_k`, with component order
`(Hbar,H,Sbar,S)`, define

\[
 P_{ijab}^{(k)} = \frac12\left(T_{ijab}^{(k)} + T_{jiba}^{(k)*}\right).
\]

Because all `lambdaHiInvN` are **real**, the operator comparison is the
real-coefficient linear system

\[
 \sum_k c_k P^{(k)}_{ijab} = \delta_{ij}\delta_{ab},
 \qquad c_k\in\mathbb R.
\]

The script stacks the **real and imaginary** parts of all components,
and uses SymPy's exact symbolic linear algebra. The native tensors are
read from the previously exported JSON, not regenerated or guessed.
Machine-precision input components are rejected.

The solver distinguishes:

1. **Unique exact embedding:** raw tensors independent, Hermitian-operator
   rank equals the number of native invariants, RGBeta singlet is in span.
2. **Non-unique exact embedding:** raw tensors independent, Hermitian rank
   drops, RGBeta singlet remains in the span. Output contains the exact
   nullspace, zero-projection tensor indices, and ONE deterministic
   pivot-column representative. That representative is **not a unique
   physical source coupling identification**.
3. **Blocked:** raw registered tensors already dependent, or the RGBeta
   singlet is outside the projected real Hermitian span, or source hashes
   or frozen blocks disagree. The error specifies `class/scalar`, exact
   ranks and zero projected indices where relevant.

For valid operator embeddings, the conditional scalar mass shift is

\[
 \Delta A_i = -\frac{v^2}{2}\lambda_{Hi}\,I_{d_i}
\]

under the existing declared `V=-L` and `H^0=v/sqrt(2)` conventions.
The all-16 frozen scalar block audit is read-only. No numerical neutrino
mass matrices, production inputs or EFT thresholds are changed.

**This is not full physical Matchete/RGBeta source-to-Weyl certification.**
The native four-spinor phases, MF convention, RGBeta scalar-action factors,
T3 mixing portal, optional couplings and EFT thresholds remain unresolved.
Physical model certification remains `0/16`.

## Files (complete replacements)

- `Numerical/lo_certification/direct_uv/native_portal_dictionary.py`
- `scripts/RunLODirectUVPortalDictionary.py`
- `tests/test_lo_direct_uv_native_portal_dictionary.py`
- `Numerical/lo_certification/direct_uv/NATIVE_PORTAL_DICTIONARY.md`

The native Wolfram exporter is unchanged. **Do not rerun Matchete.**

## PowerShell

From the project root after extracting the ZIP:

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline

.\.venv\Scripts\python.exe -m pytest `
  tests\test_lo_direct_uv_native_portal_dictionary.py `
  tests\test_lo_direct_uv_native_portal_export_guard.py -q

.\.venv\Scripts\python.exe `
  scripts\RunLODirectUVPortalDictionary.py `
  --project-root . `
  --native-json output\general_observable\lo_certification_direct_uv\native_portal_tensor_01\native_five_class_portal_tensors.json `
  --output-dir output\general_observable\lo_certification_direct_uv\native_portal_dictionary_02
```

The output is
`output/general_observable/lo_certification_direct_uv/native_portal_dictionary_02/native_portal_operator_dictionary.json`.
A new output directory avoids accidental overwrite of prior reports.

Local synthetic tests cover exact unique and non-unique embeddings,
nullspace witnesses, raw-dependence rejection, missing-singlet span,
source/hash guards, synthetic all-16 scalar blocks, and assembly-guard
regressions. They **do not** stand in for running the user's actual
native five-class JSON (not present in this packaging environment).

## EFT boundary (unchanged)

UV input at `mu0=10^7 GeV` → integrate **all** BSM fermions together at
`muF=10^5 GeV`, preserving **both BSM scalars and all SM fermions** in EFT1 →
integrate **both BSM scalars** at `muS=10^3 GeV` → SMEFT.
This batch implements none of those matching or running steps.
