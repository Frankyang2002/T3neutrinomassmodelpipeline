# Workstream A: combined UV audit v1

## Scope and provenance

Read-only GitHub inspection at `a43acc3e1357cfdfad1f94cbd1b051ce20e4abe8` (2026-10-31 07:09 UTC).
This test works against existing `T3Fields.wl`, `T3Topology.wl`,
`LagrangianBuilder.wl`, and the independently catalogued scalar quartics.
No production source files are changed. No GitHub modification, matching,
RGEs, or numerics.

## Tests

- `quartics`: 16 models, opt-in diagonal fermion masses, six extra Yukawa
  matrices across the model catalogue, and all 50 independent complex extra
  scalar quartics. Additional scalar bilinears are *disabled*.
- `full`: adds all eight extra complex quadratic scalar bilinear structures.
  At B(-1) and C(-1), a holomorphic scalar mixing term produces a non-diagonal
  heavy scalar mass matrix. The existing builder is expected to report
  `RequiresMixedScalarMassDiagonalization`, not `Success`.

An expected B/C mass-basis block is *not* a Matchete full-UV pass. The audit
returns exit code 0 to distinguish an anticipated physical-basis limitation
from an unexpected code failure, but records the exact status separately.
This script checks the builder's summary and full ordinary `CheckLagrangian`;
it is not the detailed non-mass diagnostic audit. For detailed checks on
B/C(-1), see `Lagrangian/RunGeneralQuadraticMatcheteAudit.wl`.

This test does not certify UV operator independence/normalizations, exact Z2
registration in Matchete, or physical completeness. It does not test the
impact on matching or RGEs and must not activate these experimental flags
for Workstreams B, C or D.

## Execution

From the project root, use a separate PowerShell session. After installing
this script under `Lagrangian/`, set:

```powershell
$env:T3_EXPERIMENTAL_FLAVOR_MASSES = '1'
$env:T3_EXPERIMENTAL_GENERAL_FERMIONS = '1'
$env:T3_EXPERIMENTAL_GENERAL_QUARTICS = '1'
$env:T3_EXPERIMENTAL_GENERAL_QUADRATICS = '0'
```

Then execute `RunT3CombinedUVAudit.wl OUTPUT.json CLASS alpha quartics`.
For `full`, set only `T3_EXPERIMENTAL_GENERAL_QUADRATICS='1'` before invoking.
Supported `alpha` codes: `m4`, `m3`, `m2`, `m1`, `p0`, `p1`, `p2` (per class).

Successful `quartics` reports: `PassedCombinedUVCheckLagrangian`.
Successful `full` reports outside B/C(-1): `PassedCombinedUVCheckLagrangian`.
B/C(-1) `full` reports should show `ExpectedHeavyScalarMassBasisBlock`;
this must never be interpreted as permission to run `Match`.

The script exports per-model JSON, `*.free.inputform.txt`, and
`*.uv.inputform.txt` reports. After testing, deactivate the four
experimental environment variables in the session.
