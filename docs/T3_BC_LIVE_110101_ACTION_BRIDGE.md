# B/C(−1): `110101` exact action bridge gated by live evidence

## Purpose and verified prerequisites

The user-local 2026-10-31 Wolfram/Matchete and Python runs reported:

- Installed RGBeta v1.2.0: B(−1) 64/64 and C(−1) 64/64 direct/conjugate native `Lam` component coefficients verified with **zero** mismatch (`native_lam_coeff_v4.json`).
- Installed Matchete: both native `110101` interactions individually pass `CheckLagrangian` and their sum passes for B(−1) and C(−1) (`matchete_110101_v2.json`).
- Earlier Stage22 native Matchete tensors and live RGBeta SU(2) covariance exports give two independent tensor bases in ordered species `H, Hdag, S1dag, S2dag`.

This follow-up checks a precise **action identity in the exported common SU(2) component basis**, **only after** verifying both live reports. The audit is read-only and does **not** rerun Wolfram, use Matchete `Match[]`, or modify an RGE.

## Exact map

For Matchete component tensors `M1,M2`, RGBeta invariant tensors `R1,R2`, and ordered native coupling vectors `h,g`, the exact tensor identity is

\[
(R_1,R_2)=(M_1,M_2) A,\quad
A=\begin{pmatrix}1/2&1\\-\sqrt3/2&0\end{pmatrix},
\]

and the **source-convention signed** action identity is

\[
\mathcal L_M^{110101}=-V_R^{110101},\qquad
h=-A g.
\]

The installed RGBeta `Lam` audit separately confirmed its ordered coefficient `6/24=1/4`, which is **internal** to its `Lam` real-scalar tensor; it is **not** a further factor on `g` in the physical potential. The bridge checks all 16 ordered SU(2) components for each of the two independent `g` couplings, both the direct and conjugate expressions: **64 coefficient identities/model**. Wrong sign, spurious `/4`, and swapped Matchete contractions are negative controls.

## Project paths

- New: `scripts/AuditT3BCLive110101ActionBridge.py`
- New: `tests/test_t3_bc_live_110101_action_bridge.py`
- New: `docs/T3_BC_LIVE_110101_ACTION_BRIDGE.md`
- Updated full replacement: `docs/T3_GENERAL_SCALAR_MATCHING_RESOLUTION_LOG.md`

Existing dependencies (from earlier user-installed audit packages):

- `scripts/AuditT3BC110101MassBasisMap.py`
- `scripts/AuditT3BCNativeQuarticMassBasis.py`
- `scripts/VerifyT3BCNativeLam110101Coefficients.py`
- `scripts/VerifyT3BCMatchete110101Action.py`

Default inputs are the saved live B/C v2 Matchete and v4 RGBeta reports, the native Stage22 Matchete/RGBeta tensor records, and the unchanged source `Lagrangian/interactions/GeneralScalarOperators.wl`. No new Mathematica execution is required.

## PowerShell commands

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "C:\path\to\T3_BC_live_110101_action_bridge_v1.zip" -DestinationPath . -Force
python -X utf8 -m pytest -q tests/test_t3_bc_live_110101_action_bridge.py
python -u -X utf8 -m scripts.AuditT3BCLive110101ActionBridge --run-id bc_live_action_bridge_v1
```

The new report will be written to
`output/audits/t3_bc_110101_live_action_bridge/bc_live_action_bridge_v1/summary.json`.
The audit refuses to overwrite an existing report. Any mismatch terminates without writing a successful summary.

## Evidence boundary

The saved Matchete live JSON contains the actual selected `PlusHc` terms and successful `CheckLagrangian` return values but **not the independently expanded numeric components of the live `CG[...]` symbols**. The Stage22 native component tensors were exported separately by Matchete's invariant-basis exporter, which uses the same source tensor builder. Therefore the bridge establishes **exact agreement in the exported source/native component basis, backed by two separately passing live package gates**. It does **not** directly subtract the two fully evaluated runtime actions in one common native Wolfram session.

The following remain false: independent live Matchete CG component-to-action equality; direct live cross-tool action subtraction; mixed-doublet mass-basis `HeavyMassBasis`/`Match[]`; full sequential EFT matching; physical neutrino-mass certification. Do not use this report to set existing production certification flags to true. The next genuinely independent target is extracting the **live Matchete registered CG components** and comparing those to its exported Stage22 tensors, before treating the signed map as a physically certified cross-tool map.
