# B/C Majorana flavour-action batch v1

## What was verified

The repository source `Lagrangian/model/T3Fields.wl` defines `NewFermion` as self-conjugate for zero hypercharge and odd SU(2) dimension (B(-1) singlet, C(-1) triplet). `Lagrangian/interactions/T3Topology.wl` writes the first vertex with `CConj[F]`, the second with `F`. The preceding native chiral probe passes B/C abstract identities but does **not** measure the absolute Matchete charge-conjugation orientation.

This batch independently evaluates a generic complex **two-component Majorana action** under `chi_old=U chi_new`:

`M_new=U.T M_old U`, `y1_new=y1_old U`, `y2_new=y2_old U`.

Both chiral sectors of the four-component Majorana spinor rotate with `U` and `U*` as appropriate; `CConj[F]=F`. The test preserves the mass bilinear and the two chiral Yukawa terms, and includes an intentionally wrong `y1 U*` negative control. 64 randomized complex tests are default.

## Scope and limitations

The same-U result is **conditional on choosing a shared left-Weyl Yukawa convention** and does not establish an absolute Matchete↔RGBeta coupling map, group-theory tensor signs or normalization. Therefore it is not enough by itself to edit `Numerical/experimental/T3BCRunningOneScaleTriangle.py` or certify C5. Neither source nor production outputs are modified.

The batch can consume `output/studies/t3_bc_majorana_chiral_rotation_v1/native_matchete.json` to check the observed native six-check evidence and preserve its uncertified status.

## Run

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$env:USERPROFILE\Downloads\t3_bc_majorana_action_batch_v1.zip" -DestinationPath . -Force
& .\scripts\RunT3BCMajoranaActionBatch.ps1
```

Output: `output/studies/t3_bc_majorana_action_batch_v1/summary.json`.
