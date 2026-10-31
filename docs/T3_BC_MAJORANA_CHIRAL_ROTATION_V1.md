# B/C(-1) native Majorana chiral probe

Runs fresh native Matchete model builds and abstract CConj/projector identities for B(-1) and C(-1), plus an independent numerical 12x12 chiral transformation witness. The abstract identities do **not** identify the installed concrete charge-conjugation sign, RGBeta epsilon orientation, or the physical cross-package Yukawa dictionary. Production C5 and RGEs remain untouched.

Run from repository root:

```powershell
& .\scripts\RunT3BCMajoranaChiralRotation.ps1
python -m unittest tests.test_t3_bc_majorana_chiral_rotation_v1
```

Outputs under `output/studies/t3_bc_majorana_chiral_rotation_v1/`: `summary.json`, `matrix_witness.json`, `native_matchete.json` (if native Wolfram export succeeds), `native_matchete.log`.

If a native identity fails, retain the JSON and log; do not change the production Yukawa map based on the matrix witness alone.
