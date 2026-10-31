# B/C restricted triangle: controlled heavy-Yukawa prescription comparison

## Scope

At the same running UV endpoint (1e5 GeV), re-evaluate the **archived Matchete 480-term restricted direct one-loop ordered coefficient** for each of the three Takagi fermion eigenmasses. All masses, scalar U, quartics and y2 are identical between prescriptions. Only y1 changes:

- production baseline: `y1 @ conj(V)`;
- same-U Weyl candidate: `y1 @ V`;
- both: `y2 @ V` and `V.T @ MF @ V = diag(m_i)`.

Use the stored full running report and native direct-source SHA256 to reject accidental mismatches. Recompute original C5 and require agreement with archived C5 to relative `2e-8`.

Outputs: `output/studies/t3_bc_convention_matching_batch_v1/{B_m1_conventions.json,C_m1_conventions.json,summary.json}`.

## Interpretation

Same-U is supported by previous Weyl action and RGBeta basis covariance, but the native four-spinor Matchete to RGBeta coupling/phase dictionary is still not certified. Both calculated spectra are **conditional diagnostics**, not physical neutrino-mass predictions. This does not perform full finite UV→EFT1→EFT2 matching, add independent mixed-quadratic native one-loop insertions, or alter the production triangle.

## Run

From repository root with the previous native/action/UV result files present:

```powershell
& .\scripts\RunT3BCConventionMatchingBatch.ps1
```

All files in this ZIP are new. Existing production and Workstream A files are untouched.

## Fix: real Takagi benchmark

The running UV input can produce a real (or effectively real) Takagi basis
`V`, in which case `y1 @ V == y1 @ conj(V)`. The comparison records
`rotation_prescriptions_numerically_distinguishable=false`, with zero C5
and spectrum difference, rather than throwing an error. This does **not**
prove the two field-theory prescriptions equivalent for generic complex
Majorana mass matrices. Test the generic complex case separately.
