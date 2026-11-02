# v6 native Matchete self-triplet Ward fix (rough numerical benchmark)

This fixes the four v6 **runtime** Weinberg SU(2) projection failures in
`A_-2`, `D_+0`, `E_-2`, and `E_+0` by correcting a verified inconsistency
between the v6 scalar rotation and the Matchete source Lagrangian. This is a
candidate remedy; the actual all-16 repository run must confirm whether the
four projection residuals are now within tolerance.

The original Matchete source `Lagrangian/RunT3GeneralSelfMassBasis.wl`
explicitly defines, for each identical SU(2) component index `k`,

    S[k] = (T3ScalarMinus[k] + i*T3ScalarPlus[k])/sqrt(2).

No SU(2) spherical-index intertwiner `W` acts on this field definition.
The previous v6 implementation incorrectly used the **separate** spherical
prototype `T3SelfBilinearMassBasis.real_intertwiner(dim)` to convert these
Matchete components. The error was invisible for singlets (dim=1) and
manifested in the four special triplet models (dim=3).

The corrected v6 retains the native `Minus` and `Plus` scalar field definitions
at **both** muF and muS. It does not rephase the already matched LLSS
boundary when a running quadratic coupling has a different sign or phase.
Instead it builds a quadratic mass matrix in the fixed named scalar basis.
For the declared Matchete Lagrangian convention

    L_int = +mu (S^dag)^2 + h.c.

and `S=(Minus+i Plus)/sqrt(2)`, the quadratic potential Hessian is

    M2 = [[(m2-2 Re(mu)) I_d,   -2 Im(mu) I_d],
          [  -2 Im(mu) I_d,    (m2+2 Re(mu)) I_d]].

Its eigenvalues are `m2 +/- 2|mu|`, each repeated `d` times. For complex
`mu`, off-diagonal mixing **is not discarded**: the existing joint-scalar
bubble diagonalizes the full real mass matrix along with all other scalar
components at one shared scalar threshold muS. For negative real mu,
`Minus/Plus` labels remain fixed rather than being silently swapped.

The script `RunT3JointScalarAll16V6` and its CLI are unchanged. Only the
full replacement files below are included:

    Numerical/joint_scalar_matching/special_mass_v6.py
    tests/test_t3_joint_scalar_special_mass_v6.py
    Numerical/joint_scalar_matching/README_NATIVE_SELF_WARD_V6_FIX.md

No other files should be overwritten. The tests include a native-index
mapping check for dim=3, the exact quadratic potential for positive,
negative and complex mu, and a fixed-field-basis eigenvalue check. Local
regression: **58/58 tests passed**.

PowerShell, from the existing Windows repository:

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
$zip = Get-ChildItem "$HOME\Downloads" -Filter 'T3_joint_scalar_v6_native_ward_fix*.zip' -File |
  Sort-Object LastWriteTime -Descending | Select-Object -First 1
if (-not $zip) { throw 'Download the native Ward fix ZIP first.' }
Expand-Archive -LiteralPath $zip.FullName -DestinationPath . -Force
python -m unittest discover -s tests -p 'test_t3_joint_scalar_*py' -v
python -m scripts.RunT3JointScalarAll16V6 `
  --allow-exploratory-map `
  --allow-exploratory-special-bases `
  --output-dir 'output\general_observable\lo_joint_scalar_matching\exploratory_all16_v6_native_fix'
```

The runner deliberately refuses to overwrite an existing summary, so use
this **new** output directory. Check the final line `Conditional C5: X/16;
errors=Y` and the `summary.json` model statuses. No physical absolute
GammaCC/Weyl normalization, certified Matchete-RGBeta quartic alias, or
complete one-loop C5 is claimed.

This fix is narrowly targeted to the observed v6 projection failures;
any remaining Ward failure must be investigated, **not suppressed**.
