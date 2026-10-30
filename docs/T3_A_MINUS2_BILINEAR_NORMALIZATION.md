# A(-2): native scalar bilinear normalization and independent C5 loop-shape test

This is a **focused check**, not an implementation of the full general-scalar
physical neutrino-mass matching. It reuses the 16-model symbolic run without
modifying existing RGEs, fitted configurations or neutrino outputs.

## Field content and operator

A(-2) has `S1 ~ (1,-1)`, `S2 ~ (3,0)`, `F ~ (2,-1/2)` under
`SU(2)_L x U(1)_Y`. The additional allowed quadratic is

`(S2dag)^T C (S2dag) + h.c.`

with complex coefficient of mass dimension two.

In the opt-in Matchete Lagrangian the operator is defined by `PlusHc`:

`L_M contains +mu_M (S2dag)^T T_M (S2dag) + h.c.`

The new Wolfram helper asks the **installed Matchete** for the actual
symmetrized quadratic invariant tensor, and verifies
`T_M = k_native I3` exactly. It does not assume `k_native = 1`.

RGBeta's existing experimental convention is the scalar **potential**

`V_R contains +mu_R (S2dag)^T C_spherical (S2dag) + h.c.`

where the spherical triplet metric `C_spherical` has elements
`C[m,n] = (-1)^(1-m) delta[m,-n]`. The independent real-basis
intertwiner obeys `W.T C_spherical W = I3`. Relative to the normalized
native invariant, the corresponding mass-parameter relation is

**`mu_R = -k_native * mu_M`**

(the minus sign is from `L = kinetic - V`). Because no extra `1/2` is
included, the real scalar squared masses are

`m_±² = mS2Sq ± 2 abs(mu_R)`.

The script tests the complete quadratic potential on random real field
vectors and compares with the analytically diagonal form. This is an
independent test of the **quadratic metric and its coefficient**, not an
independent certification of the full Matchete/RGBeta triplet generator,
Yukawa, or quartic vertex intertwiners.

## One-loop C5 mass-response test

The current repository already contains a conditional self-bilinear
spectral C5 continuation (`T3SelfBilinearMassDependentKernel.py`). Here,
the saved one-scale Matchete A(-2) coefficient is independently compared,
**by mass-dependent ratios only**, against

`J(MF²,MS1²,MS2²) = integral_0^infinity dt t / ((t+MF²)(t+MS1²)(t+MS2²))`.

This confirms whether the existing one-scale kernel has the correct
three-propagator scalar-mass response needed for normal/anomalous spectral
splitting. It deliberately does **not** claim the absolute C5 normalization,
the new HH quartic sign or mapping, or sequential EFT matching.

## Files (new; no replacements)

- `Lagrangian/RunA2NativeBilinearTensor.wl`: extract native A(-2) bilinear CG tensor.
- `scripts/VerifyA2BilinearNormalization.py`: inspect saved symbolic UV/EFT1 betas,
  derive and check the mass map, independently test the one-scale C5 mass dependence.
- `tests/test_a2_bilinear_normalization.py`: offline unit tests plus one saved-source test.
- `docs/T3_A_MINUS2_BILINEAR_NORMALIZATION.md`: this document.

## PowerShell

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath 'C:\path\to\T3_Aminus2_bilinear_normalization.zip' -DestinationPath . -Force

python -X utf8 -m pytest -q tests/test_a2_bilinear_normalization.py
python -u -X utf8 -m scripts.VerifyA2BilinearNormalization --symbolic-run-id all16_v1 --run-id a2_norm_v1
```

New output:

`output/studies/general_scalar_extended/a2_bilinear_normalization/a2_norm_v1/certificate.json`

and `native_matchete_bilinear.json` plus the Wolfram log.
The script fails closed if any native metric, sign, existing RGE, source
response, or mass-eigenvalue check fails. Use a **new** run ID to retry;
it does not overwrite prior results.

## What still remains after a successful run

- Derive the physical Matchete<->RGBeta **full component** map (not just the
  quadratic invariant), including `SU(2)` generator and interaction tensors.
- Certify normalization/conjugation of the independent `H H S1dag S2dag`
  quartic responsible for A(-2)'s anomalous C5 channel.
- Compute one-loop C5 with **all interactions** in a consistent scalar
  mass basis, then perform sequential matching/running in the same scheme.
- Only then connect the additional complex `mu2` to the production
  numerical state and fit parameters.

The status field explicitly says which restricted checks passed; the
report always sets `independent_one_loop_mass_eigenstate_matching_performed`
and `complete_sequential_EFT_matching_performed` to `false`.
