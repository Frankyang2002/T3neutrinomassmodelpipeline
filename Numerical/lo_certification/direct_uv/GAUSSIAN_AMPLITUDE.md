# Batch: all-heavy-fermion Gaussian LLSS and four-spinor witness

## Exact ownership and inputs

All code is additive under `Numerical/lo_certification/direct_uv/` and
`scripts/RunLODirectUVAmplitudeBridge.py`; tests live under
`tests/test_lo_direct_uv_gaussian_amplitude.py`. Only the directory
`output/general_observable/lo_certification_direct_uv/<new_run>/` is written.
Existing `Numerical/direct_one_loop/core.py`, Matchete/RGBeta, frozen
benchmarks and all EFT codes remain unchanged.

The runner consumes the *existing provisional* component inputs at
`output/general_observable/lo_direct_one_loop/provisional_inputs/<MODEL>.json`.
It checks source definitions in:

- `RGE/running/rgbeta/T3RGBetaModel.wl`
- `Lagrangian/interactions/T3Topology.wl`
- `Lagrangian/model/T3Fields.wl`

It also verifies the original benchmark byte hash and the UV input scale.
The input components **remain provisional**: valid algebra is not proof of
physical Matchete / RGBeta / Weyl normalisation.

## Chosen 4-component/Weyl convention

In the explicit Weyl basis,

```
gamma^0 = [[0,I],[I,0]]
gamma^2 = [[0,sigma_2],[-sigma_2,0]]
C = +i gamma^2 gamma^0 = diag(epsilon,-epsilon)
epsilon = [[0,1],[-1,0]]
F=(chi_alpha, -epsilon bar(xi))
CConj(F)=(xi_alpha, -epsilon bar(chi))
```

Independent **exterior/Grassmann algebra** proves:

\[
\bar L F^c=(\ell\epsilon\chi)^\dagger,\qquad
\bar L F=(\ell\epsilon\xi)^\dagger,
\]
\[
\bar F F=(\chi\epsilon\xi)+(\chi\epsilon\xi)^\dagger.
\]

For Majorana `F=F^c` set `xi=chi`. The factor
`chi epsilon chi = 2 chi_0 chi_1` is verified; the 4-spinor Majorana mass
carries the conventional **1/2**. A wrong sign for `C` is rejected by a
negative test. This is a **chosen numerical gamma basis**, NOT a measured
matrix representation of the installed Matchete `GammaCC` or RGBeta `eps`.

The source Matchete vertices are `y1 Bar[L] CConj[F] S1` and
`y2 Bar[L] F S2dag`, both with H.c. In the declared basis their H.c.
contains `conj(y1) (l epsilon chi) S1dag` and
`conj(y2) (l epsilon xi) S2`. The direct one-loop engine defines
`L_Y=-nu psi (Y phi + Z phi*) +h.c.`; therefore the *candidate* corresponding
input tensors carry the opposite coefficient sign, plus native CG conjugation.
**This conditional identification is not a source-certified numerical input.**

In this source-H.c. orientation, with `h_i=conjugate(y_i_M)`, the
rephasing invariant is

\[
h_1 h_2 M_F^*\lambda_{T3}
= y_{1M}^* y_{2M}^* M_F^*\lambda_{T3M}.
\]

This differs in apparent conjugation from the earlier
`spinor_mass_dictionary.py` *conditional re-labelled action* because the
latter starts from a different choice of holomorphic scalar/Yukawa names.
Neither statement fixes the installed package-to-package dictionary.

## Complete tree Gaussian source tensor

For ALL left-handed heavy Weyl components, without separating heavy families
into different thresholds, use the **existing numerical core** convention:

\[
\mathcal L_{\mathrm{UV}}\supset
-\tfrac12\psi^T M_F\psi-\psi^T J+\mathrm{h.c.},\quad
J_i=\sum_{a,s}\nu_a(Y_{ais}\phi_s+Z_{ais}\phi_s^*).
\]

Solving all heavy equations together gives

\[
\psi=-M_F^{-1}J,\qquad
\mathcal L_{\rm EFT1,tree}=\tfrac12J^TM_F^{-1}J+\mathrm{h.c.}
\]

and four scalar orientation sectors

\[
K_{abst}^{\sigma\tau}=V_{ais}^\sigma (M_F^{-1})_{ij}
 V_{bjt}^\tau,\qquad V^+=Y,\quad V^-=Z.
\]

The **full coefficient** `K` is saved, with `1/2` occurring in the
Lagrangian, not multiplied into saved tensor coefficients. Every ordered
external scalar orientation `YY`, `YZ`, `ZY`, `ZZ`, every source vertex pair,
and all heavy eigenstates are retained. Contributions are **summed**, not
filtered by size. For unequal source pairs `A x B`, both `AB` and `BA`
are present exactly once each in `K`; there is no spurious factor 2.

The independent Takagi/eigenstate construction is checked against the direct
matrix inverse. The computation checks the Schur complement, exchange
symmetry, U(1)_em charge conservation and field-rephasing covariance. For
nearly zero tensors, reconstruction uses the scale of the uncontracted
amplitude to avoid meaningless ratios from cancellation.

## EFT threshold caveat: do not compare different scales

The frozen benchmark is evaluated at `mu_0 = 1e7 GeV`. Existing independent
fermion-boundary source releases are evaluated at `mu_F = 1e5 GeV`.
This code **does not run** renormalisable UV couplings down to `mu_F`.
The new `K` is an exact **algebraic tree source tensor at the supplied UV
input scale**, NOT a physically matched `mu_F` Wilson coefficient.
The optional `eft1_release_gate` records the scale mismatch and unverified
native spinor normalisation, without equating the tensors. No heavy SM
fermion is eliminated.

Required full EFT chronology (not executed by this batch):

1. UV at `mu_F`: eliminate **ALL heavy BSM fermions together**, retain both
   BSM scalars and every SM fermion; include full hard matching.
2. EFT1 run from `mu_F` to `mu_S`.
3. At `mu_S`: eliminate **BOTH heavy BSM scalar multiplets together**.
4. SMEFT, with the physical neutrino `C5` convention and running.

### Explicit non-claims

- No actual native Matchete four-point amplitude or installed `GammaCC` sign
  was measured by this Python batch.
- No physical RGBeta/MATCHETE coupling dictionary was measured.
- No hard one-loop `LLHH` fermion matching was included; this may be nonzero.
- No RG running, scalar-loop matching or `C5` is generated.
- Extra source-allowed Yukawas are not assigned invented input coefficients.
- Unit tests use **synthetic all-16 fixtures**; their numerical amplitudes
  are not physical benchmarks.

## Run (PowerShell from repository root)

```powershell
Expand-Archive -LiteralPath "$HOME\Downloads\T3_LO_direct_UV_gaussian_amplitude_batch.zip" -DestinationPath .
python -m unittest discover -s tests -p "test_lo_direct_uv_*.py" -v
python scripts/RunLODirectUVAmplitudeBridge.py --project-root . `
  --output-dir output/general_observable/lo_certification_direct_uv/gaussian_amplitude_01
```

Outputs: `all_16_gaussian_amplitude_status.json`,
`all_16_gaussian_amplitude_summary.csv`,
`models/<MODEL>.json` (source-pair and full `LLSS` tensors; real/imag encoding).
To avoid large tensors, add `--summary-only` (norms still retained).
Use a **new** `--output-dir` if any output already exists; no overwrite.
