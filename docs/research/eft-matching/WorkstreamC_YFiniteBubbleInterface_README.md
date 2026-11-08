# Workstream C — conditional finite Y scalar-bubble interface (v12)

## Scope

This file **does not certify the physical finite Y threshold**. It assembles only the independently checked scalar-quartic two-propagator subset at strict dimension five, with `muY = MY` and hard-region `MX/MY -> 0`. It is a **callable, source-anchored candidate contribution** which can later be combined with full EFT matching, inherited corrections, and running.

Earlier audits established (for B(-1), C(-1)):

- v5: source-derived `(KYY + KXY)` UV kernel, independent scalar-bubble MSbar counterterm, and `KXY` native formal finite coefficient.
- v6: exact `KXX` component equality in EFT1 and EFT2; light-XX UV and finite terms cancel in the EFT difference.
- v8: exact PL/PR Hermitian-conjugate scalar/flavour polynomials, with the ordered `y1*y2` term appearing *once*.
- v11.2: the live installed-Matchete test returned `ExactPlusEquality` for native PL `NCM[Transp[l], GammaCC PL, l]` versus the Warsaw `Q_llHH` operator, positive sign.

The new interface **checks** these successful report files, hashes the validated constant ledger, and rereads its 16 PL and 16 PR exact symbolic monomials. On a missing or changed required gate it fails without producing a partial physical coefficient.

## Exact operator convention

Define the **native ordered** scalar-bubble `PL` operator coefficient at `muY=MY` by

\[
A_{ij} = \frac{1}{16\pi^2 M_F} \left[
 a_{11}\,y_{1i}^*y_{1j}^*
 +a_{12}\,y_{1i}^*y_{2j}^*
 +a_{22}\,y_{2i}^*y_{2j}^*
\right],
\]

where each `a_ab` is the *exact formal finite scalar-bubble coefficient* assembled from the B/C native `PL` ledger, **excluding** the `hbar/MF` common factor. The `a_ab` include the real `1/sqrt(3)` factor for C(-1) and the actual native mixing-spurion and quartic-coupling conjugation.

For the installed Warsaw model the flavour sum is `A_ij Q_ij + h.c.` (no explicit `1/2`). Given the verified positive-sign native PL operator identity, the conventional `L5 = (1/2) C5_ij Q_ij + h.c.` normalization implies

\[
C_{5,\mathrm{bubble}}^{Y,\mathrm{conditional}} = A + A^T.
\]

**Transpose, not Hermitian conjugate.** The native PR polynomial is the Hermitian-conjugate operator and is **not added again**. With `L = [y1*  y2*]` as a `3x2` matrix, the exact coupling-polynomial interface is

\[
C_{5,\mathrm{bubble}}^{Y,\mathrm{conditional}} =
\frac{1}{16\pi^2 M_F}\,
L \begin{pmatrix}2a_{11}&a_{12}\\a_{12}&2a_{22}\end{pmatrix}L^T.
\]

This is a **partial contribution** only. The outputs explicitly set `physical_C5=null` and `full_physical_matching_certified=false`.

## Project files and commands

Put `scripts/BuildT3BCYFiniteBubbleInterface.experimental.py` at exactly that repository-relative path. It imports the **existing** v8 helper and checks existing reports; it never modifies Matchete, RGBeta, or any production files.

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline

python scripts\BuildT3BCYFiniteBubbleInterface.experimental.py --self-test

python scripts\BuildT3BCYFiniteBubbleInterface.experimental.py `
  --output output\audits\t3_bc_y_conditional_finite_interface_v1\summary.json
```

Required existing reports:

- `output/audits/t3_bc_y_scalar_bubble_msbar_v1/summary.json`
- `output/audits/t3_bc_y_eft2_light_xx_overlap_v1/summary.json`
- `output/audits/t3_bc_y_flavour_symmetry_v1/summary.json`
- `output/audits/t3_bc_y_warsaw_direct_identity_v11/summary_live.json`
- `output/audits/t3_bc_y_renormalization_ledger_v1/{B_m1,C_m1}_y_llhh_constant_neutral_flavour.json`

If an earlier report is absent or stale, rerun **only the existing relevant audit**. Do not fabricate a report or bypass the gate.

### Optional conditional numerical evaluation

Supply your actual Yukawa vectors, native source quartic coupling values and the fermion mass. The mixing U is read from the already verified benchmark parameter JSON rather than guessed. Complex numbers use explicit `[real, imag]` arrays:

```json
{
  "model": "B_-1",
  "MF_GeV": 1800.0,
  "y1": [[0.01,0.0],[0.02,0.0],[0.03,0.0]],
  "y2": [[0.02,0.0],[0.03,0.0],[0.01,0.0]],
  "quartics": {
    "lambdaT3": [0.1,0.0],
    "lambdaT3ExtraQ020020I1": [0.05,0.0],
    "lambdaT3ExtraQ020200I1": [-0.03,0.0]
  }
}
```

**The example is synthetic, not a B/C numerical benchmark.** Use your actual MF, yukawas and coupling values before interpreting a numerical candidate. Run:

```powershell
python scripts\BuildT3BCYFiniteBubbleInterface.experimental.py `
  --numeric-input path\to\your\real_B_minus1_inputs.json `
  --output output\audits\t3_bc_y_conditional_finite_interface_v1\numeric_B.json
```

This returns a complex symmetric `3x3` **conditional scalar-bubble-only** coefficient in `GeV^-1`. It does not return a physical neutrino mass, even if inputs are numerical. `--require-physical` fails by design.

## Unresolved physics gate

To certify complete finite Y matching one still needs the **full** EFT1/EFT2 renormalized-operator and amplitude subtraction (including non-bubble contributions, any induced matching of field/operators, and inherited fermion-hard terms), and a consistent prescription for Y-inherited `LLXX`/`LLHX` contributions, RG evolution, X threshold, and soft-region higher-order operators. The independent canonical scalar EFT versus original fermion-first `EOMEFT` equivalence must be established for every operator needed by this full prescription, not only the checked bubble subset.

Mass corrections in the full XY bubble are large for `MX^2/MY^2 ~ 0.76`; this diagnostic deliberately does not substitute the full-mass `0.13...` reference factor for the strict dimension-five hard coefficient `1`.
