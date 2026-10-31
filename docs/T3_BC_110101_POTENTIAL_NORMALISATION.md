# B(-1)/C(-1) `110101`: source-reconstructed potential normalization

## Purpose and exact scope

This audit advances the previous **native SU(2) tensor-geometry** and **mixed-scalar
mass-basis** audits to an action-level *reconstruction* for the two independent
`H Hdag S1dag S2dag` quartic invariants. It uses live-exported RGBeta quartic
**metadata**, the saved native Matchete quartic components, the Matchete source
`PlusHc` convention, and the standard RGBeta generic Lagrangian convention.

**It does not directly run the installed RGBeta `Lam` on test fields and does not
run native Matchete `Match[]`.** The result is a testable source-convention map,
not a full physical matching certificate. Neither production Lagrangians nor
neutrino fit results are modified.

## Source-backed derivation

Use the ordered component field list `(H, Hdag, S1dag, S2dag)`. The previously
validated native tensors satisfy

\[
 (R_1,R_2)=(M_1,M_2)A, \qquad
 A=\begin{pmatrix}1/2&1\\-\sqrt3/2&0\end{pmatrix}.
\]

The native RGBeta metadata for each of these two couplings contains:

- `native_coupling_function_inputform = 6*Matrix[...][##1] &`;
- exactly 24 permutations in `native_unique_arrangements_inputform`;
- `self_conjugate = false`, so the complex conjugate is included.

RGBeta's `Package/FieldsAndCouplings.m`, `AddQuartic` and `Lam`, construct the
symmetric real-scalar quartic coefficient by using `24` times
`(1/Sqrt[2])^4=1/4`. The internal normalization is therefore `6`.
RGBeta's generic scalar interaction is conventionally

\[
 \mathcal L^{\rm RGB}_4=-\frac1{4!}\Lambda_{abcd}
 \Phi_a\Phi_b\Phi_c\Phi_d\quad\Longleftrightarrow\quad
 V^{\rm RGB}_4=+\frac1{4!}\Lambda_{abcd}\Phi_a\Phi_b\Phi_c\Phi_d.
\]

This sign and the real-scalar normalization appear in Thomsen,
*Introducing RGBeta*, Eur. Phys. J. C **81** (2021) 408, Eq. (2.1).
When the complex fields are \(\phi=(R+iI)/\sqrt2\), a fixed ordered invariant
appears with coefficient

\[
 \frac{6}{4!}(\sqrt2)^4=1.
\]

In the real-component reconstruction the factor `(sqrt2)^4` converts the
complex field product to the unnormalized real combinations used by the real
quartic tensor. The 24 arrangements symmetrize the real coefficient; their
average is unchanged for a commuting quartic monomial. Neither `6` nor `24`
multiplies the **physical** complex quartic coupling.

The source `Lagrangian/interactions/GeneralScalarOperators.wl` defines the
Matchete extra interaction with the expression `PlusHc[contribution]` and
`contribution = coupling[] product * CG[...]`. Hence

\[
 V^{\rm RGB}=g_1R_1\phi^4+g_2R_2\phi^4+\mathrm{h.c.},
 \qquad
 \mathcal L^{\rm Mat}=h_1M_1\phi^4+h_2M_2\phi^4+\mathrm{h.c.}
\]

in the declared common component convention. **Within these source
conventions**, matching \(V^{\rm RGB}=-\mathcal L^{\rm Mat}\) gives

\[
 \boxed{\begin{pmatrix}h_1\\h_2\end{pmatrix}
 =-A\begin{pmatrix}g_1\\g_2\end{pmatrix}},
 \qquad h_1=-g_1/2-g_2,\quad h_2=\sqrt3\,g_1/2.
\]

The map is a *full 2x2 mixing of couplings*, not one minus sign per symbol.

## Verification performed by the new script

- Reads independently saved native Matchete `110101` tensors, live exported
  RGBeta invariant components, and the actual B/C saved scalar mass bases.
- Re-derives the exact two-dimensional tensor map using unrounded components.
- Rejects a changed RGBeta real-tensor coupling factor, missing permutation,
  absent Hermitian conjugate, or changed Matchete `PlusHc` source.
- Numerically compares the RGBeta real-tensor reconstructed potential with
  the ordinary complex-invariant potential for independent complex couplings.
- Compares \(V^{\rm RGB}+\mathcal L^{\rm Mat}\), **including h.c.**, both before
  and after rotation to the mixed-scalar mass basis.
- Outputs per-model residuals and input hashes, but leaves physical matching
  and runtime action-level certification flags `false`.

## Install and execute (PowerShell)

From the repository root:

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath 'C:\path\to\T3_BC_110101_potential_normalisation.zip' -DestinationPath . -Force
python -X utf8 -m pytest -q tests/test_t3_bc_110101_potential_map.py
python -u -X utf8 -m scripts.AuditT3BC110101PotentialMap --all --run-id bc_110101_potential_v1
```

Dependencies are the already installed `numpy`, `sympy`, the prior audit
scripts `AuditT3BCNativeQuarticMassBasis.py`, `AuditT3BC110101MassBasisMap.py`,
and existing Stage22 and B/C mass-basis JSON files. A new run ID is required
if an output folder already exists.

Outputs:

`output/audits/t3_bc_110101_potential_map/bc_110101_potential_v1/`

Files: `B_m1.json`, `C_m1.json`, `summary.json`.

## Limits and next gate

The reconstructed sign and normalizations use RGBeta's **published generic
Lagrangian** plus the previously saved native RGBeta record and Matchete
source. They do **not** independently evaluate the installed RGBeta
`Lam[a,b,c,d]` and its `sDelS`/real-scalar projector, and they do **not**
validate a newly built mass-eigenstate Matchete Lagrangian.

The next independent step should evaluate one complex coupling at a time in
*installed* RGBeta's `Lam` (real components) and compare against the native
Matchete `PlusHc` interaction at explicit real component values. Only then
should these two couplings' full action map be certified and fed into
numerical matching. The other B/C quartic blocks, Yukawa vertices, native
mass-eigenstate `CheckLagrangian`/`Match[]`, and sequential Wilson thresholds
remain further gates.

## References and provenance

- Anders Eller Thomsen, *Introducing RGBeta*, Eur. Phys. J. C **81**, 408
  (2021), Eq. 2.1: https://link.springer.com/article/10.1140/epjc/s10052-021-09142-4
- RGBeta `Package/FieldsAndCouplings.m`, inspected read-only:
  https://github.com/aethomsen/RGBeta/blob/master/Package/FieldsAndCouplings.m
- Repository files (read-only inspection):
  `Lagrangian/interactions/GeneralScalarOperators.wl`,
  `output/general_t3_stage22/rgbeta_native_invariants.json`,
  `output/general_t3_stage22/stage22_110101_wolfram_native_covariance.json`,
  `output/general_t3_stage22/matchete_tensors.json`,
  `output/general_t3_mass_basis/{B_m1,C_m1}_final_mass_basis.json`.

The local test fixture used when this package was authored contains the same
saved B/C tensor components, with a minimal Matchete source-string fixture and
reproduced native 24-permutation metadata. This is **not** a successful run
against the user's installed Wolfram/RGBeta instance.
