# B/C(-1) quartic coupling normalization: source evidence and next check

## Current verified result

The local B/C tensor-geometry audit returned

| Model | 110101 rank | Max rotated-tensor residual |
| --- | ---: | ---: |
| B(-1) | 2 | `2.220e-16` |
| C(-1) | 2 | `1.110e-16` |

This resolves the invariant-tensor geometry in the saved mixed-scalar basis,
not the physical coupling map. The source-based exact tensor identity in the
ordered `(H,Hdag,S1dag,S2dag)` field convention is

\[
(R_1,R_2)=(M_1,M_2)A, \qquad
A=\begin{pmatrix}1/2&1\\-\sqrt{3}/2&0\end{pmatrix}.
\]

`M_i` are saved Matchete invariant tensors and `R_i` are saved *live*
RGBeta invariant tensors. Neither the RGBeta potential normalization nor the
relative Lagrangian sign is contained in this identity.

## Direct source evidence inspected

Upstream source: [RGBeta `Package/FieldsAndCouplings.m`](https://github.com/aethomsen/RGBeta/blob/master/Package/FieldsAndCouplings.m), branch `master`.

- `AddQuartic` assigns an internal prefactor `24` and multiplies it by
  `1/Sqrt[2]` for each charged/complex scalar occurrence. Thus the native
  internal tensor coupling is multiplied by `6` for four complex fields;
  this agrees with the exported `6*Matrix[...][##1] &` strings.
- The source creates the invariant tensor through the user-supplied
  `GroupInvariant` and uses `CouplingPermutations` and `AvgPermutations`.
- `Lam` includes the conjugate interaction when `SelfConjugate -> False`.
- The [RGBeta paper](https://link.springer.com/article/10.1140/epjc/s10052-021-09142-4) explicitly states that it does not assume quartic normalization factors automatically; model-defined invariant prefactors matter.

The factor `6` is an **internal real-scalar tensor conversion factor**, not
proof that the physical coefficient in Matchete's Lagrangian is six times
the RGBeta input coupling. Such an identification would be unjustified.

## The remaining normalization/sign proof

Write the RGBeta potential and Matchete interaction in one **common field
component convention**, for independent complex couplings `g_1,g_2` and
`h_1,h_2`:

\[
V_{\rm RGB} = \sum_{i=1}^2 c_{R,i}(g_i)\, R_i\,\Phi^4 + {\rm h.c.},
\qquad
\mathcal{L}_{\rm Mat}=\sum_{a=1}^2 c_{M,a}(h_a)\,M_a\,\Phi^4+{\rm h.c.}.
\]

The independently verified geometric transformation gives the tensor relation
but not yet the physical coefficient relation. To establish that relation:

1. **RGBeta action-level test:** derive the scalar potential from the
   installed RGBeta version's internal `Lam` and its real-to-complex scalar
   convention, including `1/4!`, `UniqueArrangements`, `SelfConjugate`, and
   both coupling and conjugate terms. Test two independent complex couplings
   one at a time. Avoid interpreting raw `6*Matrix` as the potential factor.
2. **Matchete action-level test:** evaluate the actual native Lagrangian
   vertex in exactly the same ordered SU(2) component convention, including
   `PlusHc` and the sign of the interaction in the action.
3. **Field test:** evaluate both scalar potentials at generic complex Higgs,
   S1, S2 field components and require numerical equality for several
   independent assignments; also compare exact coefficient tensors.
4. **Mass-basis test:** transport the resulting *physically normalized*
   coefficients through the saved B/C scalar mixing matrix and reproduce
   the same action before and after rotation.
5. **Fail closed:** refuse physical certification if any scalar-field,
   conjugation, normalization, or sign convention remains unidentified.

Only after this test can the exact `110101` geometry matrix `A` be used as a
physical parameter map in numerical matching. Full Yukawa-vertex mapping,
remaining quartic blocks, `HeavyMassBasis` validation, one-loop matching,
and sequential EFT thresholds remain further tasks.

## Source provenance / status

- Local `bc_110101_map_v1` terminal result supplied by the user: PASS for B
  and C; saved `summary.json` not separately read here.
- Upstream RGBeta `FieldsAndCouplings.m` inspected read-only: internal
  normalization and H.c. control flow confirmed.
- Absolute physical prefactor and Matchete relative sign: **NOT VERIFIED**.
- Production Lagrangians and physics code: **UNCHANGED**.

The next calculation should be an action-level normalization comparison, **not**
another rank or random-tensor-geometry audit.
