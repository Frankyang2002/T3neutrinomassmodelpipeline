# B2 batch 3 — literal RGBeta SU(2) Ward and charge witness

## Scope and provenance

- **Checked source:** `RGE/running/rgbeta/T3RGBetaModel.wl`, unchanged Git blob `7af7a4ddf99a7e38210543a21b3b1cda203eaf02` on GitHub `main` commit `e2a5def5326cc452f668ef593bd81ae13fdff30b` (2026-11-03 Sydney). All other five B2-tracked source blob IDs are inherited from the previously reviewed B2 source inventory. The Python runner checks **all six** tracked paths against local `HEAD`, the index, canonical worktree contents and a clean Git diff. It does **not** require the older historical HEAD when unchanged source blobs have been carried into a newer commit.
- **Native executed evidence already in the repository:** `output/general_observable/physical_normalisation/b2_rgbeta_native_01/formal_invariants.json`, Git blob `34bf2f1d7b325cab0d5162b235ea2df0e8034c78`. The snapshot shows five native formal classes and explicitly reports no numerical RGBeta component axes or physical action normalisations.
- **This batch:** a standalone, exact SymPy SU(2) Ward witness and a fresh-kernel extended registry probe. The Python symbolic tests were executed here; Wolfram/installed RGBeta **was not run here**. No production file is changed, and no RGBeta operator is automatically replaced.

## 1. A concrete diagnostic, not an asserted package bug

The source registers the topology mixing quartic with ordered fields

\[
\lambda_{T3}: \quad (H_i, H_j, S_{1,s}, S_{2,t}^*)
\]

and defines, for `(dS1,dS2)=(2,2)`,

```wl
T3RGBetaMixInvariant[2, 2] :=
  (delS2[SU2L, a, #1, #2] * delS2[SU2L, a, #3, #4] &);
```

The previously executed native RGBeta *formal* export simplifies this expression to

\[
C_{ij s\bar t}^{(\delta\delta)}
=\frac12(\delta_{it}\delta_{js}+\delta_{jt}\delta_{is}).
\tag{1}
\]

This is a statement about the **formal expression** exported by the package. It does not say that `del` means a Kronecker delta on four *numerically specified* field component axes. An unexamined conjugate-index or pseudoreality map inside RGBeta could alter its interpretation.

### Literal-index test

Choose SU(2) fundamentals of weights `(+1/2,-1/2)` with generators

\[
T_1=\tfrac12\sigma_1,\quad T_2=\tfrac12\sigma_2,
\quad T_3=\tfrac12\sigma_3,\quad
\epsilon=\begin{pmatrix}0&1\\-1&0\end{pmatrix}.
\]

A barred scalar is the conjugate representation, with generator
\(T_a^{\bar{2}}=-T_a^*\). For a candidate coefficient tensor \(C\) and ordered reps \(R_k\), the invariant condition is

\[
(\delta_a C)_{i_1\ldots i_n}
=\sum_{k=1}^n\sum_{j} C_{i_1\ldots j\ldots i_n}
(T_a^{R_k})_{j i_k}=0.
\tag{2}
\]

The exact component sum \(\sum_i |\delta_a C_i|^2\), for `reps=(2,2,2,bar2)`, gives:

| Tensor on `(H,H,S1,Bar[S2])` | \(T_1\) | \(T_2\) | \(T_3\) |
|---|---:|---:|---:|
| Literal exported \(\delta\delta\), Eq. (1) | 3 | 0 | 3 |
| Chosen SU(2) singlet \(\delta\epsilon\), Eq. (3) | 0 | 0 | 0 |

The convention-chosen valid tensor is

\[
C_{ij s\bar t}^{(\delta\epsilon)}=
\frac12\left(\delta_{it}\epsilon_{js}+\delta_{jt}\epsilon_{is}\right).
\tag{3}
\]

Both tensors are symmetric under \(i\leftrightarrow j\) (the two Higgs fields). The delta–epsilon tensor is a nonzero SU(2) singlet; its **phase and amplitude relative to RGBeta's actual operator remain unmeasured**.

### Independent charge test

Setting both Higgs indices to the neutral component \(i=j=1\) (zero-based), Eq. (1) supports only \((s,t)=(1,1)\), i.e. equal doublet weights. Yet the hypercharges satisfy \(Y_{S_2}=Y_{S_1}+1\), hence \(Q_{S_2}-Q_{S_1}=1\): this does **not** preserve the electromagnetic charge in scalar mixing under literal indices. Eq. (3) supports \((s,t)=(0,1)\) and gives equal charges. The Python runner repeats this exactly for B/C at all three alphas `{-3,-1,+1}`: **six conditional model witnesses**.

This finding concerns the source's **literal symbolic transcription**. It is **not yet evidence that RGBeta's abstract `delS2` or `Bar` implementation produces the literal tensor (1) when interpreted on physical component fields**. Do not silently modify `T3RGBetaModel.wl` on this basis alone.

## 2. Yukawa counterpart, B class

The RGBeta source defines `T3RGBetaYukawaInvariant[2,1]` using `eps` and registers `y1` on `{Bar[S1], l, F}`, with the heavy `F` an SU(2) singlet in class B. Under literal `Bar[S1]∈bar2`, `l∈2`, this is the `bar2⊗2` invariant, and the expected contraction is \(\delta_{ab}\), whereas the registered expression is \(\epsilon_{ab}\).

Using Eq. (2):

| B-class \(y_1\) tensor on `(Bar[S1],l)` | \(T_1\) | \(T_2\) | \(T_3\) |
|---|---:|---:|---:|
| Literal `eps` | 2 | 0 | 2 |
| Chosen `delta` | 0 | 0 | 0 |

The same `eps` in `y2` on `{S2,l,F}` contracts two **unbarred** fundamentals, and passes (0,0,0). Converting `Bar[S1]` into a pseudoreal doublet basis using an epsilon changes the coordinate interpretation. A **direct native index/basis witness is required** before identifying a source defect or signed coupling conversion. For class C, `delS2` similarly has abstract triplet and barred-index semantics which this diagnostic does not assign numerical values.

## 3. What gets computed and what does not

`b2_rgbeta_ward_witness.py` constructs tensors explicitly from Pauli matrices, epsilon and Kronecker delta. The results are **algebraic identities in a declared common basis**, independent of the native RGBeta numerical component conventions. They can diagnose tension with a formal source expression; they do not calibrate the package's generator basis.

The optional guarded Matchete comparator reads *actual existing* `component_support_audit.json` and `matchete_generator_analysis.json`. For the six B/C models it solves

\[
C^{(M)}_{ij s\bar t}=r_M C^{(\delta\epsilon)}_{ij s\bar t}
\quad\text{and reports}\quad
\|C^{(M)}-r_M C^{(\delta\epsilon)}\|^2.
\]

`r_M` is a ratio to a **chosen invariant**, NOT a physical \(\lambda_{T3}^{(RGBeta)}/\lambda_{T3}^{(Matchete)}\) conversion. It cannot be promoted to a package-to-package ratio without the native RGBeta component operator and action normalisation. Six representative checks do not cover A/D/E triplet component bases.

The new `RunPhysicalB2RGBetaRegistryAudit.wl` gets and serialises the live `y1`, `y2`, `MF`, and `lambdaT3` registry metadata for one representative of each of A–E. This can distinguish a missing/changed native registry from a source-comment inference. `Coupling` and `CouplingBar` functions, if present, are reported *as formal metadata*: their internal numeric prefactors do **not** establish action factors.

## 4. Certification and next decisive native evidence

- **Exact synthetic Ward tests:** passed in the authored workspace. The runner will independently reproduce them in the user's local checkout.
- **Native formal registry evidence:** previously captured; guarded against changed fields/status/expressions when supplied locally.
- **Native numeric RGBeta `eps`/`delS2` coordinates, index raising map and scalar action factors:** **not established**.
- **B1 spinor/`GammaCC` physical phase and mass dictionary:** independent B1 work, still required for full signed Yukawas and `MF`.
- **Physical cross-package conversion** for \(y_1,y_2,M_F,\lambda_{T3}\): **NULL**.
- **Portal sectors:** the prior 10/10 Hermitian singlet **operator embeddings** remain the existing evidence. Four sectors have real null directions; no full RGBeta action multiplier measured in this batch.
- **Physical source certification:** `false`, **0/16**; safe-to-integrate physical dictionary **no**. A warning here is NOT an automatic proof of a defect.

**Next decisive calculation:** inspect the installed RGBeta group-index interpretation and generated action-level vertex for `y1` and `lambdaT3` in class B (especially how `Bar[S1]` and `Bar[S2]` are represented), then repeat C and the triplet classes. If the generated action really uses the literal tensors of Eqs. (1)–(2), these are gauge-invariance failures that must be resolved in the RGBeta source before claiming any physical normalisation. If it maps barred indices pseudoreally, record that map explicitly and redo the Ward/conversion calculation with that map. In either case **do not apply a guessed source fix**.

**Outside scope:** EFT, RG running, direct one-loop mass evaluation or modification of frozen production inputs. The project's future EFT ordering remains UV → integrate out all fermions → EFT1 with both scalars → integrate out both scalars → SMEFT.
