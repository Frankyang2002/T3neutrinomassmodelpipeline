# All-16 native mass Ward and conditional spinor-spurion dictionary

## Scope

Additive, read-only source-convention audit for all 16 ordinary T3 models. It
consumes existing local `RGE/running/rgbeta/T3RGBetaModel.wl`,
`Lagrangian/interactions/T3Topology.wl`, and `Lagrangian/model/T3Fields.wl`.
It does **not** run Wolfram, change any Lagrangian, edit the provisional UV
benchmark, or recalculate the already completed 16/16 UV masses.

Run from the repository root:

```powershell
python -m unittest discover -s tests -p "test_lo_direct_uv_*.py" -v
python scripts/RunLODirectUVSpinorMassDictionary.py `
  --project-root . `
  --output-dir output/general_observable/lo_certification_direct_uv/spinor_mass_01
```

Output: `output/general_observable/lo_certification_direct_uv/spinor_mass_01/all_16_spinor_mass_dictionary.json`.
Existing report files are never overwritten. The output directory must be
inside the existing dedicated direct-UV certification tree. The test suite has
43 tests (32 existing plus 11 new) in the isolated replica of the previously
inspected direct-loop Python modules.

## New source finding: SU(2) fermion mass contraction

The **current** RGBeta `T3RGBetaAddMasses` source registers:

- Dirac `dF=1`: scalar invariant `1`;
- Dirac `dF=2`: `eps[T3RGBetaRep[dF],#1,#2]` between `F` and `Fc`;
- Dirac `dF=3`: `del[T3RGBetaRep[dF],#1,#2]`;
- self-conjugate `dF=1,3`: `F,F` with `1` or `del`.

Independent **exact SU(2) Ward** identities verify
`T_a^T C + C T_a = 0` for singlet, fundamental doublet, and Cartesian
adjoint triplet. The two-left-fundamental delta `C=I_2` fails generators
1 and 3 and is an explicit negative control. The Cartesian triplet basis is
a *declared representation convention*, not a measured native RGBeta
numerical-basis transform.

**Version warning:** the earlier
`Lagrangian/RunT3A0DiracMassWard.py` source guard expected the old doublet
`del` mass invariant and may reject the current corrected epsilon form. This
batch deliberately does not edit that script; any optional prior A(0) native
mass report with an old SHA or delta is labelled stale/incompatible, never
reused as a current physical certificate.

## Source vertex orientation and conditional phase transport

The native source definitions agree on the operator *families* but do **not**
by themselves fix their common two-component spinor phase:

| Source | y1 | y2 |
|---|---|---|
| Matchete T3Topology | `NCM[Bar[l],CConj[F]] S1` | `NCM[Bar[l],F] Bar[S2]` |
| RGBeta T3 model | `{Bar[S1],l,F}` | `{S2,l,Fc}` (Dirac); `F` (Majorana) |

The specific model-level all-left-handed **candidate action** is defined as

\[
{\cal L}_\mathrm{candidate}\supset
 y_1\nu F^c S_1+y_2\nu F S_2^*+M_F F F^c+
 \lambda_{T3}HH S_1 S_2^*+\mathrm{h.c.}
\]

with the conventional factor `1/2` multiplying the self-conjugate Majorana
mass bilinear. Under `F -> exp(i a) F`, `Fc -> exp(i b) Fc`,
`S1 -> exp(i c) S1`, `S2 -> exp(i d) S2`, the spurion rephasings are:

\[
\begin{split}
q(y_1)&=(0,1,1,0),\\
q(y_2)&=(1,0,0,-1),\\
q(M_F)&=(1,1,0,0),\\
q(\lambda_{T3})&=(0,0,1,-1).
\end{split}
\]

For a Majorana multiplet set `b=a`; the fermion mass charge becomes `2a`.
Hence the conditional invariant

\[
\boxed{\ I_\mathrm{candidate}=y_1y_2M_F^*\lambda_{T3}^*\ }
\]

has exactly zero charge for the four field rephasings in both cases. The
incorrect nonconjugated `lambdaT3` assignment is rejected algebraically. For
Dirac cases, assigning both Yukawas to the same left-Weyl component is also
rejected. The six additional source-allowed Yukawa couplings `y1Conj` and
`y2Conj` are **independent operators**, not Hermitian-conjugate duplicates;
we explicitly retain their distinct rephasing vectors without assigning
missing numerical values.

This is a conditional source-operator **rephasing statement**, **not** a proof
that the physical one-loop matrix is proportional to this particular complex
product in the installed 4-spinor and RGBeta conventions. A native amplitude
projection must determine complex conjugations, phases, possible `1/2`
normalisations, and the detailed SU(2) tensor contraction. No existing UV
mass is changed by this audit.

## Certification status and next native gate

- Algebraic all-16 source mass Ward and candidate phase-spurion identities:
  **checked** (subject to local source-contract guard).
- Installed Matchete `GammaCC` 4x4 sign relative to the numerical Weyl basis:
  **not measured**.
- Installed RGBeta fundamental epsilon component orientation and full
  `F/Fc` correspondence: **not measured**.
- Signed physical dictionary across `y1,y2,MF,lambdaT3`: **not certified**.
- Scalar portal RGBeta/Matchete physical sign: **not certified**.
- All 16 original full-UV numerical results: **unchanged; still provisional**.

Do not promote a source-operator Ward identity or rephasing invariant to a
physical source-verified coupling dictionary. The next step is a **native
amplitude-level projection** for one representative Dirac doublet, one
neutral Majorana state, then the remaining representation classes, with a
common installed Matchete/RGBeta Weyl component basis and the explicit
`MF` mass-insertion sign and symmetry factor.

The EFT order remains **UV -> integrate out ALL heavy fermions at muF ->
EFT1 with both scalars and all SM fermions -> integrate out BOTH scalars at
muS -> SMEFT**.
