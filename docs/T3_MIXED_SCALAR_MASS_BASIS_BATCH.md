# T3 mixed scalar mass-basis batch — experimental B(-1)/C(-1)

## Source and scope

Inspected repository `main` definitions:

- `Lagrangian/model/T3ModelCatalog.wl`: for alpha = -1, the two doublets have
  `Y(S1)=-1/2` and `Y(S2)=+1/2`.
- `RGE/running/rgbeta/T3RGBetaModel.wl`: 10 baseline scalar quartic slots
  (including SM lambdaH), 6 pre-existing structural quartics for each model,
  and 2 Yukawa scalar legs (`Bar @ S1` in `y1`, `S2` in `y2`).
- `RGE/general/generated/{B,C}_-1_extra_quartics.wl`: a seventh, independent
  RGBeta candidate `lambdaT3GeneralQ110101I1`, alongside the pre-existing
  `lambdaHHdagS1barS2barCross` contraction of the same field content.
- `Numerical/experimental/GeneralT3ScalarQuadraticMasses.py`: existing
  eight-real-mode Hessian in a specified SU(2) spherical convention.
- `scripts/RunGeneralT3QuadraticEFT1.py`: recorded Stage EFT1 data at
  `output/general_t3_quadratic_eft1/{B,C}m1_eft1_quadratic_running.json`.

No production files are overwritten. No GitHub changes. The two new Python
files operate on saved numerical data, so they do not need Wolfram for execution.

## Exact mixing convention

Set `epsilon=[[0,1],[-1,0]]`, `tilde(S2)=epsilon S2*`. For each weak-isospin
component define `Phi=(S1,tilde(S2))` and 

\[
  V_2=\Phi^\dagger
  \begin{pmatrix}m_1^2&\mu^2\\(\mu^2)^*&m_2^2\end{pmatrix}\Phi,
  \qquad \Phi=U\begin{pmatrix}X\\Y\end{pmatrix}.
\]

All `X,Y` are SU(2) **doublets of hypercharge -1/2**. By construction
`U^dagger A U=diag(m_-^2,m_+^2)` with

\[
 m_\pm^2=\frac{m_1^2+m_2^2}{2}\pm
 \sqrt{\left(\frac{m_1^2-m_2^2}{2}\right)^2+|\mu^2|^2}.
\]

The inverse species substitutions are

\[
 S_1=U_{11}X+U_{12}Y, \quad
 S_2=-\epsilon(U_{21}^*X^*+U_{22}^*Y^*).
\]

Here `U` indices are 1-based in the equations and 0-based in the JSON.
`mu2` is the source coefficient of `S1dag epsilon S2dag+h.c.`; there is NO
extra factor 1/2 for this mixed term. If both squared masses are identical and
`mu2=0`, the eigenbasis is nonunique and reported as degenerate.

## What is transported, and what is not

For **each of the 17 source quartic symbols** in B(-1) and C(-1) we record
one separate multilinear tensor transformation. The basis holds 10 ordinary
couplings, 6 pre-existing structural analogues and the 1 additional RGBeta
candidate. Both original Yukawa scalar legs are transformed too.

For a source interaction `lambda T[i1,...,in] f1[i1]...fn[in]`, every rotated
term stores an **ordered list of new fields and its full 2x2 SU(2) index map**
`M1,...,Mn`. The interaction is

\[
 \lambda\,T_{i_1\ldots i_n}
 M^{(1)}_{i_1j_1}\cdots M^{(n)}_{i_nj_n}
 f'_1{}_{j_1}\cdots f'_n{}_{j_n}.
\]

The source `T` is **preserved rather than replaced by an assumed common
CG**. This is why two different contractions with identical field lists
stay independent. The transformation is valid for **any** supplied original
SU(2) component tensor; the supplied 2x2 matrices can be used when a physical
Matchete/RGBeta common-basis tensor has been certified. The report separately
stores the corresponding source numerical RGBeta couplings. It does not
pretend that these are already physical mass-basis independent coupling
constants.

Only the scalar factor in the Yukawa vertices is transported. The fermion
and lepton SU(2) tensor factors and flavor matrices are unchanged symbols;
they cannot be numerically evaluated from scalar-only EFT1 endpoints because
F has already been integrated out there.

The **whole quadratic mass tensor**, including diagonal terms and the
holomorphic mixing term plus its hermitian conjugate, is independently
transformed. It is checked against the diagonal `U^dagger A U` and against
the existing eight-real-mode Hessian spectrum (each complex-doublet mass
squared repeated four times as a real eigenvalue).

## Install (all files in one ZIP)

Expand the ZIP into `C:\thesistemp\T3neutrinomassmodelpipeline`. The files go to:

- `Numerical/experimental/GeneralT3MixedScalarMassBasis.py` (new)
- `scripts/RunGeneralT3MixedScalarMassBasis.py` (new)
- `tests/python/test_general_t3_mixed_scalar_mass_basis.py` (new)
- `docs/T3_MIXED_SCALAR_MASS_BASIS_BATCH.md` (new)

### Commands (PowerShell, from the repository root)

```powershell
python -m pytest -q tests\python\test_general_t3_mixed_scalar_mass_basis.py
python -m scripts.RunGeneralT3MixedScalarMassBasis --model B_-1
python -m scripts.RunGeneralT3MixedScalarMassBasis --model C_-1
python -m scripts.RunGeneralT3MixedScalarMassBasis --all

# Optional: check the boundary immediately below the fermion threshold
python -m scripts.RunGeneralT3MixedScalarMassBasis --all --which initial
```

If the existing EFT1 inputs are absent, first run the already-installed
`python -m scripts.RunGeneralT3QuadraticEFT1 --all --mu2-fraction 0.01`.

Outputs (one for each model and endpoint):

- `output/general_t3_mass_basis/B_m1_final_mass_basis.json`
- `output/general_t3_mass_basis/C_m1_final_mass_basis.json`

Each report includes the mass eigenvalues, unitary rotation, source couplings,
**all** vertex-index transport terms, the exact quadratic mass tensor check,
SU(2) covariance residual, tensor-identity residuals, provenance hashes and
explicit uncertified matching flags. The printed eigenmasses are *running*
scalar mass parameters, not physical pole masses.

## Verified versus blocked

- **Verified offline:** analytic eigenvalue identity; unitary and inverse map;
  positive/negative mixing phases; SU(2) gauge covariance; mass potential
  equality for arbitrary field configurations; consistency with the existing
  B/C quadratic Hessian; formal quartic tensor transport for 17/17 source
  vertices for each model; Yukawa scalar-leg transport; source field-inventory
  consistency; fail-closed corrupt input cases.
- **Source verified:** field charge assignments, 17 RGBeta scalar couplings,
  2 original Yukawa scalar legs, 7 additional quartic structures with one
  two-dimensional contraction sector, Stage EFT1 report contract.
- **Not certified:** Matchete ↔ RGBeta signed component and coupling prefactors,
  physical field component CGs, native Matchete mass-eigenfield definitions,
  heavy-mass-basis matching, transformed Wilson coefficient matching, and
  low-energy neutrino masses. The code does **not** claim Matchete
  `HeavyMassBasis` succeeds on a redefined Lagrangian.

## Next physical task

Build certified Matchete component tensors and a physical signed coupling map;
install the two mass-eigen heavy doublets into the Matchete Lagrangian with all
transformed vertices and verify `CheckLagrangian`/`HeavyMassBasis` there. Only
then calculate scalar threshold Wilson coefficients with the two nondegenerate
mass eigenstates. **Do not** substitute the present tensor-routing report for
an actual Matchete mass-basis model or loop-matching computation.
