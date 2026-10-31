# Workstream A: exact scalar identities and conditional UV normalization

**Project scope:** the 16 ordinary T3 models with the assumed exact BSM-odd
\(Z_2\) assignment. This stage reads existing reports and native CG data. It does
not alter the UV production builder, RGEs, matching, or numerical workstreams.

## Entry point

```powershell
python Lagrangian\RunT3UVExactEquivalence.py --self-test
python Lagrangian\RunT3UVExactEquivalence.py --all
```

Requirements: the Python environment used by the prior Workstream A stages,
including NumPy and SymPy. The prior successful `RunT3UVFinalEvidence.py --all`
and six-model self-bilinear eigenbasis reports must exist. The repository must
also contain:

- `output/general_t3_stage22/matchete_tensors.json` (Matchete native tensors)
- `output/general_t3_stage22/rgbeta_native_invariants.json` (native RGBeta source metadata)
- `output/general_t3_stage22/stage23_quartic_prefactor_candidates.json` (source-derived, conditional normalization maps)
- `Reports/output/workstream_a/final_uv_evidence/final_uv_evidence.json`
- `Reports/output/workstream_a/six_self_eigenbasis/{A,D,E}_{m2,p0}.json`

Outputs (new, never overwritten without `--overwrite`):

- `Reports/output/workstream_a/exact_uv_equivalence/exact_uv_equivalence.json`
- `Reports/output/workstream_a/exact_uv_equivalence/exact_uv_equivalence.md`

## Exact scalar mass and free-kinetic relations

For six real SU(2) singlet/triplet neutral-scalar sectors, use the *native
Matchete Cartesian* self-bilinear tensor \(\kappa\delta_{ij}\). In the source
Lagrangian the extra bilinear is **plus** \(\lambda_{\rm extra}\,\kappa (S^*)^2\)
plus its Hermitian conjugate. Consequently, in the *potential* the coefficient
is \(b=-\kappa\lambda_{\rm extra}\). The convention is:

\[
 b=\rho p^2,\quad p^*p=1,\quad \rho=|b|,\qquad
 S_i=\frac{p}{\sqrt2}(R_i+iI_i),\quad
 m_R^2=m_S^2+2\rho,\quad m_I^2=m_S^2-2\rho.
\]

The independent symbolic calculation confirms equality of the full quadratic
potential and ordinary derivative kinetic polynomial for dimensions 1 and 3.
It obtains \(\kappa\) from each actual prior Matchete report. A vanishing
bilinear is the continuous limiting case (arbitrary rotation phase). For the
six numerical benchmark points a positive mass spectrum is assumed; the
symbolic identity does not assume non-tachyonic parameters.

For B(−1) and C(−1), the source potential mass matrix is

\[
\Phi=\begin{pmatrix}S_1\\ \epsilon S_2^*\end{pmatrix},\quad
 \epsilon=\begin{pmatrix}0&1\\-1&0\end{pmatrix},\quad
 \mathcal M^2=U\operatorname{diag}(M_X^2,M_Y^2)U^\dagger,
\quad \Phi=U\begin{pmatrix}X\\Y\end{pmatrix}.
\]

The symbolic calculation uses a general SU(2) flavour rotation up to field
rephasings, with a nontrivial complex phase and an exact rational parameter
for its mixing angle. It proves \(U^\dagger U=I\), the componentwise
quadratic-potential identity, epsilon invertibility and free-kinetic identity.
These are **analytic parameter-map identities**, not substitution into the
actual Matchete post-rotation `.wl` Lagrangian.

## Native extra-quartic transport

Every *actual saved* Stage22 extra quartic is parsed to exact integer/rational/
radical components, rejecting floating-point data. The eight exceptional
models are checked via exact rational/radical component contractions for two
nontrivial unit phases and scalar assignments per source tensor. This yields
**92 exact witnesses** for 46 exceptional-model extra-quartic operator slots,
if the expected source inventory is present. The remaining four native extra
quartics belong to unmixed models and are validated against the source
inventory but require no exceptional rotation.

The universal tensor pullback is a consequence of distributivity,
\(T_{i_1\cdots i_n}\prod_a M^{(a)}_{i_a j_a}\chi^{(a)}_{j_a}\).
The finite exact witnesses check the source components and implemented
transformation maps. They do **not** prove that every term in the *actual*
post-rotation Wolfram Lagrangian has been symbolically identified with its
source term. The earlier Matchete acceptance and coupling-inventory checks
remain separate evidence.

## Conditional Matchete–RGBeta normalization

The saved Stage23 analysis provides 46 blocks for 50 extra-quartic slots.
For each block, this stage checks symbol registration against BOTH Stage22
catalogues, non-singularity of its exact source-derived \(R\) matrix, and

\[
 P=-R^{-1},\qquad RP=PR=-I.
\]

This is the **candidate** map under Stage23's stated relative-potential sign
and component-basis assumptions. The script records each exact candidate map
and its geometry class, including the 2×2 independent-invariant mixing blocks.
The source's oriented quartic prefactor `1/4` is independently checked as
metadata. A mathematical inverse check cannot certify that the *physical*
Matchete and RGBeta fields use the same SU(2) component convention. Stage23
explicitly identifies several triplet and conjugate-doublet mappings as
conditional. Those conditions are **not promoted to physical certificates**.

## Remaining after this stage

1. Compare the complete *actual Matchete expressions* for all source UV terms,
   including baseline scalar quartics, topology and Yukawa CGs, under a single
   constrained symbolic eigenmass/phase substitution. (This stage proves
   scalar quadratic/kinetic identities and extra-quartic tensor witnesses.)
2. Certify one physical Matchete↔RGBeta component and coupling dictionary,
   and check the four-scalar signed action for all 50 slots. Inverting saved
   candidate matrices is not such a proof.
3. If strict package-level symmetry registration is desired, register and
   test \(Z_2\) explicitly in Matchete; prior 16/16 tests establish parity of
   assembled UV terms, not discrete-charge registration.

**Interpretation of successful status:**
`ExactScalarIdentityAndConditionalNormalizationPassedNotFullPhysicalCertificate`.
This is deliberately a *conditional* Workstream A progress milestone, not
completion of physical UV equivalence or any EFT matching.
