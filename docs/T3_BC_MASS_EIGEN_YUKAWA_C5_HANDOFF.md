# B/C(−1): native mass-eigenstate Yukawa → one-scale C5 handoff

## What was already in the repository

**Not new:** `Lagrangian/RunGeneralT3MassEigenMatchete.wl` and
`scripts/RunGeneralT3MassEigenMatchete.py` already build and validate the
full mixed-scalar `X/Y` UV mass-eigenstate Lagrangian. They include both
native Yukawa CG tensors, 17 scalar quartics and the `HeavyMassBasis`
check. The saved `B_m1_final_matchete.json` and
`C_m1_final_matchete.json` report `HeavyMassBasis=True`, both Yukawa
interactions present and physical matching not performed.

`output/general_t3_mass_eigen_matchete/one_loop_probe/` already contains
Matchete one-loop *direct all-heavy* ordered Weinberg coefficients. The
existing `scripts.RunT3QuadraticEndpointC5` uses the saved expression,
frozen reference Yukawas and the candidate EFT1 scalar endpoint.

## New code in this package

`scripts/RunT3BCMassEigenYukawaC5Handoff.py` is a **read-only audit and
numerical handoff**, not a new symbolic Lagrangian or physical matcher.
It requires the user's already-passed `bc_live_registered_cg_v1` gate,
the existing native heavy-basis and one-loop Matchete certificates, the
saved scalar mass parameters, and the conditional one-scale neutrino mass
matrices. It verifies that their inputs and numerical masses agree and
saves the actual 3×N mass-eigenstate Yukawa matrices for both X/Y.

For the verified project convention,

\[
\Phi=\begin{pmatrix}S_1\\\epsilon S_2^*\end{pmatrix}
=U\begin{pmatrix}X\\Y\end{pmatrix},\quad
\epsilon=\begin{pmatrix}0&1\\-1&0\end{pmatrix}.
\]

The source's two **different fermion structures** are retained: `y1`
contains `CConj[NewFermion]`, while `y2` contains `NewFermion`.
Native Matchete CG tensor projections have factors

\[
\alpha_1=1,\qquad
\alpha_2^{B(-1)}=-1,\qquad
\alpha_2^{C(-1)}=+1.
\]

Therefore the coefficients *in the native mass-eigenstate CG bases* are

\[
Y^{(1)}_{X,Y}=\alpha_1\,U_{1,X/Y}\,y_1,\qquad
Y^{(2)}_{X,Y}=-\alpha_2\,U_{2,X/Y}\,y_2.
\]

The minus in the second expression arises from
`Bar[S2] = -epsilon × (U21 X + U22 Y)`; it is not a free sign choice.
Complex `U` values are retained with their phases, and both source scalar
rows remain normalized and orthogonal. **Do not combine y1 and y2 into
one Yukawa vertex**: they carry distinct spinor charge-conjugation rules.

### Source benchmark (saved repository evidence)

| Model | mX [GeV] | mY [GeV] | y1X/y1 | y1Y/y1 | y2X/y2 | y2Y/y2 |
|---|---:|---:|---:|---:|---:|---:|
| B(−1) | 1035.583 | 1188.529 | −0.046900 | +0.998900 | +0.998900 | +0.046900 |
| C(−1) | 1032.613 | 1185.998 | −0.046505 | +0.998918 | −0.998918 | −0.046505 |

The saved **conditional** Takagi masses in eV are:

- B(−1): `[0.02734103876, 0.02072521457, 0.00361669005]`
- C(−1): `[0.05353064456, 0.02939104094, 0.00561172448]`

These are already computed one-scale results. This package does not
claim they are new predictions. The default run rechecks the saved 3×3
complex mass matrices, their Takagi singular values, source file digests,
scalar masses and Yukawa input provenance. The `--recompute-c5` option
**re-evaluates** the existing large Matchete coefficient with the
existing `numerical_c5` evaluator and compares the resulting neutrino
mass matrices, without regenerating the Matchete expression.

### Install and run (PowerShell, from repository root)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "C:\path\to\T3_BC_mass_eigen_yukawa_c5_handoff_v1.zip" -DestinationPath . -Force
python -X utf8 -m pytest -q tests/test_t3_bc_mass_eigen_yukawa_c5_handoff.py
python -u -X utf8 -m scripts.RunT3BCMassEigenYukawaC5Handoff --all --run-id bc_yukawa_c5_v1
```

For a **new independent numerical replay from the saved one-loop
expression** (slower than the metadata/consistency check):

```powershell
python -u -X utf8 -m scripts.RunT3BCMassEigenYukawaC5Handoff --all --run-id bc_yukawa_c5_replay_v1 --recompute-c5
```

Do not reuse a run ID; the script refuses to overwrite existing run
folders. New output goes to
`output/audits/t3_bc_mass_eigen_yukawa_c5/<run-id>/`:
`B_m1_yukawa_c5.json`, `C_m1_yukawa_c5.json`, `summary.json`.

Inputs include current `output/general_t3_mass_eigen_matchete` reports,
`output/general_t3_mass_basis` masses, the recent registered CG summary,
`output/studies/general_scalar_extended/quadratic_endpoint_c5` and the
frozen fitted-reference Yukawas. Source hashes are checked.

### Scope / blockers

**Verified if this run passes:** completeness of the saved native Matchete
heavy-basis/Yukawa certificates, field/mixing signs and eigenmasses,
flavour-dependent Y1/Y2 mass-eigenstate vertices, saved candidate
one-scale C5/mnu numerical consistency, and optionally actual replay
of the previously extracted one-loop coefficient.

**Not certified:** physical running of Yukawa couplings on the full candidate
background, complete Matchete/RGBeta conventions across *all* scalar
couplings, independent direct live subtraction of full actions, fermion-first
EFT1 Wilson threshold matching, separate X/Y decoupling, finite scalar
threshold contributions, or a physical 16-model neutrino fit.

In particular, the saved scalar state is at 1000 GeV, below both X/Y
masses in these two benchmark points. That endpoint is not a valid
replacement for evaluation at each separate scalar threshold in a
sequential EFT scheme.

No production or GitHub files are changed.
