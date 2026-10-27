# B/C(-1) — native Matchete mass eigenfields and matching preparation

This experimental batch follows the successful `RunGeneralT3MixedScalarMassBasis --all`
outputs. It does **not** modify production UV/EFT states, RGBeta, or threshold matching.

## Physics and representation

For `B_-1` and `C_-1` the original scalar charges are
`Y(S1)=-1/2`, `Y(S2)=+1/2`. The numerical convention is

\[
\epsilon=\begin{pmatrix}0&1\\-1&0\end{pmatrix},\quad
\Phi=\begin{pmatrix}S_1\\\epsilon S_2^*\end{pmatrix}
=U\begin{pmatrix}X\\Y\end{pmatrix},\qquad
U^\dagger A U=\mathrm{diag}(m_X^2,m_Y^2).
\]

`X` and `Y` are **distinct heavy complex SU(2) doublets**, both with `Y=-1/2`.
The Matchete heavy free Lagrangian contains only their diagonal masses `MSX`,
`MSY`. The original off-diagonal quadratic term is represented in the numerical
matrix `A`; it is **not** reinserted into the diagonal Matchete Lagrangian.

The Wolfram script loads the existing ordinary model and all seven extra Matchete
quartic candidates, constructs two native pseudoreal SU(2) bridge tensors, and
rewrites all original **eleven baseline interaction candidates** and **seven extra
quartic candidates**. The resulting Matchete Lagrangian contains the SM, heavy
T3 fermion, heavy X and Y, both Yukawas, and all quartic structures. This is a
**UV interaction-basis verification**: the provided numerical rotation comes
from the scalar-only EFT1 endpoint; the Yukawa vertices are not claimed to be
active in that EFT1.

The script fails if original scalar fields remain, any source vertex is absent,
`CheckLagrangian` fails, or the detailed `HeavyMassBasis` check is not `True`.

## Key limitations

- Both native Matchete bridge Clebsches are individually checked to be
  proportional to the antisymmetric epsilon tensor, aligned to the fixed
  numerical epsilon component convention and verified unitary. This is a
  component-order check *inside Matchete*, not an independent physical
  Matchete-to-RGBeta source-tensor map or quartic normalization certificate.
- The numerical unitary matrix is used as an interaction-coordinate choice.
  Matchete's free scalar masses `MSX`, `MSY` are independent symbolic couplings;
  they are not inferred to reproduce a perturbative EFT1 mass RG trajectory.
- No `Match[...]`, scalar threshold Wilson matching, or physical neutrino-mass
  prediction is performed. The matching flags remain false.
- The scalar-only EFT1 sample at 1000 GeV is below **both** benchmark scalar
  thresholds (~1030 and ~1185 GeV). Proper sequential matching must first run
  parameters to the *heavier* matching scale, integrate out Y, run the remaining
  EFT, and match X at its own threshold. Using the 1000 GeV endpoint directly
  as either threshold boundary is not justified.
- This installation contains no Wolfram kernel. Native `CheckLagrangian`
  success can only be established by running on your Windows installation.
- If the full numerical source has changed since the saved mass-basis JSON,
  preparation fails rather than quietly reusing stale coefficients.

## Install (PowerShell, repository root)

```powershell
Expand-Archive -LiteralPath "C:\path\to\T3_Native_Mass_Eigen_Matchete_Batch.zip" -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_mass_eigen_matchete.py
python -m scripts.RunGeneralT3MassEigenMatchete --all --prepare-only
python -m scripts.RunGeneralT3MassEigenMatchete --model B_-1
python -m scripts.RunGeneralT3MassEigenMatchete --model C_-1
python -m scripts.RunGeneralT3MassEigenMatchete --all
```

Output directory: `output/general_t3_mass_eigen_matchete/`.

Each model gets a `*_parameters.json`, a Wolfram `*_matchete.json`, and a full
`*_matchete.log`. Summary JSON includes `physical_matching_performed: false`.
A non-zero exit from Wolfram or any failed heavy-basis diagnostic is an actual
blocker; use the log, rather than assuming `HeavyMassBasis` passed.

## Source status

Checked against the public read-only GitHub repository on 2026-10-27:
`Lagrangian/model/{T3Fields,LagrangianBuilder,SU2Invariants}.wl`,
`Lagrangian/interactions/{T3Topology,GeneralScalarOperators}.wl`, and
`RGE/general/generated/general_rgbeta_quartic_registry.json`.
The two prior mixed-basis Python modules and saved JSON schemas were inspected
from the delivered mass-basis ZIP. GitHub was not modified.
