# Workstream A — final conditional UV evidence

This is an opt-in **read-only** evidence collector. It does **not** modify the
production Lagrangian, run `Match[]`, extend RGEs or perform EFT matching.
All results are per the ordinary 16-model T3 catalogue, exact assumed odd-BSM
`Z_2`, and the repository's current gauge/conjugation conventions.

## Installation

Install the new file under `Lagrangian/RunT3UVFinalEvidence.py` and this
explanation under `docs/WORKSTREAM_A_FINAL_UV_EVIDENCE.md`. Existing files are
not replaced by the package.

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
python Lagrangian\RunT3UVFinalEvidence.py --self-test
python Lagrangian\RunT3UVFinalEvidence.py --all
```

If a previous output already exists, back it up first and use `--overwrite`.
Default output directory:
`Reports/output/workstream_a/final_uv_evidence/`.

The runner requires the passing prior
`Reports/output/workstream_a/uv_closure_stage/uv_closure_stage.json`, all six
passing self-bilinear eigenbasis reports, the two passing B/C `MF` eigenbasis
reports and the B/C prepared mass matrices. It also requires the **existing**
`output/general_t3_stage22/matchete_tensors.json`. No Wolfram kernel is run.

## Independent fermionic census

The BSM content has `F`, `S1`, `S2` color singlets and odd under the exact
`Z_2`. With standard lepton doublet `L` of hypercharge `-1/2`, the gauge
singlets of the modeled renormalizable `Lbar`–`F/Fc`–`S/Sdag` Yukawa class
obey

\[
\frac12 + Y(F\text{ or }F^c)+Y(S\text{ or }S^\dagger)=0,
\qquad 2\otimes d_F\otimes d_S\supset 1.
\]

All four scalar orientations and two fermion-charge orientations are enumerated.
The `S1` and `S2` dimensions satisfy `dS=dF±1`, hence a unique SU(2) singlet
exists for each hypercharge-allowed orientation. The right-handed SM charged
lepton singlet cannot couple because `dF != dS`. One colored SM quark cannot
make a color singlet with two colorless BSM fields. `FFH` cannot make an SU(2)
singlet, because `jF ⊗ jF` contains only integer isospins but the Higgs has
`j=1/2`. A single BSM field mixing with an SM fermion is odd under `Z_2`.

For `B(-1)` and `C(-1)`, `F` is a neutral real SU(2) representation and can be
Majorana, so `F^c` and `F` are not counted as separate fields/operators.
For the resulting chiral Yukawa basis:

- Generic models: two independent Yukawa coupling matrices.
- `A/D/E` at `alpha=-2` or `alpha=0`: three matrices.
- Total over 16 models: **38** matrices; **six** extra relative to topology.

These are *coupling matrix types*, not the number of independent complex
entries in a generic flavour matrix. The implemented three heavy-flavour
masses describe an already-diagonal mass basis, not arbitrary independent
flavour rotations. A full cross-package Lorentz/flavour/spurion basis proof
remains outside this script.

## Native SU(2) tensor and normalization evidence

The runner reads the existing Stage22 **actual Matchete invariant tensors**:
50 independently-labelled extra-quartic slots across the 16 models.
It requires each tensor to have the correct ordered SU(2) axes, expected
canonical species, source prefactor `1`, nonzero numerical components and
linearly independent CG vectors for every repeated structure. It cross-checks
model-by-model counts with the independent `Sym^n(SU(2))` scalar basis and the
actual UV Matchete-builder counts in the prior Z2 reports.

For the eight exceptional mass bases, it evaluates each actual native extra
quartic tensor before and after its full field substitution and checks the
complex interaction plus its Hermitian conjugate. This is a **numerical
component-level equality** in the chosen Matchete-native representation;
not a test of RGBeta normalization, and not a symbolic identity imposed in
the existing Mathematica eigenbasis driver. Baseline quartic, topology and
Yukawa CG **source tensors** are not processed in this stage.

## Mass parameter maps

For the six neutral self-bilinears, if the full scalar potential is
`V2=m² S†S + b(S*)² + b* S²` in the native Cartesian basis with
`b = -kappa lambda_extra` (the original Lagrangian contains the negative of
this potential), the exact field replacement

\[
S=e^{i\arg(b)/2}\frac{R+iI}{\sqrt2},\qquad
m_R^2=m^2+2|b|,\quad m_I^2=m^2-2|b|
\]

diagonalizes `V2`. `kappa` is the model's **extracted** native bilinear tensor
factor (from its saved passing report). For `b=0`, the phase is arbitrary.

For `B/C(-1)`, writing `Phi=(S1,epsilon S2*)`, the Hermitian mass-squared
matrix is `[[m1²,mu],[mu*,m2²]]`. Its exact eigenvalues are

\[
\lambda_\pm=\frac{m_1^2+m_2^2}{2}\pm
\sqrt{\left(\frac{m_1^2-m_2^2}{2}\right)^2+|\mu|^2}.
\]

The earlier numerical B/C module additionally checks the epsilon and unitary
field rotation. This stage symbolically proves the generic eigenvalue formula
and re-reads existing B/C parameter/eigenbasis evidence; it does not
constrain new Matchete mass spurions to these formulas.

## Interpretation and remaining work

A successful status is
`ConditionalUVOperatorEvidencePassedNotFullCompleteness`, **not** a fully
certified UV Lagrangian. Outstanding: symbolic Matchete mass-spurion map;
actual rotated baseline and Yukawa/topology native CG normalizations;
Matchete↔RGBeta dictionary; explicit Matchete discrete-Z2 registration; and
full general Dirac/Majorana flavour convention documentation. Workstreams
B, C and D are unchanged.
