# Workstream C v16 — exact formal common-scale envelope (B(-1)/C(-1))

Paths in repository:

- `scripts/AssembleT3BCFormalFiniteScaleEnvelope.experimental.py`
- `scripts/RunT3BCFormalFiniteEnvelope.ps1`
- `docs/WorkstreamC_FormalFiniteEnvelope_v16_README.md`

**Purpose**: make direct numerical and symbolic use of the **verified** v15 neutral-component identity `P_F + P_Y + P_X = 0`, while explicitly refusing to identify a raw amplitude sum with a physical matched Wilson coefficient. This is a finite **formal source combination**, not renormalized matching. No source files, RGEs, UV Lagrangian or GitHub content are modified.

For each of three independent ordered Yukawa structures (y1*y1, y1*y2, y2*y2), strip the common `hbar/MF` from the native PL/Warsaw coefficient. Let `K_a` denote the formal scale-independent raw additive terms, `P_a` the UV-pole residues (equal algebraically to the corresponding log coefficients), where `a=F,Y,X`. The **source identity** verified in v15 is

`P_F + P_Y + P_X = 0`.

Therefore, at a common scale, the formal raw sum has no explicit renormalization-scale logarithm:

`K_formal = K_F + K_Y + K_X + P_F ln(M_X^2/M_F^2) + P_Y ln(M_X^2/M_Y^2)`.

Here `K_X=0` **only as a raw standalone constant**, not as a physical threshold. All scalar and fermion denominators, signs, SU(2) epsilons, H.c. conjugations and operator normalization are parsed from the original saved B/C source by the existing project-owned parsers. The conversion to the 1/2-normalized conventional Weinberg tensor is `C5_formal=A_PL+A_PL.T`; do **not** add PR a second time. The common loop factor `hbar/MF` is conventionally `1/(16 pi^2 MF)`. The optional numeric mode provides a symmetric matrix for a **formal, nonphysical** diagnostic only.

**Not certified:** Whether `Gamma_high - Gamma_low + counterterms` yields this sum, subtraction scheme, any finite EOM/field/parameter shifts, complete source algebraic equivalence in native Matchete, running, physical `C5`, `m_nu` or phenomenology. The v15 pole identity is *necessary consistency evidence* but is **not** a counterterm calculation. The independent component-equivalence audit previously established the 25-sector HHXY residual zero in scalar components, not full native Matchete equality; do not repeat it.

## PowerShell

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$HOME\Downloads\T3_WorkstreamC_FormalFiniteEnvelope_v16.zip" -DestinationPath . -Force
& .\scripts\RunT3BCFormalFiniteEnvelope.ps1
```

Optional source-anchored numerical matrix (one B/C model or a B/C manifest, same input schema as v14):

```powershell
python scripts\AssembleT3BCFormalFiniteScaleEnvelope.experimental.py `
  --numeric-input path\to\benchmark_numeric_inputs.json `
  --output output\audits\t3_bc_formal_finite_envelope_v16\numeric.json
```

Output: `output/audits/t3_bc_formal_finite_envelope_v16/summary.json`. The script fails closed for missing raw source, modified SHA from v14, v15 pole cancellation failure, unexpected chirality/CG or failed conjugation. `--require-physical` explicitly fails.
