# Native extended C5 coefficient reduction — source-derived, conditional numerical bridge

This package adds new files only. It **does not alter** the production Matchete Lagrangian, RGBeta model, UV numerical solver, or hierarchical matching code.

## Exact source findings

The input is the previously exported set of six 3-heavy-flavour native Matchete reports and 18 coefficient expressions. The driver independently validates the three exported sectors for each model, including their input-name and lepton/heavy-index patterns.

Define the three-mass loop function using squared masses:

\[
I(a,b,c)=\frac{a b\ln(a/b)-a c\ln(a/c)-b c\ln(a/b)+b c\ln(a/c)}{(a-b)(a-c)(b-c)}.
\]

With \(\hbar=(16\pi^2)^{-1}\), and **under the explicit per-generation identification** `Coupling[MF,{},0] = Coupling[MF,{Index[r,NFlavor]},0] = M_r`, the exported ordered \(A_5^{pq}\) is exactly:

\[
A_5^{pq} = k\hbar\sum_{r=1}^{3}M_r I(M_r^2,M_{S1}^2,M_{S2}^2)
\left[\lambda_{T3}\overline{y_{1,pr}}\overline{y_{2,qr}}+
\overline{\lambda_{\rm extra}}\,
\overline{y_{a,pr}}\overline{y_{b,qr}}\right].
\]

`k` is \(-2/\sqrt3\) for A/D, and \(+2\sqrt{2/3}\) for E. For alpha=0, `(a,b)=(y1Conj,y2)` and `lambda_extra=lambdaT3ExtraQ021010I1`. For alpha=-2, `(a,b)=(y1,y2Conj)` and `lambda_extra=lambdaT3ExtraQ020101I1`.

The project reports symmetric \(C_5^{pq}=A_5^{pq}+A_5^{qp}\). **The native file itself does not have the symmetric flavour partner**, and must not be fed into a 3x3 matrix as though it were already symmetrized. The numerical evaluator implements the explicit p/q symmetrization.

**Interpretation**: These results establish a **conditional algebraic reduction** of the uploaded native one-step expressions. They do not prove the unindexed/indexed mass identification at the level of a native flavour-basis certificate, an absolute Matchete↔RGBeta coupling map, the complete general scalar quadratic sector, or hierarchical matching. They do not authorize using RGBeta or fitted-reference numbers as though they are already Matchete-basis values.

## Project paths

- `Numerical/experimental/T3ExtendedC5ExpressionReduction.py` — restricted InputForm parser, 6-model exact symbolic proof, optional conditional 3x3 numerical evaluator.
- `scripts/RunT3ExtendedC5ExpressionReduction.py` — CLI and JSON outputs.
- `scripts/TestT3ExtendedC5ExpressionReduction.py` — real native artifact validation and negative controls.
- `Numerical/configs/extended_c5_matchete_toy_A_p0.json` — deliberately **synthetic**, not a fitted or physical benchmark. Uses explicit Matchete-coordinate declaration for demonstration.
- `README_T3ExtendedC5ExpressionReduction.md` — these instructions.

## PowerShell commands (local repository only)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
$zip = "$HOME\Downloads\T3C5ExpressionReduction_Local.zip"
Expand-Archive -LiteralPath $zip -DestinationPath . -Force

# Numeric, parser and sign/index/loop negative controls, using your own real native files:
python scripts\RunT3ExtendedC5ExpressionReduction.py --self-test
python scripts\TestT3ExtendedC5ExpressionReduction.py

# Algebraically verify all six native matched expressions (no fresh Wolfram needed):
python scripts\RunT3ExtendedC5ExpressionReduction.py --all

# Optional: one deliberately synthetic, conditional Matchete-basis 3x3 example:
python scripts\RunT3ExtendedC5ExpressionReduction.py --all --input Numerical\configs\extended_c5_matchete_toy_A_p0.json
```

The source verification accepts `--zip PATH` to read the exact uploaded `ExtendedC5Expressions.zip` instead of the locally saved native export directory.

Outputs go to `Reports/output/workstream_a/extended_c5_expression_reduction/summary.json` and, when `--input` is used, `A_p0_conditional_c5.json` (or corresponding model). The report will always label a numerical result as **illustrative/conditional**.

### About the example

The sample input is not derived from an RGBeta trajectory or a neutrino-oscillation fit. It merely tests the Matchete-source formula, positive loop kernel, three flavour generations, charged-lepton basis rotation, and Takagi observable interface. To evaluate real parameters, provide a new JSON with independent Matchete-basis couplings; do not use `--input` to disguise a still-unknown coupling dictionary.

### Physical blockers retained

1. Native proof that the unindexed `MF` display placeholder can be identified with each heavy-generation `MF[r]` under the export's flavour sum; independent physical mass-basis and phase validation.
2. Absolute Matchete–RGBeta Yukawa, mass and extra-quartic dictionary, including triplet basis for A/D/E.
3. Complete hierarchical matching and general scalar quadratic effects; the uploaded expressions cover a simultaneous-heavy one-step truncation only.
4. Only after (1)–(3), physical neutrino masses and PMNS predictions from the full UV trajectories.
