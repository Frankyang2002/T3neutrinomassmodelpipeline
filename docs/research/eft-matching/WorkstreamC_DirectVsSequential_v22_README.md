# Workstream C v22 — Direct all-heavy vs sequential full-mass finite matching

## Purpose

Compare two **distinct Matchete matching paths** for archived B(-1) and C(-1) at order `y_i y_j lambda`:

1. **Direct**: simultaneous full one-loop `Match[UV,5,1]` with F, X and Y heavy. Source: `output/general_t3_mass_eigen_matchete/one_loop_probe/<stem>_final_one_loop_probe_weinberg_ordered_coefficient.inputform.txt`.
2. **Sequential**: F-hard, Y, X source operators (v14–v17), combined into the exact three-mass XX/XY/YY triangle with the tested scalar-bubble hard/soft subtraction (v19).

The direct source is an **ordered** lepton-flavour coefficient A, not symmetric C5. Its matching report explicitly states `C5=A+A.T` in the one-half-normalized Weinberg convention. The F/Y/X source coefficients use the independent v11.2 Warsaw operator identity and v16/v17 projection.

This resolves a major practical question: whether the *finite mass dependence, signs, 1/2 factors and normalization* from these two Matchete paths agree for B/C. It is stronger than comparing their UV residues but is **not an external Feynman-rule derivation**: both paths originate in Matchete, and frozen high-scale couplings are not a running solution.

## Exact implementation

- Strictly source-checks the direct matched EFT metadata and prior v19/v21 model scope. No GitHub writes.
- Parses all **480** additive monomials of the direct A (B and C separately). An arithmetic-only Python AST evaluator rejects unknown symbols, function calls or malformed Yukawas. `mpmath` 65-digit arithmetic protects against large cancellations between direct rational-mass terms.
- Direct term inventory: `y1_i y1_j` = 144, `y1_i y2_j` = 192, `y2_i y2_j` = 144. Quartic inventory: `lambdaT3` = 192 and the two T3 extra LLHH quartics = 144 each.
- Includes the native outer `/3` for C(-1) if and only if the matched expression has that exact structure.
- Independently recomputes W_XX, W_YY, W_XY from the F/Y/X raw matching operators, using v16 source hash checks (and v17 signed exact three-mass scalar function).
- Tests four quartic scenarios: baseline-only; each extra quartic separately at a nonzero complex value; all three together. Tests each of the three different heavy-fermion masses from the fitted reference. **12 direct vs sequential comparisons per model**, each comparing the three ordered Yukawa monomial coefficients; an overall discrepancy, sign flip or factor-of-two discrepancy FAILS.
- Additionally builds the full complex-symmetric 3x3 matrix from the direct ordered coefficients for the frozen reference case, compares it with the archived v18 full-mass matrix, and exports exploratory singular masses. This validates the full flavour and Majorana convention without relabelling these results as physical.
- Writes actual signed comparison values and residuals to JSON even if the equality fails.

Source Wilson coefficients are evaluated in GeV^-1 with `hbar=1/(16 pi^2)`. Three heavy masses are used **as independent diagonal columns**: the archived direct Matchete coefficient has a single formal MF symbol, and its replacement by three diagonal masses is valid for the assumed diagonal, noninteracting heavy-flavour columns. The Yukawas are frozen from the 10^7 GeV reference; no UV/EFT1/EFT2 RGEs are run.

## Run

In PowerShell from the project root:

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$HOME\Downloads\T3_WorkstreamC_DirectVsSequential_v22.zip" -DestinationPath . -Force
& .\scripts\RunT3BCDirectVsSequential.ps1
```

Or manually:

```powershell
python scripts\DirectVsSequentialOneLoop.experimental.py --self-test
python scripts\DirectVsSequentialOneLoop.experimental.py --output output\audits\t3_bc_direct_vs_sequential_full_mass_v22\summary.json
```

## Outputs and interpretation

Output: `output/audits/t3_bc_direct_vs_sequential_full_mass_v22/summary.json`.

- `direct_vs_sequential_full_mass_agreement=true` means all 12 checks for that model fall below the relative tolerance of `2e-8`; inspect the recorded signed coefficient arrays and matrix-level residual.
- If false, the script raises a nonzero exit status **after writing the residual report**. Do not interpret a mismatch as an ultraviolet inconsistency without first checking conventions and truncation differences.
- `physical_C5=null` throughout; no full finite field/EOM and parameter matching, RG flow, general nonzero mixed quadratic or physically matched running couplings are certified.

## Assumptions / caveats

This comparison is limited to the archived two-Yukawa B/C(-1) UV model with one baseline and two additional direct quartics, no independent extra mixed quadratic turned on, unbroken electroweak phase, leading dimension-five and one loop. It tests absolute matching normalization **between two independent paths using the same Matchete engine**, not an external first-principles sign convention for the complete UV quantum effective action. There is no claim that bare or MSbar-renormalized EFT Wilson coefficients can be combined indiscriminately, or that the benchmark Yukawas were run to the threshold.

No production source modules, UV Lagrangian, beta functions or Github files are modified. Only this new experimental script, PowerShell runner and documentation are included.
