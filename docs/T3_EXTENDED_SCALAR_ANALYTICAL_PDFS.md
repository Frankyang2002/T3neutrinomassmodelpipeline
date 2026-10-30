# Extended-scalar analytical PDFs (all 16 T3 assignments)

This is a **read-only report generator**, not a new physics derivation. The output is a standalone extended-scalar addendum in the same three analytical categories as the previous reports. It DOES NOT replace the old PDFs or rewrite the old Matchete Lagrangian. It does not run Matchete, Wolfram, RGBeta or any neutrino fits.

## Sources used

- `RGE/general/generated/general_rgbeta_quartic_registry.json`: model-dependent additional quartic operator inventory and candidate RGBeta SU(2) tensor contractions; 50 extra coupling slots = 38 newly named + 12 structurally reused, across 16 models.
- `docs/general_scalar_trial_operator_manifest.json`: independently catalogued SU(2) operator multiplicities and eight allowed extra holomorphic scalar quadratic structures.
- `output/general_t3_quadratic_eft1/{A_p0,...}_{uv,eft1}_quadratic_live_rgbeta.json`: 16 x 2 saved **candidate** one-loop beta exports. The builder validates source status, model labels, loop convention, key completeness and registry counts. UV and scalar-EFT equations are reported separately.

The report PDFs now display couplings in the same paper-style notation as the original reports. The topology coupling remains $\lambda_5$, previously used B/C symbols remain $\lambda_7$--$\lambda_{12}$, and newly named independent quartics have globally stable labels $\lambda_{13},\lambda_{14},\ldots$. The exact 1:1 source-name-to-paper-symbol dictionary is written to `Reports/output/analytical/extended_scalar/manifest.json` as `display_notation`. Each PDF includes a notation key with field content and model assignments. No raw `lambdaT3GeneralQ...` names appear in the report body.

The numerical coefficients and Mathematica InputForm structure of all saved RGE expressions are **unchanged**. RGE listings display the paper symbols in inline LaTeX math mode instead of raw implementation identifiers; this is a notation-only operation, not an independent algebraic simplification or a Matchete/RGBeta mapping certification.

## Output categories

| New PDF (project-relative path) | Content |
|---|---|
| `Reports/output/analytical/extended_scalar/Lagrangian/UV_extended_scalar.pdf` | Operator count per model; every additional quartic monomial, newly numbered vs established scalar couplings; allowed extra quadratic monomials; the unchanged original potential is referred to, not silently regenerated |
| `Reports/output/analytical/extended_scalar/GroupFactors/GF_UV_extended_scalar.pdf` | Representation/Casimir context, model-by-model extra invariant tensors and multiplicities, and **source-exported** beta(lambdaT3) expressions. **Not** an independent certified set of new group factors |
| `Reports/output/analytical/extended_scalar/RGE/RGE_UV_EFT1_extended_scalar.pdf` | Full candidate scalar-quartic and scalar-mass beta functions, including original lambdaT3, for UV and EFT1, each model. Pass `--rge-scope all` to include all exported Yukawa, fermion and gauge equations too |

The original PDFs remain at `Reports/output/analytical/{Lagrangian,RGE,GroupFactors}/`.

## PowerShell, from repository root

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath .\T3_ExtendedScalar_PaperNotation_TeXFix.zip -DestinationPath . -Force
python -m pytest -q tests\test_build_t3_extended_scalar_analytical_pdfs.py
python -m scripts.BuildT3ExtendedScalarAnalyticalPDFs --check
python -u -m scripts.BuildT3ExtendedScalarAnalyticalPDFs
```

The generator requires `pdflatex` (MiKTeX or TeX Live) to be on `PATH`. To make only `.tex`, use `--tex-only`. To include every exported running coupling rather than scalar sector only, use:

```powershell
python -u -m scripts.BuildT3ExtendedScalarAnalyticalPDFs --rge-scope all
```

Only the new `Reports/output/analytical/extended_scalar/` folder is written. A `manifest.json` records SHA-256 values for every imported registry, manifest and beta export. The generator fails instead of filling in missing or changed source equations.

## Physics status and limitations

- The operator **existence** and inventory come from the repo's scalar operator enumeration and Matchete audits.
- The additional quartic beta expressions are already computed candidate native RGBeta exports; this script **transcribes** them and does not claim an independent group-theory proof.
- The correspondence between Matchete scalar-potential coefficients and RGBeta coefficients is **not** generally normalised/certified for the extended operator basis.
- The additional holomorphic mass operators were catalogued and their candidate beta functions saved, but complete UV Matchete mass-eigenstate and sequential threshold matching is not implemented.
- The group-factor addendum therefore distinguishes SU(2) tensor contractions and RGBeta source coefficients from independently verified group factors. **Do not cite it as a full normalised scalar group-factor derivation.**
- No general scalar C12/C5 Wilson matching is inferred from the renormalisable UV/EFT1 beta exports.

For the thesis, this is suitable as a **supplementary technical appendix / exploratory extension** until the general scalar coefficient conventions and matching are certified.

### Optional: preserve the historical report pages in the new PDFs

To copy the old analytical UV TeX reports into the new documents and append the scalar extension, use:

```powershell
python -u -m scripts.BuildT3ExtendedScalarAnalyticalPDFs --append-original
```

This builds separate `*_with_original.pdf` files under the same `extended_scalar/` directory. The old TeX/PDF files remain unchanged. Because the original reports are very large, use the standalone addenda for a faster build. If the old TeX has compilation issues, the standalone reports do not depend on it.

## Notation-only revision

This package replaces the full `scripts/BuildT3ExtendedScalarAnalyticalPDFs.py` and its unit tests, not the underlying source data or original analytical reports. It does not rename RGBeta symbols in calculation files. After extraction, regenerate the new extended-scalar reports with the commands above. The metadata `manifest.json` keeps the exact original identifiers for audit and reproducibility. No additional Wolfram calculation is needed.

The earlier `Reports/output/analytical/extended_scalar/previews/` PDFs used raw identifiers and are superseded. Open the regenerated full PDFs in the three category folders instead. The display numbering is deterministic from the current registry and must be regenerated if that registry changes.

## TeX conjugation fix (30 October 2026)

The previous paper-notation generator could emit an invalid double superscript
(e.g. `\mu_{S_2}^2^{*}`) for a source expression such as
`Bar[muT3GeneralS2Sq]`. This triggered the `Double superscript` error while
compiling `RGE_UV_EFT1_extended_scalar.tex`. The corrected formatter wraps
any decorated symbol that already contains a superscript before applying the
complex-conjugation star, producing `(\mu_{S_2}^2)^{*}` in the PDF.
It also handles conjugated baseline couplings carrying superscripts,
including `lambdaH1`, `lambdaS1Adj`, and `lambda12Cross`.

This is a **typesetting-only** correction. It does not alter native beta
expressions, physics inputs, or any historical report. It is safe to re-run
the same output command; the existing generated TeX/PDF files in the
`extended_scalar` output folder will be replaced by corrected reports.

Validation: unit tests cover conjugated scalar bilinears and decorated
quartics, and all three PDFs were compiled against a *synthetic* 16-model
fixture containing the problematic mass-conjugation forms. The user's full
native beta exports still need to be compiled locally.

## Mathematical RGE display (revision)

The RGE PDF converts saved `report_betas` Mathematica InputForm expressions into mathematical LaTeX using SymPy's Mathematica parser and a custom printer for the repository's existing paper notation. Matrix `Dot`, `Tr`, `Trans`, and `Bar` are rendered as matrix multiplication, traces, transposes, and complex conjugation, respectively. Additive terms are printed on separate, page-breakable lines. Long terms are width-constrained. Coefficients are not intentionally rescaled. Parsing may reorder algebraically commutative factors/terms for display. Exact saved raw expressions and input hashes are unchanged; the renderer does not validate physical RGBeta–Matchete normalisation.

**Important:** These candidate source exports contain `report_betas` (InputForm) rather than the `report_beta_latex` field consumed by the historical `Reports/RGEComparison.py` code. The revision provides the missing symbolic-to-LaTeX conversion, while preserving the established coupling naming conventions.

Requires `sympy`, already used by the existing `Reports/RGEComparison.py` reporting module.
