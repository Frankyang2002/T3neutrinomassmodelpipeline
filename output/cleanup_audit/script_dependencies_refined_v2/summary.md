# Refined T3 script dependency audit (read-only)

- Checkout HEAD: `fe0f511104cce8d7d914d114e3fd884eb47c6a95`
- Script modules: **98**
- Python files parsed: **587**
- Text files scanned by first pass: **9879**
- Originally in broad `RETAIN_STATIC_CODE_OR_TEST_REFERENCE`: **77**
- From that group now needing *manual* review rather than presumed code/test retention: **11**

## All scripts by refined status

| Status | Count |
|---|---:|
| `RETAIN_DECLARED_WORKFLOW` | 11 |
| `RETAIN_OTHER_PYTHON_DEPENDENCY` | 24 |
| `RETAIN_TEST_EXECUTION_OR_IMPORT` | 31 |
| `RETAIN_TRANSITIVE_ACTIVE_DEPENDENCY` | 11 |
| `REVIEW_DOC_OR_CONFIG_ONLY` | 4 |
| `REVIEW_PYTHON_TEXT_MENTION_ONLY` | 1 |
| `REVIEW_SHELL_ENTRYPOINT` | 6 |
| `REVIEW_TEST_TEXT_MENTION_ONLY` | 10 |

## Breakdown of the original 77

| Refined status | Count |
|---|---:|
| `RETAIN_OTHER_PYTHON_DEPENDENCY` | 24 |
| `RETAIN_TEST_EXECUTION_OR_IMPORT` | 31 |
| `RETAIN_TRANSITIVE_ACTIVE_DEPENDENCY` | 11 |
| `REVIEW_PYTHON_TEXT_MENTION_ONLY` | 1 |
| `REVIEW_TEST_TEXT_MENTION_ONLY` | 10 |

## Review candidates from the original 77

These are **not** approved for archival. See `inventory.csv` and `strong_evidence.csv`.

| Script | New status | Code text | Test text | Docs | Saved result |
|---|---|---:|---:|---:|---:|
| `GenerateGeneralScalarWLManifest` | `REVIEW_TEST_TEXT_MENTION_ONLY` | 0 | 2 | 1 | 0 |
| `InspectGeneralScalarMatching` | `REVIEW_TEST_TEXT_MENTION_ONLY` | 0 | 2 | 4 | 0 |
| `ProbeT3ReducedUVScalarBoundaries` | `REVIEW_TEST_TEXT_MENTION_ONLY` | 0 | 1 | 3 | 0 |
| `ReportGeneralScalarTrial` | `REVIEW_TEST_TEXT_MENTION_ONLY` | 0 | 3 | 6 | 0 |
| `run_comparison_study` | `REVIEW_TEST_TEXT_MENTION_ONLY` | 0 | 1 | 5 | 0 |
| `run_regression` | `REVIEW_TEST_TEXT_MENTION_ONLY` | 0 | 1 | 6 | 0 |
| `RunGeneralT3QuadraticEFT1` | `REVIEW_PYTHON_TEXT_MENTION_ONLY` | 1 | 1 | 5 | 0 |
| `RunGeneralT3QuadraticMatcheteAudit` | `REVIEW_TEST_TEXT_MENTION_ONLY` | 0 | 1 | 16 | 0 |
| `RunGeneralT3Stage13` | `REVIEW_TEST_TEXT_MENTION_ONLY` | 0 | 1 | 2 | 0 |
| `RunGeneralT3Stage14` | `REVIEW_TEST_TEXT_MENTION_ONLY` | 0 | 2 | 2 | 0 |
| `RunGeneralT3Stage15` | `REVIEW_TEST_TEXT_MENTION_ONLY` | 0 | 3 | 3 | 0 |

## Declared-workflow dependency closure

Only AST imports and literal execution calls contribute to the closure:

`AnalyzeT3BCMatchingConventions`, `BuildGeneralScalar48FinalReports`, `BuildReducedScalarExistingInteractiveReports`, `BuildReducedScalarLambdaStyleReports`, `BuildT3ReducedScalarFinalResults`, `FitT3BPlus1GeneralThreeStudies`, `RunGeneralT3FrozenScalarC5`, `RunGeneralT3NumericalC5`, `RunGeneralT3ReferenceC5`, `RunT3BPlus1GeneralEFT1`, `RunT3GeneralScalar48Campaign`, `RunT3GeneralScalarBCExploratoryFit`, `RunT3GeneralScalarBCPilot`, `RunT3ReducedScalarSobolComparison`, `RunT3ZeroBoundaryLowEnergyC5`, `RunT3ZeroExtraBoundaryBenchmark`, `Stage30CoupledA0UVRunning`, `Stage31CoupledUVRunning`, `Stage33FermionThresholdHandoff`, `Stage34IntermediateEFTRunning`, `Stage35ScalarThresholdHandoff`, `Stage36FinalWeinbergRunning`

## Evidence interpretation

- `PYTHON_IMPORT` and `TEST_IMPORT`: syntax-level imports (including relative imports in scripts/).
- `*_LITERAL_MODULE_LOAD`: literal importlib/runpy/__import__ module loads.
- `*_LITERAL_EXECUTION`: literal script/module arguments of known execution APIs.
- `CODE_REFERENCE` and `TEST_REFERENCE` from the first pass are **text mentions**, not evidence of an actual import.
- `SHELL_REFERENCE` remains a manual executable entrypoint review, not proof of inactivity.
- Saved-result references are preserved as provenance; do not rewrite historical JSON to remove them.
- This tool does not move, delete, rename or execute any script. No `--apply` mode exists.
