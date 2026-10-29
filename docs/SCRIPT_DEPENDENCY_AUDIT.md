# Historical script dependency audit (read-only)

**Scope:** identification of direct and indirect dependencies before any possible archiving of old numerical or symbolic scripts. This tool does **not** delete, rename, archive, import, or execute project physics scripts. It only reads repository text and writes an audit report.

At the inspected GitHub commit `fe0f511104cce8d7d914d114e3fd884eb47c6a95` (29 October 2026), `scripts/` contains 99 Python files, including `__init__.py`. The tool inventories all non-`__init__` top-level Python script modules **from your own local checkout**, not from a hard-coded snapshot.

## Installation

The package adds these two files, preserving project-relative paths:

- `tools/AuditScriptDependencies.py`
- `docs/SCRIPT_DEPENDENCY_AUDIT.md`

From the local repository root in PowerShell:

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
$zip = Join-Path $HOME 'Downloads\T3_Script_Dependency_Audit.zip'
if (!(Test-Path $zip)) { throw "ZIP not found: $zip" }
if ((Test-Path .\tools\AuditScriptDependencies.py) -or (Test-Path .\docs\SCRIPT_DEPENDENCY_AUDIT.md)) { throw 'Audit file(s) already exist. Review before overwriting.' }
Expand-Archive -LiteralPath $zip -DestinationPath . -Force
python -X utf8 -m tools.AuditScriptDependencies
Get-Content -Encoding UTF8 .\output\cleanup_audit\script_dependencies_v1\summary.md | Select-Object -First 55
```

The script does not require Wolfram, Matchete, RGBeta, a live internet connection, or numerical reruns. It does not run any `scripts/` module.

## Generated local files

- `output/cleanup_audit/script_dependencies_v1/summary.md`: human-readable classification and import-dependency list.
- `output/cleanup_audit/script_dependencies_v1/inventory.csv`: every local `scripts/*.py` module (except `__init__.py`), classification, source hash and reference counts.
- `output/cleanup_audit/script_dependencies_v1/references.csv`: evidence by source path, line number and snippet.
- `output/cleanup_audit/script_dependencies_v1/audit.json`: full machine-readable metadata including skipped files and parse/read errors.

## Classification meaning

- `RETAIN_DECLARED_WORKFLOW`: a documented active or currently required numerical/reporting entrypoint. This is a conservative preservation policy, not a guarantee that the script's numerical results are physically certified.
- `RETAIN_STATIC_CODE_OR_TEST_REFERENCE`: source imports the script or it is referenced in Python/tests. A test referencing an old stage runner counts as a dependency until separately reviewed.
- `REVIEW_DOCUMENTED_ENTRYPOINT`: referenced in documentation, shell scripts, or configuration; check whether the invocation remains relevant.
- `REVIEW_FROZEN_PROVENANCE`: appears in saved numerical output; retain the reproducibility reference unless an external archive and durable path are arranged.
- `REVIEW_NO_STATIC_REFERENCE`: no reference found in the scanned text. **This is not permission to delete**; external invocation or dynamically constructed names can be missed.

The AST parser identifies ordinary `import scripts.X`, `from scripts import X`, `from scripts.X import ...`, relative script imports, and literal `importlib.import_module('scripts.X')` or `runpy.run_module('scripts.X')` expressions. It also scans source, docs, shell commands, and saved outputs for literal script names. A text mention is less conclusive than an import. The scanner skips `.git`, virtual environments, cache directories, and its own generated report directory.

Input text files exceeding the default 8 MiB limit are reported as *skipped*, not assumed clean; increase the limit using `--max-text-mb N` when necessary.

## Do not retire yet

In particular, `BuildT3ReducedScalarFinalResults.py` and `BuildReducedScalarExistingInteractiveReports.py` import `RunT3GeneralScalar48Campaign.py`. The script's historical name therefore does **not** make it obsolete. The original figure/report builders, matching evaluators, threshold-running scripts and source-hash-verified fit records may have similar dependencies.

Review `references.csv` and the direct-import relationships before proposing any historical script archive. The frozen benchmark inputs and report generators must remain reproducible. Once a shortlist is reviewed against the local checkout, any archival should use an **external archive with original relative paths and SHA-256 manifests**, followed by `python -X utf8 -m pytest -q --tb=line`. This package performs no such archival.
