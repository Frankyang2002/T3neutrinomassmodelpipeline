# Script dependency audit — refinement of the 77 broad-reference modules

**Purpose:** reclassify the first-pass `RETAIN_STATIC_CODE_OR_TEST_REFERENCE` scripts without mistaking textual mentions for imports. This is a **read-only local audit**. No scripts, tests, results, PDFs, numerical fit records, or GitHub content are modified. No archival action is supplied.

## Provenance and scope

The first-pass `tools/AuditScriptDependencies.py` reported **98** script modules: 11 declared workflows, 77 with import/test/code mentions, and 10 documented-only entrypoints on the user's checkout at `fe0f511104cce8d7d914d114e3fd884eb47c6a95` (29 October 2026). The new audit rescans the **current local checkout** rather than relying on these counts or an old cached CSV. It uses the installed first-pass scanner to cover project text and independently parses all Python files to find concrete executable dependencies.

The refinement distinguishes:

1. **Syntactic Python imports:** `import scripts.X`, `from scripts import X`, `from scripts.X import ...`, and single-level relative imports inside `scripts/`.
2. **Literal dynamic imports:** literal `importlib.import_module('scripts.X')`, `__import__('scripts.X')`, or `runpy.run_module('scripts.X')`.
3. **Literal execution calls:** known `subprocess`, `os.system`, or `runpy.run_path` calls with an actual literal script/module target.
4. **Tests:** actual test imports/execution versus tests that merely mention a name in source strings or comments.
5. **Weak evidence:** Python comments/string mentions, shell wrappers, documentation, configuration, and frozen-output provenance. These are kept in separate columns; they are not automatically promoted to executable dependencies.

The closure of explicitly preserved workflows is computed by following confirmed script-to-script AST/import/execution edges, including transitive edges. Other executable dependencies and test imports remain protected by their own categories.

**Limitations:** constructed module names, reflection, subprocess commands assembled from variables, externally invoked entrypoints, and stale documentation may not be resolvable statically. `REVIEW` means **manual investigation**, *never* safe for deletion. Source metadata in historical JSON must not be edited to make a script appear unused. The earlier two full-study convenience wrappers remain potential review candidates, not automatic archive targets.

## Install (PowerShell)

Extract the ZIP at the project root, preserving exactly these new paths:

- `tools/RefineScriptDependencyAudit.py`
- `tests/python/test_t3_refined_script_dependencies.py`
- `docs/SCRIPT_DEPENDENCY_REFINEMENT.md`

The original `tools/AuditScriptDependencies.py` must remain installed. The new script imports its read-only `examine()` function; it does **not** replace or modify it.

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
$zip = Join-Path $HOME 'Downloads\T3_Script_Dependency_Refinement_v2.zip'
if (!(Test-Path $zip)) { throw "ZIP not found: $zip" }
$targets = @('tools/RefineScriptDependencyAudit.py', 'tests/python/test_t3_refined_script_dependencies.py', 'docs/SCRIPT_DEPENDENCY_REFINEMENT.md')
foreach ($target in $targets) { if (Test-Path $target) { throw "Already exists: $target" } }
Expand-Archive -LiteralPath $zip -DestinationPath . -Force
python -X utf8 -m pytest -q tests/python/test_t3_refined_script_dependencies.py --tb=line
python -X utf8 -m tools.RefineScriptDependencyAudit
Get-Content -Encoding UTF8 .\output\cleanup_audit\script_dependencies_refined_v2\summary.md | Select-Object -First 100
```

## Output files

Under `output/cleanup_audit/script_dependencies_refined_v2/`:

- `summary.md`: counts and named review candidates from the original 77, plus the active workflow transitive dependency closure.
- `inventory.csv`: all 98 script modules with prior and refined classes, exact hard versus weak evidence counts, and SHA-256 fingerprints.
- `strong_evidence.csv`: each resolvable import or literal execution, with referring source file, line and source excerpt.
- `audit.json`: all above information and first-pass weak references for additional review.

Run the full Python suite **only if production files are subsequently changed or scripts are archived**; no physics computation is performed by the audit.

Do not archive the 77-script group from classification totals alone. The next decision should be made from actual names in `summary.md` and the references in `strong_evidence.csv`, while preserving tested entrypoints, live imports, archived derivations, fitting provenance, and reproducibility.
