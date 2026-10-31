"""Read-only static / optional test audit for a local T3 repository.
Run from repository root: python scripts/maintenance/Audit-T3Pipeline.py
Optional: python scripts/maintenance/Audit-T3Pipeline.py --run-tests
Writes output/audit/t3_pipeline_audit.json. Never modifies source files.
"""
from __future__ import annotations
import argparse
import ast
import json
from pathlib import Path
import re
import subprocess
import sys
from datetime import datetime, timezone

EXCLUDE={'.git','.venv','venv','__pycache__','.pytest_cache','.mypy_cache','node_modules','output','Reports','dist','build'}
KEYS=('solve_ivp','beta_values_to_derivative','evaluate_inputform_expression','Takagi','takagi','PMNS','pmns','MNSParameters','C5','c5','make_transfer','frozen','matching','wilson')

def scan(root:Path):
    rows=[]; imports={}; errs=[]
    for p in sorted(root.rglob('*.py')):
        relative=p.relative_to(root)
        if any(x in EXCLUDE for x in relative.parts):continue
        try:src=p.read_text(encoding='utf-8-sig'); tree=ast.parse(src,filename=str(relative))
        except (OSError,SyntaxError,UnicodeError) as exc:
            errs.append({'path':str(relative),'error':str(exc)});continue
        funcs=[n.name for n in ast.walk(tree) if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef))]
        references={k:len(re.findall(r'\b'+re.escape(k)+r'\b',src,re.IGNORECASE)) for k in KEYS}
        references={k:v for k,v in references.items() if v}
        imported=[]
        for n in ast.walk(tree):
            if isinstance(n,ast.Import): imported.extend(x.name for x in n.names)
            if isinstance(n,ast.ImportFrom) and n.module:imported.append(n.module)
        imports[str(relative)]=imported
        if references or ('main' in funcs and relative.parts[0]=='scripts'):
            rows.append({'path':str(relative),'functions':funcs[:100],'keywords':references,'imports':imported})
    return rows,imports,errs

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root',type=Path,default=Path.cwd())
    ap.add_argument('--run-tests',action='store_true',help='Run pytest (may take time, requires dependencies).')
    ap.add_argument('--pytest-args',nargs='*',default=['-q'])
    a=ap.parse_args();root=a.root.resolve()
    if not (root/'Numerical').is_dir() or not (root/'scripts').is_dir():
        ap.error('Point --root to the repository root containing Numerical/ and scripts/')
    rows,imports,errs=scan(root)
    tests=sorted(str(p.relative_to(root)) for p in (root/'tests').rglob('test_*.py')) if (root/'tests').exists() else []
    tests_result=None
    if a.run_tests:
        proc=subprocess.run([sys.executable,'-m','pytest',*a.pytest_args],cwd=root,text=True,capture_output=True)
        tests_result={'returncode':proc.returncode,'stdout_tail':proc.stdout[-18000:],'stderr_tail':proc.stderr[-6000:]}
    report={'timestamp_utc':datetime.now(timezone.utc).isoformat(),'root':str(root),
            'scope':'static code inventory; optional local pytest execution; not a physics certification',
            'scanned_relevant_python_files':rows,'parse_errors':errs,
            'test_files':tests,'test_execution':tests_result,
            'manual_checks_required':['Identify actual current production CLI and trace complete call chain',
               'Confirm coupling basis and parameter inventory per model',
               'Confirm one RG-evolved parameter state feeds mass calculation',
               'Confirm neutrino masses and charged-lepton-corrected PMNS output',
               'Check matching appropriateness for requested scale hierarchy',
               'Run numerical known-limit comparisons']}
    out=root/'output'/'audit'/'t3_pipeline_audit.json'
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,indent=2),encoding='utf-8')
    print('Audit:',out);print('Relevant Python modules:',len(rows));print('Tests found:',len(tests));print('Parse errors:',len(errs))
    if tests_result is not None:print('pytest exit:',tests_result['returncode'])
    print('This report is an inventory, not a completion certificate.')
if __name__=='__main__':main()
