"""Read-only coupled-IVP baseline and A(-4) fermion-boundary audit.

Intended project path: scripts/AuditT3CoupledIVPBaseline.py
Standard library only. Does not import or run numerical solvers.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path

ROOTS = (
    ('A4', 'output/studies/general_scalar_extended/experimental_simultaneous_eft1_seven_term/simultaneous_A_m4_v2'),
    ('nonmixed13', 'output/studies/general_scalar_extended/experimental_simultaneous_eft1_general_nonmixed/simultaneous_nonmixed13_v1'),
    ('mixedBC', 'output/studies/general_scalar_extended/experimental_simultaneous_eft1_mixed_BC/simultaneous_mixed_BC_v1'),
)
REQUIRED = ('1111', '0111', '1001', '1110')
RUNNERS = (
    'scripts/RunT3SimultaneousEFT1A4.py',
    'scripts/RunT3SimultaneousEFT1All16.py',
    'scripts/RunT3SimultaneousEFT1MixedBC.py',
    'scripts/RunT3GeneralSequentialFitAblation.py',
    'Numerical/running/FermionThresholdBoundary.py',
)

def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def finite_check(obj, loc='$', *, depth=0):
    """Find nonfinite JSON numbers, including nonstandard NaN/Infinity JSON literals."""
    if depth > 120:
        return [loc + ': nested depth exceeded']
    if isinstance(obj, float) and not math.isfinite(obj):
        return [loc + ': nonfinite float']
    if isinstance(obj, dict):
        return [e for k,v in obj.items() for e in finite_check(v, f'{loc}.{k}', depth=depth+1)]
    if isinstance(obj, list):
        return [e for i,v in enumerate(obj) for e in finite_check(v, f'{loc}[{i}]', depth=depth+1)]
    return []


def case_map(data):
    if not isinstance(data, dict):
        return None, 'top level not a JSON object'
    cases = data.get('cases')
    if isinstance(cases, dict):
        return cases, 'cases dictionary'
    return None, 'missing/non-dictionary cases'


def audit_one(project: Path, category: str, path: Path):
    entry = dict(group=category, model_dir=path.parent.name, path=str(path.relative_to(project)),
                 exists=path.is_file(), sha256=None, model=None, status=None,
                 case_count=0, missing_required=list(REQUIRED), nonfinite_count=None,
                 nonfinite_examples=[], issues=[], provenance={})
    if not path.is_file():
        entry['issues'].append('MISSING_FILE')
        return entry
    entry['sha256'] = sha256(path)
    try:
        with path.open('r', encoding='utf-8') as f:
            data = json.load(f)
    except (ValueError, UnicodeError, OSError) as exc:
        entry['issues'].append('INVALID_JSON: ' + str(exc))
        return entry
    if not isinstance(data, dict):
        entry['issues'].append('INVALID_TOP_LEVEL')
        return entry
    entry['model'] = data.get('model')
    entry['status'] = data.get('status')
    cases, structure = case_map(data)
    entry['case_structure'] = structure
    if entry['status'] != 'Computed':
        entry['issues'].append('STATUS_NOT_COMPUTED')
    if not isinstance(entry['model'], str):
        entry['issues'].append('MODEL_ID_MISSING')
    if cases is None:
        entry['issues'].append('CASES_NOT_FOUND')
    else:
        entry['case_count'] = len(cases)
        entry['missing_required'] = [code for code in REQUIRED if code not in cases]
        if entry['missing_required']:
            entry['issues'].append('MISSING_REQUIRED_CASES')
        for code in REQUIRED:
            if code in cases and not isinstance(cases[code], dict):
                entry['issues'].append(f'CASE_{code}_NOT_OBJECT')
            elif code in cases:
                record = cases[code]
                if record.get('case', code) != code:
                    entry['issues'].append(f'CASE_{code}_CODE_MISMATCH')
                if record.get('full_physical_one_loop_matching') is True:
                    entry['issues'].append(f'CASE_{code}_CLAIMS_FULL_HARD_MATCHING')
                if record.get('matching_approximation') is None:
                    entry['issues'].append(f'CASE_{code}_MATCHING_LABEL_ABSENT')
                if record.get('combined_ode_nfev') is None:
                    entry['issues'].append(f'CASE_{code}_ODE_NFEV_UNRECORDED')
    bad = finite_check(data)
    entry['nonfinite_count'] = len(bad)
    entry['nonfinite_examples'] = bad[:8]
    if bad:
        entry['issues'].append('NONFINITE_NUMBERS')
    for key in ('format','matching_approximation','full_sequential_general_C5_certified',
                'full_general_C5_matchete_path','full_general_C5_sha256','fermion_tree',
                'matching_kind','fitted_uv_sha256','C5_source_path','C5_source_sha256',
                'running_background','hard_matching_residual_independently_calculated'):
        if key in data:
            value = data[key]
            entry['provenance'][key] = value if isinstance(value, (str, int, float, bool, type(None))) else '<structured value>'
    for k in ('metadata','provenance','inputs','source'):
        if k in data:
            entry['provenance'][k] = '<present; inspect original JSON>'
    return entry


def source_audit(project: Path):
    result = {'files': {}, 'A4_boundary': {}}
    for relative in RUNNERS:
        p = project / relative
        result['files'][relative] = {'exists':p.is_file(), 'sha256':sha256(p) if p.is_file() else None}
    p = project / RUNNERS[0]
    if not p.is_file():
        result['A4_boundary']['status'] = 'runner missing; cannot audit local source'
        return result
    src = p.read_text(encoding='utf-8', errors='replace')
    checks = {
        'boundary_created_by_nonmixed_boundary': 'original.nonmixed_boundary(model,uv_end' in src,
        'boundary_raw_passed_to_coupled_ivp': 'np.asarray(tree.raw,complex)*complex(spinor_factor)' in src,
        'fermion_threshold_eft_state': 'build_intermediate_scalar_boundary(uv_end' in src,
        'matchete_source_path_returned': "'full_general_C5_matchete_path':str(matchete_path)" in src,
        'matchete_source_hash_returned': "'full_general_C5_sha256':original.hash_file(matchete_path)" in src,
        'hard_one_loop_marked_absent': "'full_physical_one_loop_matching':False" in src,
        'new_matching_residual_not_calculated': "data['hard_matching_residual_independently_calculated']=False" in src,
    }
    result['A4_boundary'] = {'status': 'static source-text indicators only; no physical matching validation',
                             'checks': checks,
                             'unresolved': ['GammaCC-to-Weyl sign/symmetry factor',
                                            'renormalised UV-hard minus EFT1-soft finite matching',
                                            'one-loop source operator completeness']}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    project = args.project.expanduser().resolve()
    out = args.output.expanduser().resolve()
    if not project.is_dir():
        parser.error(f'project directory missing: {project}')
    if out == project or not str(out).startswith(str(project) + str(Path('/'))):
        parser.error('output must be a NEW directory inside the project')
    if out.exists():
        parser.error(f'output already exists; refusing overwrite: {out}')
    rows = []
    for group, rel in ROOTS:
        directory = project / rel
        discovered = sorted(directory.glob('*/cases.json')) if directory.is_dir() else []
        for path in discovered:
            rows.append(audit_one(project, group, path))
        if not discovered:
            rows.append({'group':group,'model_dir':'<none>','path':rel+'/*/cases.json',
                         'exists':False,'issues':['NO_CASES_DISCOVERED'], 'missing_required':list(REQUIRED)})
    counts = {group: sum(r['exists'] for r in rows if r['group'] == group) for group,_ in ROOTS}
    model_ids = [r.get('model') for r in rows if r.get('model')]
    if len(model_ids) != len(set(model_ids)):
        for r in rows:
            if r.get('model') in {m for m in model_ids if model_ids.count(m)>1}:
                r['issues'].append('DUPLICATE_MODEL_ID')
    report = {
        'created_utc':datetime.now(timezone.utc).isoformat(),
        'project':str(project),
        'expected_groups':{'A4':1, 'nonmixed13':13, 'mixedBC':2},
        'observed_counts':counts,
        'required_cases':list(REQUIRED),
        'all_16_files_found':counts=={'A4':1,'nonmixed13':13,'mixedBC':2},
        'all_16_required_cases_pass': counts=={'A4':1,'nonmixed13':13,'mixedBC':2} and all(not r['issues'] for r in rows),
        'note':'Static JSON/source audit only; no claimed physical certification or numerical rerun.',
        'models':rows,
        'source':source_audit(project),
    }
    out.mkdir(parents=True)
    (out/'audit.json').write_text(json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False), encoding='utf-8')
    with (out/'models.csv').open('w', encoding='utf-8-sig', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['group','model_dir','model','status','exists','case_count','missing_required','nonfinite_count','issues','path','sha256'])
        writer.writeheader()
        for r in rows:
            writer.writerow({k: (';'.join(v) if isinstance(v,list) else r.get(k,'')) for k,v in [(key,r.get(key,'')) for key in writer.fieldnames]})
    summary = [f'Expected 16 model files; observed groups: {counts}',
               f'All files found: {report["all_16_files_found"]}',
               f'All required-case and static checks pass: {report["all_16_required_cases_pass"]}',
               '', 'Model findings:']
    summary += [f'- {r.get("model") or r["model_dir"]}: {", ".join(r["issues"]) or "static checks passed"}' for r in rows]
    summary += ['', 'A(-4) local fermion-boundary static indicators:',
                json.dumps(report['source']['A4_boundary'], indent=2),
                '', 'NOTE: static validation is not proof of one-loop matching or Wilson normalisation.']
    (out/'SUMMARY.txt').write_text('\n'.join(summary)+'\n', encoding='utf-8')
    print('\n'.join(summary[:5]))
    print(f'Wrote {out / "audit.json"}, {out / "models.csv"}, {out / "SUMMARY.txt"}')
    return 0 if report['all_16_required_cases_pass'] else 2

if __name__ == '__main__':
    raise SystemExit(main())
