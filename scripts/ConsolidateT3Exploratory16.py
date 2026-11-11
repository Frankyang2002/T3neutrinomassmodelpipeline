"""Consolidate the completed exploratory 16-model T3 U/E/W/S benchmark.

Reads existing results only; never changes or re-runs the physics. Exact
source run folders are required so mixed/nonmixed outputs cannot be confused.
No claim of independently verified RGBeta/Matchete normalization is made.

Usage: python -X utf8 -m scripts.ConsolidateT3Exploratory16
       python -X utf8 -m scripts.ConsolidateT3Exploratory16 --out output/my_report
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
STUDIES = ROOT / 'output/studies/general_scalar_extended'
SOURCES = {
    'nonmixed_initial': STUDIES / 'experimental_simultaneous_eft1_nonmixed_fullq' / 'exploratory_nonmixed_fullq_v1',
    'nonmixed_recovered': STUDIES / 'experimental_simultaneous_eft1_nonmixed_fullq' / 'exploratory_nonmixed_fullq_v2',
    'mixed': STUDIES / 'experimental_simultaneous_eft1_mixed_BC' / 'exploratory_mixed_fullq_v3',
}
EXPECTED = (
    'A_+0', 'A_-2', 'A_-4', 'B_+1', 'B_-1', 'B_-3',
    'C_+1', 'C_-1', 'C_-3', 'D_+0', 'D_+2', 'D_-2',
    'E_+0', 'E_+2', 'E_-2', 'E_-4',
)
MIXED = {'B_-1', 'C_-1'}
RECOVERED = {'A_+0', 'A_-2', 'D_+0', 'D_-2', 'E_+0', 'E_-2'}
CASES = {''.join(map(str, (u,e,w,s))) for u in (0,1) for e in (0,1)
         for w in (0,1) for s in (0,1)}


class ConsolidationError(ValueError):
    pass


def _read(path: Path) -> dict:
    if not path.is_file():
        raise ConsolidationError(f'Missing source: {path}')
    try:
        data = json.loads(path.read_text(encoding='utf-8-sig'))
    except (ValueError, OSError) as exc:
        raise ConsolidationError(f'Invalid JSON: {path}: {exc}') from exc
    if not isinstance(data, dict):
        raise ConsolidationError(f'Expected JSON object: {path}')
    return data


def _hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _selected_source(model: str) -> str:
    if model in MIXED:
        return 'mixed'
    if model in RECOVERED:
        return 'nonmixed_recovered'
    return 'nonmixed_initial'


def _number(v: Any) -> float | None:
    if isinstance(v, bool) or not isinstance(v, (int, float)):
        return None
    f = float(v)
    return f if math.isfinite(f) else None


def _metric(case: dict, *keys: str) -> float | None:
    for key in keys:
        value = _number(case.get(key))
        if value is not None:
            return value
    return None


def collect(sources: dict[str, Path] | None = None) -> dict:
    sources = SOURCES if sources is None else sources
    for key in ('mixed', 'nonmixed_recovered', 'nonmixed_initial'):
        if key not in sources:
            raise ConsolidationError(f'Missing source category {key}')
    all_records: dict[str, dict] = {}
    extra = []
    for label, folder in sources.items():
        summary = _read(folder / 'summary.json')
        rows = summary.get('model_rows')
        if not isinstance(rows, list):
            raise ConsolidationError(f'{folder}: summary lacks model_rows')
        reported = {}
        for entry in rows:
            if not isinstance(entry, dict) or not isinstance(entry.get('model'), str):
                raise ConsolidationError(f'{folder}: malformed model_rows entry')
            model = entry['model']
            if model in reported:
                raise ConsolidationError(f'{folder}: duplicate {model} in summary')
            reported[model] = entry
        if label == 'mixed':
            desired = MIXED
        elif label == 'nonmixed_recovered':
            desired = RECOVERED
        elif label == 'nonmixed_initial':
            desired = set(EXPECTED) - MIXED - RECOVERED
        else:
            continue
        for model in sorted(desired):
            meta = reported.get(model)
            if meta is None or meta.get('status') != 'Computed':
                raise ConsolidationError(f'{folder}: {model} not Computed in summary')
        for filepath in folder.glob('*/cases.json'):
            record = _read(filepath)
            name = record.get('model')
            if name not in EXPECTED:
                extra.append(str(filepath))
                continue
            if name not in desired:
                continue
            if name in all_records:
                raise ConsolidationError(f'Model repeated: {name}')
            if record.get('status') != 'Computed':
                raise ConsolidationError(f'{name}: case file not Computed')
            cases = record.get('cases')
            if not isinstance(cases, dict) or set(cases) != CASES:
                raise ConsolidationError(f'{name}: expected 16 U/E/W/S cases, got {sorted(cases) if isinstance(cases,dict) else type(cases)}')
            for code, case in cases.items():
                if not isinstance(case, dict) or case.get('case') != code:
                    raise ConsolidationError(f'{name}/{code}: malformed case record')
                for k, digit in zip(('uv_on','eft1_on','wilson_log_on','smeft_on'), code):
                    if k in case and bool(case[k]) != bool(int(digit)):
                        raise ConsolidationError(f'{name}/{code}: wrong {k}')
            all_records[name] = {
                'model': name,
                'status': 'Computed',
                'source_category': label,
                'source_file': str(filepath.resolve()),
                'source_sha256': _hash(filepath),
                'cases': cases,
            }
    missing = set(EXPECTED) - set(all_records)
    if missing:
        raise ConsolidationError(f'Missing computed case files: {sorted(missing)}')
    rows: list[dict] = []
    for model in EXPECTED:
        record = all_records[model]
        for code in sorted(CASES):
            c = record['cases'][code]
            rows.append({
                'model': model, 'case': code,
                'uv_on': int(code[0]), 'eft1_on': int(code[1]),
                'wilson_log_on': int(code[2]), 'smeft_on': int(code[3]),
                'mixed': int(model in MIXED),
                'source_category': record['source_category'],
                'eft1_nfev': _metric(c,'eft1_nfev','combined_ode_nfev'),
                'simultaneous_vs_fixed_log_rel': _metric(
                    c, 'simultaneous_vs_fixed_log_rel',
                    'simultaneous_vs_original_fixed_log_rel'),
                'simultaneous_vs_frozen_7term_rel': _metric(
                    c,'simultaneous_vs_frozen_7term_rel'),
                'rgbeta_endpoint_residual': _metric(
                    c, 'separate_rgbeta_endpoint_max_block_rel'),
                'c5_log_frobenius': _metric(c,'c5_log_frobenius'),
            })
    records_public = [{k:v for k,v in all_records[m].items() if k!='cases'}
                      for m in EXPECTED]
    report = {
        'format':'t3-exploratory-16-model-consolidation-v1',
        'status':'ConsolidatedExploratoryResults',
        'complete':True,
        'model_count':len(records_public),
        'case_count':len(rows),
        'expected_models':list(EXPECTED),
        'source_runs':{k:str(v.resolve()) for k,v in sources.items()},
        'models':records_public,
        'assumptions': [
            'Matchete extra-quartic tensor to RGBeta mapping uses provisional invariant order',
            'The six recovered ordinary models freeze unsupported complex quadratic betas',
            'Nonmixed newly introduced quartics start at zero at the fermion threshold; their UV running is not included',
            'All fermions integrated out at muF; both scalar multiplets integrated together at muS',
            'Not a certified physical scalar normalization or independent hard one-loop matching',
        ],
        'unexpected_case_files':extra,
    }
    return {'report':report,'rows':rows}


def write_outputs(data: dict, out: Path) -> None:
    out = out.resolve()
    if out.exists():
        raise ConsolidationError(f'Output already exists (no overwrite): {out}')
    out.mkdir(parents=True)
    with (out/'all16_cases.csv').open('w',newline='',encoding='utf-8') as fp:
        w=csv.DictWriter(fp,fieldnames=list(data['rows'][0]))
        w.writeheader()
        w.writerows(data['rows'])
    (out/'all16_manifest.json').write_text(
        json.dumps(data['report'],indent=2,allow_nan=False)+'\n',encoding='utf-8')
    focus=[r for r in data['rows'] if r['case']=='1111']
    md=['# T3 exploratory simultaneous EFT1: 16-model overview','',
        'Status: 16 computed models; 256 U/E/W/S results. Provisional quartic normalization.','',
        '| Model | Joint nfev (1111) | Simultaneous vs fixed LL (1111) | Source |',
        '|---|---:|---:|---|']
    for r in focus:
        num=r['simultaneous_vs_fixed_log_rel']
        md.append(f"| {r['model']} | {r['eft1_nfev'] or '—'} | "
                  f"{num:.4%}" if num is not None else
                  f"| {r['model']} | {r['eft1_nfev'] or '—'} | —")
        md[-1] += f" | {r['source_category']} |"
    md += ['', '**Exploratory qualifications:**'] + [f'- {a}' for a in data['report']['assumptions']]
    (out/'all16_overview.md').write_text('\n'.join(md)+'\n',encoding='utf-8')


def main(argv=None):
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out',type=Path,default=ROOT/'output/studies/general_scalar_extended/exploratory_full_quartic_all16_consolidated_v1')
    args=ap.parse_args(argv)
    data=collect()
    write_outputs(data,args.out)
    print(f"PASS: {data['report']['model_count']} models and {data['report']['case_count']} cases consolidated")
    print(f'Output: {args.out.resolve()}')
    return 0


if __name__=='__main__':
    raise SystemExit(main())
