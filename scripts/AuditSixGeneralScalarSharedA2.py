"""Read-only preflight for six-model general-scalar, shared A(-2) Yukawa study.

Does not run matching, RGEs, neutrino observables, or fits. Refuses to
substitute reduced-scalar inputs/results for general-scalar records.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'output/studies/general_scalar_extended/bulk_48/bulk48_v1'
MODELS = ('A_-2', 'D_-2', 'B_-1', 'C_-1', 'B_-3', 'C_-3')
SCENARIO = 'large_uv_frozen_0p5'
SLUGS = {'A_-2':'A_m2','D_-2':'D_m2','B_-1':'B_m1','C_-1':'C_m1','B_-3':'B_m3','C_-3':'C_m3'}
OUTPUT = ROOT / 'output/studies/general_scalar_extended/six_model_shared_a2_audit'
BASE_STATUS = ('CandidateConstrainedAccepted', 'ExistingConditionalAccepted')


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_matrix(obj: dict, name: str) -> dict:
    if 'base_state' not in obj or 'ordinary' not in obj['base_state']:
        raise ValueError('Missing base_state.ordinary in A(-2) fitted_config')
    matrix = obj['base_state']['ordinary'][name]
    if not isinstance(matrix, dict) or 'real' not in matrix:
        raise ValueError('Missing Yukawa matrix: '+name)
    for channel in ('real', 'imag'):
        rows = matrix.get(channel)
        if not isinstance(rows, list) or len(rows) != 3 or any(not isinstance(row, list) or len(row) != 3 for row in rows):
            raise ValueError(f'{name}.{channel} must be 3x3')
        for row in rows:
            if any(not isinstance(v, (int,float)) for v in row):
                raise ValueError(f'{name}.{channel} must be numeric')
    return matrix


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scenario', choices=('large_uv_frozen_0p5', 'small_uv_frozen_5e-9'), default=SCENARIO)
    parser.add_argument('--output-dir', type=Path, default=OUTPUT)
    args = parser.parse_args()
    src = SOURCE / 'campaign_frozen.csv'
    if not src.is_file():
        raise FileNotFoundError(f'General-scalar campaign CSV missing: {src}')
    with src.open(encoding='utf-8-sig', newline='') as f:
        rows = list(csv.DictReader(f))
    if len(rows) != 48 or len({(r['model'],r['study']) for r in rows}) != 48:
        raise ValueError('Unexpected general-scalar campaign inventory; expecting 48 unique cases')
    reference_path = SOURCE / args.scenario / 'A_m2' / 'fitted_config.json'
    if not reference_path.is_file():
        raise FileNotFoundError(f'General-scalar A(-2) reference missing: {reference_path}')
    reference = json.loads(reference_path.read_text(encoding='utf-8'))
    y1, y2 = read_matrix(reference,'y1'), read_matrix(reference,'y2')
    rep = reference.get('representation', {})
    if (rep.get('d_s1'),rep.get('d_s2'),rep.get('d_f'),rep.get('alpha')) != (1,3,2,-2):
        raise ValueError(f'Unexpected A(-2) representation in {reference_path}')
    selected = {(r['model'],r['study']):r for r in rows}
    report = []
    for model in MODELS:
        r = selected.get((model,args.scenario))
        if r is None: raise ValueError(f'Missing campaign record: {model} {args.scenario}')
        target_path = SOURCE / args.scenario / SLUGS[model] / 'fitted_config.json'
        # A target's *individually fitted* config is not a shared-Yukawa result.
        status = r['status']
        report.append({
            'model':model,'scenario':args.scenario,'campaign_status':status,
            'campaign_matching_scheme':r['c5_source_kind'],
            'candidate_extra_quartic_count':'',
            'target_fitted_config_present':target_path.is_file(),
            'campaign_chi2_individual':r['chi2'],
            'campaign_delta_m21_sq_ev2_individual':r['delta_m21_sq_ev2'],
            'campaign_delta_m3l_sq_ev2_individual':r['delta_m3l_sq_ev2'],
            'shared_A2_C5':'NOT_CALCULATED',
            'shared_A2_masses_and_mixings':'NOT_CALCULATED',
            'explanation':('Existing individually fitted/candidate result; NOT shared-Yukawa prediction' if status in BASE_STATUS
                           else 'No accepted individual endpoint; NOT a shared-Yukawa prediction'),
        })
    out = args.output_dir
    out.mkdir(parents=True,exist_ok=True)
    fields = list(report[0])
    with (out/'preflight.csv').open('w',encoding='utf-8',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=fields)
        writer.writeheader(); writer.writerows(report)
    (out/'reference_A_m2_yukawas.json').write_text(json.dumps({
        'source_kind':'general-scalar fitted_config',
        'scenario':args.scenario,
        'source_path':reference_path.relative_to(ROOT).as_posix(),
        'source_sha256':sha256(reference_path),
        'y1':y1,'y2':y2,
        'note':'Do not treat as a cross-model result; feed into existing general-scalar evaluator with each target model configuration.'
    },indent=2)+'\n',encoding='utf-8')
    (out/'preflight_provenance.json').write_text(json.dumps({
        'scenario':args.scenario,'campaign_csv':src.relative_to(ROOT).as_posix(),
        'campaign_csv_sha256':sha256(src),
        'reference_config_sha256':sha256(reference_path),
        'models':list(MODELS),
        'actual_shared_input_neutrino_evaluations':0,
        'notes':['No reduced-scalar data used','No refit performed','No new matching or RG integration performed',
                 'Failed individual fits do not prove shared-point failure or representation exclusion',
                 'Nonuniform C5 matching schemes must be retained when evaluating comparisons']
    },indent=2)+'\n',encoding='utf-8')
    print('General-scalar six-model shared-A(-2) PRE-FLIGHT only')
    print('scenario:',args.scenario,'reference SHA256:',sha256(reference_path))
    for r in report:
        print(f"{r['model']:6} {r['campaign_status']:31} {r['campaign_matching_scheme']}")
    print('Generated:',out.relative_to(ROOT))
    print('No C5 or neutrino observables calculated. Do not interpret preflight rows as predictions.')

if __name__=='__main__':
    main()
