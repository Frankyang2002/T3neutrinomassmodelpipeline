"""Run an existing T3 RGBeta -> C5 -> SM Weinberg benchmark without source edits.

This is a CONDITIONAL one-step sudden-decoupling calculation. It does not
certify the Matchete/RGBeta field map, nor the sequential EFT matching.
Run from the project root:
 python scripts/maintenance/RunBatch2Numerics.py --model A_p0
Optional: --scenario path/to/scenario.json --study path/to/study.json
"""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import traceback

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / 'output' / 'audit' / 'batch2'
SIX = ('A_p0', 'A_m2', 'D_p0', 'D_m2', 'E_p0', 'E_m2')
SCENARIO_FORMAT = 't3-conditional-rgbeta-matchete-bridge-scenario-v1'
STUDY_FORMAT = 't3-extended-c5-conditional-low-energy-v1'


def json_read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def discover(fmt, explicit):
    if explicit:
        path = Path(explicit).resolve()
        if json_read(path).get('format') != fmt:
            raise ValueError(f'{path} is not a {fmt} configuration')
        return path
    # Search configurations, not old run-output reports; prefer user-defined
    # configs over generated artifacts. Never silently manufacture a mapping.
    for root in (ROOT / 'configs', ROOT / 'Numerical' / 'benchmarks', ROOT / 'Reports'):
        if not root.is_dir():
            continue
        for path in sorted(root.rglob('*.json')):
            if 'output' in path.relative_to(ROOT).parts:
                continue
            try:
                if json_read(path).get('format') == fmt:
                    return path
            except (OSError, ValueError, UnicodeError, TypeError):
                continue
    return None


def run_tests():
    paths = [
        'tests/python/test_general_t3_eft1_hard_pole_audit.py',
        'tests/python/test_general_t3_quartic_action_identity.py',
        'tests/python/test_general_t3_triplet_intertwiner.py',
        'tests/python/test_general_t3_triplet_native_structure.py',
        'tests/python/test_general_t3_coupled_uv_running.py',
    ]
    env = os.environ.copy()
    env['PYTHONUTF8'] = '1'  # Addresses cp932 synthetic-fixture failure, NOT physics.
    rows = []
    for rel in paths:
        p = ROOT / rel
        if not p.exists():
            rows.append({'path': rel, 'status': 'missing'})
            continue
        command = [sys.executable, '-m', 'pytest', '-q', rel]
        result = subprocess.run(command, cwd=ROOT, env=env,
                                capture_output=True, text=True, encoding='utf-8', errors='replace')
        log = OUTPUT / (p.stem + '.log')
        log.write_text(result.stdout + '\n' + result.stderr, encoding='utf-8')
        rows.append({'path':rel,'exit_code':result.returncode,'log':str(log.relative_to(ROOT)),
                     'last_lines':result.stdout.splitlines()[-10:]})
    return rows


def run_numerics(model, scenario_path, study_path):
    from Numerical.experimental import T3RGAssistedExtendedC5 as BR
    from Numerical.experimental import T3ExtendedLowEnergyStudy as LE
    scenario = BR.read_scenario(json_read(scenario_path))
    study = LE.load_study(json_read(study_path))
    if not study['allow_unverified_sudden_decoupling']:
        raise ValueError('Study config must explicitly opt into unverified sudden decoupling')
    # Reuse original full coupled RGBeta, exact native C5, and Weinberg solvers.
    result = LE.run_model(model, scenario, study, root=ROOT)
    rows = []
    for scan in result['matching_scale_scan']:
        for phase in scan['phase_scenarios']:
            obs = phase['low_energy_observables']
            rows.append({
                'model': model,
                'mu_match_gev': scan['mu_match_gev'],
                'mu_low_gev': scan['mu_low_gev'],
                'phase_rad': phase['extra_relative_phase_radians'],
                'masses_ev': obs['masses_ev'],
                'delta_m21_sq_ev2': obs['delta_m21_sq_ev2'],
                'delta_m31_sq_ev2': obs['delta_m31_sq_ev2'],
                'pmns_abs': obs['pmns_abs'],
                'takagi_residual': obs['takagi_residual'],
                'c5_linearity_residual': phase['c5_linearity_residual'],
            })
    return result, rows


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--model', choices=SIX, default='A_p0')
    ap.add_argument('--scenario', type=Path)
    ap.add_argument('--study', type=Path)
    ap.add_argument('--skip-tests', action='store_true')
    args = ap.parse_args()
    os.chdir(ROOT)
    OUTPUT.mkdir(parents=True, exist_ok=True)
    report = {
        'model':args.model,
        'calculation_type':'CONDITIONAL: 1-loop RGBeta UV + one-step C5 + SM Weinberg RGE',
        'physical_prediction_certified':False,
        'sequential_eft1_matching_performed':False,
        'tests':[], 'numerical_status':'not_run', 'physics_blockers':[
            'Absolute Matchete/RGBeta coupling and heavy-flavor basis mapping not independently certified',
            'Separate fermion and scalar threshold matching is absent in this selected path'
        ]
    }
    if not args.skip_tests:
        report['tests'] = run_tests()
    try:
        scenario = discover(SCENARIO_FORMAT, args.scenario)
        study = discover(STUDY_FORMAT, args.study)
        report['scenario_file'] = str(scenario.relative_to(ROOT)) if scenario else None
        report['study_file'] = str(study.relative_to(ROOT)) if study else None
        if not scenario or not study:
            missing = [name for name,p in [('scenario',scenario),('study',study)] if not p]
            raise FileNotFoundError('Missing existing '+', '.join(missing)+' JSON; pass --scenario and --study explicitly. No conversion assumptions will be invented.')
        full, rows = run_numerics(args.model,scenario,study)
        (OUTPUT / 'conditional_full_result.json').write_text(json.dumps(full,indent=2,allow_nan=False)+'\n',encoding='utf-8')
        (OUTPUT / 'conditional_observables.json').write_text(json.dumps(rows,indent=2,allow_nan=False)+'\n',encoding='utf-8')
        report['numerical_status'] = 'conditional_result_generated'
        report['points'] = len(rows)
        print('NUMERICAL RESULT:',OUTPUT/'conditional_observables.json')
        for row in rows:
            print(f"{row['model']} mu={row['mu_match_gev']:.6g} GeV: m_i={row['masses_ev']} eV, Dm21^2={row['delta_m21_sq_ev2']:.6g} eV^2")
    except Exception as exc:
        report['numerical_status']='blocked'
        report['error'] = str(exc)
        report['traceback'] = traceback.format_exc()
        print('NUMERICAL BLOCKER:',str(exc),file=sys.stderr)
    out=OUTPUT/'batch2_result.json'
    out.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('REPORT:',out)
    return 0 if report['numerical_status']=='conditional_result_generated' else 2

if __name__=='__main__':
    raise SystemExit(main())
