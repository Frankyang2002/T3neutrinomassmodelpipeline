#!/usr/bin/env python3
"""Run all 16 native Matchete one-loop F-only threshold drivers (raw evidence).

Repository destination: scripts/RunLOSpinorHardThreshold.py

Never touches existing reference results, Matchete source files or GitHub.
All heavy F integrated together at one mu_F; both scalar multiplets remain.
One fresh Wolfram kernel per model. Raw one-loop sectors are NOT scheme-
renormalized, phase-normalized, or suitable as physical EFT1 Wilson inputs.
"""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--source-release', type=Path,
                        default=Path('output/general_observable/lo_fermion_boundary/releases/source_021dcc1e8fe8ba0f'))
    parser.add_argument('--output', type=Path,
                        default=Path('output/general_observable/lo_certification_spinor'))
    parser.add_argument('--run-id', default='batch01',
                        help='Unique run identifier; existing output files are never overwritten')
    parser.add_argument('--models', nargs='*', default=None,
                        help='Optional restricted native probe; default is all 16')
    parser.add_argument('--wolframscript', default='wolframscript')
    parser.add_argument('--phase-probe', action='store_true',
                        help='Also run existing native Matchete A_p0 abstract CConj probe in its own kernel')
    parser.add_argument('--preflight', action='store_true',
                        help='Check all prerequisites/hashes and print commands; no native runs or writes')
    parser.add_argument('--allow-partial', action='store_true',
                        help='Exit 0 with explicit incomplete batch status if native model(s) fail')
    args = parser.parse_args(argv)
    root = args.repo_root.resolve()
    sys.path.insert(0, str(root))
    from Numerical.lo_certification.spinor.boundary import (
        BoundaryError, MODELS, certify_source_release, checked_write, file_sha, tag)
    from Numerical.lo_certification.spinor.hard_threshold import (
        FORMAT, preflight, model_command, inspect_raw, blocked_record)

    def relative(path):
        return path.resolve() if path.is_absolute() else (root / path).resolve()

    try:
        require_root = (root / 'output/general_observable/lo_certification_spinor').resolve()
        release = relative(args.source_release)
        output = relative(args.output)
        if not output.is_relative_to(require_root):
            raise BoundaryError('Refusing any output outside exclusive spinor worker directory')
        if not re.fullmatch(r'[A-Za-z0-9_-]{1,40}', args.run_id):
            raise BoundaryError('run-id must be 1..40 ASCII alphanumerics, underscores or hyphens')
        selected = list(MODELS) if args.models is None or not args.models else args.models
        if len(selected) != len(set(selected)) or any(m not in MODELS for m in selected):
            raise BoundaryError('Invalid or repeated model in --models')
        manifest = certify_source_release(root, release, None, strict_references=True,
                                         strict_counts=True)
        source_stamp = manifest['source_release_summary_sha256'][:16]
        base = output / f'native_hard_{source_stamp}_{args.run_id}'
        plans = [preflight(root, m) for m in selected]
        if args.phase_probe:
            from Numerical.lo_certification.spinor.hard_threshold import PHASE_DRIVER
            phase_driver = root / PHASE_DRIVER
            if not phase_driver.is_file():
                raise BoundaryError(f'Native abstract GammaCC probe driver missing: {phase_driver}')
        else:
            phase_driver = None
        commands = {
            m: model_command(root, m, base / tag(m), args.wolframscript)
            for m in selected
        }
        common = {
            'format': FORMAT,
            'source_release_summary_sha256': manifest['source_release_summary_sha256'],
            'source_snapshot_sha256': manifest['source_snapshot_sha256'],
            'mu_F_gev': manifest['mu_F_gev'],
            'models_requested': selected,
            'model_plans': plans,
            'single_simultaneous_fermion_threshold': True,
            'both_scalar_multiplets_retained_at_F_threshold': True,
            'SM_fermions_remain_active': True,
            'scalar_threshold_mu_S_performed': False,
            'native_GammaCC_to_master_Weyl_phase': None,
            'physical_master_wilson_boundary': None,
            'hard_one_loop_MSbar_certified': False,
            'complete_leading_one_loop_neutrino_mass_certified': False,
            'no_other_worker_outputs_modified': True,
            'native_phase_probe_requested': args.phase_probe,
        }
        if args.preflight:
            print(json.dumps({**common, 'status': 'PREFLIGHT_ONLY',
                              'commands': [commands[m] for m in selected]},
                             indent=2))
            return 0
        if base.exists():
            raise BoundaryError(f'Immutable run directory already exists: {base}. Use a new --run-id.')
        if shutil.which(args.wolframscript) is None and not Path(args.wolframscript).is_file():
            raise BoundaryError('wolframscript unavailable: cannot execute native Matchete matching')
        # Create a new directory exclusively inside our owned output tree.
        base.mkdir(parents=True, exist_ok=False)
        records = []
        from Numerical.lo_certification.spinor.hard_threshold import (
            native_environment, phase_environment)
        for m in selected:
            place = base / tag(m)
            place.mkdir(exist_ok=False)
            cmd = commands[m]
            print(f'NATIVE F THRESHOLD {m}: starting one-loop Matchete, fresh process',
                  file=sys.stderr, flush=True)
            try:
                # Capturing all output avoids writes to another worker's paths.
                proc = subprocess.run(cmd, cwd=root,
                                      env=native_environment(m, os.environ),
                                      stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                      encoding='utf-8', errors='replace', check=False)
                log = proc.stdout
                code = proc.returncode
            except (OSError, ValueError) as exc:
                log, code = f'Cannot execute native Wolfram: {exc}\n', 127
            # Log is an immutable evidence artifact, even for a failed model.
            (place / 'wolfram.log').open('x', encoding='utf-8').write(log)
            try:
                if code != 0:
                    raise BoundaryError(f'Native executable returned {code}; see wolfram.log')
                rec = inspect_raw(m, place)
                rec['model'] = m
                rec['native_exit_code'] = code
                rec['native_log_sha256'] = file_sha(place / 'wolfram.log')
            except (BoundaryError, OSError, ValueError, KeyError) as exc:
                rec = blocked_record(m, str(exc), code)
                rec['native_log_sha256'] = file_sha(place / 'wolfram.log')
            checked_write(place / 'inspection.json', rec)
            records.append(rec)
            print(f'NATIVE F THRESHOLD {m}: {rec["status"]}',
                  file=sys.stderr, flush=True)
        phase_result = None
        if args.phase_probe:
            from Numerical.lo_certification.spinor.hard_threshold import inspect_native_phase
            pdir = base / 'phase_probe'
            pdir.mkdir(exist_ok=False)
            phase_output = pdir / 'native_A_p0.json'
            # Pre-existing native A_p0 probe requires quadratic flag off.
            # Isolate it from both other drivers in a new kernel.
            phase_env = phase_environment(os.environ)
            try:
                phase_proc = subprocess.run(
                    [args.wolframscript, '-file', str(phase_driver), str(phase_output)],
                    cwd=root, env=phase_env, stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT, encoding='utf-8',
                    errors='replace', check=False)
                phase_code, phase_log = phase_proc.returncode, phase_proc.stdout
            except (OSError, ValueError) as exc:
                phase_code, phase_log = 127, f'Native phase probe failed to start: {exc}\n'
            (pdir / 'wolfram.log').open('x', encoding='utf-8').write(phase_log)
            try:
                if phase_code != 0:
                    raise BoundaryError(f'Native abstract CConj returned {phase_code}')
                phase_result = inspect_native_phase(phase_output)
                phase_result['native_exit_code'] = phase_code
            except (BoundaryError, OSError, ValueError) as exc:
                phase_result = {
                    'status': 'BLOCKED', 'reason': str(exc),
                    'native_exit_code': phase_code,
                    'GammaCC_to_master_Weyl_factor': None,
                    'physical_master_wilson_boundary': None,
                }
            phase_result['native_log_sha256'] = file_sha(pdir/'wolfram.log')
            checked_write(pdir/'inspection.json', phase_result)
        complete = sum(r['status'] != 'BLOCKED' for r in records) == len(selected)
        summary = {**common, 'model_records': records,
                   'native_abstract_phase_probe_result': phase_result,
                   'n_models_with_raw_native_one_loop_evidence': sum(
                       r['status'] != 'BLOCKED' for r in records),
                   'n_models_physically_normalized': 0,
                   'n_models_with_MSbar_hard_sector_certified': 0,
                   'run_directory': str(base.relative_to(root)),
                   'status': ('RAW_F_FIRST_ONE_LOOP_EVIDENCE_COLLECTED_NOT_PHYSICAL'
                              if complete else 'INCOMPLETE_NATIVE_F_THRESHOLD_BATCH')}
        checked_write(base / 'summary.json', summary)
        print(json.dumps(summary, indent=2, sort_keys=True))
        phase_pass = (not args.phase_probe or
                      phase_result['status'] != 'BLOCKED')
        return 0 if (complete and phase_pass) or args.allow_partial else 2
    except (OSError, ValueError, KeyError, AssertionError) as exc:
        print(f'BLOCKED: {type(exc).__name__}: {exc}', file=sys.stderr)
        return 2


if __name__ == '__main__':
    raise SystemExit(main())
