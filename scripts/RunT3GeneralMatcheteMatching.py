"""Generate source-derived full-general one-scale Matchete C5 for the 16 T3 models.

Self-bilinear A/D/E (-2/0) models are FIRST rotated to two self-conjugate
scalar mass eigenfields, retaining the full extra-Yukawa and quartic basis.
This pilot uses real positive self-bilinears. Complex phase-dependent general
matching, sequential Wilson matching, and physical basis certification are
NOT claimed. B/C(-1) mixed-doublet masses use their separate existing code.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
STANDARD = ROOT / 'Lagrangian' / 'RunModel.wl'
SELF_RUNNER = ROOT / 'Lagrangian' / 'RunT3GeneralSelfMassBasis.wl'
MODELS = {
    'A': (-4, -2, 0),
    'B': (-3, -1, 1),
    'C': (-3, -1, 1),
    'D': (-2, 0, 2),
    'E': (-4, -2, 0, 2),
}
MIXED = frozenset(('B_-1', 'C_-1'))
SELF = frozenset(('A_-2', 'A_+0', 'D_-2', 'D_+0', 'E_-2', 'E_+0'))
ALL = tuple(f'{family}_{alpha:+d}' for family, alphas in MODELS.items()
            for alpha in alphas)
SCHEME_SELF = 'FullGeneralSelfBilinearRealMassEigenfieldsOneScaleCandidate'
SCHEME_NORMAL = 'FullGeneralUnsplitMatcheteOneScaleCandidate'


def tag(model: str) -> str:
    return model.replace('+', 'p').replace('-', 'm')


def source_spec(model: str) -> tuple[Path, str]:
    if model in MIXED:
        raise ValueError('B/C(-1) needs separate mixed-scalar one-scale matching')
    return (SELF_RUNNER, SCHEME_SELF) if model in SELF else (STANDARD, SCHEME_NORMAL)


def _write(path: Path, record: dict) -> None:
    path.write_text(json.dumps(record, indent=2) + '\n', encoding='utf-8')


def run_model(model: str, root: Path, executable: str, *, reuse: bool) -> dict:
    if model not in ALL:
        raise ValueError('Unknown T3 model: ' + model)
    output = Path(root) / tag(model)
    output.mkdir(parents=True, exist_ok=True)
    logfile = output / 'general_matching.log'
    recordfile = output / 'general_matching_status.json'
    if model in MIXED:
        result = {'model': model, 'status': 'BlockedMixedScalarMassBasis',
                  'reason': 'Use the existing B/C(-1) mixed scalar one-scale C5'}
        _write(recordfile, result)
        return result
    runner, scheme = source_spec(model)
    if reuse and recordfile.is_file():
        prev = json.loads(recordfile.read_text(encoding='utf-8'))
        coeff = output / 'c5_coefficient.txt'
        if (prev.get('model') == model and prev.get('status') == 'Matched' and
                prev.get('matching_scheme') == scheme and coeff.is_file() and
                coeff.stat().st_size > 0):
            return prev
    family, alpha = model.split('_')
    alpha_cli = 'm'+alpha[1:] if alpha.startswith('-') else 'p'+alpha[1:]
    candidate = output / 'c5_coefficient.txt'
    candidate.unlink(missing_ok=True)  # never report stale C5 as newly generated
    if model in SELF:
        command = [executable, '-file', str(runner), str(output), family, alpha_cli]
    else:
        command = [executable, '-file', str(runner), str(output),
                   '5', '1', family, alpha_cli]
    env = os.environ.copy()
    env.update({'T3_EXPERIMENTAL_GENERAL_FERMIONS': '1',
                'T3_EXPERIMENTAL_GENERAL_QUADRATICS': '1',
                'T3_EXPERIMENTAL_GENERAL_QUARTICS': '1'})
    try:
        p = subprocess.run(command, cwd=ROOT, env=env,
                           capture_output=True, text=True, check=False)
        logfile.write_text(p.stdout + '\n' + p.stderr,
                           encoding='utf-8', errors='replace')
        if model in SELF:
            model_ok = 'GENERAL SELF-BILINEAR MASS-BASIS C5: SUCCESS' in p.stdout
            metadata = output / 'self_mass_basis_metadata.json'
            try:
                recorded = json.loads(metadata.read_text(encoding='utf-8'))
                metadata_ok = (recorded.get('model') == model and
                               recorded.get('self_bilinear_nonzero_by_mass_splitting') is True and
                               recorded.get('complex_bilinear_phase_supported') is False)
            except (ValueError, OSError, KeyError):
                metadata_ok = False
            model_ok = model_ok and metadata_ok
        else:
            model_ok = f'Model: T3-{family}, alpha = {int(alpha)}' in p.stdout
        success = (p.returncode == 0 and model_ok and candidate.is_file() and
                   candidate.stat().st_size > 0)
        result = {'model': model, 'status': 'Matched' if success else 'Failed',
                  'exit_code': p.returncode, 'model_selection_verified': model_ok,
                  'wolfram_alpha_argument': alpha_cli,
                  'matching_scheme': scheme,
                  'self_bilinear_in_mass_eigenvalues': model in SELF,
                  'complex_bilinear_phase_supported': model not in SELF,
                  'c5_coefficient_file': str(candidate) if candidate.is_file() else None,
                  'log': str(logfile)}
    except OSError as exc:
        result = {'model': model, 'status': 'Failed',
                  'matching_scheme': scheme, 'error': str(exc), 'log': str(logfile)}
    _write(recordfile, result)
    return result


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--model', choices=ALL)
    ap.add_argument('--output-dir', type=Path,
                    default=ROOT / 'output' / 'general_matchete_matching')
    ap.add_argument('--wolframscript', default='wolframscript')
    ap.add_argument('--reuse-success', action='store_true')
    ap.add_argument('--list', action='store_true')
    args = ap.parse_args(argv)
    selected = (args.model,) if args.model else ALL
    if args.list:
        for m in selected:
            print(m, '(mixed-scalar separate)' if m in MIXED else
                  '(self-bilinear mass eigenbasis)' if m in SELF else '')
        return 0
    if not STANDARD.is_file() or not SELF_RUNNER.is_file():
        ap.error('Missing general Matchete source: ' + str(SELF_RUNNER))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for i, model in enumerate(selected, 1):
        result = run_model(model, args.output_dir, args.wolframscript,
                           reuse=args.reuse_success)
        rows.append(result)
        print(f'[{i}/{len(selected)}] {model}: {result["status"]}', flush=True)
        _write(args.output_dir / 'summary.json', {
            'models': rows,
            'matched_count': sum(r['status'] == 'Matched' for r in rows),
            'blocked_count': sum(r['status'].startswith('Blocked') for r in rows),
            'general_numerical_c5_implemented': False,
            'self_bilinear_real_mass_basis_candidate': True})
    return int(any(r['status'] == 'Failed' for r in rows))


if __name__ == '__main__':
    sys.exit(main())
