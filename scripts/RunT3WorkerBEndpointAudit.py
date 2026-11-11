"""Read-only A_-4/D_-2 Worker B EFT1@muS quadratic endpoint audit.

Uses the repository's validated illustrative running.json reader and the installed
Worker B quadratic Hessian. Does NOT generate missing source HHSS, invent quartics,
or claim full 16-model scalar matching. Writes a fresh report only if --out given.

python -X utf8 -m scripts.RunT3WorkerBEndpointAudit --models A_-4,D_-2
"""
from __future__ import annotations
import argparse
import json
import math
from hashlib import sha256
from pathlib import Path
from typing import Any, Mapping
import numpy as np

from Numerical.running.T3ScalarRunningSourceBridge import (
    source_quadratic_tensor, load_saved_running_state,
)

SUPPORTED = ('A_-4', 'D_-2')


def as_complex(x: Any) -> complex:
    if isinstance(x, Mapping):
        if set(x) != {'real', 'imag'}:
            raise ValueError('complex running coupling requires exactly real and imag')
        x = complex(float(x['real']), float(x['imag']))
    z = complex(x)
    if not np.isfinite([z.real, z.imag]).all():
        raise ValueError('nonfinite running coupling')
    return z


def analytic_spectrum(model: str, state: Mapping[str, Any]) -> np.ndarray:
    """Independent quadratic eigenvalues in fixed real S1,S2 ordering."""
    if model not in SUPPORTED:
        raise ValueError('unsupported model')
    m1 = as_complex(state['mS1Sq'])
    m2 = as_complex(state['mS2Sq'])
    if abs(m1.imag) > 1e-12 or abs(m2.imag) > 1e-12:
        raise ValueError('mass squared must be real')
    if model == 'A_-4':
        return np.sort(np.r_[np.repeat(m1.real, 2), np.repeat(m2.real, 6)])
    if 'muT3GeneralS2Sq' not in state:
        raise ValueError('D_-2 bilinear absent')
    mu = as_complex(state['muT3GeneralS2Sq'])
    return np.sort(np.r_[np.repeat(m1.real, 6), m2.real-2*abs(mu), m2.real+2*abs(mu)])


def audit_state(model: str, mu_s: float, state: Mapping[str, Any],
                *, reference: np.ndarray | None = None) -> dict:
    """Compare independent analytic and fixed-source Hessian diagonalization."""
    if not math.isfinite(mu_s) or mu_s <= 0:
        raise ValueError('invalid scalar matching scale')
    m = source_quadratic_tensor(model, state)
    heavy = m[4:,4:]
    actual, u = np.linalg.eigh(heavy)
    analytic = analytic_spectrum(model, state)
    err = float(np.max(np.abs(actual - analytic)))
    scale = max(1.0, float(np.max(np.abs(analytic))))
    orth = float(np.linalg.norm(u.T@u-np.eye(8)))
    diagonal = float(np.linalg.norm(u.T@heavy@u-np.diag(actual))/scale)
    trace_err = float(abs(np.trace(heavy)-np.sum(analytic))/scale)
    info = {
        'model': model, 'mu_s_gev': float(mu_s),
        'source_mass_eigenvalues_gev2': actual.tolist(),
        'analytic_mass_eigenvalues_gev2': analytic.tolist(),
        'max_abs_analytic_residual_gev2': err,
        'max_rel_analytic_residual': err/scale,
        'orthogonality_residual': orth,
        'diagonalization_relative_residual': diagonal,
        'trace_relative_residual': trace_err,
        'tachyonic': bool(np.min(actual) < 0),
        'source_quadratic': ('mS1Sq,mS2Sq,muT3GeneralS2Sq' if model == 'D_-2'
                             else 'mS1Sq,mS2Sq'),
        'source_quartic_evaluated': False,
        'source_hhss_evaluated': False,
        'scope': 'running quadratic endpoint only; no independent scalar potential certification',
    }
    info['internal_pass'] = bool(max(err/scale, orth, diagonal, trace_err)<1e-9)
    if reference is not None:
        ref = np.sort(np.asarray(reference, float))
        if ref.shape != (8,) or not np.isfinite(ref).all():
            raise ValueError('historical reference must be 8 mass-squared eigenvalues')
        ref_err = float(np.max(np.abs(actual-ref)))
        info['historical_mass_eigenvalues_gev2'] = ref.tolist()
        info['historical_max_abs_residual_gev2'] = ref_err
        info['historical_max_rel_residual'] = ref_err/max(1.,float(np.max(np.abs(ref))))
    return info



def _extract_endpoint(path: Path, model: str) -> tuple[float, dict, str]:
    """Accept ONLY documented running.json shape; no inference from case outputs."""
    doc = json.loads(path.read_text(encoding='utf-8-sig'))
    if not isinstance(doc, dict) or doc.get('model') != model:
        raise ValueError(f'{path}: model mismatch')
    scales = doc.get('scales_gev')
    state = doc.get('eft1_at_scalar_threshold')
    if not isinstance(scales, dict) or not isinstance(state, dict):
        raise ValueError(f'{path}: missing explicit EFT1 scalar-threshold state or scales')
    mu = scales.get('S')
    if not isinstance(mu, (float, int)) or not math.isfinite(mu) or mu <= 0:
        raise ValueError(f'{path}: invalid scalar threshold')
    needed = {'mS1Sq','mS2Sq'} | ({'muT3GeneralS2Sq'} if model == 'D_-2' else set())
    if not needed <= set(state):
        raise ValueError(f'{path}: missing required parameters {sorted(needed-set(state))}')
    return float(mu), state, sha256(path.read_bytes()).hexdigest()


def _discover(project: Path, model: str) -> list[Path]:
    """Read-only search with strict schema checks, not arbitrary JSON substitution."""
    root = project/'output'
    if not root.is_dir():
        return []
    found=[]
    for path in sorted(root.rglob('running.json')):
        try:
            _extract_endpoint(path, model)
        except (ValueError, OSError, TypeError, json.JSONDecodeError):
            continue
        found.append(path)
    return found


def _explicit_sources(items: list[str], project: Path) -> dict[str, Path]:
    result={}
    for item in items:
        if '=' not in item:
            raise ValueError('--state-file expects MODEL=PATH')
        model,path=item.split('=',1)
        if model not in SUPPORTED or model in result or not path.strip():
            raise ValueError('invalid or duplicate --state-file model')
        raw=Path(path)
        result[model]=(raw if raw.is_absolute() else project/raw).resolve()
    return result


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--project', type=Path, default=Path(__file__).resolve().parents[1])
    ap.add_argument('--models', default='A_-4,D_-2')
    ap.add_argument('--out', type=Path, help='New JSON report; existing file never overwritten')
    ap.add_argument('--reference-npz', type=Path,
                    help='Optional NPZ containing A_-4 and D_-2 arrays, each shape (8,) GeV^2')
    ap.add_argument('--find-inputs', action='store_true',
                    help='List schema-compatible saved running.json files; do not calculate')
    ap.add_argument('--state-file', action='append', default=[], metavar='MODEL=PATH',
                    help='Explicit saved running.json per model; may be repeated')
    args=ap.parse_args(argv)
    project=args.project.resolve()
    models=[m.strip() for m in args.models.split(',') if m.strip()]
    if not models or len(models)!=len(set(models)) or not set(models)<=set(SUPPORTED):
        ap.error('models must be unique subset of A_-4,D_-2')
    if args.find_inputs:
        for model in models:
            paths=_discover(project,model)
            print(f'{model}: {len(paths)} schema-compatible running.json candidates')
            for path in paths:
                print(f'  {path.relative_to(project) if path.is_relative_to(project) else path}')
        return 0
    explicit=_explicit_sources(args.state_file,project)
    if not set(explicit)<=set(models):
        ap.error('--state-file contains a model not selected by --models')
    references={}
    if args.reference_npz is not None:
        with np.load(args.reference_npz,allow_pickle=False) as data:
            for model in models:
                if model not in data:
                    raise ValueError(f'historical reference missing {model}')
                references[model]=np.asarray(data[model],float).copy()
    results={}; provenance={}
    for model in models:
        if model in explicit:
            path=explicit[model]
        else:
            candidates=_discover(project,model)
            if len(candidates)!=1:
                print(f'{model}: found {len(candidates)} schema-compatible running.json inputs.',
                      'Run --find-inputs then provide --state-file MODEL=PATH.',flush=True)
                return 3
            path=candidates[0]
        try:
            mu,state,digest=_extract_endpoint(path,model)
        except (FileNotFoundError, ValueError, OSError) as exc:
            print(f'{model}: invalid saved running input: {exc}',flush=True)
            return 3
        results[model]=audit_state(model,mu,state,reference=references.get(model))
        provenance[model]={'path':str(path),'sha256':digest}
    report={'schema':'t3-worker-b-running-quadratic-endpoint-audit-v2',
            'status':'PassedInternalQuadraticChecks' if all(x['internal_pass'] for x in results.values())
                else 'FailedInternalQuadraticChecks',
            'project_root':str(project), 'source_files':provenance,'results':results,
            'historical_spectra_compared':args.reference_npz is not None,
            'physically_complete_scalar_potential_verified':False,
            'joint_scalar_matching_performed':False}
    print(json.dumps(report,indent=2,allow_nan=False))
    if args.out is not None:
        path=args.out.resolve()
        if path.exists(): raise FileExistsError(f'refusing to overwrite {path}')
        path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    return 0 if report['status']=='PassedInternalQuadraticChecks' else 2

if __name__=='__main__':
    raise SystemExit(main())
