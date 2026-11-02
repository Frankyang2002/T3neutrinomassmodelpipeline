"""Generic all-16 Matchete -> symmetric polynomial Wilson-coordinate bridge.

Read-only source inputs. Output is not physical matching normalization or EFT1
Wilson running; no gamma matrices are invented.

  python -X utf8 -m scripts.RunT3GeneralWilsonPolynomialBatch --all
"""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

from Numerical.running.GeneralT3WilsonEvolution import (
    MODELS, tag, verify_batch, read_json, save_json, digest,
)
from Numerical.running.GeneralT3WilsonPolynomialBatch import (
    FORMAT, STATUS, build_polynomial_bridge, master_weinberg_polynomial_check,
)

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = Path('output/general_observable/ordered_wilson_tensor_v3')
DEFAULT_OUTPUT = Path('output/general_observable/general_wilson_polynomial_batch_v1')


def absolute(path: Path) -> Path:
    return path.resolve() if path.is_absolute() else (ROOT / path).resolve()


def verified_source_path(model: str, payload: dict, input_root: Path) -> Path:
    directory = input_root / tag(model)
    if model in ('B_-1','C_-1'):
        path = directory / 'source_flavour.json'
        if not path.is_file():
            raise FileNotFoundError(f'{model}: missing native matched source_flavour.json')
    else:
        raw = payload.get('source_boundary')
        if not isinstance(raw, str) or not raw:
            raise ValueError(f'{model}: missing original Matchete boundary path')
        path = absolute(Path(raw))
        if not path.is_file():
            raise FileNotFoundError(f'{model}: original source boundary missing: {path}')
        if digest(path) != payload.get('source_boundary_sha256'):
            raise ValueError(f'{model}: original source boundary SHA-256 mismatch')
    return path


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    selector = ap.add_mutually_exclusive_group(required=True)
    selector.add_argument('--all', action='store_true')
    selector.add_argument('--models', help='Comma-separated subset (for diagnostics)')
    ap.add_argument('--boundary-root', type=Path, default=DEFAULT_INPUT)
    ap.add_argument('--output-root', type=Path, default=DEFAULT_OUTPUT)
    ap.add_argument('--skip-master-check', action='store_true',
                    help='Fixture-only: omit one-generation master Weinberg check')
    args = ap.parse_args(argv)
    if args.all:
        chosen = MODELS
    else:
        chosen = tuple(x.strip() for x in args.models.split(',') if x.strip())
        if not chosen or len(set(chosen)) != len(chosen) or not set(chosen) <= set(MODELS):
            ap.error('Unknown/duplicate model selection')
    src = absolute(args.boundary_root)
    out = absolute(args.output_root)
    if out == src or src in out.parents or out in src.parents:
        ap.error('Output and input must be distinct non-nested directories')
    try:
        checked = verify_batch(src)  # require full verified v3 input even for subsets
        master = None if args.skip_master_check else master_weinberg_polynomial_check()
    except (OSError, ValueError, KeyError, RuntimeError, ImportError) as exc:
        print(f'BLOCKED: {type(exc).__name__}: {exc}', file=sys.stderr)
        return 2
    reports = []
    for model in chosen:
        try:
            original_file = src / tag(model) / 'ordered_tensor.json'
            payload = read_json(original_file)
            source_file = verified_source_path(model, payload, src)
            source = read_json(source_file)
            bridge = build_polynomial_bridge(model, payload, source)
            bridge.update({
                'ordered_boundary_sha256':digest(original_file),
                'operator_source_sha256':digest(source_file),
                'operator_source_file':str(source_file),
            })
            path = out / tag(model) / 'polynomial_bridge.json'
            save_json(path, bridge)
            row = {'model':model, 'status':bridge['status'],
                   'source_PL_term_count':bridge['source_PL_term_count'],
                   'independent_symmetric_coordinate_count':bridge['independent_symmetric_coordinate_count'],
                   'polynomial_relative_error':bridge['polynomial_identity']['relative_monomial_difference'],
                   'flavour_permutation_histogram':bridge['flavour_permutation_histogram'],
                   'ordered_boundary_sha256':checked[model]['tensor_sha256'],
                   'output_file':str(path)}
        except (OSError, ValueError, TypeError, KeyError, IndexError, RuntimeError) as exc:
            row = {'model':model, 'status':'Blocked', 'reason':f'{type(exc).__name__}: {exc}'}
            save_json(out / tag(model) / 'error.json', row)
        reports.append(row)
        print(f'{model}: {row["status"]}', flush=True)
        save_json(out / 'summary.json', {
            'format':FORMAT,
            'status':('PolynomialBatchVerified_NoSpinorNormalization'
                      if len(reports)==len(chosen) and all(r['status']==STATUS for r in reports)
                      else 'PartialOrBlocked'),
            'models_requested':len(chosen),
            'models_verified':sum(r['status']==STATUS for r in reports),
            'models_blocked':sum(r['status']=='Blocked' for r in reports),
            'full_16_model_run':args.all,
            'master_weinberg_internal_polynomial_check':master,
            'matchete_to_master_spinor_normalization_certified':False,
            'physical_eft1_wilson_rge_computed':False,
            'results':reports,
        })
    n = sum(r['status']==STATUS for r in reports)
    print(f'Polynomial coordinate bridges: {n}/{len(chosen)}; '
          'physical Matchete spinor normalization and EFT1 running: NOT CERTIFIED')
    return 0 if n==len(chosen) else 2

if __name__ == '__main__':
    raise SystemExit(main())
