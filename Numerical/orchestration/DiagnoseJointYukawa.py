"""Diagnose exceptions hidden by RunJointYukawaClusters' penalty fallback.

Run from the repository root using python -m Numerical.orchestration.DiagnoseJointYukawa.
This script does not change any fitted files.
"""
from __future__ import annotations

import argparse
import json
import traceback
from pathlib import Path
import numpy as np
from Numerical.orchestration import RunJointYukawaClusters as joint
from Numerical.orchestration import RunSmallLambdaRefit as base


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--artifact-map', type=Path)
    parser.add_argument('--cross-probes', type=int, default=3)
    args = parser.parse_args()
    models = joint.load_models()
    mapping = base.load(args.artifact_map) if args.artifact_map else None
    failures = []
    successful_self = 0
    successful_cross = 0
    for i, model in enumerate(models):
        name = model['name']
        print(f'\n[{i+1}/{len(models)}] {name}', flush=True)
        try:
            artifacts = base.locate_artifacts(name, mapping)
            parameters, evaluator = base.get_evaluator(model['config'], artifacts)
        except Exception as exc:
            print('  SETUP FAILURE:', type(exc).__name__, str(exc), flush=True)
            failures.append({'model': name, 'probe': 'setup', 'error': traceback.format_exc()})
            continue
        for source_index in [i] + [j for j in range(len(models)) if j != i][:max(0,args.cross_probes)]:
            source = models[source_index]
            tag = 'SELF' if source_index == i else 'CROSS'
            try:
                result = base.evaluate_point(evaluator, parameters, source['y'])
                residual = np.asarray(result.whitened_residual, dtype=float)
                if residual.shape != (5,) or not np.all(np.isfinite(residual)):
                    raise ValueError(f'Invalid residual: shape={residual.shape}, values={residual}')
                chi2 = float(np.dot(residual, residual))
                print(f'  {tag:5} source={source["name"]}: chi2={chi2:.6g}', flush=True)
                if source_index == i: successful_self += 1
                else: successful_cross += 1
            except Exception as exc:
                err = traceback.format_exc()
                print(f'  {tag:5} source={source["name"]}: {type(exc).__name__}: {exc}', flush=True)
                failures.append({'model': name, 'source': source['name'], 'probe': tag,
                                 'error': err})
    out = base.STUDY / 'joint_clusters' / 'evaluation_diagnostics.json'
    base.save(out, {'successful_self': successful_self, 'successful_cross': successful_cross,
                    'failures': failures})
    print(f'\nSuccessful self tests: {successful_self}/{len(models)}; cross tests: {successful_cross}', flush=True)
    print(f'Full exception traces saved: {out}', flush=True)
    if failures:
        print('First exception traceback:\n' + failures[0]['error'], flush=True)


if __name__ == '__main__':
    main()
