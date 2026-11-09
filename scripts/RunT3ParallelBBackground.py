"""All-16 Workstream B CLI; read-only on physics source/B1 cache/GitHub.

PowerShell: python -X utf8 scripts/RunT3ParallelBBackground.py --all
"""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from Numerical.parallel_t3.background.contracts import MODELS, DEFAULT_BATCH1
from Numerical.parallel_t3.background.runner import Options, run_all


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    group = ap.add_mutually_exclusive_group(required=True)
    group.add_argument('--all', action='store_true')
    group.add_argument('--model', choices=MODELS)
    ap.add_argument('--run-id', default='B_all16_qualitative_v1')
    ap.add_argument('--batch1', default=DEFAULT_BATCH1)
    ap.add_argument('--mu-uv', type=float, default=1e6)
    ap.add_argument('--mu-f', type=float, default=1e5)
    ap.add_argument('--mu-s', type=float, default=1e3)
    ap.add_argument('--uv-samples', type=int, default=17)
    ap.add_argument('--eft1-samples', type=int, default=33)
    ap.add_argument('--rtol', type=float, default=1e-8)
    ap.add_argument('--atol', type=float, default=1e-11)
    ap.add_argument('--hard-matching-root', default=None,
                    help='Optional, independently computed same-run hard C5 package')
    ap.add_argument('--force', action='store_true',
                    help='Overwrite ONLY this B run-id output files; not symbolic sources')
    a = ap.parse_args(argv)
    options = Options(root=ROOT, run_id=a.run_id, batch1=a.batch1,
                      mu_uv=a.mu_uv, mu_f=a.mu_f, mu_s=a.mu_s,
                      uv_samples=a.uv_samples, eft1_samples=a.eft1_samples,
                      rtol=a.rtol, atol=a.atol,
                      hard_matching_root=a.hard_matching_root, force=a.force)
    result = run_all(options, MODELS if a.all else (a.model,))
    print('Candidate renormalisable runs:',result['background_success_count'],
          '/',len(result['model_results']))
    print('Source-exact formal F-tree runs:',result['formal_f_tree_success_count'],
          '/',len(result['model_results']))
    # A partially blocked batch is a VALID audited output. Nonzero exit when
    # no models could evolve at all (usually environment/provenance failure).
    return 0 if result['background_success_count'] else 2


if __name__ == '__main__':
    raise SystemExit(main())
