"""Read-only source evidence bundle for the T3 fermion one-loop matching audit.

Collects only relevant A(-4)/D(-2) Matchete source, Laurent audits, and
associated F/S inputs. Does not run Matchete, alter matching, or certify C5.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import zipfile

REQUIRED = (
    'Numerical/parallel_t3/complete/a4_direct_hard_soft.py',
    'Numerical/parallel_t3/complete/hard_llhh_source.py',
    'Numerical/parallel_t3/wilson/general_batch2_laurent.py',
    'Numerical/parallel_t3/wilson/general_batch2_source.py',
    'Numerical/parallel_t3/complete/two_model_hard_general_batch10.py',
    'Numerical/running/GeneralT3DynamicFermionWilson.py',
    'Lagrangian/RunT3GeneralFThreshold.wl',
)
SOURCE_PATTERNS = (
    'output/parallel_t3/full_general_batch2/**/A_-4/**/*',
    'output/parallel_t3/full_general_batch2/**/D_-2/**/*',
    'output/parallel_t3/**/A_-4/**/fermion_matching_input.json',
    'output/parallel_t3/**/D_-2/**/fermion_matching_input.json',
    'output/parallel_t3/**/A_-4/**/scalar_matching_input.json',
    'output/parallel_t3/**/D_-2/**/scalar_matching_input.json',
    'output/parallel_t3/**/A_-4/**/manifest.json',
    'output/parallel_t3/**/D_-2/**/manifest.json',
    'output/eft_sequential_audit/T3_D_alpha_m2/data/*one_loop_inputform.txt',
    'output/full16_symbolic/T3_A_alpha_m4/data/*one_loop_inputform.txt',
    'output/full16_symbolic/T3_D_alpha_m2/data/*one_loop_inputform.txt',
)
# Keep the exact raw matching text and audit files; WXF is binary and needed
# for full Laurent verification, so include it if present.
ALLOW_SUFFIXES = {'.json', '.txt', '.wxf', '.wl', '.py', '.m'}


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo', type=Path, default=Path(__file__).resolve().parents[1])
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--max-file-mib', type=float, default=25.)
    args = p.parse_args()
    root = args.repo.resolve()
    dest = args.output.resolve()
    if not root.is_dir() or dest.exists():
        p.error('Repository must exist and output ZIP must be a new path')
    if args.max_file_mib <= 0:
        p.error('--max-file-mib must be positive')
    selected: dict[str, Path] = {}
    missing = []
    for name in REQUIRED:
        path = root / name
        if path.is_file():
            selected[path.relative_to(root).as_posix()] = path
        else:
            missing.append(name)
    for pattern in SOURCE_PATTERNS:
        for path in root.glob(pattern):
            if path.is_file() and path.suffix.lower() in ALLOW_SUFFIXES:
                selected[path.relative_to(root).as_posix()] = path
    cap = int(args.max_file_mib * 1024**2)
    rows=[]
    for name, path in sorted(selected.items()):
        size = path.stat().st_size
        rows.append({'path':name, 'size_bytes':size,
                     'sha256':hashlib.sha256(path.read_bytes()).hexdigest() if size <= cap else None,
                     'included':size <= cap})
    manifest = {
        'purpose':'Worker A source inspection, not finite matched Wilson',
        'created_utc':datetime.now(timezone.utc).isoformat(),
        'missing_required_sources':missing,
        'files':rows,
        'verified_physical_one_loop_boundary':False,
    }
    dest.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(dest,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for row in rows:
            if row['included']:
                z.write(selected[row['path']],arcname=row['path'])
        z.writestr('WORKER_A_SOURCE_MANIFEST.json', json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({'zip':str(dest),'included_files':sum(x['included'] for x in rows),
                      'skipped_oversize_files':sum(not x['included'] for x in rows),
                      'missing_required_sources':missing},indent=2))

if __name__=='__main__':
    main()
