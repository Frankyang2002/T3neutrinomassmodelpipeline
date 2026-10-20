"""Preview-first consolidation of 16 best T3 fits; never refits or edits physics JSON.

Run: python -m tools.ConsolidateFittedBenchmarks [--apply] [--archive-old]
--apply copies and validates 16 canonical models, refusing conflicts.
--archive-old moves old fit-stage directories to a local quarantine only after verification.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parents[1]
FULL = ROOT / 'Reports' / 'output' / 'full'
STAGES = ('fit_feasibility', 'fit_feasibility_retry', 'fit_feasibility_focused')
CANONICAL = 'fitted_benchmarks'

def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()

def inventory(scenario: str):
    candidates = {}
    for stage in STAGES:
        base = FULL / stage / scenario
        for summary in sorted(base.glob('T3_*/feasibility_summary.json')):
            source = summary.parent
            config = source / 'fitted_candidate_config.json'
            if not config.is_file():
                raise RuntimeError(f'Missing config: {config}')
            obj = json.loads(summary.read_text(encoding='utf-8-sig'))
            chi2 = float(obj['best']['chi2'])
            import math
            if not math.isfinite(chi2):
                continue
            candidates.setdefault(source.name, []).append((chi2, stage, source))
    if len(candidates) != 16:
        raise RuntimeError(f'Expected exactly 16 fitted models; found {len(candidates)}: {sorted(candidates)}')
    return {model: min(opts, key=lambda x: (x[0], STAGES.index(x[1])))
            for model, opts in sorted(candidates.items())}

def verify_destination(scenario, winners):
    root = FULL / CANONICAL / scenario
    manifest = root / 'selection_manifest.json'
    if not manifest.is_file():
        raise RuntimeError(f'Missing consolidated manifest: {manifest}')
    rows = json.loads(manifest.read_text(encoding='utf-8'))
    if set(rows) != set(winners):
        raise RuntimeError('Consolidated manifest does not contain exactly the 16 selected models')
    for model, (chi2, stage, origin) in winners.items():
        dest = root / model
        row = rows[model]
        if row['fit_stage'] != stage or row['chi2'] != chi2:
            raise RuntimeError(f'Manifest selection changed for {model}')
        for filename in ('feasibility_summary.json', 'fitted_candidate_config.json'):
            src, dst = origin / filename, dest / filename
            if not src.is_file() or not dst.is_file() or digest(src) != digest(dst):
                raise RuntimeError(f'Content verification failed: {model}/{filename}')
    return root

def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--scenario', default='smallY_smallL')
    p.add_argument('--apply', action='store_true')
    p.add_argument('--archive-old', action='store_true', help='Requires --apply; archive fit-stage directories after validation')
    args = p.parse_args(argv)
    if args.archive_old and not args.apply:
        p.error('--archive-old requires --apply')
    if not (ROOT / 'pipeline.py').is_file():
        p.error('Run from the correct T3 repository (pipeline.py is missing)')
    winners = inventory(args.scenario)
    destroot = FULL / CANONICAL / args.scenario
    print('Canonical destination:', destroot)
    for model, (chi2, stage, _) in winners.items():
        print(f' {model}: {stage}, chi2={chi2:.9g}')
    print(f'Selected {len(winners)} models; original numerical configs will be copied byte-for-byte.')
    if not args.apply:
        print('PREVIEW: no changes. Use --apply to copy and verify.')
        return 0
    destroot.mkdir(parents=True, exist_ok=True)
    manifest = {}
    for model, (chi2, stage, origin) in winners.items():
        dest = destroot / model
        dest.mkdir(parents=True, exist_ok=True)
        hashes = {}
        for name in ('feasibility_summary.json', 'fitted_candidate_config.json'):
            src, dst = origin / name, dest / name
            if dst.exists() and digest(dst) != digest(src):
                raise RuntimeError(f'Refusing to overwrite changed canonical file: {dst}')
            if not dst.exists():
                shutil.copy2(src, dst)
            if digest(src) != digest(dst):
                raise RuntimeError(f'Copy verification failed: {dst}')
            hashes[name] = digest(dst)
        manifest[model] = {'fit_stage': stage, 'chi2': chi2,
                           'source_folder': str(origin.relative_to(ROOT)).replace('\\', '/'),
                           'sha256': hashes}
    manifest_path = destroot / 'selection_manifest.json'
    if manifest_path.exists():
        old = json.loads(manifest_path.read_text(encoding='utf-8'))
        if old != manifest:
            raise RuntimeError('Existing canonical manifest differs; refusing to overwrite')
    else:
        manifest_path.write_text(json.dumps(manifest, indent=2) + '\n', encoding='utf-8')
    verify_destination(args.scenario, winners)
    print('VERIFIED: 16 selected models, matching SHA256 checksums; all physics JSON unchanged.')
    if args.archive_old:
        archive = ROOT / '.cleanup_quarantine' / ('fitted_stages_' + datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ'))
        for stage in STAGES:
            source = FULL / stage
            if source.exists():
                target = archive / stage
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.move(str(source), str(target))
                print('ARCHIVED', source.relative_to(ROOT), '->', target.relative_to(ROOT))
        print('Old fitting-stage folders moved, not deleted:', archive)
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
