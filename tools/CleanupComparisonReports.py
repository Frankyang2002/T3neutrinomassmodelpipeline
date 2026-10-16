"""Conservative, reversible cleanup of generated comparison report graphics.

Run from the repository root:
    python -m tools.CleanupComparisonReports
    python -m tools.CleanupComparisonReports --apply
    python -m tools.CleanupComparisonReports --restore

Only direct generated figures and dashboards under model/scenario report trees
are candidates. Source numerical data, TeX, PDFs and all other files remain.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / 'Reports' / 'output' / 'full' / 'comparison'
DATA = ROOT / 'output' / 'full' / 'comparison'
QUARANTINE = ROOT / '.cleanup_quarantine' / 'comparison_graphics'
MANIFEST = QUARANTINE / 'manifest.json'
CSV = ROOT / 'comparison_cleanup_manifest.csv'
SCENARIOS = {'smallY_smallL','smallY_largeL','largeY_smallL','largeY_largeL'}
FIGURES = {'uv_coupling_running.png','intermediate_direct_weinberg_running.png',
           'c5_threshold_contributions.png','c5_running.png',
           'neutrino_mass_splitting_running.png','neutrino_mixing_running.png'}

def digest(p: Path) -> str:
    h=hashlib.sha256()
    with p.open('rb') as stream:
        for chunk in iter(lambda:stream.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()

def candidates() -> list[dict]:
    rows=[]
    for scenario in sorted(SCENARIOS):
        result=REPORT/scenario
        source=DATA/scenario
        if not result.is_dir() or not source.is_dir():continue
        successful=[]
        for p in source.rglob('running_diagnostics.json'):
            try:
                if json.loads(p.read_text(encoding='utf-8-sig')).get('status')=='Success':successful.append(p)
            except (OSError,ValueError,TypeError):pass
        if not successful:continue
        # Dashboard is generated from saved numerical diagnostics.
        dashboard=result/'interactive_comparison.html'
        if dashboard.is_file():rows.append({'path':dashboard.relative_to(ROOT).as_posix(),'reason':'dashboard_from_saved_diagnostics'})
        for p in result.rglob('*.png'):
            if p.name not in FIGURES or 'figures' not in p.parts:continue
            # Match exact model key, present in output diagnostics directory structure.
            matched=any(part.startswith('T3_dS1_') and part in p.parts for d in successful for part in d.parts)
            if matched:rows.append({'path':p.relative_to(ROOT).as_posix(),'reason':'selected_figure_from_saved_diagnostics'})
    return sorted(rows,key=lambda x:x['path'])

def write_csv(rows):
    with CSV.open('w',encoding='utf-8',newline='') as stream:
        w=csv.DictWriter(stream,fieldnames=['path','reason']);w.writeheader();w.writerows(rows)

def apply(rows):
    if MANIFEST.exists():raise RuntimeError('Existing quarantine manifest: restore before another apply')
    staged=[]
    QUARANTINE.mkdir(parents=True,exist_ok=True)
    try:
        for row in rows:
            src=ROOT/row['path'];dest=QUARANTINE/row['path'];dest.parent.mkdir(parents=True,exist_ok=True)
            sha=digest(src);shutil.move(str(src),str(dest));staged.append({**row,'sha256':sha})
        MANIFEST.write_text(json.dumps(staged,indent=2)+'\n',encoding='utf-8')
    except Exception:
        for entry in reversed(staged):
            src=QUARANTINE/entry['path'];dest=ROOT/entry['path'];dest.parent.mkdir(parents=True,exist_ok=True)
            if src.exists():shutil.move(str(src),str(dest))
        raise
    print(f'Quarantined {len(staged)} files: {QUARANTINE}')

def restore():
    if not MANIFEST.exists():raise FileNotFoundError(f'No manifest at {MANIFEST}')
    rows=json.loads(MANIFEST.read_text(encoding='utf-8'))
    for row in rows:
        src=QUARANTINE/row['path'];dst=ROOT/row['path']
        if dst.exists():raise FileExistsError(f'Restore would overwrite: {dst}')
        if not src.is_file() or digest(src)!=row['sha256']:raise RuntimeError(f'Quarantine altered: {src}')
    for row in rows:
        src=QUARANTINE/row['path'];dst=ROOT/row['path'];dst.parent.mkdir(parents=True,exist_ok=True)
        shutil.move(str(src),str(dst))
    MANIFEST.unlink()
    print(f'Restored {len(rows)} files')

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    g=p.add_mutually_exclusive_group();g.add_argument('--apply',action='store_true');g.add_argument('--restore',action='store_true')
    args=p.parse_args(argv)
    if args.restore:restore();return 0
    rows=candidates();write_csv(rows)
    print(f'{len(rows)} candidates; manifest preview: {CSV}')
    print('Source JSON/CSV, TeX, PDFs and unknown files are NOT selected.')
    if args.apply:apply(rows)
    else:print('Dry run only. Use --apply after reviewing CSV.')
    return 0
if __name__=='__main__':raise SystemExit(main())
