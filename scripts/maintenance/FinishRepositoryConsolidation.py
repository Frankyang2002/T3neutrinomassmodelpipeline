"""Preview/apply safe root README consolidation and inventory residual clutter.

Default mode never changes source. --apply copies original README into an
external hash-checked archive and installs README_REPLACEMENT.md. Residual files
are reported, NOT automatically removed: dynamic Wolfram/Python dependencies
make filename-only deletion unsafe. The intended project root is inferred from
this script's path, never GitHub.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
import shutil
import sys

ROOT = Path(__file__).resolve().parents[2]
CANDIDATE_DIRS = ('Lagrangian','Numerical','RGE','scripts','tools','docs','tests')
PATTERNS = ('*.backup-*','*.bak','*.orig','*.tmp','*.pyc','*~')

def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=ROOT)
    parser.add_argument('--apply',action='store_true')
    parser.add_argument('--restore',type=Path,help='Archive directory from --apply')
    args=parser.parse_args()
    root=args.root.resolve()
    if not (root/'pipeline.py').is_file() or not (root/'Numerical').is_dir():
        parser.error('Expected T3 repository root')
    archive_parent=root.parent/'T3_repository_archives'
    old=root/'README.md'
    new=Path(__file__).resolve().parents[2]/'README_REPLACEMENT.md'
    if args.restore:
        archive=args.restore.resolve()
        if archive.parent!=archive_parent.resolve():parser.error('Unknown archive location')
        manifest=json.loads((archive/'manifest.json').read_text(encoding='utf-8'))
        saved=archive/'README.md'
        if manifest.get('format')!='t3-readme-consolidation-v1' or digest(saved)!=manifest['old_sha256']:
            raise RuntimeError('Invalid README archive')
        if digest(old)!=manifest['new_sha256']:
            raise RuntimeError('Current README differs from installed version; refusing overwrite')
        shutil.copy2(saved,old)
        print('RESTORED',old)
        return
    if not old.is_file() or not new.is_file():
        raise FileNotFoundError('README.md or README_REPLACEMENT.md missing')
    candidates=[]
    for directory in CANDIDATE_DIRS:
        target=root/directory
        if not target.is_dir():continue
        for pattern in PATTERNS:
            for file in target.rglob(pattern):
                if file.is_file():candidates.append(file.relative_to(root).as_posix())
    # Root files also include obvious temporary logs.
    for pattern in ('*.log','*.bak','*.tmp','*~'):
        candidates += [p.relative_to(root).as_posix() for p in root.glob(pattern) if p.is_file()]
    report={'format':'t3-consolidation-preview-v1',
            'root':str(root),'replace_root_readme':digest(old)!=digest(new),
            'README_original_sha256':digest(old),'README_replacement_sha256':digest(new),
            'residual_candidates_NOT_AUTOMATICALLY_REMOVED':sorted(set(candidates)),
            'next_step':'Review references before archiving individual residual candidates; the script never deletes them.'}
    preview=root/'consolidation_preview.json'
    preview.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print('Preview:',preview)
    print('Residual file candidates (retained):',len(report['residual_candidates_NOT_AUTOMATICALLY_REMOVED']))
    if not args.apply:return
    if old.read_bytes()==new.read_bytes():
        print('README already consolidated; no changes')
        return
    archive_parent.mkdir(parents=True,exist_ok=True)
    stamp=datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S_%f_UTC')
    archive=archive_parent/('readme_consolidation_'+stamp)
    archive.mkdir()
    shutil.copy2(old,archive/'README.md')
    if digest(archive/'README.md')!=digest(old):raise RuntimeError('Archive checksum failure')
    manifest={'format':'t3-readme-consolidation-v1','root':str(root),
              'old_sha256':digest(old),'new_sha256':digest(new)}
    (archive/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    shutil.copy2(new,old)
    if digest(old)!=manifest['new_sha256']:raise RuntimeError('Installed README checksum mismatch')
    print('README UPDATED. Original archived:',archive)
    print('Other source files were untouched. Test your existing pipeline before further cleanup.')

if __name__=='__main__':
    try:main()
    except Exception as e:
        print('ERROR:',e,file=sys.stderr)
        raise
