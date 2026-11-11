"""Read-only collection of native scalar-soft EFT1 matching evidence.

Searches saved Matchete InputForm exports for T3EFT1W12 (or equivalent
retained EFT1 Wilson) and emits a fresh evidence ZIP. No physics recomputation.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import re
import zipfile

MARKERS = ('T3EFT1W12', 'EFT1W12')
SOFT_NAMES = ('one_loop_inputform', 'matching_inputform', 'after_S', 'after_scalar')
LIMIT = 12 * 1024 * 1024

def sha(path: Path) -> str:
    h=hashlib.sha256()
    with path.open('rb') as f:
        for part in iter(lambda:f.read(1024*1024), b''):
            h.update(part)
    return h.hexdigest()

def collect(root: Path, destination: Path, max_files: int=40) -> dict:
    root=root.resolve(); destination=destination.resolve()
    if not root.is_dir(): raise FileNotFoundError(root)
    if destination.exists(): raise FileExistsError(destination)
    matches=[]; scanned=0; oversize=0
    # Target specific Matchete text dumps, not arbitrary output, to prevent
    # sweeping up historical results or large numerical trajectory files.
    for p in root.glob('output/**/*.txt'):
        name=p.name.lower()
        if not (('inputform' in name or 'soft' in name) and ('loop' in name or 'scalar' in name or 'after_s' in name)):
            continue
        if not p.is_file():continue
        scanned+=1
        size=p.stat().st_size
        if size>LIMIT:
            oversize+=1;continue
        raw=p.read_text(encoding='utf-8-sig',errors='replace')
        if not any(m in raw for m in MARKERS):continue
        relative=p.relative_to(root).as_posix()
        matches.append({'relative_path':relative,'bytes':size,'sha256':sha(p),
                        'w12_occurrences':sum(raw.count(m) for m in MARKERS),
                        'has_pole':('ϵ' in raw or 'ε' in raw),
                        'has_scalar_mass':('MS1' in raw or 'MS2' in raw),
                        'has_HH':('Field[H,' in raw or 'Field[H, ' in raw),
                        'model_hint': 'A_-4' if re.search(r'A[_-](?:m4|\-4)|T3_A_alpha_m4',relative,re.I) else
                                      'D_-2' if re.search(r'D[_-](?:m2|\-2)|T3_D_alpha_m2',relative,re.I) else 'unspecified'})
    matches.sort(key=lambda x:(x['model_hint']=='unspecified',-x['has_scalar_mass'],-x['has_pole'],x['relative_path']))
    selected=matches[:max_files]
    manifest={'schema':'worker-a-scalar-soft-evidence-v1',
              'searched_root':str(root),'candidate_text_files_scanned':scanned,
              'oversize_skipped':oversize,'matching_source_count':len(matches),
              'selected_count':len(selected),'selected':selected,
              'matching_completed':False,
              'note':'Raw source evidence only; no UV-minus-EFT1 subtraction or Wilson normalisation is certified.'}
    destination.parent.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(destination,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as archive:
        archive.writestr('WORKER_A_SCALAR_SOFT_MANIFEST.json',json.dumps(manifest,indent=2)+'\n')
        for x in selected:archive.write(root/x['relative_path'],x['relative_path'])
    return manifest

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo',type=Path,default=Path('.'))
    p.add_argument('--output',required=True,type=Path)
    p.add_argument('--max-files',type=int,default=40)
    args=p.parse_args()
    if args.max_files<1:raise ValueError('max-files must be positive')
    report=collect(args.repo,args.output,args.max_files)
    print(json.dumps({k:report[k] for k in ('candidate_text_files_scanned','matching_source_count','selected_count')},indent=2))
    print('Wrote',args.output.resolve())
    if not report['selected_count']:
        print('No scalar-soft source was found. Missing evidence must not be interpreted as zero soft contribution.')
if __name__=='__main__':main()
