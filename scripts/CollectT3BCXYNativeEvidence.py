"""Make a compact, source-hashed evidence ZIP for Y-only and native X followup.

Read-only collector. Missing optional files are explicitly listed; it does
not certify physical threshold matching or parse Matchete input as code.

python -m scripts.CollectT3BCXYNativeEvidence --all
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from zipfile import ZIP_DEFLATED,ZipFile

ROOT=Path(__file__).resolve().parents[1]
DIR=ROOT/'output/general_t3_mass_eigen_matchete'
YDIR=DIR/'eft1_scalar_Y_only_canonical_probe'
XDIR=DIR/'eft1_scalar_xy_native_followup_v2'
DEFAULT=ROOT/'output/audits/t3_bc_xy_native_followup/xy_native_matching_evidence_v2.zip'
MODELS=('B_-1','C_-1')


def collect(models:tuple[str,...],output:Path=DEFAULT)->dict:
    chosen=[]
    missing=[]
    required=[]
    for model in models:
        if model not in MODELS:raise ValueError('invalid model')
        tag=model.replace('+','p').replace('-','m')
        ystem=f'{tag}_final_y_only_canonical_match'
        xstem=f'{tag}_final_xy_native_match_x'
        cstem=f'{tag}_final_xy_native_check_x'
        # Status and raw first Y-only amplitudes: present from previous success.
        required += [YDIR/(ystem+'.json'),
                     YDIR/(ystem+'.provenance.json'),
                     XDIR/(xstem+'_x_followup.json')]
        for d,stem in [(YDIR,ystem),(XDIR,xstem),(XDIR,cstem)]:
            chosen.extend(sorted(d.glob(stem+'*')))
        for more in (DIR/'eft1_fermion_first_tree').glob(
              f'{tag}_final_eft1_fermion_first_tree_eft1_transition_eom.inputform.txt'):
            chosen.append(more)
        chosen += [DIR/(tag+'_final_parameters.json')]
    for d in (ROOT/'output/studies/general_scalar_extended/continuous_full_scalar/continuous_v1').glob('*/summary.csv'):
        if d.parent.name in ('bc_eft2_sequential_v1','bc_scalar_bubble_candidate_v1'):
            chosen.append(d)
    for p in required:
        if not p.is_file():missing.append(str(p.relative_to(ROOT)))
    if missing:
        raise RuntimeError('Required prior Y/X matching files missing: '+str(missing))
    selected=[];seen=set()
    for p in chosen:
        if p.is_file() and p not in seen:
            selected.append(p);seen.add(p)
    provenance={'format':'t3-bc-xy-evidence-archive-v1',
        'models':list(models),
        'physical_sequential_matching_certified':False,
        'source_equivalence_certified':False,
        'files':{str(p.relative_to(ROOT)).replace('\\','/'):{
            'sha256':sha256(p.read_bytes()).hexdigest(),
            'size_bytes':p.stat().st_size} for p in sorted(selected)}}
    output.parent.mkdir(parents=True,exist_ok=True)
    if output.exists():raise FileExistsError(f'Archive exists; choose a new --output: {output}')
    with ZipFile(output,'w',compression=ZIP_DEFLATED,compresslevel=6) as z:
        for p in sorted(selected):
            z.write(p,arcname=str(p.relative_to(ROOT)).replace('\\','/'))
        z.writestr('MANIFEST.json',json.dumps(provenance,indent=2)+'\n')
    return {'archive':str(output),'file_count':len(selected),
            'size_bytes':output.stat().st_size,
            'physical_sequential_matching_certified':False}


def main()->int:
    arg=argparse.ArgumentParser(description=__doc__)
    m=arg.add_mutually_exclusive_group(required=True)
    m.add_argument('--all',action='store_true')
    m.add_argument('--model',choices=MODELS)
    arg.add_argument('--output',type=Path,default=DEFAULT)
    args=arg.parse_args()
    print(json.dumps(collect(MODELS if args.all else (args.model,),args.output.resolve()),indent=2))
    return 0

if __name__=='__main__':raise SystemExit(main())
