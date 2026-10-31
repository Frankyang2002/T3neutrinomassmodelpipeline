"""Collect available Y/X Matchete evidence, including incomplete stage failures.

No requirement that both X matching runs succeed. This is a diagnostic
archive, never a physical matching certificate.

python -m scripts.CollectT3BCXFreshEvidence --all
"""
from __future__ import annotations
import argparse
from hashlib import sha256
import json
from pathlib import Path
from zipfile import ZipFile,ZIP_DEFLATED

ROOT=Path(__file__).resolve().parents[1]
MAT=ROOT/'output/general_t3_mass_eigen_matchete'
V2=MAT/'eft1_scalar_xy_native_followup_v2'
V3=MAT/'eft1_scalar_xy_fresh_kernel_v3'
Y=MAT/'eft1_scalar_Y_only_canonical_probe'
DEFAULT=ROOT/'output/audits/t3_bc_xy_native_followup/xy_fresh_native_evidence_v3.zip'
MODELS=('B_-1','C_-1')


def collect(models:tuple[str,...],destination:Path=DEFAULT)->dict:
    if not models or any(model not in MODELS for model in models):
        raise ValueError('invalid model list')
    paths=[];status={};missing=[]
    for model in models:
        tag=model.replace('+','p').replace('-','m')
        stems=[(Y,f'{tag}_final_y_only_canonical_match'),
               (V2,f'{tag}_final_xy_native_diagnose_x'),
               (V2,f'{tag}_final_xy_native_match_x'),
               (V3,f'{tag}_final_x_from_saved_y_check'),
               (V3,f'{tag}_final_x_from_saved_y_match')]
        for directory,stem in stems:paths.extend(directory.glob(stem+'*'))
        params=MAT/f'{tag}_final_parameters.json'
        paths.append(params)
        paths.extend((MAT/'eft1_fermion_first_tree').glob(
             f'{tag}_final_eft1_fermion_first_tree_eft1_transition_eom.inputform.txt'))
        xpaths=[V3/f'{tag}_final_x_from_saved_y_match.provenance.json',
                V3/f'{tag}_final_x_from_saved_y_check.provenance.json',
                V2/f'{tag}_final_xy_native_diagnose_x.provenance.json']
        summary={}
        for path in xpaths:
            if path.is_file():
                try:
                    info=json.loads(path.read_text(encoding='utf-8'))
                    summary[path.parent.name+'/'+path.name]={
                        'status':info.get('status','unknown'),
                        'fresh_x_matchete_matching_completed':info.get('fresh_x_matchete_matching_completed',False),
                        'physical_matching_certified':False}
                except (ValueError,OSError):summary[path.name]={'unreadable':True}
            else:missing.append(str(path.relative_to(ROOT)).replace('\\','/'))
        status[model]=summary
    for p in [ROOT/'Lagrangian/RunT3BCXFromSavedY.experimental.wl',
              ROOT/'scripts/RunT3BCXFreshKernel.py']:
        if p.is_file():paths.append(p)
    included=sorted(set(p for p in paths if p.is_file()),key=lambda x:str(x))
    if not included:
        raise RuntimeError('No available evidence files')
    files={str(p.relative_to(ROOT)).replace('\\','/'): {
        'sha256':sha256(p.read_bytes()).hexdigest(), 'size_bytes':p.stat().st_size}
        for p in included}
    manifest={'format':'t3-bc-partial-xy-native-evidence-v3',
        'models':list(models), 'status_by_model':status,
        'unavailable_reports':missing,
        'y_one_loop_inherited_in_x_source':False,
        'f_first_eomeft_equivalence_certified':False,
        'running_boundary_substituted':False,
        'physical_sequential_matching_certified':False,
        'files':files}
    destination.parent.mkdir(parents=True,exist_ok=True)
    if destination.exists():raise FileExistsError(f'evidence archive exists: {destination}')
    with ZipFile(destination,'w',compression=ZIP_DEFLATED,compresslevel=6) as z:
        for p in included:z.write(p,arcname=str(p.relative_to(ROOT)).replace('\\','/'))
        z.writestr('MANIFEST.json',json.dumps(manifest,indent=2)+'\n')
    return {'archive':str(destination),'file_count':len(included),
        'missing_reports':len(missing),
        'physical_sequential_matching_certified':False}


def main()->int:
    p=argparse.ArgumentParser(description=__doc__)
    g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--all',action='store_true')
    g.add_argument('--model',choices=MODELS)
    p.add_argument('--output',type=Path,default=DEFAULT)
    a=p.parse_args()
    print(json.dumps(collect(MODELS if a.all else (a.model,),a.output.resolve()),indent=2))
    return 0

if __name__=='__main__':raise SystemExit(main())
