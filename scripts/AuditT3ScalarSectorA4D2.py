"""Read-only A(-4)/D(-2) scalar-sector coverage and convention gate.

Project path: scripts/AuditT3ScalarSectorA4D2.py
No RGE changes, no matching injection, no fitting, no replacement of existing outputs.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import sys
import zipfile

MODELS = ('A_-4', 'D_-2')
EXPECTED = {'A_-4': (0, 0), 'D_-2': (5, 1)}
SCHEMA = 't3-scalar-source-coverage-a4-d2-v1'
REQUIRED_SCRIPTS = (
    'scripts/RunT3SimultaneousEFT1A4.py',
    'scripts/RunT3SimultaneousEFT1All16.py',
    'Numerical/running/GeneralT3SequentialLeadingLog.py',
    'Numerical/parallel_t3/complete/two_model_scalar_coordinate_batch16.py',
)
SOURCE_SEED = 'output/general_matchete_matching/A_m4/data/eft1_after_F_S1_S2_scalar_quartic_seed.json'

def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()

def read_json(path: Path) -> dict:
    x=json.loads(path.read_text(encoding='utf-8-sig'))
    if not isinstance(x,dict):
        raise ValueError(f'Expected JSON object: {path}')
    return x

def sector(counts):
    if len(counts)!=6 or sum(counts)!=4 or any(type(x)!=int or x<0 for x in counts):
        raise ValueError(f'Invalid six-field quartic counts {counts!r}')
    h=sum(counts[:2]); s1=sum(counts[2:4]); s2=sum(counts[4:6])
    return 'HHSS' if h==2 and (s1+s2)==2 else ('H4' if h==4 else 'scalar_only' if h==0 else 'other')

def analyze_model(model, doc):
    if doc.get('model')!=model: raise ValueError(f'{model}: registry identity mismatch')
    rows=list(doc.get('new_rgbeta_symbols',[]))+list(doc.get('preexisting_rgbeta_symbols',[]))
    if len(rows)!=doc.get('total_additional_quartics'): raise ValueError(f'{model}: inconsistent quartic count')
    seen=set(); out=[]
    for item in rows:
        name=item.get('name')
        if not isinstance(name,str) or name in seen: raise ValueError(f'{model}: missing/duplicate symbol')
        seen.add(name)
        counts=item.get('counts')
        typ=sector(counts)
        out.append({'name':name,'sector':typ,'counts':counts,
                    'mapping_status':item.get('mapping_status','Unspecified'),
                    'matchete_basis_name':item.get('matchete_basis_name'),
                    'native_to_rgbeta_normalisation_certified':False})
    q=int(doc.get('quadratics_not_implemented',-1))
    if (len(rows),q)!=EXPECTED[model]:
        raise ValueError(f'{model}: source registry counts changed; manual review needed')
    return {'additional_quartic_count':len(rows),
            'HHSS_directions':[x['name'] for x in out if x['sector']=='HHSS'],
            'scalar_only_directions':[x['name'] for x in out if x['sector']=='scalar_only'],
            'other_directions':[x['name'] for x in out if x['sector']=='other'],
            'quartics':out,'quadratics_not_implemented':q,
            'all_additional_quartic_maps_certified':all(x['mapping_status'] in ('Certified','Verified') for x in out),
            'complete_mass_hessian_certified':False,
            'full_running_lambda_abcd_certified':False,
            'joint_scalar_one_loop_threshold_complete':False,
            'status':'REGISTRY_COVERAGE_VERIFIED_PHYSICAL_SCALAR_MATCHING_OPEN'}

def source_coverage(root):
    out={}
    for name in REQUIRED_SCRIPTS:
        path=root/name
        if not path.is_file(): raise FileNotFoundError(f'Required current source missing: {path}')
        raw=path.read_bytes(); text=raw.decode('utf-8-sig')
        out[name]={'sha256':sha(raw),'bytes':len(raw),
                   'provisional_marker':bool(re.search(r'provisional|conditional|not.*certif',text,re.I))}
    return out

def evidence_coverage(path):
    if path is None: return {'status':'NOT_PROVIDED'}
    if not path.is_file(): raise FileNotFoundError(path)
    with zipfile.ZipFile(path) as z:
        names=set(z.namelist())
        if SOURCE_SEED not in names:
            return {'status':'A4_SCALAR_SEED_NOT_PRESENT','archive_sha256':sha(path.read_bytes())}
        b=z.read(SOURCE_SEED)
        source=json.loads(b)
        fields=source.get('ScalarQuarticTerms',[])
        # This export is after F,S1,S2 have all been integrated out and is
        # therefore NOT the intermediate EFT1 quartic potential.
        return {'status':'DOWNSTREAM_SMEFT_SEED_ONLY_NOT_EFT1_POTENTIAL',
                'archive_sha256':sha(path.read_bytes()),
                'a4_downstream_seed_sha256':sha(b),
                'a4_downstream_seed_quartic_terms':len(fields),
                'integrated_fields':source.get('IntegratedFields'),
                'active_heavy_fields':source.get('ActiveHeavyFields'),
                'not_usable_as_EFT1_complete_scalar_potential':True}

def main(argv=None):
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--project',required=True,type=Path)
    ap.add_argument('--output',required=True,type=Path)
    ap.add_argument('--evidence-archive',type=Path)
    a=ap.parse_args(argv)
    root=a.project.resolve()
    reg=root/'RGE/general/generated/general_rgbeta_quartic_registry.json'
    registry=read_json(reg)
    models={m:analyze_model(m,registry['models'][m]) for m in MODELS}
    src=source_coverage(root)
    ev=evidence_coverage(a.evidence_archive)
    result={'format':SCHEMA,'registry_sha256':sha(reg.read_bytes()),
            'registry_status':registry.get('status'),
            'models':models,'source_files':src,'evidence':ev,
            'physical_scalar_sector_certified':False,
            'next_required_checks':[
                'Independent full invariant potential versus running RGBeta coordinate dictionary',
                'Polynomial reconstruction of symmetric real lambda_abcd for every active coupling',
                'Gauge-generator Ward test and real/imaginary phase orientation',
                'Running mass Hessian, including D(-2) missing quadratic direction',
                'Full joint two-scalar one-loop EFT1-to-SMEFT threshold',
            ],
            'existing_solver_modified':False}
    dst=a.output.resolve()
    if dst.exists(): raise FileExistsError('Refusing to overwrite '+str(dst))
    dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n',encoding='utf-8')
    print('Scalar registry inventory: PASS')
    for m,v in models.items():
        print(f"{m}: extras={v['additional_quartic_count']} HHSS={len(v['HHSS_directions'])} scalar-only={len(v['scalar_only_directions'])} missing-quadratic={v['quadratics_not_implemented']}")
    print('Physical scalar certification: BLOCKED (no source-inconsistent solver change)')
    print('Output:',dst)
    return 0

if __name__=='__main__': sys.exit(main())
