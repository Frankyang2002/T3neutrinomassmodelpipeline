"""Discover compatible committed source/background pairs, fail closed otherwise."""
from pathlib import Path
import argparse,json,sys,datetime
from Worker2FermionMatchingSensitivity import main as sensitivity_main

def load(p):
    try:return json.loads(p.read_text(encoding='utf-8-sig'))
    except (OSError,ValueError):return {}

def discover(repo,model):
    roots=[repo/'output/general_observable/lo_fermion_boundary',repo/'output/general_observable/lo_joint_scalar_matching']
    sources=[]; backgrounds=[]
    slug={'A_-4':'A_m4','D_-2':'D_m2'}[model]
    for p in roots[0].rglob('boundary.json') if roots[0].exists() else []:
        if p.parent.name==slug and load(p).get('model')==model:sources.append(p)
    for p in roots[1].rglob(slug+'.json') if roots[1].exists() else []:
        if p.parent.name=='numeric_backgrounds' and load(p).get('model')==model:backgrounds.append(p)
    import hashlib
    matches=[]
    for b in sources:
        h=hashlib.sha256(b.read_bytes()).hexdigest()
        for bg in backgrounds:
            meta=load(bg)
            if meta.get('source_boundary_sha256')==h and (bg.parent/meta.get('tensor_npz','missing')).is_file():
                matches.append((b,bg))
    # Stable, latest lexicographic pairing; user can pass explicit pair when preferred.
    return sorted(matches,key=lambda x:str(x[1]),reverse=True)

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--repo',type=Path,default=Path('.'))
    p.add_argument('--model',choices=('A_-4','D_-2','both'),default='both')
    p.add_argument('--assume-diagonal-ye',action='store_true')
    a=p.parse_args(); repo=a.repo.resolve()
    stamp=datetime.datetime.now().strftime('%Y%m%d_%H%M%S_%f')
    models=['A_-4','D_-2'] if a.model=='both' else [a.model]
    status=[]
    for model in models:
        options=discover(repo,model)
        if not options:
            status.append({'model':model,'status':'BLOCKED_NO_MATCHED_SOURCE_BACKGROUND_PAIR'})
            continue
        source,bg=options[0]
        out=repo/'output/audits/worker2_fermion_matching'/stamp/model
        try:
            args=['--repo',str(repo),'--model',model,'--boundary',str(source),
                  '--background',str(bg),'--output',str(out)]
            if a.assume_diagonal_ye:args.append('--assume-diagonal-ye')
            sensitivity_main(args)
            status.append({'model':model,'status':'conditional_success','output':str(out),
                           'source':str(source.relative_to(repo)),'background':str(bg.relative_to(repo))})
        except Exception as exc:
            status.append({'model':model,'status':'BLOCKED','reason':repr(exc)})
    report=repo/'output/audits/worker2_fermion_matching'/stamp/'batch_status.json'
    report.parent.mkdir(parents=True,exist_ok=True)
    report.write_text(json.dumps(status,indent=2)+'\n',encoding='utf-8')
    print(report.read_text())
    if any(r['status']!='conditional_success' for r in status):sys.exit(2)
if __name__=='__main__':main()
