"""Plot only completed source-resolved corrected LLSS 16-model ablations.

Intended path: scripts/PlotT3CorrectedEFT1Freeze.py
Runs after scripts.RunT3CorrectedEFT1Freeze for all 16 models.
The numbers are exploratory LLSS-sector outputs, not complete one-loop matching.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
SOURCE_BASE=ROOT/'output/studies/general_scalar_extended/corrected_llss_eft1_freeze'
OUT_BASE=ROOT/'Reports/output/thesis_qualitative/corrected_llss_eft1_freeze'
MODELS=('A_-4','A_-2','A_+0','B_-3','B_-1','B_+1',
        'C_-3','C_-1','C_+1','D_-2','D_+0','D_+2',
        'E_-4','E_-2','E_+0','E_+2')
STAGES=(('UV','0111'),('EFT1 (all off)','1001'),('SMEFT','1110'))
MASS=(('delta_m21_sq_ev2',r'$\Delta m_{21}^2$'),
      ('delta_m3l_sq_ev2',r'$\Delta m_{3\ell}^2$'))
ANGLES=(('sin2_theta12',r'$\sin^2\theta_{12}$'),
        ('sin2_theta13',r'$\sin^2\theta_{13}$'),
        ('sin2_theta23',r'$\sin^2\theta_{23}$'))

def load_data(source:Path):
    summary=json.loads((source/'summary.json').read_text(encoding='utf-8'))
    if summary.get('format')!='t3-independent-llss-U_E_W_S-ablation-v1':
        raise ValueError('Incorrect schema: this is NOT the corrected LLSS output')
    if summary.get('one_scale_C5_used_as_W_on_reference') is not False:
        raise ValueError('Anchored one-scale matching detected')
    if summary.get('computed')!=16 or set(summary.get('requested_models',[]))!=set(MODELS):
        raise ValueError('Need all 16 independently evaluated models')
    with (source/'all16_cases.csv').open(encoding='utf-8-sig',newline='') as f:
        rows=list(csv.DictReader(f))
    data={}
    for r in rows:
        if r['status']!='Computed': raise ValueError('Blocked model: '+r['model'])
        k=r['model'],r['case']
        if k in data:raise ValueError('Duplicate model/case: '+str(k))
        data[k]=r
    if len(data)!=256: raise ValueError('Need 16 models × 16 running configurations')
    for m in MODELS:
        for n in range(16):
            if (m,f'{n:04b}') not in data: raise ValueError('Missing case '+m+' '+f'{n:04b}')
    return data,summary

def matrix_and_long(data, obs):
    matrix=np.empty((16,3*len(obs)))
    long=[]
    for i,m in enumerate(MODELS):
        full=data[m,'1111']
        for j,(stage,off) in enumerate(STAGES):
            reference=data[m,off]
            for k,(name,_) in enumerate(obs):
                a=float(full[name]);b=float(reference[name]);
                if not math.isfinite(a) or not math.isfinite(b) or b==0:
                    raise ValueError('Undefined percentage for '+str((m,stage,name)))
                pct=100*(a-b)/b
                matrix[i,j*len(obs)+k]=pct
                long.append({'model':m,'stage_frozen':stage,'on_case':'1111',
                             'off_case':off,'observable':name,'full_value':a,
                             'frozen_value':b,'relative_change_percent':pct,
                             'scope':'independent_LLSS_plus_bubble_not_complete_matching'})
    return matrix,long

def draw(matrix,observables,path,caption):
    # Symmetric logarithmic colour range preserves the sign while avoiding
    # saturation by huge percentages; numeric annotations show exact effects.
    limit=max(1.,float(np.max(abs(matrix))))
    norm=SymLogNorm(linthresh=1.,vmin=-limit,vmax=limit,base=10)
    fig,ax=plt.subplots(figsize=(max(12,2.0*len(observables)*3),9))
    im=ax.imshow(matrix,aspect='auto',cmap='RdBu_r',norm=norm)
    names=[f'{stage}\n{label}' for stage,_ in STAGES for _,label in observables]
    ax.set_xticks(range(len(names)),names)
    ax.set_yticks(range(16),[x.replace('_','(')+')' for x in MODELS])
    ax.set_title(caption+'\nIndependent LLSS-sector calculation (not complete matching)')
    for i in range(16):
        for j in range(len(names)):
            z=matrix[i,j];s=f'{z:+.1f}%' if abs(z)>=1 else f'{z:+.3f}%'
            ax.text(j,i,s,ha='center',va='center',fontsize=7)
    for j in range(1,3):ax.axvline(j*len(observables)-.5,color='black',lw=1.5)
    for i in (3,6,9,12):ax.axhline(i-.5,color='black',lw=1.2)
    fig.colorbar(im,ax=ax,label='100 × (full − stage-off) / stage-off (%)')
    fig.tight_layout();fig.savefig(path,dpi=220,bbox_inches='tight');plt.close(fig)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run-id',required=True)
    p.add_argument('--source',type=Path)
    p.add_argument('--output',type=Path)
    args=p.parse_args()
    source=(args.source or SOURCE_BASE/args.run_id).resolve()
    out=(args.output or OUT_BASE/args.run_id).resolve()
    data,summary=load_data(source)
    out.mkdir(parents=True,exist_ok=True)
    all_rows=[]
    for obs,stem,title in ((MASS,'masssq','Mass-squared differences'),
                           (ANGLES,'mixing','Mixing angles')):
        matrix,rows=matrix_and_long(data,obs)
        all_rows.extend(rows)
        draw(matrix,obs,out/(stem+'_all_eft1_frozen.png'),title)
    path=out/'running_stage_comparisons.csv'
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(all_rows[0]));w.writeheader();w.writerows(all_rows)
    provenance={'source':str(source),'source_csv_sha256':hashlib.sha256((source/'all16_cases.csv').read_bytes()).hexdigest(),
                'cases':{'full':'1111','UV_off':'0111','all_EFT1_off':'1001','SMEFT_off':'1110'},
                'source_matching_scope':'source_resolved_LLSS_scalar_pair_only',
                'full_physical_matching_certified':False,
                'one_scale_anchored_subtraction':False,'rows':len(all_rows)}
    (out/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n',encoding='utf-8')
    print('Saved:',out,'\nAll 16 models × 3 stages × 5 oscillation observables',flush=True)

if __name__=='__main__':main()
