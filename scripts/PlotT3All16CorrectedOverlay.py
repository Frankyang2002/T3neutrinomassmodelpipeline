
#!/usr/bin/env python3
"""All-16 T3 comparison with explicit provenance and corrected D(-2) overlay.

Source: existing Reports/.../frozen_stage_observable_changes.csv, whose
UV/EFT1/SMEFT columns freeze ONE stage in otherwise independently evolved runs.
D(-2) override uses the corrected four-case E0/E1 x S0/S1 benchmark.
Only EFT1 and SMEFT mass observables are comparable as isolated switches;
D(-2) UV and mixing entries remain blank. The EFT1 D(-2) UV input is fixed
by construction, so it is NOT a regenerated full UV-stage counterfactual.

A(-4) is displayed from the original thesis table as LEGACY ONLY, not as
a newly certified corrected result. An optional --a4-scenarios path replaces
A(-4) EFT1/SMEFT masses/splittings if a genuinely matching corrected four-case
CSV exists; it is not fabricated from the approximate A(-4) single-run file.

The plots explicitly annotate provenance and mask unsupported cells.
No source files or existing plots overwritten.
"""
from __future__ import annotations
import argparse,csv,json,math,sys,subprocess
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from scripts.PlotT3OscillationHeatmaps import (
    MODEL_ORDER, MASSSQ_OBSERVABLES, MIXING_OBSERVABLES,
    read_rows, build_matrix, format_cell_value, write_matrix_csv
)

LEGACY=Path('Reports/output/thesis_qualitative/representation_running/csv/frozen_stage_observable_changes.csv')
DC=Path('output/parallel_t3/D_m2_corrected_frozen_running/D_m2_frozen_running_observables.csv')
OUT=Path('output/parallel_t3/all16_Dm2_A4_comparison')
STAGES=('UV','EFT1','SMEFT')
MASS_OBS=(('m_light_eV',r'$m_1$'),('m_middle_eV',r'$m_2$'),('m_heavy_eV',r'$m_3$'))
SPLIT=('delta_m21_sq_ev2','delta_m3l_sq_ev2')
ALL_OBS=list(MASS_OBS)+list(MASSSQ_OBSERVABLES)+list(MIXING_OBSERVABLES)

def load_four(path):
    if not path.is_file():
        raise FileNotFoundError(f'Missing actual corrected 4-scenario CSV: {path}')
    with path.open(newline='',encoding='utf-8-sig') as f:
        reader=csv.DictReader(f)
        required={'scenario','m1_eV','m2_eV','m3_eV',
                  'delta_m21_sq_eV2','delta_m31_sq_eV2'}
        if not required.issubset(reader.fieldnames or ()):
            raise ValueError(f'{path} requires columns {sorted(required)}')
        data={r['scenario']:r for r in reader}
    if set(data)!= {'E0_S0','E0_S1','E1_S0','E1_S1'}:
        raise ValueError(f'{path}: expected precisely E0_S0,E0_S1,E1_S0,E1_S1')
    for row in data.values():
        for k in required-{'scenario'}:
            v=float(row[k])
            if not math.isfinite(v):
                raise ValueError(f'{path}: nonfinite {k}')
    return data

def updates(data,model):
    # Isolate EFT1 effect at fixed SMEFT=ON, and SMEFT effect at EFT1=ON.
    pairs={'EFT1':('E1_S1','E0_S1'),
           'SMEFT':('E1_S1','E1_S0')}
    keys={'m_light_eV':'m1_eV','m_middle_eV':'m2_eV','m_heavy_eV':'m3_eV',
          'delta_m21_sq_ev2':'delta_m21_sq_eV2',
          'delta_m3l_sq_ev2':'delta_m31_sq_eV2'}
    out={}
    for stage,(yes,no) in pairs.items():
        for obs,k in keys.items():
            full=float(data[yes][k]);frozen=float(data[no][k])
            val=100*(full-frozen)/frozen if frozen!=0 else math.nan
            out[(model,stage,obs)]=val
    return out

def plot_matrix(matrix, labels, headings, provenance, filename,title,group_size,caption):
    fig_h=max(8.3, len(labels)*0.41+2.8)
    fig_w=max(11.5, len(headings)*1.08+3.1)
    fig,ax=plt.subplots(figsize=(fig_w,fig_h))
    finite=matrix[np.isfinite(matrix)]
    amplitude=max(1.,float(np.max(np.abs(finite))) if finite.size else 1.)
    masked=np.ma.masked_invalid(matrix)
    palette=plt.get_cmap('RdBu_r').copy()
    palette.set_bad('#eeeeee')
    image=ax.imshow(masked,cmap=palette,
                    norm=colors.TwoSlopeNorm(vmin=-amplitude,vcenter=0,vmax=amplitude),
                    aspect='auto')
    ax.set_xticks(range(len(headings)),headings,fontsize=9)
    ax.set_yticks(range(len(labels)),labels,fontsize=9)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            v=matrix[i,j]
            txt=('—' if not np.isfinite(v) else format_cell_value(float(v)))
            c=('white' if np.isfinite(v) and abs(v)>.55*amplitude else 'black')
            ax.text(j,i,txt,va='center',ha='center',fontsize=7.6,color=c)
    ax.set_xticks(np.arange(-.5,len(headings),1),minor=True)
    ax.set_yticks(np.arange(-.5,len(labels),1),minor=True)
    ax.grid(which='minor',color='#777777',lw=.55)
    ax.tick_params(which='minor',bottom=False,left=False)
    for row in (3,6,9,12):
        ax.axhline(row-.5,color='#333333',lw=1.7)
    for j in range(group_size,len(headings),group_size):
        ax.axvline(j-.5,color='#333333',lw=1.5)
    for i,model in enumerate(MODEL_ORDER):
        if model in ('A_-4','D_-2'):
            ax.axhspan(i-.48,i+.48,facecolor='none',edgecolor='#222222',lw=1.4)
    cbar=fig.colorbar(image,ax=ax,fraction=.024,pad=.02)
    cbar.set_label('Signed % change: (running − frozen) / frozen')
    ax.set_xlabel('Independent stage switch (not four simultaneous stage freezes)')
    ax.set_ylabel('T3 model')
    ax.set_title(title,pad=18,fontsize=12)
    fig.text(.08,.01,caption,ha='left',va='bottom',fontsize=8)
    fig.tight_layout(rect=(0,.06,1,1))
    filename.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(filename,dpi=240,bbox_inches='tight')
    plt.close(fig)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--legacy',type=Path,default=LEGACY)
    p.add_argument('--d-scenarios',type=Path,default=DC)
    p.add_argument('--a4-scenarios',type=Path,default=None,
                   help='Optional newly computed corrected A(-4) E0/E1 x S0/S1 CSV')
    p.add_argument('--output',type=Path,default=OUT)
    p.add_argument('--refresh-d',action='store_true',
                   help='Compute corrected D(-2) four cases first using existing script')
    args=p.parse_args()
    if args.refresh_d:
        command=[sys.executable,str(ROOT/'scripts'/'PlotDMinus2CorrectedFrozenRunning.py')]
        print('Regenerating D(-2) corrected four cases:', ' '.join(command),flush=True)
        subprocess.run(command,cwd=ROOT,check=True)
    legacy=read_rows(args.legacy)
    source={}
    for r in legacy:
        source[(r['model'],r['switch'],r['observable'])]=(
            float(r['relative_change_percent']),'legacy_16_stage_freeze')
    source.update({k:(v,'corrected_Dm2_four_scenario') for k,v in
                   updates(load_four(args.d_scenarios),'D_-2').items()})
    # D UV unavailable with fixed UV boundary; D mixing not measured in
    # corrected scenario script. Must not recycle incompatible legacy values.
    for stage in STAGES:
        for obs,_ in ALL_OBS:
            if stage=='UV' or obs.startswith('sin2_theta'):
                source.pop(('D_-2',stage,obs),None)
    if args.a4_scenarios is not None:
        source.update({k:(v,'corrected_Am4_four_scenario') for k,v in
                       updates(load_four(args.a4_scenarios),'A_-4').items()})
        for stage in STAGES:
            for obs,_ in ALL_OBS:
                if stage=='UV' or obs.startswith('sin2_theta'):
                    source.pop(('A_-4',stage,obs),None)
    args.output.mkdir(parents=True,exist_ok=True)
    provenance=[]
    for observable_defs,stem in (
        (list(MASSSQ_OBSERVABLES),'masssq'),
        (list(MIXING_OBSERVABLES),'mixing'),
        (list(MASS_OBS),'masses')):
        cols=[(stage,key) for stage in STAGES for key,_ in observable_defs]
        headings=[stage+'\n'+label for stage in STAGES for _,label in observable_defs]
        matrix=np.full((16,len(cols)),np.nan)
        for i,model in enumerate(MODEL_ORDER):
            for j,(stage,obs) in enumerate(cols):
                item=source.get((model,stage,obs))
                if item is not None and math.isfinite(item[0]):
                    matrix[i,j]=item[0]
                    status=item[1]
                else:
                    status='missing_noncomparable'
                provenance.append(dict(model=model,stage=stage,observable=obs,
                                       value_percent=matrix[i,j],provenance=status))
        labels=[m+(' [A4 legacy]' if m=='A_-4' and args.a4_scenarios is None
                    else ' [corrected D]' if m=='D_-2'
                    else ' [corrected A]' if m=='A_-4' else '') for m in MODEL_ORDER]
        name=args.output/f'all16_{stem}_stage_switch_provenance.png'
        plot_matrix(matrix,labels,headings,None,name,
                    'T3 all-16 stage-freeze sensitivity — provenance-labelled',
                    len(observable_defs),
                    ('D(-2): corrected EFT1/SMEFT four-case comparison; '
                     'A(-4): '+('corrected overlay' if args.a4_scenarios else 'original legacy result')+
                     '; remaining models: original legacy stage-freeze study. '
                     'Missing entries are NOT zero. Mixed methodology: compare cautiously.'))
        write_matrix_csv(args.output/f'all16_{stem}_matrix.csv',matrix,labels,headings)
        print('Saved:',name)
    with (args.output/'all16_provenance.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=['model','stage','observable','value_percent','provenance'])
        w.writeheader()
        for r in provenance:
            r=dict(r)
            if not math.isfinite(r['value_percent']):r['value_percent']=''
            w.writerow(r)
    metadata={
        'purpose':'all-16 display only, provenance shown per cell',
        'legacy_source':str(args.legacy),
        'corrected_D_source':str(args.d_scenarios),
        'corrected_A_source':str(args.a4_scenarios) if args.a4_scenarios else None,
        'A4_status':('corrected four-case overlay' if args.a4_scenarios else
                     'legacy only; corrected A(-4) four-case data NOT available'),
        'D2_status':'corrected EFT1/SMEFT switches only; UV/mixing intentionally missing',
        'other14_status':'original thesis stage-freeze results, not rerun',
        'methodology_warning':'Mixed pipeline definitions; D E/S scenarios do not regenerate UV and cannot be interpreted as a uniform 16-model refit',
        'negative_values':'possible; signed percentage differences',
        'missing_values':'NaN/blank, never zero'}
    (args.output/'all16_provenance.json').write_text(json.dumps(metadata,indent=2)+'\n',encoding='utf-8')
    print('Saved:',args.output/'all16_provenance.csv')
    print('A(-4) corrected full four-case run:', 'present' if args.a4_scenarios else 'NOT PROVIDED (legacy values labelled)')
    print('No fitted inputs or existing plots changed.')

if __name__=='__main__':
    main()
