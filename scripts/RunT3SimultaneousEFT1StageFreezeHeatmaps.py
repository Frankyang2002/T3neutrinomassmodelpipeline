"""Original thesis heatmaps for the 16-model simultaneous EFT1 stage-freeze study.

Project path: scripts/RunT3SimultaneousEFT1StageFreezeHeatmaps.py

Reads the existing fitted 16-case outputs per model, WITHOUT rerunning any RGEs.
Reuses the ORIGINAL scripts/PlotT3OscillationHeatmaps.py data parsing and matrices.
Renders with signed symmetric-logarithmic colours; the original plotter is unchanged.

Bit order: U E W S = UV renormalisable, EFT1 renormalisable,
                     EFT1 dimension-5 Wilson, SMEFT SM+Weinberg.
    FULL       1111
    UV off     0111   UV renormalisable state frozen; others active
    EFT1 off   1001   EFT1 renormalisable AND Wilson evolution frozen
    SMEFT off  1110   SM AND Weinberg evolution frozen

The finite joint scalar matching at muS remains present in every case.
This is a numerical exploratory comparison, not certified hard one-loop matching.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
STUDY = ROOT / 'output' / 'studies' / 'general_scalar_extended'
DEFAULT_A4 = STUDY / 'experimental_simultaneous_eft1_seven_term' / 'simultaneous_A_m4_v2'
DEFAULT_NONMIXED = STUDY / 'experimental_simultaneous_eft1_general_nonmixed' / 'simultaneous_nonmixed13_v1'
DEFAULT_MIXED = STUDY / 'experimental_simultaneous_eft1_mixed_BC' / 'simultaneous_mixed_BC_v1'
DEFAULT_OUTPUT = STUDY / 'reports' / 'simultaneous_16_stagefreeze_original_heatmaps_v1'
ORIGINAL_PLOTTER = ROOT / 'scripts' / 'PlotT3OscillationHeatmaps.py'
MODELS = (
    'A_-4','A_-2','A_+0',
    'B_-3','B_-1','B_+1',
    'C_-3','C_-1','C_+1',
    'D_-2','D_+0','D_+2',
    'E_-4','E_-2','E_+0','E_+2',
)
FULL = '1111'
STAGES = {'UV':'0111', 'EFT1':'1001', 'SMEFT':'1110'}
OBSERVABLES = (
    'delta_m21_sq_ev2',
    'delta_m3l_sq_ev2',
    'sin2_theta12',
    'sin2_theta13',
    'sin2_theta23',
)
LONG_FIELDS = ('model', 'switch', 'observable', 'relative_change_percent',
               'comparison_kind')
WIDE_FIELDS = ('model','stage_frozen','full_case','frozen_case','status','error',
               *[f'{o}_pct' for o in OBSERVABLES])


def tag(model: str) -> str:
    return model.replace('+', 'p').replace('-', 'm')


def read_json(path: Path) -> Any:
    with path.open('r', encoding='utf-8-sig') as f:
        return json.load(f)


def write_json(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, allow_nan=False)+'\n',encoding='utf-8')


def _read_case_file(directory: Path, model: str) -> dict:
    path = directory / tag(model) / 'cases.json'
    if not path.is_file():
        raise FileNotFoundError(f'Missing saved model cases: {path}')
    payload = read_json(path)
    if payload.get('model') != model or payload.get('status') != 'Computed':
        raise ValueError(f'{path}: not a computed case record for {model}')
    cases = payload.get('cases')
    if not isinstance(cases,dict) or len(cases) != 16:
        raise ValueError(f'{path}: expected 16 U/E/W/S cases')
    return cases


def _observable(case: dict, name: str) -> float:
    prediction = case.get('prediction')
    if not isinstance(prediction, dict):
        raise ValueError('Case lacks prediction dictionary')
    if name not in prediction:
        raise ValueError(f'Prediction lacks {name}')
    val = float(prediction[name])
    if not math.isfinite(val):
        raise ValueError(f'Nonfinite prediction {name}')
    return val


def _relative_percent(full: float, frozen: float) -> float:
    if frozen == 0.0:
        raise ValueError('Frozen observable is exactly zero; signed percent undefined')
    diff = 100.0 * (full - frozen) / frozen
    if not math.isfinite(diff):
        raise ValueError('Nonfinite signed percent')
    return diff


def _records_for_model(model: str, cases: dict) -> tuple[list[dict], list[dict], list[dict]]:
    """Produce original-plotter input and detailed diagnostics per model.

    The existing plotting script blanks omitted cells. A failed observable is
    omitted, never filled with 0. The other observables in the row survive.
    """
    if FULL not in cases:
        raise ValueError(f'{model}: missing fully running baseline {FULL}')
    for stage, code in STAGES.items():
        if code not in cases:
            raise ValueError(f'{model}: missing {stage} frozen counterpart {code}')
    long_rows, wide_rows, errors = [], [], []
    for stage, frozen_code in STAGES.items():
        result = {'model':model,'stage_frozen':stage,'full_case':FULL,
                  'frozen_case':frozen_code,'status':'VALID_EXPLORATORY','error':''}
        for obs in OBSERVABLES:
            try:
                running = _observable(cases[FULL], obs)
                frozen = _observable(cases[frozen_code], obs)
                pct = _relative_percent(running, frozen)
            except (ValueError,TypeError,KeyError) as ex:
                result[f'{obs}_pct'] = ''
                errors.append({'model':model,'stage':stage,'observable':obs,'reason':str(ex)})
            else:
                result[f'{obs}_pct'] = repr(pct)
                long_rows.append({
                    'model':model,
                    'switch':stage,
                    'observable':obs,
                    'relative_change_percent':repr(pct),
                    'comparison_kind':'independently_evolved_running_stage',
                })
        if all(result[f'{o}_pct'] == '' for o in OBSERVABLES):
            result['status']='NO_VALID_OBSERVABLES'
            result['error']='No finite nonzero-denominator observables'
        elif any(result[f'{o}_pct'] == '' for o in OBSERVABLES):
            result['status']='PARTIAL_OBSERVABLES'
            result['error']='See missing_cells.json'
        wide_rows.append(result)
    return long_rows,wide_rows,errors


def write_csv(path: Path, rows: list[dict], fieldnames: tuple[str, ...]) -> None:
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def collect(a4: Path, nonmixed: Path, mixed: Path) -> tuple[list[dict],list[dict],list[dict]]:
    all_long,all_wide,errors=[],[],[]
    sources={**{m:nonmixed for m in MODELS if m not in ('A_-4','B_-1','C_-1')},
             'A_-4':a4,'B_-1':mixed,'C_-1':mixed}
    for model in MODELS:
        cases=_read_case_file(sources[model],model)
        long_rows,wide_rows,model_errors=_records_for_model(model,cases)
        all_long.extend(long_rows)
        all_wide.extend(wide_rows)
        errors.extend(model_errors)
    return all_long,all_wide,errors



def _signed_log_heatmap(matrix, row_labels, column_labels, output_path,
                        title, stage_group_size, *, linthresh=1.0):
    """Original heatmap layout; signed logarithmic colour transform.

    Positive and negative decades are symmetric. The interval [-1%,1%]
    is linear, and missing results are masked rather than painted as zero.
    """
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import SymLogNorm
    from matplotlib.ticker import FixedLocator, FuncFormatter

    a=np.asarray(matrix,float)
    valid=np.isfinite(a)
    if not valid.any():
        raise ValueError(f'No finite cells: {output_path}')
    vmax=max(float(np.max(np.abs(a[valid]))),linthresh)
    norm=SymLogNorm(linthresh=linthresh,linscale=1.,vmin=-vmax,
                    vmax=vmax,base=10)
    cmap=plt.get_cmap('RdBu_r').copy()
    cmap.set_bad('white')
    nrows,ncols=a.shape
    fig,ax=plt.subplots(figsize=(max(8.5,.9*ncols+4),
                                  max(8.,.38*nrows+2.8)))
    img=ax.imshow(np.ma.masked_invalid(a),cmap=cmap,norm=norm,aspect='auto')
    ax.set_title(title,pad=16)
    ax.set_xticks(np.arange(ncols))
    ax.set_xticklabels(column_labels)
    ax.set_yticks(np.arange(nrows))
    ax.set_yticklabels(row_labels)
    ax.set_xlabel('Running stage and oscillation observable')
    ax.set_ylabel('T3 model')
    for i in range(nrows):
        for j in range(ncols):
            value=a[i,j]
            if not np.isfinite(value):
                continue
            magnitude=abs(value)
            if magnitude>=10000:
                label=f'{value:+.1e}'
            elif magnitude>=10:
                label=f'{value:+.1f}'
            elif magnitude>=1:
                label=f'{value:+.2f}'
            elif magnitude>=.1:
                label=f'{value:+.3f}'
            elif magnitude>=.01:
                label=f'{value:+.4f}'
            else:
                label=f'{value:+.2e}'
            r,g,b,_=cmap(norm(value))
            luminance=.2126*r+.7152*g+.0722*b
            ax.text(j,i,label,ha='center',va='center',fontsize=8.0,
                    color='white' if luminance<.48 else 'black')
    ax.set_xticks(np.arange(-.5,ncols,1),minor=True)
    ax.set_yticks(np.arange(-.5,nrows,1),minor=True)
    ax.grid(which='minor',color='white',linewidth=.6)
    ax.tick_params(which='minor',bottom=False,left=False)
    for y in (3,6,9,12):
        ax.axhline(y-.5,color='black',linewidth=1.6)
    for x in range(stage_group_size,ncols,stage_group_size):
        ax.axvline(x-.5,color='black',linewidth=1.6)
    for i,stage in enumerate(('UV off','EFT1 off','SMEFT off')):
        ax.text((i+.5)*stage_group_size-.5,-1,stage,
                ha='center',va='center',weight='bold',fontsize=10)
    top=int(math.floor(math.log10(vmax))) if vmax>=10 else 0
    decades=[10.**i for i in range(top+1)]
    ticks=sorted(set([-d for d in decades]+[0.]+decades))
    ticks=[t for t in ticks if abs(t)<=vmax]
    cbar=fig.colorbar(img,ax=ax,fraction=.035,pad=.03)
    cbar.locator=FixedLocator(ticks)
    cbar.formatter=FuncFormatter(lambda v,pos:
       '0' if v==0 else (f'{v:.0e}' if abs(v)>=1e4 else f'{v:g}'))
    cbar.update_ticks()
    cbar.set_label('Signed change (%) — symmetric log colours')
    output_path.parent.mkdir(parents=True,exist_ok=True)
    fig.tight_layout()
    fig.savefig(output_path,dpi=300,bbox_inches='tight')
    plt.close(fig)


def render_original_heatmap_layout(input_csv: Path, out: Path) -> None:
    """Reuse the original GitHub plotter's input schema, matrices and labels."""
    from scripts import PlotT3OscillationHeatmaps as original
    rows=original.read_rows(input_csv)
    mass,models,mass_labels=original.build_matrix(rows,original.MASSSQ_OBSERVABLES)
    mixing,_,mix_labels=original.build_matrix(rows,original.MIXING_OBSERVABLES)
    original.write_matrix_csv(out/'csv'/'oscillation_heatmap_masssq_table.csv',
                              mass,models,mass_labels)
    original.write_matrix_csv(out/'csv'/'oscillation_heatmap_mixing_table.csv',
                              mixing,models,mix_labels)
    _signed_log_heatmap(mass,models,mass_labels,
                        out/'plots'/'oscillation_masssq_heatmap_percent.png',
                        'Signed percentage changes in mass-squared differences',
                        len(original.MASSSQ_OBSERVABLES))
    _signed_log_heatmap(mixing,models,mix_labels,
                        out/'plots'/'oscillation_mixing_heatmap_percent.png',
                        'Signed percentage changes in mixing-angle observables',
                        len(original.MIXING_OBSERVABLES))


def self_test() -> None:
    sample={key:{'prediction':{name:float(val) for name,val in zip(OBSERVABLES,vals)}}
            for key, vals in {
                '1111': (4.,6.,0.30,0.02,0.55),
                '0111': (2.,3.,0.25,0.01,0.50),
                '1001': (1.,2.,0.20,0.01,0.50),
                '1110': (8.,12.,0.40,0.04,0.60),
            }.items()}
    long_rows,wide_rows,errors=_records_for_model('A_-4',sample)
    assert len(long_rows)==15 and len(wide_rows)==3 and not errors
    get=lambda stage,obs:next(float(r['relative_change_percent']) for r in long_rows
                               if r['switch']==stage and r['observable']==obs)
    assert math.isclose(get('UV','delta_m21_sq_ev2'),100.)
    assert math.isclose(get('EFT1','delta_m21_sq_ev2'),300.)
    assert math.isclose(get('SMEFT','delta_m21_sq_ev2'),-50.)
    # Missing observables do not become zeros.
    sample['1001']['prediction']['delta_m21_sq_ev2']=0.0
    new_rows,_,errors=_records_for_model('A_-4',sample)
    assert len(new_rows)==14 and len(errors)==1
    assert not any(r['switch']=='EFT1' and r['observable']=='delta_m21_sq_ev2'
                   for r in new_rows)
    from matplotlib.colors import SymLogNorm
    norm=SymLogNorm(linthresh=1.0,vmin=-1e5,vmax=1e5,base=10)
    assert math.isclose(float(norm(-100.)),1.-float(norm(100.)),abs_tol=1e-12)
    assert math.isclose(float(norm(0.)),.5,abs_tol=1e-12)
    print('PASS: exact stage case IDs, signed percentage and NaN/blank logic')
    print('PASS: symmetric logarithmic colour mapping')


def main(argv: list[str] | None = None) -> int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--a4',type=Path,default=DEFAULT_A4)
    p.add_argument('--nonmixed',type=Path,default=DEFAULT_NONMIXED)
    p.add_argument('--mixed',type=Path,default=DEFAULT_MIXED)
    p.add_argument('--output',type=Path,default=DEFAULT_OUTPUT)
    p.add_argument('--self-test',action='store_true')
    p.add_argument('--overwrite',action='store_true',
                   help='Explicitly overwrite only files within the selected report directory')
    args=p.parse_args(argv)
    if args.self_test:
        self_test()
        return 0
    if not ORIGINAL_PLOTTER.is_file():
        p.error(f'Existing original heatmap plotter missing: {ORIGINAL_PLOTTER}')
    for group in (args.a4,args.nonmixed,args.mixed):
        if not group.is_dir():
            p.error(f'Missing source run directory: {group}')
    out=args.output
    if out.exists() and any(out.iterdir()) and not args.overwrite:
        p.error(f'Output directory not empty (choose fresh --output or use --overwrite): {out}')
    long_rows,wide_rows,missing=collect(args.a4,args.nonmixed,args.mixed)
    out.mkdir(parents=True,exist_ok=True)
    input_csv=out/'frozen_stage_observable_changes.csv'
    write_csv(input_csv,long_rows,LONG_FIELDS)
    write_csv(out/'running_comparisons.csv',wide_rows,WIDE_FIELDS)
    write_json(out/'missing_cells.json',missing)
    write_json(out/'provenance.json',{
        'format':'t3-16-simultaneous-eft1-original-stagefreeze-heatmaps-v1',
        'sources':{'A_-4':str(args.a4),'nonmixed13':str(args.nonmixed),
                   'mixed_BC':str(args.mixed)},
        'reference_case':'1111',
        'frozen_cases':STAGES,
        'bit_order':'UV,EFT1_renormalisable,EFT1_Wilson,SMEFT',
        'eft1_off_freezes_renormalisable_and_wilson':True,
        'uv_off_freezes_uv_renormalisable':True,
        'smeft_off_freezes_sm_and_weinberg':True,
        'percentage_formula':'100*(O_1111-O_X_off)/O_X_off',
        'n_models':16,
        'n_stage_comparisons':len(wide_rows),
        'n_finite_cells':len(long_rows),
        'n_missing_cells':len(missing),
        'gammaCC_Weyl_and_full_matching_certified':False,
        'physics_status':'experimental, provisional matching and non-HHSS potential',
        'no_new_rge_integration':True,
        'original_plotter_schema':str(ORIGINAL_PLOTTER),
        'colour_norm':'SymLogNorm',
        'colour_linear_threshold_percent':1.0,
        'colour_log_base':10,
    })
    print('Using original plotter data schema with signed log colours:',
          ORIGINAL_PLOTTER,flush=True)
    render_original_heatmap_layout(input_csv,out)
    print('Computed stage comparisons:',len(wide_rows),'/48',flush=True)
    print('Valid observable cells:',len(long_rows),'/240',flush=True)
    print('Missing observable cells:',len(missing),flush=True)
    print('Figures:',out/'plots'/'oscillation_masssq_heatmap_percent.png',flush=True)
    print('        ',out/'plots'/'oscillation_mixing_heatmap_percent.png',flush=True)
    return 0


if __name__=='__main__':
    raise SystemExit(main())
