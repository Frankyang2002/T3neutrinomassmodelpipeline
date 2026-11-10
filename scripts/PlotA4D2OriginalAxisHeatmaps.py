#!/usr/bin/env python3
"""Original thesis heatmap axes for A(-4) and D(-2), with provenance.

Reads the original 16-model independently-evolved stage-ablation CSV and uses
its original plotting functions without modifying them. Produces the original
16-row plots and two-row A(-4)/D(-2) zooms, with IDENTICAL UV/EFT1/SMEFT
columns and percentage convention. These are LEGACY figures, not corrected
benchmarks.

Optionally visualises corrected D(-2) EFT1/SMEFT mass-splitting switches in a
SEPARATE two-row plot with identical columns. Unsupported A(-4), D UV, mixing
cells remain blank, never zero. Mixing cannot be computed from mass singular
values alone. The corrected D source uses E1_S1/E0_S1 (EFT1) and E1_S1/E1_S0
(SMEFT). This is not interchangeable with the original full UV-stage study.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
from pathlib import Path
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.PlotT3OscillationHeatmaps import (
    read_rows, build_matrix, plot_heatmap, write_matrix_csv,
    MASSSQ_OBSERVABLES, MIXING_OBSERVABLES, MODEL_ORDER,
    format_cell_value
)

DEFAULT_LEGACY=ROOT/'Reports/output/thesis_qualitative/representation_running/csv/frozen_stage_observable_changes.csv'
DEFAULT_D=ROOT/'output/parallel_t3/D_m2_corrected_frozen_running/D_m2_frozen_running_observables.csv'
DEFAULT_OUT=ROOT/'output/parallel_t3/A4_Dm2_original_axis_heatmaps'


def original_heatmaps(rows, out):
    indices=[MODEL_ORDER.index('A_-4'), MODEL_ORDER.index('D_-2')]
    for defs,name,title in [
        (MASSSQ_OBSERVABLES,'masssq','Signed percentage changes in mass-squared differences'),
        (MIXING_OBSERVABLES,'mixing','Signed percentage changes in mixing-angle observables'),
    ]:
        all_matrix,models,columns=build_matrix(rows,defs)
        pair=all_matrix[indices,:]
        pair_names=[models[i] for i in indices]
        # Native original plotter, including 16-row overview and two-row zoom.
        # The colour scales are set from the FULL 16-model table and are not
        # independently rescaled in the zoom: see zoom_with_fixed_scale below.
        whole=out/f'original_all16_{name}.png'
        plot_heatmap(all_matrix,models,columns,whole,title,len(defs))
        write_matrix_csv(out/f'original_all16_{name}.csv',all_matrix,models,columns)
        zoom=out/f'original_A4_Dm2_{name}.png'
        zoom_with_fixed_scale(pair,pair_names,columns,zoom,title,len(defs),
                              max_abs=float(np.nanmax(np.abs(all_matrix))),
                              footer='Original 16-model study; A(-4) and D(-2) selected; NOT newly corrected.')
        write_matrix_csv(out/f'original_A4_Dm2_{name}.csv',pair,pair_names,columns)
        print('Saved:',whole)
        print('Saved:',zoom)


def zoom_with_fixed_scale(matrix,names,columns,path,title,group_size,max_abs,footer):
    fig,ax=plt.subplots(figsize=(max(10,1.05*len(columns)+3),4.0))
    z=max(float(max_abs),1e-12)
    cmap=plt.get_cmap('RdBu_r').copy()
    cmap.set_bad('#e8e8e8')
    img=ax.imshow(np.ma.masked_invalid(matrix),aspect='auto',cmap=cmap,
                  norm=mcolors.TwoSlopeNorm(vmin=-z,vcenter=0.,vmax=z))
    ax.set_title(title,pad=16)
    ax.set_xticks(np.arange(len(columns)))
    ax.set_xticklabels(columns)
    ax.set_yticks(np.arange(len(names)))
    ax.set_yticklabels(names)
    ax.set_xlabel('Running stage and oscillation observable')
    ax.set_ylabel('T3 model')
    for i in range(len(names)):
        for j in range(len(columns)):
            v=matrix[i,j]
            if math.isfinite(v):
                ax.text(j,i,format_cell_value(float(v)),ha='center',va='center',
                        fontsize=9,color='white' if abs(v)>.55*z else 'black')
            else:
                ax.text(j,i,'N/A',ha='center',va='center',fontsize=8,color='#666666')
    ax.set_xticks(np.arange(-.5,len(columns),1),minor=True)
    ax.set_yticks(np.arange(-.5,len(names),1),minor=True)
    ax.grid(which='minor',color='black',linestyle='-',lw=.6)
    ax.tick_params(which='minor',bottom=False,left=False)
    for j in range(group_size,len(columns),group_size):
        ax.axvline(j-.5,color='black',lw=1.8)
    cbar=fig.colorbar(img,ax=ax,fraction=.035,pad=.03)
    cbar.set_label(r'Signed percentage change: $100\times (O_{\mathrm{full}}-O_{X\mathrm{-off}})/O_{X\mathrm{-off}}$')
    fig.text(.07,.018,footer,fontsize=8)
    fig.tight_layout(rect=(0,.055,1,1))
    path.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(path,dpi=300,bbox_inches='tight')
    plt.close(fig)


def corrected_d_masssq(path):
    if not path.is_file():
        print('Corrected D(-2) four-scenario CSV not found; corrected-only overlay skipped:',path)
        return None
    with path.open(newline='',encoding='utf-8-sig') as f:
        reader=csv.DictReader(f)
        required={'scenario','delta_m21_sq_eV2','delta_m31_sq_eV2'}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f'Incorrect corrected D scenario CSV columns: {path}')
        data={r['scenario']:r for r in reader}
    if set(data)!={'E0_S0','E0_S1','E1_S0','E1_S1'}:
        raise ValueError('Corrected D input must have four scenarios E0_S0, E0_S1, E1_S0, E1_S1')
    # columns: UV dm21/dm31, EFT1 dm21/dm31, SMEFT dm21/dm31
    matrix=np.full((2,6),np.nan,float)
    for j,key in enumerate(('delta_m21_sq_eV2','delta_m31_sq_eV2')):
        for stage,top,bottom in ((1,'E1_S1','E0_S1'),(2,'E1_S1','E1_S0')):
            full=float(data[top][key]);base=float(data[bottom][key])
            if not all(map(math.isfinite,[full,base])) or base==0:
                continue
            matrix[1,stage*2+j]=100*(full-base)/base
    return matrix


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--legacy',type=Path,default=DEFAULT_LEGACY)
    ap.add_argument('--corrected-d',type=Path,default=DEFAULT_D)
    ap.add_argument('--output',type=Path,default=DEFAULT_OUT)
    args=ap.parse_args()
    out=args.output
    out.mkdir(parents=True,exist_ok=True)
    rows=read_rows(args.legacy)
    original_heatmaps(rows,out)
    corr=corrected_d_masssq(args.corrected_d)
    if corr is not None:
        # Colour scale kept from the original ALL16 mass-squared plot.
        legacy_masssq,_,cols=build_matrix(rows,MASSSQ_OBSERVABLES)
        zoom_with_fixed_scale(
            corr,['A_-4 [no corrected comparison]','D_-2 [corrected]'],cols,
            out/'corrected_Dm2_masssq_original_axes.png',
            'Corrected D(-2) signed mass-squared changes (A(-4) unavailable)',2,
            float(np.nanmax(np.abs(legacy_masssq))),
            'Same axes/colour scale as original. Different fixed-UV four-scenario definition; blank = not computed.'
        )
        write_matrix_csv(out/'corrected_Dm2_masssq_original_axes.csv',corr,
                         ['A_-4 [not computed]','D_-2'],cols)
        print('Saved:',out/'corrected_Dm2_masssq_original_axes.png')
    meta={
      'original_source':str(args.legacy),
      'original_model_selection':['A_-4','D_-2'],
      'original_heatmaps_status':'exact original CSV, original stage definitions and plotter',
      'corrected_D_status':'separate fixed-UV EFT1/SMEFT four-case comparison' if corr is not None else 'not found',
      'corrected_A_status':'not available in the matching four-case framework',
      'mixing_status':'only original legacy values; corrected D mixing unavailable',
      'UV_status':'corrected D UV switch not varied',
      'blank_values':'unsupported, never numerical zero',
      'original_color_scale':'full16 original data separately for masssq and mixing'
    }
    (out/'provenance.json').write_text(json.dumps(meta,indent=2)+'\n',encoding='utf-8')
    print('Saved:',out/'provenance.json')
    print('No original source data or production files modified.')

if __name__=='__main__':
    main()
