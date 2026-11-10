#!/usr/bin/env python3
"""Run 16 T3 representations with each saved D(-2), A(-4) fitted Yukawa donor.

Uses EXISTING repo pipeline and exact original thesis heatmap plotter.
Pipeline stage-freeze definition: 100*(FULL - X_OFF)/X_OFF, where one of
UV, EFT1 (including Wilson), SMEFT is switched off per counterfactual.

This is exploratory/generalised 32-case legacy methodology. It does NOT
propagate later D(-2) triplet-only/J corrections to all 16 representation
builders; do not represent outputs as gauge-certified corrected EFT.

Writes only new output folders, preserving existing runs and fits.
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from scripts.BuildT3GeneralisedThesisTables import load_reference

DEFAULT_TABLES=Path('output/parallel_t3/thesis_generalised_tables_v1')
DEFAULT_REF=Path('output/parallel_t3/A4_Dm2_targeted_fit01')
DEFAULT_RUNNING=Path('output/parallel_t3/thesis_two_benchmarks_all16_stagefreeze_v2')
DEFAULT_FIGURES=Path('output/parallel_t3/thesis_two_benchmarks_all16_four_heatmaps_v2')
REFS=('D_-2','A_-4')
MODELS=('A_-4','A_-2','A_+0','B_-3','B_-1','B_+1',
        'C_-3','C_-1','C_+1','D_-2','D_+0','D_+2',
        'E_-4','E_-2','E_+0','E_+2')
STAGES=('UV','EFT1','SMEFT')

def run(args):
    print('EXECUTE:', ' '.join(map(str,args)), flush=True)
    return subprocess.run([str(x) for x in args],cwd=ROOT).returncode

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--tables',type=Path,default=DEFAULT_TABLES)
    ap.add_argument('--reference-results',type=Path,default=DEFAULT_REF)
    ap.add_argument('--running-output',type=Path,default=DEFAULT_RUNNING)
    ap.add_argument('--figures-output',type=Path,default=DEFAULT_FIGURES)
    ap.add_argument('--plot-only',action='store_true',
                    help='Read an existing running_comparisons.csv; skip physics')
    ap.add_argument('--wilson-solver',choices=('RK2','DOP853'),default=None)
    ap.add_argument('--wilson-steps',type=int,default=None)
    ap.add_argument('--keep-trial-files',action='store_true')
    ap.add_argument('--debug',action='store_true')
    a=ap.parse_args()
    tables=ROOT/a.tables
    rdir=ROOT/a.reference_results
    runout=ROOT/a.running_output
    figout=ROOT/a.figures_output
    for folder in (tables,rdir):
        if not folder.is_dir():ap.error('Missing input directory: '+str(folder))
    if not (tables/'results.json').is_file() or not (tables/'provenance.json').is_file():
        ap.error('32-case baseline results.json and provenance.json required')
    provenance=json.loads((tables/'provenance.json').read_text(encoding='utf-8'))
    for ref in REFS:
        found=load_reference(ref,[rdir])
        stored=provenance.get('reference_benchmarks',{}).get(ref,{})
        # Normalize Windows/Linux path slashes; avoid dependence on native C: root.
        path=str(found['source']).replace('\\','/').lower()
        recorded=str(stored.get('result_json','')).replace('\\','/').lower()
        suffix='/' + ref.lower() + '/result.json'
        if not path.endswith(suffix) or not recorded.endswith(suffix):
            ap.error('Missing/inconsistent stored donor reference for '+ref)
        if abs(float(found['original_chi2'])-float(stored.get('original_chi2',-1)))>1e-10:
            ap.error('Different fitted donor point for '+ref+'; refuse silent mixing')
        print(f'{ref} donor: {found["source"]}; chi2={found["original_chi2"]:.10g}',flush=True)
    solver=a.wilson_solver or provenance['wilson_solver']
    steps=a.wilson_steps or int(provenance['wilson_steps'])
    if solver!=provenance['wilson_solver'] or steps!=int(provenance['wilson_steps']):
        ap.error('Solver and steps must equal baseline; '
                 f'expected {provenance["wilson_solver"]}, {provenance["wilson_steps"]}')
    print('Baseline model count: 16; donor references: D_-2 and A_-4.',flush=True)
    print('Per donor: FULL from existing 32-case table, plus 16 x 3 stage freezes.',flush=True)
    print('Total fresh counterfactual cases: 96, no Yukawa refitting.',flush=True)
    if not a.plot_only:
        if runout.exists():ap.error('Running output directory already exists: '+str(runout))
        cmd=[sys.executable,str(ROOT/'scripts/RunT3GeneralisedRunningComparisons.py'),
             '--tables',str(tables),'--reference-results',str(rdir),
             '--reference','both','--output',str(runout),
             '--wilson-solver',solver,'--wilson-steps',str(steps)]
        if a.keep_trial_files:cmd.append('--keep-trial-files')
        if a.debug:cmd.append('--debug')
        code=run(cmd)
        print('Pipeline exit code:',code,flush=True)
        if not (runout/'running_comparisons.csv').is_file():
            raise RuntimeError('No completed comparison CSV. Check output logs/blocked diagnostics.')
    if figout.exists() and any(figout.iterdir()):
        ap.error('Figures output directory already exists and is nonempty: '+str(figout))
    source=runout/'running_comparisons.csv'
    if not source.is_file():ap.error('Running CSV missing: '+str(source))
    with source.open(newline='',encoding='utf-8-sig') as f:
        rows=list(csv.DictReader(f))
    valid=sum(r.get('status')=='VALID_EXPLORATORY' for r in rows)
    print(f'Available stage comparisons: {valid}/{len(rows)}',flush=True)
    cmd=[sys.executable,str(ROOT/'scripts/PlotT3GeneralisedOriginalHeatmaps.py'),
         '--running',str(runout),'--reference','both','--output',str(figout)]
    if run(cmd):raise SystemExit('Plotting failed; stage CSV retained for inspection')
    targets=[]
    for ref in REFS:
        for kind in ('masssq','mixing'):
            name='oscillation_'+kind+'_heatmap_percent.png'
            path=figout/ref/'plots'/name
            if not path.is_file():
                raise RuntimeError('Expected plot missing: '+str(path))
            targets.append(str(path))
    print('\nFOUR GRAPHS:',flush=True)
    for target in targets:print(' ',target,flush=True)
    print('Axes: 16 models, UV/EFT1/SMEFT, original observable labels.',flush=True)
    print('Undefined/failed cells blank, never zero.',flush=True)
    print('Exploratory original 32-case methodology; not the newer corrected D(-2) local run.',flush=True)
    print('No previously saved fit or production file modified.',flush=True)

if __name__=='__main__':
    main()
